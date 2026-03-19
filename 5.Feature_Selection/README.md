# 🔧 Step 5 — Feature Extraction and Selection

> **What this folder does:** Converts raw protein sequences into numerical feature vectors, selects the most informative features using a Random Forest, and trains a Support Vector Machine (SVM) classifier — all within a rigorous 5-fold cross-validation framework.

---

## 🗂️ Files in this Folder

| File | What it does |
|------|-------------|
| `custom_features.ipynb` | Defines all functions that turn a protein sequence into a feature vector |
| `feauture_selection.ipynb` | Runs the full cross-validation, feature ranking, feature selection, and SVM training |

### Output files (created when notebooks are run)

| File pattern | Contents |
|------|----------|
| `training_features_i.npz` | Feature matrix and labels for the training set at fold `i` |
| `testing_features_i.npz` | Feature matrix and labels for the test set at fold `i` |
| `validation_features_i.npz` | Feature matrix and labels for the validation set at fold `i` |

Each `.npz` file stores two arrays:
- `matrix` — shape `(n_samples, n_features)` — the numerical feature matrix
- `target` — shape `(n_samples,)` — binary labels (1 = positive, 0 = negative)

---

## 💡 Why Feature Engineering?

Unlike the deep learning model (Step 6) which learns from raw sequences, the SVM needs numerical input. The idea is to encode biologically meaningful properties of each protein sequence as numbers — so the SVM can detect patterns that distinguish signal peptides from non-signal peptides.

```
Protein sequence (text)
        │
        ▼
Feature extraction functions
        │
        ▼
38-dimensional numerical vector per protein
        │
        ▼
SVM classifier
```

---

## 🧪 Part 1 — Feature Extraction (`custom_features.ipynb`)

### Feature 1: Von Heijne Score (`vonheijne_feature`)

**What it measures:** How much the beginning of this protein looks like a canonical signal peptide, according to the PSWM built in Step 4.

**How it works:**
- Slide a 15-residue window across the first 90 amino acids of the sequence
- At each window position, sum the PSWM log-odds scores for each residue
- Return the **maximum score** across all windows

```python
# Scan only first 90 aa (SP is always at N-terminus)
end = min(90, len(seq))
for window_start in range(0, end - 15 + 1):
    window_score = sum(PSWM[pos][aa] for pos, aa in enumerate(window))
    max_score = max(max_score, window_score)
```

> ⚠️ The PSWM is always built from the **training set only** and applied to all sets. This prevents data leakage.

---

### Feature 2: Residue Composition (`residue_composition`)

**What it measures:** What fraction of the sequence is each of the 20 standard amino acids.

**How it works:** Counts each amino acid and divides by total length. Uses a fixed alphabetical order (`ACDEFGHIKLMNPQRSTVWY`) so the output vector is always comparable across sequences.

```python
aa_order = "ACDEFGHIKLMNPQRSTVWY"
counts = Counter(sequence)
frequencies = [counts.get(aa, 0) / len(sequence) for aa in aa_order]
```

**Output:** 19 features (one per amino acid — one is excluded due to collinearity; all 20 frequencies sum to 1).

---

### Feature 3: Aromatic Residue Frequency (`aromatic_residue_freq`)

**What it measures:** The fraction of aromatic residues (H, Y, F, W) in the full sequence.

**Why it matters:** Aromatic residues are often found at membrane interfaces and protein cores. Signal peptides tend to have fewer aromatics than average.

**Output:** 1 feature — a single float in [0, 1].

---

### Feature 4: N-Region Basicity (`nregion_basicity`)

**What it measures:** How basic the first 5 residues of the sequence are.

**Why it matters:** Signal peptides have a characteristic **three-region structure**:
- **N-region:** positively charged residues (K, R) — important for membrane recognition by the translocon
- **H-region:** hydrophobic core — the longest part
- **C-region:** cleavage site (AXA motif)

The N-region basicity feature specifically captures the positive charge at the very beginning of potential SPs.

**How it works:**
```python
nregion = sequence[:5]
basicity = sum(isoelectric_point[aa]**2 for aa in nregion) / 5
```

Squaring the isoelectric point amplifies the difference between strongly basic residues (Arg pI≈10.76, Lys pI≈9.74) and neutral/acidic ones.

**Output:** 1 feature.

---

### Feature 5–12: Physicochemical Scale Features (`get_scale_features`)

**What it measures:** For each of 8 biophysical scales, compute a sliding-window average along the sequence and return two summary statistics: the **maximum** and the **mean** of the resulting profile.

```python
# Biopython computes smoothed per-residue scores
profile = ProteinAnalysis(seq).protein_scale(scale, window=15, edge=1)
max_value  = max(profile)
mean_value = mean(profile)
```

**The 8 scales used:**

| Scale | Description | Why it matters for SP detection |
|-------|-------------|--------------------------------|
| `bulkiness` | Side-chain bulkiness (Zimmerman) | Large bulky residues (L, I, V) are common in the SP hydrophobic core |
| `tm_tendency` | Transmembrane insertion tendency | High TM tendency = hydrophobic stretch resembling TM helix or SP core |
| `ag` (Argos) | Argos hydrophobicity | Strong discriminator: SPs are hydrophobic, non-SPs are not |
| `mi` (Miyazawa) | Miyazawa–Jernigan hydrophobicity | Similar to Argos but uses a different energy model |
| `Flex` (Flexibility) | Backbone flexibility (Vihinen) | SPs tend to be more structured; high flexibility = atypical SP |
| `punta` | Membrane propensity (Punta-Maritan 2003) | Distinguishes TM helices from SP hydrophobic cores |
| `chou_fasman_h` | α-helix propensity (Chou-Fasman) | SP hydrophobic core forms an amphipathic helix |
| `chou_fasman_b` | β-sheet propensity (Chou-Fasman) | TM β-barrel proteins have high β-sheet propensity vs. SPs |

**Output:** 2 features per scale × 8 scales = **16 features**.

---

### Master Extractor: `get_all_features`

This function applies all the above to every sequence in a set and assembles the complete feature matrix:

| Source | Features |
|--------|:--------:|
| 8 physicochemical scales × 2 stats | 16 |
| Aromatic frequency | 1 |
| Von Heijne PSWM score | 1 |
| N-region basicity | 1 |
| Residue composition | 19 |
| **Total** | **38** |

Before extraction, ambiguous characters are cleaned: `X` is removed, `U` (selenocysteine) is replaced with `C`.

---

## 🤖 Part 2 — Feature Selection and SVM Training (`feauture_selection.ipynb`)

### Step 1: 5-Fold Cross-Validation Setup

| Fold role | Which folds |
|-----------|------------|
| Testing | Fold `i` |
| Validation | Fold `(i % 5) + 1` |
| Training | The remaining 3 folds |

The PSWM is rebuilt from the training folds at each iteration, and features are extracted using that training-specific PSWM. This is crucial — using a PSWM built from data the test set has seen would be cheating.

---

### Step 2: Baseline SVM Grid Search

Before feature selection, a quick grid search finds reasonable SVM starting parameters:

| Parameter | Grid |
|-----------|------|
| `C` | {0.1, 1.0, 10.0, 100.0} |
| `γ` | {"scale", 0.01, 0.1, 1.0} |

The best `(C, γ)` pair is selected by maximising **MCC on the test set** using the full 38-feature matrix.

---

### Step 3: Random Forest Feature Ranking

A Random Forest with **1,000 trees** is trained to rank features by their importance:

```python
rf = RandomForestClassifier(
    n_estimators=1000,
    n_jobs=-1,
    random_state=50,
    class_weight={0: 9, 1: 1}  # upweight the minority positive class
)
rf.fit(X_train, y_train)
importances = pd.Series(rf.feature_importances_, index=feature_names)
importances_sorted = importances.sort_values(ascending=False)
```

**Gini importance** (Mean Decrease Impurity) measures how much each feature reduces node impurity when it is used as a split point, averaged across all trees.

The `class_weight={0: 9, 1: 1}` setting makes the forest pay 9× more attention to correctly classifying positive (signal peptide) samples — compensating for the class imbalance where negatives outnumber positives roughly 9:1.

---

### Step 4: Incremental Feature Subset Evaluation

Rather than blindly taking the top-k features from the Random Forest ranking, the code evaluates **every possible subset size** from 2 to 38:

```python
for k in range(2, n_features + 1):
    top_k_features = importances_sorted.index[:k]
    val_MCC = train_and_evaluate(SVM, top_k_features, train, validation)
    validation_mcc_curve.append(val_MCC)

best_k = argmax(validation_mcc_curve)
```

The resulting curve shows how validation MCC improves as more features are added. The best `k` (highlighted on the plot) is the exact number that maximises generalisation — not too few (underfitting) and not too many (overfitting on noise).

---

### Step 5: Final SVM on Selected Features

With the best `k` identified:

1. Extract the top-k columns from train/validation/test matrices
2. Run a second grid search on the reduced feature set
3. Retrain the SVM with the best parameters
4. Report the **test MCC** on the held-out fold

---

### Step 6: Feature Consensus Across All Folds

After all 5 folds complete, the code counts how many folds each feature appeared in the top-k selection:

```python
features_to_use = [f for f, count in feature_counter.items() if count > 4]
```

Only features that appear in the top-k subset for **all 5 folds** are kept. This conservative requirement ensures the final feature set is universally informative, not specific to any single data split.

---

## 📊 Cross-Validation Results

| Fold | Best k | Validation MCC | Test MCC |
|:----:|:------:|:--------------:|:--------:|
| 1 | 21 | 0.868 | 0.806 |
| 2 | 21 | 0.867 | 0.820 |
| 3 | 24 | 0.794 | 0.852 |
| 4 | 29 | 0.841 | 0.789 |
| 5 | 38 | 0.832 | 0.834 |
| **Average** | — | **0.840** | **0.820** |

The consistently high MCC values across all 5 folds confirm good generalisation without overfitting.

---

## 🏆 Final 15 Features Selected (Appeared in All 5 Folds)

| # | Feature | Why it matters |
|---|---------|----------------|
| 1 | **Von Heijne score** | Best single-feature predictor — directly measures SP-like sequence motif |
| 2 | **C (Cys) frequency** | Cysteine is rare in signal peptides; high C suggests non-SP |
| 3 | **Max TM tendency** | Transmembrane-like stretch = either SP hydrophobic core or TM helix |
| 4 | **Mean helix propensity (Chou-Fasman)** | SP hydrophobic core forms an α-helix |
| 5 | **Mean Miyazawa hydrophobicity** | SP must be hydrophobic to interact with the translocon |
| 6 | **D (Asp) frequency** | Negatively charged residues disfavour SP hydrophobic core |
| 7 | **T (Thr) frequency** | Polar residues in the core weaken the hydrophobic signal |
| 8 | **R (Arg) frequency** | Arg is enriched in the N-region of SPs |
| 9 | **Max β-sheet propensity (Chou-Fasman)** | Distinguishes TM β-barrels from SP helices |
| 10 | **N (Asn) frequency** | Polar residue rare in the hydrophobic core |
| 11 | **Max flexibility** | Disordered/flexible sequences are atypical SPs |
| 12 | **Mean membrane propensity (Punta)** | Distinguishes SP cores from TM helices |
| 13 | **Mean bulkiness** | Large hydrophobic side chains (L, I, V) dominate SP cores |
| 14 | **M (Met) frequency** | Met enriched in SP cores (and at initiator Met position) |
| 15 | **Max Argos hydrophobicity** | Peak hydrophobicity — strong SP predictor |
