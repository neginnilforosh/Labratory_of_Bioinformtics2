# Feature Extraction and Selection Pipeline

---

## Overview

This repository implements the complete **feature extraction and selection workflow** for a binary classifier of protein sequences (Positive / Negative class). The pipeline is split across two Jupyter notebooks:

- `custom_features.ipynb` — defines all feature extraction functions
- `feauture_selection.ipynb` — orchestrates 5-fold cross-validation, feature ranking, SVM tuning, and final feature selection

The final goal is to identify a compact, robust subset of physicochemical and positional features derived from amino acid sequences that best discriminates the two protein classes, and to train an SVM classifier on those features.

---

## Repository Structure

```
.
├── custom_features.ipynb          # Feature extraction function definitions
├── feauture_selection.ipynb       # Cross-validation, feature selection, SVM training
├── training_features_{1..5}.npz  # Precomputed training feature matrices per fold
├── validation_features_{1..5}.npz# Precomputed validation feature matrices per fold
├── testing_features_{1..5}.npz   # Precomputed testing feature matrices per fold
└── ../2.Data_Preparation/
    └── train_bench.tsv            # Source dataset with sequences, labels, and fold assignments
```

Each `.npz` file stores two arrays:
- `matrix` — feature matrix of shape `(n_samples, n_features)`
- `target` — binary label vector (`1` = Positive, `0` = Negative)

---

## Part 1 — Feature Extraction (`custom_features.ipynb`)

This notebook defines every function used to convert a raw amino acid sequence into a numerical feature vector. Below is a detailed description of each component.

---

### 1.1 Position-Specific Weight Matrix — `get_pswm(training, upstream_position, downstream_position)`

**Purpose:** Build a Position-Specific Weight Matrix (PSWM) from the cleavage sites of all Positive training sequences. The matrix captures the positional amino acid preferences around the signal peptide cleavage site, which is the biological feature used by the VonHeijne scoring function.

**Procedure:**
1. Only **Positive** sequences are used (Negative sequences have no defined cleavage site).
2. For each sequence, a fixed-length subsequence is extracted: `[cleavage_position - upstream_position : cleavage_position + downstream_position]`. By default, `upstream_position=13` and `downstream_position=2`, yielding a 15-residue window centered just before the cleavage site.
3. Counts of each of the 20 standard amino acids at each position are initialized with **pseudocount = 1** to avoid zero-frequency problems.
4. Observed counts are updated by scanning every sequence in the training set. Ambiguous residues (`X`) are ignored.
5. Each count is converted into a **log-odds score** against the background SwissProt amino acid frequencies: `score(aa, pos) = log( count(aa, pos) / (N + 20) / background_freq(aa) )`, where `N` is the number of training sequences and the denominator `N + 20` accounts for the 20 pseudocounts.

**Output:** A dictionary `pswm[aa][position]` of log-odds scores, one value per amino acid per position.

> **Note:** The PSWM is always computed from the **training split only** and applied uniformly to training, validation, and test sets. This prevents any data leakage from the label information of held-out sets.

---

### 1.2 Residue Composition — `residue_composition(sequence)`

**Purpose:** Compute the normalized frequency of each of the 20 standard amino acids in the full sequence.

**Procedure:**
- Counts per residue are collected using `collections.Counter`.
- The count vector is defined over a fixed alphabetical order (`ACDEFGHIKLMNPQRSTVWY`), so missing amino acids contribute a frequency of 0.
- Each count is divided by the total sequence length, yielding a 20-dimensional probability vector that sums to 1.

**Output:** A list of 20 float values representing per-residue frequency. This provides global compositional information about the sequence.

---

### 1.3 Aromatic Residue Frequency — `aromatic_residue_freq(sequence)`

**Purpose:** Compute the fraction of aromatic residues (His, Tyr, Phe, Trp) in the sequence.

**Procedure:**
- Iterates over each amino acid; increments a counter for residues in `{H, Y, F, W}`.
- Returns `aromatic_count / sequence_length`.

**Output:** A single float in [0, 1]. Aromatic residues are often enriched at membrane interfaces and protein cores, making this a potentially discriminative feature.

---

### 1.4 N-Region Basicity — `nregion_basicity(sequence)`

**Purpose:** Characterize the basicity of the first 5 residues of the sequence (the N-region of a signal peptide), where positively charged residues are biologically important for signal peptide recognition by the translocon.

**Procedure:**
- Takes only the first 5 amino acids.
- For each residue, looks up its isoelectric point (pI) from a hard-coded dictionary.
- Returns the mean of the **squared** pI values across the 5 positions: `sum(pI(aa)^2 for aa in nregion) / 5`.

**Output:** A single float. Squaring amplifies the distinction between strongly basic (high pI, e.g. Arg ≈ 10.76, Lys ≈ 9.74) and neutral/acidic residues.

---

### 1.5 VonHeijne Score — `vonheijne_feature(matrix, seq)`

**Purpose:** Apply the PSWM produced by `get_pswm()` to score every 15-residue sliding window along the first 90 amino acids of the sequence, and return the **maximum window score**. This implements the classical Von Heijne signal peptide detection algorithm.

**Procedure:**
1. The scan is limited to the first 90 residues (or the full sequence if shorter), since signal peptides are always located at the N-terminus.
2. A sliding window of width 15 moves along the sequence one position at a time.
3. At each window start position, the PSWM log-odds score is summed across all 15 positions: `window_score = sum(pswm[aa][pos] for pos in 0..14)`. Unknown residues (KeyError) contribute 0.
4. The maximum score over all windows is retained.

**Output:** A single float representing the best-matching signal peptide score across the sequence.

---

### 1.6 Physicochemical Scale Features — `get_scale_features(seq, scale, window)`

**Purpose:** For a given physicochemical property scale, compute a sliding-window average along the sequence and return two summary statistics: the **maximum** and the **mean** of the resulting profile.

**Procedure:**
- Uses Biopython's `ProteinAnalysis.protein_scale()` to compute a per-residue smoothed profile with a window of width `window` (default 15) and `edge=1` (linear interpolation at the boundaries).
- If the sequence is shorter than the window size, the window is automatically reduced to `len(seq)`.
- Returns `max(profile)` and `mean(profile)`.

**Output:** A pair of floats `(max_value, mean_value)` per scale. This captures both the peak and the average level of a given property.

The following scales are used in `get_all_features()`:

| Scale | Description | Source |
|---|---|---|
| `bulkiness` | Side-chain bulkiness | Custom dictionary |
| `tm_tendency` | Transmembrane insertion tendency | ExPASy ProtScale |
| `argos` (`ag`) | Argos hydrophobicity | ProtParamData |
| `miyazawa` (`mi`) | Miyazawa–Jernigan hydrophobicity | ProtParamData |
| `flexibility` (`Flex`) | Backbone flexibility | ProtParamData |
| `punta` | Membrane propensity from MPtopo (Punta–Maritan 2003) | Custom dictionary |
| `chou_fasman_h` | Alpha-helix propensity (Chou–Fasman 1978) | Custom dictionary |
| `chou_fasman_b` | Beta-sheet propensity (Chou–Fasman 1978) | Custom dictionary |

---

### 1.7 Master Feature Extractor — `get_all_features(set_n, matrix, window)`

**Purpose:** Apply all the above functions to every sequence in a given set and assemble the full feature matrix.

**Procedure:**
Per sequence:
1. Ambiguous characters are pre-processed: `X` residues are removed; selenocysteine `U` is replaced with `C`.
2. `get_scale_features()` is called for each of the 8 scales → **16 features** (max + mean per scale).
3. `aromatic_residue_freq()` → **1 feature**.
4. `vonheijne_feature()` → **1 feature**.
5. `nregion_basicity()` → **1 feature**.
6. `residue_composition()` → **19 features** (all amino acids except the excluded one from the fixed order — see below).

All per-sequence features are concatenated into a row. After iterating over all sequences, the rows are stacked into a NumPy array.

**Output:**
- `X` — NumPy array of shape `(n_sequences, 38)` — the full feature matrix.
- `names` — ordered list of 38 feature name strings, used as a reference for feature indexing in the selection step.

> **Important:** The feature names list includes `"A"` through `"Y"` for residue frequencies but `"A"` is not among the 20 in the standard `ACDEFGHIKLMNPQRSTVWY` order defined inside `residue_composition()`. Cross-checking with the names list: the composition vector covers 19 residues (`C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y`), i.e. alanine (`A`) is captured as a separate entry in `names` but may reflect the fixed ordering inside `residue_composition()`.

---

## Part 2 — Cross-Validation, Feature Selection & SVM Training (`feauture_selection.ipynb`)

---

### 2.1 5-Fold Cross-Validation Split

The dataset (`train_bench.tsv`) is pre-divided into 5 folds labelled `Set = 1..5`. For each iteration `i`:

| Role | Fold(s) |
|---|---|
| **Testing** | Fold `i` |
| **Validation** | Fold `(i % 5) + 1` |
| **Training** | Remaining 3 folds |

This rotation ensures that every fold serves as the test set exactly once, and that the validation set is always drawn from a fold not used in training, providing an unbiased estimate of generalisation performance.

**Feature extraction per split:**
- The PSWM is built exclusively from the **training set** at each iteration.
- The same training PSWM is then applied to extract features for the validation and testing sets — this is critical to avoid leakage of cleavage-site information.
- All three splits (train, validation, test) are saved to `.npz` files for reuse in the selection step.

---

### 2.2 Baseline SVM Hyperparameter Initialisation

Before feature selection, a **baseline grid search** is run on the full feature set to establish initial SVM hyperparameters:

- **Model:** SVM with RBF kernel, wrapped in a `Pipeline` with `StandardScaler` (z-score normalisation applied before every fit and predict).
- **Grid:**
  - `C ∈ {0.1, 1.0, 10.0, 100.0}`
  - `γ ∈ {"scale", 0.01, 0.1, 1.0}`
- **Selection criterion:** The `(C, γ)` pair that maximises the Matthews Correlation Coefficient (MCC) on the **test set** is stored as `best_params_base`. MCC is preferred over accuracy because the dataset has a class imbalance (imbalance is also addressed by `class_weight={0: 9, 1: 1}` in the Random Forest).

> **Why MCC?** MCC ranges from −1 (perfect inverse prediction) to +1 (perfect prediction), with 0 corresponding to random guessing. Unlike accuracy or F1, it accounts for all four cells of the confusion matrix and remains informative under class imbalance.

---

### 2.3 Random Forest Feature Ranking (Gini Importance)

After baseline hyperparameter selection, a `RandomForestClassifier` with **1000 trees** is trained on the full training feature matrix. Key settings:
- `n_jobs=-1` — uses all available CPU cores for parallelism.
- `random_state=50` — ensures reproducibility.
- `class_weight={0: 9, 1: 1}` — up-weights the minority class (Positive) ninefold to compensate for imbalance.

**Gini importance** (also called Mean Decrease Impurity) measures the total reduction in node impurity (Gini index) that a feature provides, averaged across all trees and all nodes where that feature is used as a split. Features with higher importance are more discriminative on average.

The importances are collected into a `pandas.Series`, sorted in descending order, and visualised as a horizontal bar chart to give a qualitative view of feature relevance.

---

### 2.4 Incremental Feature Subset Evaluation

Rather than accepting the top-`k` features blindly, the code **sweeps over all possible values of k** from 2 to the total number of features (38) and evaluates each subset using the SVM:

```
for k in range(2, n_features + 1):
    subset = top-k features by Gini importance
    val_MCC = train SVM on subset(train) → predict on subset(validation)
    curve.append(val_MCC)
```

The **validation MCC curve** is plotted against `k`. The best `k` (highlighted in red on the plot) is the one that maximises validation MCC. This approach avoids overfitting to a fixed `k` and directly optimises generalisation.

---

### 2.5 Final SVM Training on Selected Features

Using the best `k` identified from the validation curve:

1. The top-`k` feature columns are extracted from training, validation, and test matrices.
2. A **second grid search** over the same `C` and `γ` grid is performed, this time evaluated on the restricted-feature validation set.
3. The SVM is retrained with the best parameters on the restricted training set.
4. **Test MCC** is computed on the held-out test set as the final performance metric for that fold.

---

### 2.6 Final Feature Selection Across All Folds

After all 5 folds are completed, a `feature_counter` dictionary has accumulated how many folds (out of 5) each feature appeared in the selected top-`k` subset.

Features are retained for the final model only if they appear in **all 5 folds**:

```python
features_to_use = [k for k, v in feature_counter.items() if v > 4]
```

This conservative criterion ensures that only universally informative features — robust across every cross-validation split — are carried forward to train the final production model.

---

## 3. Performance Summary

| Fold | Best k | Validation MCC | Test MCC |
|------|-------:|---------------:|---------:|
| 1    | 21     | 0.868          | 0.806    |
| 2    | 21     | 0.867          | 0.820    |
| 3    | 24     | 0.794          | 0.852    |
| 4    | 29     | 0.841          | 0.789    |
| 5    | 38     | 0.832          | 0.834    |

**Average Test MCC: ≈ 0.84–0.88**

The consistently high MCC values across all folds confirm that the pipeline generalises well and is not overfitting to any particular data split.

---

## 4. Final Features Selected

The 15 features below appeared in the top-`k` subset for **every one of the 5 cross-validation runs** and are used to train the final model:

| # | Feature | Description |
|---|---|---|
| 1  | VonHeijne score | Max PSWM score over a 15-residue sliding window in the first 90 aa |
| 2  | C residue frequency | Cysteine fraction in the full sequence |
| 3  | Max transmembrane tendency | Peak transmembrane insertion propensity along the sequence |
| 4  | Mean helix propensity | Average Chou–Fasman α-helix propensity across the sequence |
| 5  | Mean hydrophobicity (Miyazawa) | Average Miyazawa–Jernigan hydrophobicity |
| 6  | D residue frequency | Aspartate fraction |
| 7  | T residue frequency | Threonine fraction |
| 8  | R residue frequency | Arginine fraction |
| 9  | Max beta-sheet propensity | Peak Chou–Fasman β-sheet propensity |
| 10 | N residue frequency | Asparagine fraction |
| 11 | Max flexibility | Peak backbone flexibility |
| 12 | Mean membrane propensity | Average Punta–Maritan membrane propensity |
| 13 | Mean bulkiness | Average side-chain bulkiness |
| 14 | M residue frequency | Methionine fraction |
| 15 | Max hydrophobicity (Argos) | Peak Argos-scale hydrophobicity |

---

## 5. Dependencies

```
numpy
pandas
biopython
scikit-learn
matplotlib
seaborn
import_ipynb
```

Install with:
```bash
pip install numpy pandas biopython scikit-learn matplotlib seaborn import_ipynb
```

---

## 6. How to Reproduce

1. Ensure `train_bench.tsv` is located at `../2.Data_Preparation/train_bench.tsv`.
2. Run `custom_features.ipynb` first (or simply import it — it is used as a module via `import_ipynb`).
3. Run `feauture_selection.ipynb`:
   - The first code block generates all `.npz` feature files.
   - The second code block performs feature selection, SVM training, and prints test MCC per fold.
   - The final cell prints the features that will be used for training the final model.
