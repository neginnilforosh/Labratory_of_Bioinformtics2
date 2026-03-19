# 🧮 Step 4 — Von Heijne Method

> **What this folder does:** Implements the classical Von Heijne algorithm for signal peptide prediction using a Position-Specific Weight Matrix (PSWM). This is the simplest and most interpretable baseline model in the project.

---

## 🗂️ Files in this Folder

| File | What it does |
|------|-------------|
| `create_pswm.ipynb` | Builds the PSWM from positive training sequences and visualises it as a heatmap |
| `validation_and_testing_vonheijne.ipynb` | Scores sequences using the PSWM, tunes the classification threshold, and evaluates performance |
| `vonheijne.ipynb` | Full end-to-end pipeline combining both notebooks for the complete 5-fold cross-validation |

### Output figures (saved to `Plots/`)

| Figure | Description |
|--------|-------------|
| `Plots/Heatmap_PSWM_2.png` | Heatmap of the PSWM — rows = positions, columns = amino acids |
| `Plots/pr_curve_and_heatmap_1.png` through `_5.png` | Precision-Recall curve + PSWM heatmap for each of the 5 cross-validation folds |
| `Plots/all_pr_curves.png` | All 5 PR curves overlaid on one plot for comparison |

---

## 💡 What is the Von Heijne Method?

The Von Heijne method (1986) is a **scoring rule** that asks: *does the beginning of this protein sequence look like a signal peptide?*

It answers this question by building a statistical model of what amino acids are typically found at each position around a signal peptide cleavage site, then giving every candidate sequence a score based on how well it matches that model.

It requires **no iterative training** and is fully transparent — you can read the matrix directly and understand why a protein received a high or low score.

---

## 📐 Step 1 — Building the PSWM (`create_pswm.ipynb`)

The PSWM captures positional amino acid preferences by analysing all known signal peptides in the training set.

### What is extracted from each positive protein?

For each signal peptide, a **15-residue window** is extracted centred at the cleavage site:

```
positions: [cleavage − 13]  ...  [cleavage − 1]  |  [cleavage + 1]  [cleavage + 2]
index:           0                    12              13                  14
```

This 15-residue window spans the last 13 residues of the signal peptide plus the first 2 residues after the cleavage site.

### How the matrix is computed

**1. Count amino acids at each position with pseudocounts**

```python
# Pseudocount = 1 to avoid zero-frequency problems
count(position, aa) = observed_count + 1

# Normalise
frequency(position, aa) = count(position, aa) / (N + 20)
# N = number of training sequences; 20 = one pseudocount per amino acid
```

**2. Compute log-odds score against SwissProt background**

```python
PSWM(position, aa) = log( frequency(position, aa) / background_frequency(aa) )
```

A **positive** PSWM score at a position means that amino acid appears more often in signal peptides than in proteins generally. A **negative** score means it appears less often than expected.

### Visualisation

The heatmap `Plots/Heatmap_PSWM_2.png` shows the full PSWM. Amino acids are ordered by chemical class so patterns are easier to see:

```
Nonpolar  │ G A V L I M P
Aromatic  │ F Y W
Polar     │ S T C N Q
Negative  │ D E
Positive  │ K R H
```

Bright warm colours at a position indicate that amino acid is strongly preferred in signal peptides there. The hydrophobic region (positions −10 to −4) should show strong positive scores for L, A, V, I.

---

## 🔎 Step 2 — Scoring Sequences (`validation_and_testing_vonheijne.ipynb`)

### Sliding window scoring

For each test sequence, the algorithm slides a 15-residue window from position 1 to position 75 (the first 90 residues, since SPs are always at the N-terminus):

```python
for window_start in range(0, min(90, len(seq)) - 15 + 1):
    window = seq[window_start : window_start + 15]
    score  = sum(PSWM[pos][aa] for pos, aa in enumerate(window))
    # Handle unknown residues (key error) → contribute 0

max_score = max(all_window_scores)
```

The **maximum score** over all window positions is kept as the sequence's signal peptide score. A high maximum score means at least one 15-residue stretch looks very much like a canonical signal peptide.

### Classification threshold

The raw score is a real number. To turn it into a binary prediction (SP / non-SP), a threshold is needed:

```python
# Compute Precision-Recall curve on the VALIDATION set
precision, recall, thresholds = precision_recall_curve(val_labels, val_scores)

# Find the threshold that maximises F1 score
F1 = 2 * (precision * recall) / (precision + recall)
best_threshold = thresholds[argmax(F1)]
```

Any sequence with a score ≥ best_threshold is predicted as positive (has a signal peptide).

The PR curve and the chosen threshold are saved per fold (`Plots/pr_curve_and_heatmap_1.png` … `_5.png`).

---

## 🔁 Step 3 — 5-Fold Cross-Validation

The full pipeline is repeated 5 times:

| Iteration | Training | Validation | Testing |
|:---------:|:--------:|:----------:|:-------:|
| 1 | Folds 2, 3, 4 | Fold 5 | Fold 1 |
| 2 | Folds 3, 4, 5 | Fold 1 | Fold 2 |
| 3 | Folds 4, 5, 1 | Fold 2 | Fold 3 |
| 4 | Folds 5, 1, 2 | Fold 3 | Fold 4 |
| 5 | Folds 1, 2, 3 | Fold 4 | Fold 5 |

At each iteration:
1. A new PSWM is built from the **training folds only**
2. The optimal threshold is selected on the **validation fold**
3. Performance is measured on the **test fold**

This prevents any fold from influencing the threshold that is used to evaluate it.

---

## 📊 Performance Results

| Metric | Mean ± Std |
|--------|:----------:|
| **MCC** | **0.663 ± 0.016** |
| Precision | 0.934 ± 0.003 |
| Accuracy | 0.691 ± 0.014 |
| Sensitivity (Recall) | 0.711 ± 0.029 |

### Understanding the metrics

- **MCC (Matthews Correlation Coefficient):** The primary metric. Ranges from −1 (completely wrong) to +1 (perfect). 0 = random guessing. An MCC of 0.663 is solid for a rule-based method with no training parameters. Used because the dataset is class-imbalanced (far more negatives than positives).
- **Precision 0.934:** When the model says "this is a signal peptide", it is correct 93.4% of the time — very few false alarms.
- **Sensitivity 0.711:** The model detects 71.1% of real signal peptides — it misses about 29% (false negatives).

The high precision but lower sensitivity pattern is expected for the Von Heijne method: it is conservative (rarely fires false alarms) but misses atypical signal peptides that don't fit the canonical AXA motif perfectly.

---

## ⚖️ Strengths and Limitations

| Strengths | Limitations |
|-----------|-------------|
| Fully interpretable — you can read the matrix | Cannot learn non-linear patterns |
| No iterative training needed | Performance limited by the hand-crafted window |
| Biologically meaningful scoring | Does not use global sequence context |
| Fast — a single matrix-vector multiplication per window | Sensitive to unusual SP compositions |


