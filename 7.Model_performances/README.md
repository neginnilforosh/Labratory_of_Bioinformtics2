# Signal Peptide Prediction: Von Heijne · SVM · Deep Learning

A comparative study of three computational methods for **signal peptide (SP) classification** in protein sequences — from a classical scoring rule, through a kernel machine, to a hybrid deep neural network. Performance is evaluated on an independent benchmark set using the **Matthews Correlation Coefficient (MCC)** as the primary metric.

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Repository Structure](#2-repository-structure)
3. [Dataset](#3-dataset)
4. [Methods](#4-methods)
   - 4.1 [Von Heijne Method](#41-von-heijne-method)
   - 4.2 [Support Vector Machine (SVM)](#42-support-vector-machine-svm)
   - 4.3 [Deep Learning — SP-NN](#43-deep-learning--sp-nn)
5. [Hyperparameter Optimization](#5-hyperparameter-optimization)
   - 5.1 [SVM — Bayesian Search](#51-svm--bayesian-search)
   - 5.2 [SP-NN — Ray Tune Random Search](#52-sp-nn--ray-tune-random-search)
6. [Results](#6-results)
7. [Error Analysis](#7-error-analysis)
   - 7.1 [False Positive Analysis](#71-false-positive-analysis)
   - 7.2 [False Negative Analysis](#72-false-negative-analysis)
8. [Feature-wise Distributions](#8-feature-wise-distributions)
9. [Sequence Logos](#9-sequence-logos)
10. [Running the Code](#10-running-the-code)
11. [Dependencies](#11-dependencies)

---

## 1. Project Overview

Signal peptides are short N-terminal sequences (~15–30 amino acids) that direct nascent proteins to the secretory pathway. Their accurate prediction is a key task in computational biology. This project benchmarks three approaches of increasing complexity:

```
Sequence  ──►  Von Heijne Score  ──►  Binary label (SP / non-SP)
          ──►  Feature Engineering ──► SVM ──►  Binary label
          ──►  One-Hot Encoding   ──►  CNN+LSTM+MLP ──► Binary label
```

All models are trained on the same 5-fold cross-validation splits and evaluated on the same held-out benchmark set, enabling a fair, direct comparison.

---

## 2. Repository Structure

```
.
├── 2.Data_Preparation/
│   └── train_bench.tsv            # Master dataset (5 CV folds + Benchmark)
├── 5.Feature_Selection/
│   ├── custom_features.py         # Feature extraction utilities
│   ├── feauture_selection.ipynb   # Feature selection pipeline
│   ├── training_features_5.npz    # Pre-computed feature matrix (train)
│   ├── testing_features_5.npz     # Pre-computed feature matrix (test)
│   └── validation_features_5.npz  # Pre-computed feature matrix (validation)
├── vonHeijne/
│   ├── create_pswm.ipynb          # PSWM builder for Von Heijne scoring
│   └── validation_and_testing_vonheijne.ipynb
├── 6.SVM/
│   ├── hyperparameter_tuning.ipynb   # Bayesian search for SVM
│   ├── model_evaluation_and_plots.ipynb  # SVM + Von Heijne evaluation
│   ├── SignalPeptideSVM.pkl        # Saved SVM pipeline (after training)
│   └── benchmark_features.npz     # Benchmark features (after training)
└── 7.DL/
    ├── lstm.ipynb                  # SP-NN architecture + Ray Tune HPO
    ├── benchmark_test.ipynb        # Final DL evaluation on benchmark
    └── SignalPeptideLSTM.pt        # Saved TorchScript model (after training)
```

---

## 3. Dataset

The dataset is a tab-separated file (`train_bench.tsv`) with the following key columns:

| Column | Description |
|--------|-------------|
| `Sequence` | Full amino acid sequence |
| `Class` | `Positive` (has SP) or `Negative` (no SP) |
| `Set` | Fold assignment: `1`–`5` (cross-validation) or `Benchmark` |
| `SPStart` / `SPEnd` | Start and end positions of the signal peptide |
| `Kingdom` | Taxonomic kingdom (e.g. `Eukaryota`, `Viridiplantae`) |
| `OrganismName` | Species name |
| `HelixDomain` | Boolean flag for transmembrane helix presence |

The dataset is divided into **5 cross-validation folds** (used for training and model selection) and one **independent Benchmark set** used exclusively for final evaluation. No benchmark data is seen during training or hyperparameter tuning.

---

## 4. Methods

### 4.1 Von Heijne Method

The Von Heijne method is a classical, rule-based approach based on a **Position-Specific Weight Matrix (PSWM)**. It scores each sequence by computing a weighted sum over residues near the cleavage site, reflecting amino acid preference at each position.

**How it works:**

```
1. Build PSWM from Positive training sequences
   (window = 13 residues before + 2 after cleavage site)

2. For each test sequence:
   Slide the PSWM across the N-terminus
   → compute score at each position

3. Classify as Positive if max score > threshold
```

This method is interpretable and biologically motivated but cannot model complex non-linear relationships in the sequence.

---

### 4.2 Support Vector Machine (SVM)

The SVM operates on a **hand-crafted feature vector** of 18+ biochemical and sequence-derived descriptors computed per protein. These include:

| Feature Category | Examples |
|-----------------|---------|
| Amino acid composition | Frequencies of C, D, T, R, N, M |
| Hydrophobicity | Max/mean Miyazawa scale, Argos scale |
| Structural propensity | Chou-Fasman helix (mean), beta-sheet (max) |
| Transmembrane tendency | Max TM propensity (custom scale) |
| Flexibility | Max residue flexibility (Vihinen scale) |
| Membrane propensity | Max Punta-Maritan scale |
| Bulkiness | Mean residue bulkiness |
| Signal peptide signal | Von Heijne score (updated PSWM) |

An **RBF kernel SVM** is used with a `StandardScaler` preprocessing step, packaged into a `sklearn.Pipeline`. The Von Heijne score (feature index 17) is re-computed using the full training PSWM before fitting.

**Pipeline:**

```
Raw sequence
    │
    ▼
Feature extraction (18 descriptors)
    │
    ▼
StandardScaler (z-score normalization)
    │
    ▼
SVC (kernel=rbf, C=0.6758, γ=0.0350)
    │
    ▼
Binary prediction (SP / non-SP)
```

---

### 4.3 Deep Learning — SP-NN

SP-NN is a **hybrid CNN-LSTM-MLP** architecture that operates directly on raw sequences without manual feature engineering.

**Sequence encoding:**

All sequences are standardized to exactly **90 residues**:
- Shorter sequences are right-padded with `X` (mapped to an all-zero vector)
- Longer sequences are truncated to 90 (preserving the N-terminus)

Each position is encoded as a **21-dimensional one-hot vector** (20 standard amino acids + `X` for padding), producing a `[90 × 21]` matrix per sequence.

**Architecture:**

```
Input: [Batch, 90, 21]
       │
       ▼
┌─────────────────────────────────────────┐
│  1. CNN Block                           │
│     Conv1D(in=21, out=64, kernel=17)    │  ← kernel=17 ≈ hydrophobic core length
│     Output: [Batch, 90, 64]             │
└─────────────────────────────────────────┘
       │
       ▼
┌─────────────────────────────────────────┐
│  2. LSTM Block                          │
│     LSTM(64 → 128, layers=2)            │
│     Take last hidden state              │
│     BatchNorm1D(128)                    │
│     Output: [Batch, 128]                │
└─────────────────────────────────────────┘
       │
       ▼
┌─────────────────────────────────────────┐
│  3. MLP Head  (Hourglass topology)      │
│     Linear(128 → 256) + ReLU + Dropout  │
│     Linear(256 → 128) + ReLU + Dropout  │
│     Linear(128 →  64) + ReLU + Dropout  │
│     Linear( 64 → 1024)+ ReLU + Dropout  │
│     Linear(1024 →   1) + Sigmoid        │
└─────────────────────────────────────────┘
       │
       ▼
Output: probability ∈ [0, 1]  →  threshold 0.5  →  {SP, non-SP}
```

**Design rationale:**

- The CNN kernel size of **17** is deliberately chosen to match the typical length of the hydrophobic (H) region of a signal peptide, enabling detection of local hydrophobic motifs.
- The LSTM processes the CNN feature maps sequentially, capturing **long-range dependencies** along the sequence.
- The MLP head uses an **hourglass/flexible topology** (compress → expand) rather than a pure funnel, which the hyperparameter search found to be more expressive for final classification.
- **Gradient clipping** (`max_norm=1.0`) is applied during training to prevent the exploding gradient problem common in LSTMs.
- **Early stopping** (`patience=20`) prevents overfitting by saving the best-validation-MCC checkpoint.

**Training configuration (best trial):**

| Parameter | Value |
|-----------|-------|
| MLP hidden layers | `[256, 128, 64, 1024]` |
| LSTM hidden size | 128 |
| LSTM layers | 2 |
| Dropout | 0.498 |
| Learning rate | 2.86 × 10⁻⁴ |
| Batch size | 20 |
| Optimizer | Adam |
| Loss | BCELoss |
| Max epochs | 100 |
| Early stopping patience | 20 |

The model is exported as a **TorchScript** (`.pt`) file via JIT tracing, making it loadable without the original class definition.

---

## 5. Hyperparameter Optimization

### 5.1 SVM — Bayesian Search

Bayesian optimization (`BayesSearchCV`, 60 iterations, 5-fold CV) was used to tune two RBF kernel parameters:

| Hyperparameter | Search Range | Best Value | Description |
|---------------|-------------|-----------|-------------|
| `C` | [0.01, 100] log-uniform | **0.6758** | Regularization — trade-off between margin size and classification error |
| `γ (gamma)` | [0.01, 100] log-uniform | **0.0350** | RBF bandwidth — influence radius of single training examples |
| `kernel` | fixed: `rbf` | **rbf** | Radial Basis Function kernel |
| **Best CV MCC** | — | **0.8569** | Matthews Correlation Coefficient on validation folds |

Bayesian optimization is preferred over grid/random search because it uses a probabilistic surrogate model to propose the next configuration based on previous results, concentrating evaluations in the most promising region of the search space.

**Iteration 5 split used for SVM training:**

```
Set 1  ──► Validation (held out during tuning)
Set 2  ─┐
Set 3  ─┤──► Training  (x_training_conc)
Set 4  ─┤
Set 5  ─┘
Benchmark ──► Final evaluation only
```

---

### 5.2 SP-NN — Ray Tune Random Search

Random search with 15 trials and 5-fold cross-validation per trial was conducted using **Ray Tune**. The shared Object Store (`ray.put`) avoids redundant dataset serialization across worker processes.

**Search space:**

| Parameter | Type | Range |
|-----------|------|-------|
| `num_layers` | Categorical | {2, 3, 4, 5} |
| `hidden_sizes` | Sampled | Random subset from {1024, 512, 256, 128, 64, 32} |
| MLP topology | 50/50 | Funnel (descending) or Flexible (shuffled) |
| `dropout` | Uniform | [0.1, 0.5] |
| `lr` | Log-uniform | [1×10⁻⁴, 1×10⁻³] |
| `batch_size` | Categorical | {10, 20} |
| `num_lstm_layers` | Categorical | {1, 2} |
| `lstm_hidden_size` | Categorical | {64, 128} |

The metric optimized is the **mean MCC across 5 folds**. Ray Tune minimizes `−MCC` as its internal loss.

---

## 6. Results

All three models were evaluated on the same independent **Benchmark set** (never seen during training or tuning). The primary metric is **MCC**, which is robust to class imbalance.

### Final Benchmark Performance

| Method | MCC | Notes |
|--------|-----|-------|
| Von Heijne | **0.688** | Rule-based, no training required |
| SVM (RBF) | **0.808** | Engineered features + kernel machine |
| SP-NN (CNN+LSTM) | **0.902** | End-to-end deep learning |

```
MCC
1.00 │
0.90 │                              ████  0.902
0.80 │              ████  0.808
0.70 │  ████  0.688
0.60 │
     └───────────────────────────────────────
        Von Heijne     SVM        SP-NN (DL)
```

**Key observations:**

- **Von Heijne (0.688):** Captures the general statistical preference of residues near the cleavage site but lacks the capacity to model complex sequence patterns or longer-range interactions.

- **SVM (0.808):** A substantial improvement (+0.12 MCC) demonstrates that properly engineered biochemical features, combined with a kernel method, extract much more discriminative information than a simple position matrix.

- **SP-NN (0.902):** The deep model achieves the best performance (+0.094 over SVM), confirming that the CNN layer captures local hydrophobic motifs while the LSTM captures the sequential context across the N-terminus, together yielding a richer representation than any handcrafted feature set.

---

## 7. Error Analysis

Both the SVM (`model_evaluation_and_plots.ipynb`) and SP-NN (`benchmark_test.ipynb`) are subjected to the same structured error analysis to understand *why* misclassifications occur.

Each benchmark entry is assigned one of four labels:

| Label | Meaning |
|-------|---------|
| **TP** | True Positive — SP correctly detected |
| **TN** | True Negative — Non-SP correctly rejected |
| **FP** | False Positive — Non-SP incorrectly classified as SP |
| **FN** | False Negative — SP missed by the classifier |

---

### 7.1 False Positive Analysis

**Transmembrane (helix) domain proteins:**

A key hypothesis is that transmembrane (TM) proteins are over-represented among false positives because their N-terminal hydrophobic TM helix structurally resembles a signal peptide hydrophobic core.

The analysis computes the FPR separately for:
- All negative-class proteins (baseline FPR)
- Negative-class proteins with a known transmembrane domain (`HelixDomain == True`)

```
Example bar chart output (model_evaluation/HD_prediction_2.png):

  % of group
  100 │
   80 │  ████ no helix     ████ no helix
   60 │  ████              ████
   40 │  ████ helix        ████ helix
   20 │  ████              ████
    0 └──────────────────────────────
          FP              TN + FP
```

**Kingdom-level FPR pie chart** (`pieplot_kingdom_FPR.png`):

FPR is computed per kingdom (Eukaryota, Viridiplantae, Bacteria, Archaea) to test whether plant proteins show inflated FPR due to **chloroplast transit peptides** — N-terminal targeting sequences structurally similar to SPs.

**Viridiplantae deep-dive:**

Within plant entries, FPR is further split between *Arabidopsis thaliana* (the dominant plant species in the dataset) and all other plant species, to test whether apparent plant-FPR is driven by data imbalance rather than biology.

**FPR comparison summary (SVM):**

```
FPR (entire benchmark)          ~X.XX
FPR (transmembrane only)        ~X.XX  ← expected to be higher
FPR Von Heijne (all)            ~X.XX
FPR Von Heijne (TM only)        ~X.XX

Von Heijne FPR / SVM FPR  ≈  reported as % increase in-notebook
```

*(Exact values are printed at runtime — see cell output in `model_evaluation_and_plots.ipynb`)*

---

### 7.2 False Negative Analysis

**Signal peptide length distribution** (`Signalength_distribution.png`):

Histogram comparing SP length (`SPEnd − SPStart`) for TP vs FN sequences. Unusually short or long signal peptides may be harder for the model to score correctly.

```
Counts
 ▲
 │  ██ TP
 │  ██ ██
 │  ██ ██ FN
 │  ██ ██ ██
 └─────────────────► SP length (aa)
   5  10  15  20  25  30
```

**Hydrophobicity distribution** (`Hydrophobicity_distribution.png`, `Boxplot_hydrophobicity.png`):

The Miyazawa hydrophobicity scale is applied over the SP region. FN sequences may have an atypically low hydrophobicity in their H-region, making them harder to distinguish from non-SPs.

```
Mean hydrophobicity per group (barplot):
  TP  ─── typical SP hydrophobicity
  FN  ─── similar or slightly lower
  → "virtually no difference" noted in notebook comment
```

**Amino acid frequency comparison** (`AA_fequencies_3.png`, `figure9.png`):

Relative amino acid frequencies are compared across TP, FN, FP, and all Positives/Negatives. Residues are grouped and color-coded by chemical class:

| Color | Category | Residues |
|-------|----------|---------|
| 🔵 Blue | Nonpolar | G A V L I M P |
| 🟠 Orange | Aromatic | F Y W |
| 🟢 Green | Polar | S T N Q C |
| 🔴 Red | Positive charge | K R H |
| 🟣 Purple | Negative charge | D E |

Selected residues for spotlight analysis: **C, D, T, R, N, M** — residues whose frequencies diverge most between FP/FN and TP/TN.

**SP-region amino acid frequencies** (`AA_fequencies_3.png`):

The same analysis restricted to the SP subsequence only (`Sequence[SPStart−1 : SPEnd]`), isolating compositional signals within the cleavage region.

**Taxonomic composition** (`Taxa_composition_FP_FN.png`):

Three pie charts side-by-side:

```
  Total Positives        False Negatives       True Positives
  ┌──────────┐           ┌──────────┐          ┌──────────┐
  │  Euk 60% │           │  Euk ??% │          │  Euk ??% │
  │  Bac 25% │           │  Bac ??% │          │  Bac ??% │
  │  Pla 10% │           │  Pla ??% │          │  Pla ??% │
  │  Arc  5% │           │  Arc ??% │          │  Arc ??% │
  └──────────┘           └──────────┘          └──────────┘
```

Enrichment of plant sequences in FN would support the hypothesis that transit peptides generate atypical SP-like signals that confuse the classifier.

---

## 8. Feature-wise Distributions

For the SVM, the 15 most informative features (identified during feature selection) are visualised as **KDE plots** and **box plots** comparing their distributions across the four prediction groups (TP, TN, FP, FN). This pinpoints which biochemical properties drive misclassification.

| Feature | Scale/Source | Plot type |
|---------|-------------|-----------|
| Max transmembrane tendency | Custom TM scale | KDE + Boxplot |
| Chou-Fasman helix propensity (mean) | Chou & Fasman 1978 | Boxplot |
| Chou-Fasman β-sheet propensity (max) | Chou & Fasman 1978 | Boxplot |
| Residue flexibility (max) | Vihinen (Flex) | Boxplot |
| Membrane propensity (max) | Punta-Maritan 2003 | Boxplot |
| Mean residue bulkiness | Zimmerman 1968 | Boxplot |
| Hydrophobicity (max, Argos) | Argos scale (`ag`) | Boxplot |

**Key transmembrane propensity plot** (`KDE_TranmembranP.png`):

```
Density
  ▲
  │     TP / SP        FP          TN          FN
  │    ████           ████        ████
  │   ██████          ████      ████████
  │  ████████        ██████    ██████████
  └──────────────────────────────────────► Max TM tendency
    -4  -2   0   2   4
```

False positives (non-SP predicted as SP) are expected to show elevated transmembrane tendency scores, overlapping with the TP distribution — confirming that TM helices are the primary confusion factor.

---

## 9. Sequence Logos

Sequence logos are generated using **logomaker** for all four prediction groups (TP, TN, FP, FN). Two metrics are used:

**Information content (bits)** — highlights strongly conserved positions:
```
Bits
2.0 │  [TP / FN logos: y-axis 0–2.0 bits]
    │  Position weight proportional to conservation
0   └─────────────────────────────────────
    -13   -10   -5   cleavage  +4

1.0 │  [TN / FP logos: y-axis 0–1.0 bits]
    │  Lower conservation expected for non-SP sequences
0   └─────────────────────────────────────
```

**Probability** — shows residue frequency at each position without information-content normalization.

Positive sequences are aligned **around the cleavage site** (`SPEnd − 13` to `SPEnd + 4`), capturing the characteristic AXA cleavage motif. Negative sequences are aligned from position `1:18` of the N-terminus.

**Saved figures:**

| File | Content |
|------|---------|
| `figure12A-1.png` / `Sequence_Logo-1.png` | True Positives |
| `figure12A-2.png` / `Sequence_Logo-2.png` | False Negatives |
| `figure12A-3.png` / `Sequence_Logo-3.png` | True Negatives |
| `figure12A-4.png` / `Sequence_Logo-4.png` | False Positives |

---

## 10. Running the Code

Run the notebooks in this exact order:

### Step 1 — `lstm.ipynb`
Defines the SP-NN architecture, utility functions, and runs Ray Tune hyperparameter search. The best config is printed at the end and hardcoded into `benchmark_test.ipynb`.

```bash
# Resource note: set gpu=0 if no GPU is available
resources_per_trial={"cpu": 4, "gpu": 0}
```

### Step 2 — `hyperparameter_tuning.ipynb`
Trains the SVM pipeline with Bayesian search. Outputs:
- `SignalPeptideSVM.pkl` — trained SVM model
- `benchmark_features.npz` — encoded benchmark feature matrix

### Step 3 — `benchmark_test.ipynb`
Retrains the SP-NN with the best config on all 4 folds, exports `SignalPeptideLSTM.pt`, and runs the full deep learning error analysis.

### Step 4 — `model_evaluation_and_plots.ipynb`
Loads the SVM model and benchmark features, runs the Von Heijne baseline, and produces all comparison plots.

---

## 11. Dependencies

```
Python >= 3.9

# Deep Learning
torch >= 2.0
ray[tune] >= 2.0

# Classical ML
scikit-learn >= 1.2
scikit-optimize          # BayesSearchCV

# Bioinformatics
biopython                # ProteinAnalysis, ProtParamData

# Visualization
matplotlib
seaborn
logomaker

# Utilities
numpy
pandas
importnb                 # Import .ipynb as modules
pickle                   # Standard library
```

Install all at once:

```bash
pip install torch scikit-learn scikit-optimize biopython \
            matplotlib seaborn logomaker numpy pandas \
            "ray[tune]" importnb
```

---

## Conclusion

This study demonstrates a clear and consistent performance hierarchy:

```
Von Heijne  ──(+0.12 MCC)──►  SVM  ──(+0.09 MCC)──►  SP-NN
   0.688                      0.808                    0.902
```

Each step up in model complexity yields meaningful gains: feature engineering over a simple PSWM, and learned representations over hand-crafted features. The MCC of **0.902** for the CNN+LSTM model confirms that deep learning can reliably learn both local motifs (hydrophobic core) and global sequence context from raw protein sequences — without any domain-expert feature engineering.

The error analyses consistently point to **transmembrane domain proteins** as the primary source of false positives, and to **atypical signal peptides** (unusual length, composition, or taxonomy) as the main source of false negatives — providing a roadmap for future improvements.
