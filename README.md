<div align="center">

# 🧬 Signal Peptide Prediction
## Laboratory of Bioinformatics II — Module 2

</div>


## 🗺️ Quick Navigation

| # | Folder | What it does | README |
|---|--------|-------------|--------|
| 1 | [`1.Data_Collection/`](#1-data-collection) | Downloads proteins from UniProt API | [→ README](1.Data_Collection/README.md) |
| 2 | [`2.Data_Preparation/`](#2-data-preparation) | Removes redundancy, splits train/benchmark | [→ README](2.Data_Preparation/README.md) |
| 3 | [`3.Data_Analysis/`](#3-data-analysis) | Explores dataset properties with plots | [→ README](3.Data_Analysis/README.md) |
| 4 | [`4.vonHeijne/`](#4-von-heijne-method) | Classical PSWM-based baseline method | [→ README](4.vonHeijne/README.md) |
| 5 | [`5.Feature_Selection/`](#5-feature-extraction--selection) | Builds 38 features, selects best 15, trains SVM | [→ README](5.Feature_Selection/README.md) |
| 6 | [`6.Deep_Learning/`](#6-deep-learning--sp-nn) | Hybrid CNN+LSTM neural network | [→ README](6.Deep_Learning/README.md) |
| 7 | [`7.Model_performances/`](#7-model-performance-comparison) | Final benchmark evaluation and error analysis | [→ README](7.Model_performances/README.md) |

---

## 🔬 What Is This Project About?

### The biological problem

Every cell is a miniature city, and proteins are its workers. To reach the right compartment — the membrane, the extracellular space, the endoplasmic reticulum — many proteins carry a **molecular address label** at their very tip: the **signal peptide**.

A signal peptide is a short stretch of ~15–30 amino acids at the N-terminus of a protein. It acts as a zip code, routing the protein through the **secretory pathway** (ER → Golgi → membrane or extracellular space). Once the protein arrives, a dedicated enzyme (signal peptidase) clips the signal peptide off — the protein carries the label only for the journey, not the destination.

```
Signal peptide structure:

[++++ n-region ++++][━━━━━━ h-region (hydrophobic core) ━━━━━━][AXA c-region]
 Positively charged   Interacts with membrane translocon         Cleavage site
 (K, R enriched)      (L, A, V, I enriched)                     (small neutral aa)
                                                                     ↑
                                                              Cleaved here
```

Knowing whether a protein has a signal peptide — and where the cleavage site is — is crucial for:
- Predicting subcellular localisation
- Identifying drug targets (secreted and membrane proteins are highly druggable)
- Designing vaccines (exposed bacterial surface proteins)
- Functional annotation of unannotated genomes

Experimental determination is slow and expensive. **Computational prediction** is the scalable alternative.

### The computational challenge

Signal peptide prediction is a **binary sequence classification** problem:

```
Input:  Protein amino acid sequence  (e.g. MLPGLALLLLAAWTARALEV...)
Output: Positive (has signal peptide) or Negative (does not)
```

The challenge is that the signal peptide hydrophobic core closely resembles transmembrane helices — a major source of false positives. And some signal peptides are atypical in length or composition — a major source of false negatives.

---

## 🎯 Project Goals

This project was developed as part of the **Laboratory of Bioinformatics II** course. The stated objectives were:

1. **Retrieve** a high-quality, non-redundant dataset from UniProtKB
2. **Preprocess** the data for cross-validation and independent benchmarking
3. **Analyse** the dataset statistically and visually
4. **Implement** the classical Von Heijne (1986) motif-based method
5. **Engineer features** and train a Support Vector Machine (SVM) classifier
6. **Evaluate** both methods using 5-fold cross-validation and a blind benchmark set
7. **Extend** the project with a hybrid deep learning architecture (SP-NN) — optional extension
8. **Analyse errors** biologically to understand model failure modes

---

## 📐 The Full Pipeline

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        COMPLETE PROJECT PIPELINE                            │
└─────────────────────────────────────────────────────────────────────────────┘

  UniProt API
      │
      ▼
┌─────────────────────┐
│  1. DATA COLLECTION │  get_dataset_pos.ipynb   → positive_dataset.tsv/.fasta
│                     │  get_dataset_neg.ipynb   → negative_dataset.tsv/.fasta
│  2,932 positives    │
│  20,615 negatives   │
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  2. DATA PREP       │  MMseqs2 clustering (30% identity, 40% coverage)
│                     │  → 1,093 non-redundant positives
│  Redundancy removal │  → 8,934 non-redundant negatives
│  80/20 split        │  prepare_datasets.ipynb → train_bench.tsv
│  5-fold CV          │
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  3. DATA ANALYSIS   │  plots.ipynb
│                     │  Length distributions · AA composition vs SwissProt
│  Exploratory plots  │  Taxonomy · Sequence logos at cleavage site
└────────┬────────────┘
         │
         ├──────────────────────────────────────┐
         │                                      │
         ▼                                      ▼
┌─────────────────────┐              ┌──────────────────────┐
│  4. VON HEIJNE      │              │  5. FEATURE SELECT.  │
│                     │              │                      │
│  PSWM from training │              │  38 features per seq │
│  Sliding window     │              │  Random Forest rank  │
│  F1-optimised thr.  │              │  Best 15 features    │
│  MCC = 0.663±0.016  │              │  SVM (RBF kernel)    │
└────────┬────────────┘              └──────────┬───────────┘
         │                                      │
         │                          ┌──────────────────────┐
         │                          │  6. DEEP LEARNING    │
         │                          │                      │
         │                          │  One-hot encoding    │
         │                          │  CNN + LSTM + MLP    │
         │                          │  Ray Tune HPO        │
         │                          └──────────┬───────────┘
         │                                      │
         └──────────────┬───────────────────────┘
                        │
                        ▼
              ┌─────────────────────┐
              │  7. BENCHMARKING    │  hyperparameter_tuning.ipynb
              │                     │  model_evaluation_and_plots.ipynb
              │  Independent        │  benchmark_test.ipynb
              │  benchmark set      │
              │  (never seen)       │  → model_evaluation/
              │                     │  → model_evaluation_DL/
              └─────────────────────┘
```

---

## 📊 Final Results at a Glance

All three models are evaluated on the same independent benchmark set (never seen during training or hyperparameter tuning). The primary metric is **MCC** (Matthews Correlation Coefficient) — chosen because the dataset is class-imbalanced and MCC accounts for all four cells of the confusion matrix.

| Method | MCC | Precision | Sensitivity | Accuracy |
|--------|:---:|:---------:|:-----------:|:--------:|
| Von Heijne (PSWM) | 0.688 ± 0.016 | 0.934 | 0.711 | 0.691 |
| SVM (RBF kernel) | 0.808 | ~0.92 | ~0.93 | ~0.93 |
| **SP-NN (CNN+LSTM)** | **0.902** | **~0.97** | **~0.97** | **~0.97** |

```
MCC
 1.0 │
 0.9 │                                        ████████  0.902
 0.8 │               ████████  0.808
 0.7 │  ████████  0.688
 0.6 │
 0.5 │
     └──────────────────────────────────────────────────────
          Von Heijne          SVM              SP-NN
         (Rule-based)    (Feature ML)     (Deep Learning)

         +0.12 MCC over VH ──┘    +0.09 MCC over SVM ──┘
```

> **MCC interpretation:** 0.0 = random guessing · 1.0 = perfect · −1.0 = perfectly wrong

---

## 📂 Folder-by-Folder Guide

### 1. Data Collection
**`1.Data_Collection/`** · [Detailed README →](1.Data_Collection/README.md)

Downloads protein sequences directly from the **UniProt REST API** (Release 2025_03). Two separate queries are designed using the UniProt Advanced Search interface to ensure clean, high-quality, non-overlapping datasets.

**Positive set** — proteins *with* an experimentally confirmed signal peptide:
- Reviewed (Swiss-Prot), existence at protein level, eukaryotic, ≥40 aa, non-fragment, experimentally annotated SP (`ft_signal_exp:`)
- Extra filter: SP must be ≥14 residues and have a defined cleavage site

**Negative set** — proteins *without* a signal peptide:
- Same quality filters + explicitly no signal peptide annotation (`NOT ft_signal:`)
- Must be experimentally localised to a purely intracellular compartment (cytoplasm, nucleus, mitochondrion, plastid, peroxisome, or cellular membrane)
- Records whether a transmembrane helix starts within the first 90 residues (`HelixDomain` flag)

| Dataset | Proteins collected |
|---------|:-----------------:|
| Positive (SP) | 2,932 |
| Negative (no SP) | 20,615 |
| Negative with N-terminal TM helix | 1,384 |

**Key notebooks:** `get_dataset_pos.ipynb` · `get_dataset_neg.ipynb` · `output_recap.ipynb`

---

### 2. Data Preparation
**`2.Data_Preparation/`** · [Detailed README →](2.Data_Preparation/README.md)

Removes near-duplicate sequences and builds the master dataset file used by every downstream step.

**Why redundancy reduction?** If similar sequences appear in both training and test sets, the model appears to perform better than it truly does (data leakage). MMseqs2 clusters all sequences at **30% identity / 40% coverage** — two sequences in the same cluster cannot appear on opposite sides of the train/test split.

**Splitting strategy:**
- 80% of non-redundant sequences → **training set** (further split into 5 folds for cross-validation)
- 20% of non-redundant sequences → **benchmark set** (held out completely until final evaluation)

| Subset | Negatives | Positives | Total |
|--------|:---------:|:---------:|:-----:|
| Training (Folds 1–5) | 7,147 | 874 | 8,021 |
| Benchmark | 1,787 | 219 | 2,006 |
| **Total** | **8,934** | **1,093** | **10,027** |

The output `train_bench.tsv` is the single master file used by all models. Its key columns: `EntryID`, `Sequence`, `Class` (Positive/Negative), `Set` (1–5 or Benchmark), `SPStart`, `SPEnd`, `Kingdom`, `HelixDomain`.

**Key notebooks:** `recover_ids_clusters.ipynb` · `prepare_datasets.ipynb`

---

### 3. Data Analysis
**`3.Data_Analysis/`** · [Detailed README →](3.Data_Analysis/README.md)

Explores the dataset visually before any modelling — confirming the data is biologically sensible and the train/benchmark split is consistent.

**Analyses performed** (all figures saved automatically when `plots.ipynb` runs):

| Analysis | What it checks | Output folder |
|----------|---------------|---------------|
| Protein length distribution | Are positive and negative proteins the same size? | `Sequence_lengths_comparison/` |
| Signal peptide length | Do SPs fall in the expected 15–30 aa range? | `Sequence_lengths_comparison/` |
| Amino acid composition | Are SPs enriched in hydrophobic residues vs SwissProt? | `3.AA_Comparison/` |
| Taxonomic composition | Is the dataset representative? Training ≈ benchmark? | `4.Taxonomy_classification/` |
| Sequence logos | Is the AXA cleavage motif visible at the correct positions? | `5.SequenceLogo/` |

**Key insight confirmed:** Signal peptides are strongly enriched in hydrophobic residues (L, A, V, I) — the basis for all prediction methods. The AXA motif at positions −3/−1 relative to the cleavage site is clearly visible in the sequence logos.

**Key notebook:** `plots.ipynb`

---

### 4. Von Heijne Method
**`4.vonHeijne/`** · [Detailed README →](4.vonHeijne/README.md)

Implements the classical Von Heijne (1986) algorithm — the earliest and most interpretable method for signal peptide detection. Serves as the **rule-based baseline** against which machine learning methods are compared.

**How it works in three steps:**

1. **Build a PSWM** from training sequences — count amino acid frequencies at each of the 15 positions around the cleavage site (window: −13 to +2), normalise with pseudocounts, compute log-odds against SwissProt background
2. **Score sequences** by sliding the 15-residue window across the first 90 aa of each test protein and returning the maximum window score
3. **Select threshold** by maximising F1 on the validation fold, then apply to the test fold

**5-fold cross-validation results:**

| Metric | Mean ± Std |
|--------|:----------:|
| **MCC** | **0.663 ± 0.016** |
| Precision | 0.934 ± 0.003 |
| Sensitivity | 0.711 ± 0.029 |
| Accuracy | 0.691 ± 0.014 |

The PSWM heatmap and per-fold PR curves are saved to `Plots/`.

**Key notebooks:** `create_pswm.ipynb` · `validation_and_testing_vonheijne.ipynb` · `vonheijne.ipynb`

---

### 5. Feature Extraction & Selection
**`5.Feature_Selection/`** · [Detailed README →](5.Feature_Selection/README.md)

Converts raw sequences into a **38-dimensional numerical feature vector** per protein, then uses a Random Forest to rank features by importance and selects the most robust subset to train the final SVM.

**The 38 features** come from 5 categories:

| Category | Features | Count |
|----------|---------|:-----:|
| Physicochemical scales (max + mean) | Hydrophobicity, TM tendency, helix/β-sheet propensity, flexibility, membrane propensity, bulkiness, Argos hydrophobicity | 16 |
| Sequence motif | Von Heijne PSWM score | 1 |
| N-region | Basicity of first 5 residues (squared isoelectric points) | 1 |
| Composition | Aromatic residue frequency | 1 |
| Residue composition | Frequency of 19 standard amino acids | 19 |

**Feature selection pipeline:**
1. Train a **Random Forest** (1,000 trees, class-balanced) on the full 38 features → rank by Gini importance
2. Sweep over all possible subset sizes k = 2 to 38 → evaluate SVM validation MCC for each k
3. Pick the k that maximises validation MCC (not fixed — found independently per fold)
4. Count how many folds each feature appears in the best top-k subset
5. **Keep only features that appear in all 5 folds** → ensures universally informative features

**Final 15 features selected (appeared in all 5 folds):**
Von Heijne score · Cys frequency · Max TM tendency · Mean helix propensity · Mean Miyazawa hydrophobicity · Asp frequency · Thr frequency · Arg frequency · Max β-sheet propensity · Asn frequency · Max flexibility · Mean membrane propensity · Mean bulkiness · Met frequency · Max Argos hydrophobicity

**Cross-validation SVM performance (15 features):** Average Test MCC ≈ **0.82–0.84**

**Key notebooks:** `custom_features.ipynb` · `feauture_selection.ipynb`

---

### 6. Deep Learning — SP-NN
**`6.Deep_Learning/`** · [Detailed README →](6.Deep_Learning/README.md)

Implements SP-NN, a hybrid **CNN + LSTM + MLP** neural network that processes raw sequences directly — no feature engineering required. This was the optional deep learning extension of the project.

**Sequence encoding:**
- Standardise all sequences to **90 residues** (truncate or pad with `X`)
- One-hot encode each position as a **21-dimensional binary vector** (20 amino acids + `X`)
- Each sequence becomes a `[90 × 21]` matrix — the direct network input

**Architecture:**

```
[90 × 21 one-hot matrix]
        │
   Conv1D(kernel=17, 64 maps)   ← detects the ~17 aa hydrophobic core in one pass
        │
   LSTM × 2 (hidden=128)        ← models sequential context across the N-terminus
   → last hidden state           ← compresses 90 positions into one vector
        │
   BatchNorm1D
        │
   MLP: 128→256→128→64→1024→1   ← hourglass topology; Sigmoid output
        │
   P(SP) ∈ [0,1]
```

**Training stability mechanisms:**
- Gradient clipping (max_norm=1.0) — prevents LSTM exploding gradients
- Early stopping (patience=20) — saves best-validation-MCC checkpoint

**Hyperparameter optimisation:** Ray Tune random search, 15 trials × 5-fold CV. Best config: `hidden_sizes=[256,128,64,1024]`, `lstm_hidden=128`, `layers=2`, `dropout=0.498`, `lr=2.86e-4`, `batch=20`.

**Benchmark MCC: 0.9019**

The trained model is exported as a **TorchScript `.pt` file** — loadable anywhere without the original source code.

**Key notebooks:** `lstm.ipynb` (architecture + HPO) · `benchmark_test.ipynb` (final training + evaluation)

---

### 7. Model Performance Comparison
**`7.Model_performances/`** · [Detailed README →](7.Model_performances/README.md)

The final evaluation stage. All three models are tested on the **independent benchmark set** and subjected to a deep biological analysis of their errors.

**SVM fine-tuning:** The SVM is retrained using `BayesSearchCV` (60 iterations, 5-fold CV), finding optimal `C=0.6758` and `γ=0.0350`. Final benchmark MCC: **0.808**.

**Error analysis for both SVM and SP-NN** — figures saved to `model_evaluation/` (SVM) and `model_evaluation_DL/` (SP-NN):

**False Positive (FP) analysis — Why does the model call non-SPs as SPs?**
- TM helix proteins: their N-terminal hydrophobic helix mimics a signal peptide hydrophobic core
- Plant proteins (Viridiplantae): chloroplast transit peptides have a similar structure to SPs

**False Negative (FN) analysis — Why does the model miss real SPs?**
- Atypically short or long signal peptides (outside the canonical 15–30 aa range)
- Low hydrophobicity in the H-region (weaker than the canonical core)
- Enrichment in polar/charged residues (S, T, N, C, R, D) instead of hydrophobic ones

**Sequence logos** (`figure12A-1.png` to `_4.png`) visually confirm:
- TP sequences show the canonical AXA motif and hydrophobic core
- FN sequences show a weaker, attenuated version of the same motif
- FP sequences show a hydrophobic spike in the core region — these are TM helices

**Key notebooks:** `hyperparameter_tuning.ipynb` · `model_evaluation_and_plots.ipynb` · `benchmark_test.ipynb`

---

## 🧮 Evaluation Methodology

### Why MCC over accuracy?

The dataset has ~8× more negatives than positives. A model that always predicts "Negative" would achieve >88% accuracy but MCC = 0. MCC is the only standard metric that accounts for all four cells of the confusion matrix and remains meaningful under class imbalance.

```
           Predicted Positive   Predicted Negative
Actual Positive  │ TP (correct)  │  FN (missed SP)  │
Actual Negative  │ FP (false alarm)│  TN (correct)   │

MCC = (TP×TN − FP×FN) / √((TP+FP)(TP+FN)(TN+FP)(TN+FN))
Range: −1 (perfectly wrong) → 0 (random) → +1 (perfect)
```

### 5-fold cross-validation

Used during training and model selection. In each of 5 iterations:
- **3 folds** → training
- **1 fold** → validation (threshold / hyperparameter selection)
- **1 fold** → testing

Final cross-validation performance = mean ± standard error across 5 test folds.

### Independent benchmark

The 20% holdout set is used **once**, at the very end, to report unbiased final performance. It is never touched during training, cross-validation, or hyperparameter tuning.

## Conclusion

This project systematically compared three generations of signal peptide prediction methods using a high-quality, non-redundant dataset derived from UniProtKB/Swiss-Prot (2025 release):

- **Von Heijne (1986)** rule-based approach using a position-specific weight matrix (PSWM)  
- **SVM** classifier with carefully engineered and selected biochemical + sequence features  
- **SP-NN** — a modern hybrid CNN+LSTM deep learning model operating directly on one-hot encoded sequences

**Key quantitative findings** (independent benchmark set):

| Method              | MCC    | Precision | Sensitivity | Notes                              |
|---------------------|--------|-----------|-------------|------------------------------------|
| Von Heijne (PSWM)   | 0.688  | 0.934     | 0.711       | Strong baseline, high precision    |
| SVM (15 features)   | 0.808  | ~0.92     | ~0.93       | Clear improvement via features     |
| SP-NN (CNN+LSTM)    | **0.902** | ~0.97  | ~0.97       | Best overall performance           |

The deep learning model achieved the highest Matthews Correlation Coefficient (**+0.21** over the classical method, **+0.094** over the feature-based SVM), confirming that modern sequence modeling architectures capture patterns that are difficult to hand-craft.

**Main biological lessons learned from error analysis:**

- The dominant source of **false positives** in all models remains **N-terminal transmembrane helices**, whose hydrophobic character strongly resembles the signal peptide h-region.
- **False negatives** are most frequent among atypically short (<15 aa) or unusually long (>35 aa) signal peptides, and among sequences with reduced hydrophobicity or atypical amino acid composition in the core region.
- Plant chloroplast transit peptides and fungal signal peptides cause only minor performance degradation, suggesting good generalization across eukaryotic kingdoms.

Taken together, the results demonstrate clear progress in signal peptide prediction accuracy over the past four decades — from simple position-specific scoring matrices, through kernel methods with domain knowledge, to end-to-end deep learning. Nevertheless, transmembrane helix confusion remains a fundamental limitation shared by all current approaches.

This work provides a reproducible pipeline, well-documented feature extraction code, a trained TorchScript model, and detailed biological interpretation of errors — resources that can serve as a foundation for future improvements in secretory pathway prediction, protein localization annotation, and synthetic biology applications.

## 📁 Repository Structure (TO BE DELETED)

```
.
├── README.md                          ← You are here
│
├── 1.Data_Collection/
│   ├── README.md
│   ├── get_dataset_pos.ipynb          ← Positive dataset download
│   ├── get_dataset_neg.ipynb          ← Negative dataset download
│   └── output_recap.ipynb             ← Collection summary
│
├── 2.Data_Preparation/
│   ├── README.md
│   ├── recover_ids_clusters.ipynb     ← Parse MMseqs2 cluster output
│   ├── prepare_datasets.ipynb         ← Build train_bench.tsv
│   └── output_recap.ipynb
│
├── 3.Data_Analysis/
│   ├── README.md
│   ├── plots.ipynb                    ← All exploratory figures
│   ├── Sequence_lengths_comparison/   ← Length distribution figures
│   ├── 3.AA_Comparison/               ← Amino acid frequency figures
│   ├── 4.Taxonomy_classification/     ← Kingdom/species figures
│   └── 5.SequenceLogo/               ← Cleavage site logos
│
├── 4.vonHeijne/
│   ├── README.md
│   ├── create_pswm.ipynb              ← PSWM construction + heatmap
│   ├── validation_and_testing_vonheijne.ipynb
│   ├── vonheijne.ipynb                ← Full end-to-end pipeline
│   └── Plots/                         ← PR curves + PSWM heatmaps
│
├── 5.Feature_Selection/
│   ├── README.md
│   ├── custom_features.ipynb          ← All 38 feature extraction functions
│   ├── feauture_selection.ipynb       ← RF ranking + SVM training
│   └── *_features_*.npz              ← Pre-computed feature matrices
│
├── 6.Deep_Learning/
│   ├── README.md
│   ├── lstm.ipynb                     ← SP-NN architecture + Ray Tune HPO
│   ├── benchmark_test.ipynb           ← Final training + DL evaluation
│   └── SignalPeptideLSTM.pt           ← Exported TorchScript model
│
└── 7.Model_performances/
    ├── README.md
    ├── hyperparameter_tuning.ipynb    ← SVM Bayesian optimisation
    ├── model_evaluation_and_plots.ipynb ← SVM + Von Heijne analysis
    ├── benchmark_test.ipynb           ← SP-NN benchmark evaluation
    ├── SignalPeptideSVM.pkl           ← Saved SVM pipeline
    ├── model_evaluation/             ← SVM error analysis figures
    └── model_evaluation_DL/          ← SP-NN error analysis figures
```


