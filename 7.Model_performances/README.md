<div align="center">

# 🧬 Signal Peptide Prediction
### A Comparative Study: Von Heijne · SVM · Deep Learning

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python&logoColor=white)](https://python.org)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0%2B-EE4C2C?logo=pytorch&logoColor=white)](https://pytorch.org)
[![scikit-learn](https://img.shields.io/badge/scikit--learn-1.2%2B-F7931E?logo=scikit-learn&logoColor=white)](https://scikit-learn.org)

*From handcrafted rules to end-to-end deep learning — a rigorous benchmark on protein signal peptide classification*

</div>

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Dataset](#-dataset)
- [Methods](#-methods)
- [Hyperparameter Optimization](#-hyperparameter-optimization)
- [Results](#-results)
- [Error Analysis — SVM](#-error-analysis--svm-model_evaluation)
- [Error Analysis — Deep Learning](#-error-analysis--deep-learning-model_evaluation_dl)
- [Feature Distributions](#-feature-wise-distributions-svm)
- [Sequence Logos](#-sequence-logos)
- [Running the Code](#-running-the-code)
- [Dependencies](#-dependencies)

---

## 🔬 Overview

Signal peptides are short N-terminal sequences (~15–30 amino acids) that route nascent proteins into the secretory pathway. Accurate computational prediction is a fundamental challenge in proteomics and drug discovery.

This project benchmarks **three approaches of increasing complexity** on the same dataset and held-out benchmark, enabling a fair head-to-head comparison:

```
Protein Sequence (N-terminus)
         │
         ├─── Von Heijne Score  (PSWM)  ─────────────────► Binary Label
         ├─── Feature Engineering  →  RBF SVM ───────────► Binary Label
         └─── One-Hot Encoding  →  CNN → LSTM → MLP ──────► Binary Label
```

All output figures are saved automatically when the notebooks are run:
- **`model_evaluation/`** — SVM and Von Heijne analysis plots
- **`model_evaluation_DL/`** — SP-NN deep learning analysis plots

---

## 📊 Dataset

The dataset (`train_bench.tsv`) is partitioned into **5 cross-validation folds** and one **independent Benchmark set** withheld from all training and tuning:

| Column | Description |
|--------|-------------|
| `Sequence` | Full amino acid sequence |
| `Class` | `Positive` (has SP) / `Negative` (no SP) |
| `Set` | Fold: `1`–`5` or `Benchmark` |
| `SPStart` / `SPEnd` | Signal peptide cleavage boundaries |
| `Kingdom` | Eukaryota · Bacteria · Viridiplantae · Archaea |
| `HelixDomain` | `True` if sequence contains a transmembrane helix |

> ⚠️ The benchmark is **never used during training or hyperparameter search** — it serves exclusively as the final unbiased evaluation set.

---

## 🧪 Methods

### 1. Von Heijne Method

A classical **Position-Specific Weight Matrix (PSWM)** scoring approach. The matrix is built from positive training sequences over a window of ±13/+2 residues around the cleavage site. A test sequence is classified as a signal peptide if its maximum PSWM score exceeds a threshold.

**Strengths:** Interpretable, biologically motivated, no iterative training.  
**Limitations:** Cannot model non-linear interactions or long-range dependencies.

---

### 2. Support Vector Machine (SVM)

An **RBF-kernel SVM** (`sklearn.Pipeline: StandardScaler → SVC`) trained on 18 hand-crafted biochemical descriptors per sequence:

| Feature | Scale | Aggregation |
|---------|-------|-------------|
| Von Heijne score | PSWM | Max |
| Hydrophobicity | Miyazawa (`mi`) | Mean over SP region |
| Transmembrane tendency | Custom TM scale | Max |
| Chou-Fasman helix propensity | Chou & Fasman 1978 | Mean |
| Chou-Fasman β-sheet propensity | Chou & Fasman 1978 | Max |
| Residue flexibility | Vihinen (`Flex`) | Max |
| Membrane propensity | Punta-Maritan 2003 | Max |
| Mean bulkiness | Zimmerman 1968 | Mean |
| Hydrophobicity (Argos) | Argos (`ag`) | Max |
| Amino acid frequencies | — | C, D, T, R, N, M |

The Von Heijne feature (index 17) is recomputed using the full training PSWM before fitting to ensure consistency between training and benchmark encoding.

---

### 3. SP-NN — Hybrid Deep Learning

SP-NN processes raw sequences directly — **no feature engineering required**.

**Encoding:** All sequences standardized to **90 residues** (padded with `X` or truncated from C-terminus), then encoded as a **21-dimensional one-hot vector** per position → `[90 × 21]` matrix.

**Architecture:**

```
Input [Batch, 90, 21]
    │
    ▼
Conv1D(in=21, out=64, kernel=17, padding='same')   ← kernel=17 ≈ SP hydrophobic core length
    │
    ▼
LSTM(64 → 128, layers=2)  →  last hidden state  →  BatchNorm1D(128)
    │
    ▼
MLP: Linear(128→256) → ReLU → Dropout
     Linear(256→128) → ReLU → Dropout
     Linear(128→ 64) → ReLU → Dropout
     Linear( 64→1024)→ ReLU → Dropout    ← Hourglass: compress then expand
     Linear(1024→  1)→ Sigmoid
    │
    ▼
P(SP) ∈ [0,1]  →  threshold 0.5  →  {Positive, Negative}
```

| Design choice | Rationale |
|--------------|-----------|
| Conv1D kernel=**17** | Matches the typical length of the signal peptide hydrophobic core |
| LSTM × 2 | Captures long-range sequential dependencies across the N-terminus |
| Hourglass MLP `[256→128→64→1024]` | Compress-then-expand topology found optimal by hyperparameter search |
| Gradient clipping (max_norm=1.0) | Prevents exploding gradients inherent to LSTM training |
| Early stopping (patience=20) | Saves best-validation-MCC checkpoint, prevents overfitting |

**Best configuration (from Ray Tune):**

```python
config = {
    'hidden_sizes':     [256, 128, 64, 1024],
    'lstm_hidden_size': 128,
    'num_lstm_layers':  2,
    'dropout':          0.498,
    'lr':               2.86e-4,
    'batch_size':       20,
}
```

The model is exported as a **TorchScript** `.pt` file — loadable anywhere without the original class definition.

---

## ⚙️ Hyperparameter Optimization

### SVM — Bayesian Search (`BayesSearchCV`)

60 iterations with 5-fold CV. Bayesian optimization builds a probabilistic surrogate of the objective and proposes the next configuration from the most promising region — far more efficient than grid or random search.

**Optimal hyperparameters:**

| Hyperparameter | Search Range | Best Value | Interpretation |
|---------------|-------------|-----------|----------------|
| `C` | [0.01, 100] log-uniform | **0.6758** | Soft margin — noisy feature space |
| `γ (gamma)` | [0.01, 100] log-uniform | **0.0350** | Wide RBF — avoids tight overfitting |
| **Best CV MCC** | — | **0.8569** | Validation performance |

> A low `C` (0.68) with a low `γ` (0.035) indicates the SVM needed a wide, soft-margin classifier. The feature space is separable but noisy, and tight boundaries would overfit.

---

### SP-NN — Ray Tune Random Search

15 trials × 5-fold CV each. Ray's shared Object Store (`ray.put`) avoids memory explosion when distributing the one-hot encoded dataset across worker processes.

| Parameter | Range | Best |
|-----------|-------|------|
| `num_layers` | {2, 3, 4, 5} | 4 |
| MLP topology | Funnel / Flexible | Flexible (hourglass) |
| `dropout` | [0.1, 0.5] uniform | 0.498 |
| `lr` | [1e-4, 1e-3] log-uniform | 2.86e-4 |
| `lstm_hidden_size` | {64, 128} | 128 |
| `num_lstm_layers` | {1, 2} | 2 |

---

## 🏆 Results

All models evaluated on the **independent Benchmark set** (never seen during training or tuning). Primary metric: **MCC** — robust to class imbalance.

| Method | MCC | Notes |
|--------|:---:|-------|
| Von Heijne | 0.688 | Rule-based, no training required |
| SVM (RBF) | 0.808 | +0.12 over Von Heijne |
| **SP-NN (CNN+LSTM)** | **0.902** | **+0.09 over SVM** |

**Interpretation:**

- 🔵 **Von Heijne (0.688):** Captures cleavage-site residue preferences but cannot model non-linear interactions or sequence context beyond the immediate cleavage window.
- 🟢 **SVM (0.808):** Properly engineered biochemical features and the RBF kernel together handle the non-linear separability in feature space that a linear PSWM score cannot.
- 🟡 **SP-NN (0.902):** End-to-end learning avoids feature engineering entirely. The CNN captures local hydrophobic motifs; the LSTM models the full sequential context of the N-terminus — together learning a richer representation than any handcrafted descriptor set.

---

## 🔍 Error Analysis — SVM (`model_evaluation/`)

Each benchmark prediction is labelled **TP / TN / FP / FN** and subjected to systematic biological analysis. Results are saved in the `model_evaluation/` directory.

---

### Transmembrane Domain & False Positives

**Hypothesis:** Transmembrane (TM) proteins are over-represented among false positives because their N-terminal hydrophobic TM helix structurally resembles a signal peptide hydrophobic core.

<div align="center">

![Helix Domain Predictions](model_evaluation/HelixDomain_predictions.png)

*Percentage of sequences with/without a transmembrane helix domain across all four prediction groups (TP, TN, FP, FN). False Positives are strongly enriched for HelixDomain=True, confirming the TM-helix confusion hypothesis.*

</div>

<div align="center">

![HD Prediction 2](model_evaluation/HD_prediction_2.png)

*Comparison of transmembrane helix prevalence between False Positives only (FP) and all negative-class sequences (TN+FP). The FP group has a disproportionately higher fraction of TM proteins, showing the SVM is specifically confused by TM helices.*

</div>

---

### Kingdom-level False Positive Rate

<div align="center">

![FPR per Kingdom](model_evaluation/pieplot_kingdom_FPR.png)

*False Positive Rate (FPR) broken down by taxonomic kingdom. Plant sequences (Viridiplantae) are expected to show elevated FPR due to **chloroplast transit peptides** — N-terminal targeting sequences that are structurally similar to signal peptides and can fool the classifier.*

</div>

---

### Taxonomic Composition of Misclassifications

<div align="center">

![Taxa Composition](model_evaluation/Taxa_composition_FP_FN.png)

*Three side-by-side pie charts: total positive sequences, False Negatives, and True Positives, broken down by kingdom. Enrichment of a particular kingdom in the FN chart indicates that kingdom's signal peptides are systematically harder to detect.*

</div>

---

### Signal Peptide Length Distribution

<div align="center">

![SP Length Distribution](model_evaluation/Signalength_distribution.png)

*Signal peptide length (SPEnd − SPStart) for correctly classified True Positives vs missed False Negatives. FN sequences tend to have atypically short or long signal peptides — deviating from the canonical 15–30 aa range that the model has learned to recognise.*

</div>

---

### Hydrophobicity Analysis

<div align="center">

![Hydrophobicity Distribution](model_evaluation/Hydrophobicity_distribution.png)

*Distribution of mean hydrophobicity (Miyazawa scale) over the signal peptide region for TP vs FN sequences. FN sequences show a leftward shift — a weaker hydrophobic core is the most common reason signal peptides are missed.*

</div>

<div align="center">

![Hydrophobicity Boxplot](model_evaluation/Boxplot_hydrophobicity.png)

*Boxplot version of the same comparison. The median hydrophobicity of FN sequences is visibly lower than TP, reinforcing the finding that atypically hydrophilic signal peptides drive most false negatives.*

</div>

<div align="center">

![Mean Hydrophobicity](model_evaluation/Mean_hydrophobicity.png)

*Mean hydrophobicity per prediction group. Despite FN sequences having lower median hydrophobicity, the overall means between TP and FN are close — suggesting hydrophobicity alone is not sufficient to explain all false negatives.*

</div>

---

### Amino Acid Frequency Analysis

<div align="center">

![AA Frequencies SP](model_evaluation/AA_fequencies_3.png)

*Relative amino acid frequencies **within the signal peptide region** (SPStart:SPEnd) for TP vs FN sequences. Residues are grouped by chemical class. FN signal peptides show enrichment in polar and charged residues (S, T, N, C, R, D) at the expense of canonical hydrophobic residues (L, A, V, I).*

</div>

<div align="center">

![AA Frequency All Groups](model_evaluation/figure9.png)

*Extended 4-group comparison: FP vs all Negatives vs FN vs all Positives. This reveals which amino acid compositional differences are specific to misclassified sequences versus simply reflecting class-level differences between SPs and non-SPs.*

</div>

---

## 🔍 Error Analysis — Deep Learning (`model_evaluation_DL/`)

The same systematic analysis is repeated for the SP-NN model, enabling direct comparison of error patterns between the SVM and the deep learning approach.

---

### Transmembrane Domain & False Positives

<div align="center">

![HD Prediction DL 1](model_evaluation_DL/HD_prediction_1.png)

*Percentage of sequences with/without a transmembrane helix across all prediction groups for SP-NN. The deep model also shows FP enrichment for TM proteins, but to a lesser degree than the SVM — the LSTM's sequential context helps distinguish some TM helices from true SPs.*

</div>

<div align="center">

![HD Prediction DL 2](model_evaluation_DL/HD_prediction_2.png)

*FP vs TN+FP helix domain comparison for SP-NN. Compared to the SVM, the gap between FP and TN+FP is smaller, indicating the DL model handles TM proteins better.*

</div>

---

### Kingdom-level False Positive Rate

<div align="center">

![FPR Kingdom DL](model_evaluation_DL/FPR_Kingdom.png)

*Kingdom-level FPR for the SP-NN. Comparing this with the SVM equivalent reveals whether the deep model disproportionately improves on specific taxonomic groups (e.g. plant sequences with transit peptides).*

</div>

---

### Taxonomic Composition of Misclassifications

<div align="center">

![Taxonomy DL](model_evaluation_DL/Pieplot_species.png)

*Taxonomic breakdown of Total Positives, False Negatives, and True Positives for the SP-NN. A more uniform FN composition compared to the SVM would indicate the DL model generalises better across kingdoms.*

</div>

---

### Signal Peptide Length Distribution

<div align="center">

![SP Length DL](model_evaluation_DL/Signalength_distribution.png)

*SP length distribution for TP vs FN in the SP-NN. The DL model's FN distribution is narrower than the SVM's — the CNN kernel (size 17) is specifically tuned to the canonical SP hydrophobic core length, reducing length-related false negatives.*

</div>

---

### Hydrophobicity Analysis

<div align="center">

![Hydrophobicity DL](model_evaluation_DL/Hydrophobicity_SP.png)

*Hydrophobicity distribution over the SP region for the SP-NN predictions. The DL model shows a cleaner separation between TP and FN hydrophobicity distributions compared to the SVM.*

</div>

<div align="center">

![Hydrophobicity Boxplot DL](model_evaluation_DL/Boxplot_Hydrophobicity_SP.png)

*Boxplot comparison of SP hydrophobicity for TP vs FN in the SP-NN. The reduced overlap relative to the SVM confirms that the deep model is more robust to low-hydrophobicity signal peptides.*

</div>

<div align="center">

![Mean Hydrophobicity DL](model_evaluation_DL/Mean_hydrophobicity.png)

*Mean hydrophobicity per prediction group for the SP-NN.*

</div>

---

### Amino Acid Frequency Analysis

<div align="center">

![AA Frequencies DL](model_evaluation_DL/AA_frequencies.png)

*Amino acid frequency comparison (TP vs FN) for SP-NN predictions. Residues coloured by chemical class: 🔵 nonpolar (G,A,V,L,I,M,P) · 🟠 aromatic (F,Y,W) · 🟢 polar (S,T,N,Q,C) · 🔴 positive (K,R,H) · 🟣 negative (D,E).*

</div>

<div align="center">

![AA Frequencies SP DL](model_evaluation_DL/AA_frequencies_SP.png)

*AA frequencies restricted to the **SP sub-sequence only** for the SP-NN, isolating compositional signals within the cleavage region rather than the full protein.*

</div>

<div align="center">

![AA Frequency Distribution DL](model_evaluation_DL/AA_frequency_Distribution.png)

*4-group comparison for the SP-NN: FP vs Negatives vs FN vs Positives. Cross-referencing with the SVM equivalent (`figure9.png`) reveals which compositional biases are model-independent versus specific to one classifier.*

</div>

---

## 📈 Feature-wise Distributions (SVM)

For the SVM, the 15 most informative features are visualised as KDE plots and boxplots comparing distributions across all four prediction groups (TP, TN, FP, FN). These plots identify exactly which biochemical properties drive misclassification.

---

### Transmembrane Propensity

<div align="center">

![KDE TM Propensity](model_evaluation/KDE_TranmembranP.png)

*Kernel Density Estimation of Max Transmembrane Propensity (custom TM scale). TP and FP sequences share a high-propensity peak — confirming TM helices mimic signal peptide hydrophobic cores. TN sequences cluster at low propensity. FN sequences occupy an intermediate zone.*

</div>

<div align="center">

![Boxplot TM](model_evaluation/Boxplot_TRTE.png)

*Boxplot version of transmembrane propensity across all four groups. The FP median is close to TP, confirming the TM-confusion hypothesis. The wider FP distribution reflects the structural diversity of TM proteins.*

</div>

---

### Chou-Fasman Propensities

<div align="center">

![Chou-Fasman Helix](model_evaluation/Boxplot_Chou_Fasman.png)

*Mean Chou-Fasman **helix** propensity across groups. Signal peptides (TP) show higher helical tendency than non-SPs (TN), consistent with the amphipathic helix structure of the hydrophobic core. FN sequences have reduced helix propensity — atypically non-helical signal peptides.*

</div>

<div align="center">

![Chou-Fasman Beta](model_evaluation/Chou_Fasman_BProp.png)

*Max Chou-Fasman **β-sheet** propensity. FP sequences (TM proteins) show elevated β-sheet propensity compared to TP — TM helices have a distinct secondary structure signature that partially distinguishes them from true SPs.*

</div>

---

### Flexibility, Membrane Propensity & Bulkiness

<div align="center">

![Residue Flexibility](model_evaluation/Boxplot_ResidueFlex.png)

*Max residue flexibility (Vihinen scale). FN sequences tend toward higher flexibility — a disordered, flexible N-terminus correlates with weaker signal peptide structure and increases misclassification risk.*

</div>

<div align="center">

![Punta-Maritan](model_evaluation/Boxplot_PuntaMaritan.png)

*Max membrane propensity (Punta-Maritan 2003 scale). FP sequences (TM proteins) show elevated membrane propensity compared to TP sequences, providing a discriminative signal the SVM can exploit to reduce the FPR on TM proteins.*

</div>

<div align="center">

![Bulkiness](model_evaluation/boxplot_Bulkiness.png)

*Mean residue bulkiness (Zimmerman 1968). TP sequences show slightly higher mean bulkiness than non-SPs — consistent with the preference for large hydrophobic residues (L, I, V, F) in signal peptide hydrophobic cores.*

</div>

---

### Hydrophobicity (Argos)

<div align="center">

![Argos Hydrophobicity](model_evaluation/Boxplot_ArgosMax.png)

*Max residue hydrophobicity (Argos scale). This is one of the strongest individual discriminators: TP sequences have substantially higher Argos hydrophobicity than TN. FN sequences overlap with the lower end of the TP distribution, explaining why weak hydrophobic cores lead to missed detections.*

</div>

---

## 🔤 Sequence Logos

Sequence logos (generated with **logomaker**) align all sequences around the cleavage site (positions −13 to +4) and visualise positional residue conservation as **information content** in bits.

---

### SVM — `model_evaluation/`

<div align="center">

| True Positives | False Negatives |
|:--------------:|:---------------:|
| ![TP Logo SVM](model_evaluation/figure12A-1.png) | ![FN Logo SVM](model_evaluation/figure12A-2.png) |
| **Strong AXA cleavage motif** at positions −3/−1 and a dominant hydrophobic core (L/A/V) from −10 to −4. This is the canonical signal peptide signature. | **Attenuated motif signal** — lower bit scores throughout. The hydrophobic core is weaker and the AXA motif less conserved. These are atypical SPs the model fails to detect. |

| True Negatives | False Positives |
|:--------------:|:---------------:|
| ![TN Logo SVM](model_evaluation/figure12A-3.png) | ![FP Logo SVM](model_evaluation/figure12A-4.png) |
| No structured motif — broad, unspecific residue distribution. Non-SP sequences lack both the hydrophobic core and the AXA cleavage signal. | Elevated hydrophobic signal at positions −10 to −4, mimicking the SP hydrophobic core. These are predominantly transmembrane proteins whose N-terminal TM helix confuses the classifier. |

*Sequence logos for SVM — Information Content (bits). Cleavage site at position 0.*

</div>

---

### Deep Learning — `model_evaluation_DL/`

<div align="center">

| True Positives | False Negatives |
|:--------------:|:---------------:|
| ![TP Logo DL](model_evaluation_DL/figure12A-1.png) | ![FN Logo DL](model_evaluation_DL/figure12A-2.png) |
| SP-NN TP logo shows the same canonical AXA motif and hydrophobic core, confirming the model has learned the correct biological signal. | DL FN sequences show a weaker but still partially structured motif — the SP-NN misses fewer edge cases than the SVM, as reflected in higher overall sensitivity. |

| True Negatives | False Positives |
|:--------------:|:---------------:|
| ![TN Logo DL](model_evaluation_DL/figure12A-3.png) | ![FP Logo DL](model_evaluation_DL/figure12A-4.png) |
| TN logo for SP-NN is similarly unstructured — the model correctly rejects sequences with no SP-like pattern. | DL FP sequences still show hydrophobic enrichment in the core region, but the peak is narrower than the SVM equivalent — the LSTM provides enough sequential context to reject more TM proteins. |

*Sequence logos for SP-NN — Information Content (bits). Cleavage site at position 0.*

</div>

---

## 🚀 Running the Code

Run notebooks in this exact order:

```
1. lstm.ipynb
   └── Defines SP_NN architecture, utility functions, runs Ray Tune (15 trials × 5-fold CV)
       Prints best config at the end → hardcode into benchmark_test.ipynb
           ↓
2. hyperparameter_tuning.ipynb
   └── Trains SVM with BayesSearchCV (60 iter)
       Saves: SignalPeptideSVM.pkl  +  benchmark_features.npz
           ↓
3. benchmark_test.ipynb
   └── Retrains SP-NN on folds 1–4, validates on fold 5, evaluates on benchmark
       Saves: SignalPeptideLSTM.pt
       Generates: model_evaluation_DL/ (all DL plots)
           ↓
4. model_evaluation_and_plots.ipynb
   └── Loads SVM + Von Heijne, evaluates on benchmark
       Generates: model_evaluation/ (all SVM + Von Heijne plots)
```

### Prerequisites

```bash
pip install torch scikit-learn scikit-optimize biopython \
            matplotlib seaborn logomaker numpy pandas \
            "ray[tune]" importnb
```

> 💡 **No GPU?** In `lstm.ipynb`, change `resources_per_trial={"cpu": 4, "gpu": 0}`.  
> ⚠️ **SVM feature shape:** Ensure `training` in `hyperparameter_tuning.ipynb` contains only Sets 2–5 (not Set 1) to match `x_training_conc`.

---

## 📦 Dependencies

| Package | Purpose |
|---------|---------|
| `torch >= 2.0` | SP-NN training and TorchScript export |
| `ray[tune] >= 2.0` | Distributed hyperparameter search |
| `scikit-learn >= 1.2` | SVM, metrics, preprocessing |
| `scikit-optimize` | `BayesSearchCV` for SVM tuning |
| `biopython` | `ProteinAnalysis`, `ProtParamData` scales |
| `matplotlib / seaborn` | All visualisations |
| `logomaker` | Sequence logo generation |
| `numpy / pandas` | Numerical and tabular data |
| `importnb` | Import `.ipynb` as Python modules |

---

## 📌 Conclusion

<div align="center">

| | Von Heijne | SVM | SP-NN |
|:---|:---:|:---:|:---:|
| **MCC** | 0.688 | 0.808 | **0.902** |
| Feature engineering | Manual (PSWM) | Manual (18 descriptors) | **None** |
| Training required | ✗ | ✓ | ✓ |
| Sequence representation | Cleavage window | Full protein | N-terminus (90 aa) |
| Architecture | Rule-based | Kernel machine | CNN + LSTM + MLP |

</div>

Each step up in model complexity yields a meaningful, consistent gain. The **MCC of 0.902** for SP-NN confirms that a CNN+LSTM architecture can learn both local hydrophobic motifs and sequential N-terminal context directly from raw sequences — outperforming both handcrafted rules and classical ML without any domain-expert feature engineering.

Error analysis consistently identifies **transmembrane-domain proteins** as the primary source of false positives across both models, and **atypically polar or short signal peptides** as the main source of false negatives — providing a concrete roadmap for future model improvements.

---

<div align="center">
<sub>Built with PyTorch · scikit-learn · Ray Tune · Biopython · logomaker</sub>
</div>
