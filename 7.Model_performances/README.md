# 🏆 Step 7 — Model Performance Comparison

> **What this folder does:** Brings all three models together — Von Heijne, SVM, and SP-NN — evaluates them on the same independent benchmark set, and performs a deep biological analysis of each model's errors to understand *why* misclassifications happen.

---

## 🗂️ Files in this Folder

| File | What it does |
|------|-------------|
| `hyperparameter_tuning.ipynb` | Fine-tunes the SVM using Bayesian search and saves the final trained model |
| `model_evaluation_and_plots.ipynb` | Evaluates the SVM and Von Heijne method; generates all SVM error analysis figures |
| `benchmark_test.ipynb` | Trains the final SP-NN with the optimal config; evaluates it and generates all DL error analysis figures |
| `lstm.ipynb` | Defines SP-NN architecture and runs Ray Tune hyperparameter search (same as Step 6) |

### Where figures are saved

| Folder | Model | Contents |
|--------|-------|----------|
| `model_evaluation/` | SVM + Von Heijne | All SVM and Von Heijne analysis figures |
| `model_evaluation_DL/` | SP-NN | All deep learning analysis figures |

---

## 🧪 The Three Models Compared

| Model | Strategy | Feature engineering? | Training? |
|-------|----------|:--------------------:|:---------:|
| Von Heijne | Score sequences against a PSWM | No — uses cleavage-site statistics | No — rule-based |
| SVM (RBF) | Classify 15 hand-crafted biochemical features | Yes — manual | Yes |
| SP-NN (CNN+LSTM) | Learn directly from one-hot sequences | No — end-to-end | Yes |

---

## ⚙️ SVM Hyperparameter Tuning (`hyperparameter_tuning.ipynb`)

The SVM from Step 5 used a simple grid search. Here, a more thorough **Bayesian optimisation** is applied using `BayesSearchCV` (60 iterations, 5-fold CV):

**Why Bayesian over grid search?** Bayesian optimisation builds a statistical model (surrogate) of the objective function as it evaluates configurations. Instead of testing random points or an exhaustive grid, it intelligently focuses on the most promising region of the search space — finding better parameters with fewer evaluations.

```python
bayes = BayesSearchCV(
    estimator = Pipeline([("scaler", StandardScaler()), ("svm", SVC())]),
    search_spaces = {
        "svm__kernel": Categorical(["rbf"]),
        "svm__C":     Real(0.01, 100, prior="log-uniform"),
        "svm__gamma": Real(0.01, 100, prior="log-uniform"),
    },
    scoring = make_scorer(matthews_corrcoef),
    cv = 5,
    n_iter = 60,
)
```

### Optimal hyperparameters found

| Parameter | Search range | Best value | What it means |
|-----------|-------------|:----------:|---------------|
| `C` | [0.01, 100] | **0.6758** | Regularisation — controls margin width. Low C = wide, tolerant margin |
| `γ (gamma)` | [0.01, 100] | **0.0350** | RBF kernel width. Low γ = each training point has a wide influence radius |
| **Best CV MCC** | — | **0.8569** | Mean MCC across the 5 validation folds |

> A low `C` and low `γ` together indicate the model needed a soft, wide-radius decision boundary — the feature space is separable but noisy, and tight boundaries would memorise the training data instead of generalising.

The final SVM is saved as `SignalPeptideSVM.pkl` for use in evaluation.

---

## 📊 Final Results on the Independent Benchmark Set

All three models are evaluated on the benchmark set that was **never used during training or hyperparameter tuning**.

### Performance table

| Method | MCC | Notes |
|--------|:---:|-------|
| Von Heijne | 0.688 | Rule-based, no training required |
| SVM (RBF) | 0.808 | +0.12 over Von Heijne |
| **SP-NN (CNN+LSTM)** | **0.902** | **+0.09 over SVM** |

### What MCC means

MCC (Matthews Correlation Coefficient) ranges from **−1** (completely wrong) to **+1** (perfect). A score of **0** means the model is no better than random guessing. MCC is used instead of accuracy because the benchmark set has more negatives than positives — accuracy can be misleadingly high even with a poor model.

### Why each model performs as it does

- **Von Heijne (0.688):** Captures the AXA cleavage motif and hydrophobic core via the PSWM, but cannot model anything beyond the 15-residue cleavage window. It has no way to account for N-region basicity, full-sequence hydrophobicity, or taxonomic variation.

- **SVM (0.808):** The 15 hand-crafted biochemical features encode much richer information than a PSWM alone — global hydrophobicity, transmembrane tendency, flexibility, and amino acid composition all contribute. The RBF kernel handles non-linear separability that the linear PSWM cannot.

- **SP-NN (0.902):** End-to-end learning from raw sequences avoids the need for feature engineering entirely. The CNN kernel (size 17) directly detects the hydrophobic core in one operation; the LSTM captures the entire sequential context from the N-terminus; together they learn a richer representation than any handcrafted set could provide.

---

## 🔍 Error Analysis — SVM (`model_evaluation/`)

Each benchmark entry is labelled TP, TN, FP, or FN and the following analyses identify the biological causes of misclassification.

### Transmembrane domain proteins and False Positives

**Hypothesis:** Transmembrane (TM) proteins are over-represented among false positives because their N-terminal hydrophobic TM helix closely resembles a signal peptide hydrophobic core — confusing the model.

| Figure | What it shows |
|--------|--------------|
| `model_evaluation/HelixDomain_predictions.png` | Percentage of sequences with/without a TM helix across TP, TN, FP, FN. FP is strongly enriched for `HelixDomain=True` |
| `model_evaluation/HD_prediction_2.png` | Direct comparison: FP only vs all Negatives (TN+FP). FP has a higher TM-helix fraction, confirming the hypothesis |

### Kingdom-level False Positive Rate

| Figure | What it shows |
|--------|--------------|
| `model_evaluation/pieplot_kingdom_FPR.png` | FPR broken down by taxonomic kingdom. Plant sequences (Viridiplantae) are expected to show elevated FPR due to chloroplast **transit peptides** — N-terminal targeting sequences structurally similar to signal peptides |

### Taxonomic composition of misclassifications

| Figure | What it shows |
|--------|--------------|
| `model_evaluation/Taxa_composition_FP_FN.png` | Three pie charts side by side: total Positives, False Negatives, and True Positives — broken down by kingdom. Enrichment of a kingdom in the FN chart means its signal peptides are systematically harder to detect |

### Signal peptide length distribution

| Figure | What it shows |
|--------|--------------|
| `model_evaluation/Signalength_distribution.png` | SP length (SPEnd − SPStart) for TP vs FN. FN sequences tend to have atypically short or long SPs — outside the 15–30 aa canonical range the model has learned |

### Hydrophobicity analysis

| Figure | What it shows |
|--------|--------------|
| `model_evaluation/Hydrophobicity_distribution.png` | Distribution of mean hydrophobicity (Miyazawa scale) over the SP region — TP vs FN. FN sequences shift left (weaker hydrophobic core) |
| `model_evaluation/Boxplot_hydrophobicity.png` | Boxplot version — the FN median is lower than TP |
| `model_evaluation/Mean_hydrophobicity.png` | Mean hydrophobicity bar chart per group |

### Amino acid frequency analysis

| Figure | What it shows |
|--------|--------------|
| `model_evaluation/AA_fequencies_3.png` | AA frequencies **within the SP region** (SPStart to SPEnd) — TP vs FN. FN SPs are enriched in polar/charged residues (S, T, N, C, R, D) at the expense of hydrophobic ones (L, A, V, I) |
| `model_evaluation/figure9.png` | 4-group comparison: FP / all Negatives / FN / all Positives |

### Feature-wise distributions (SVM-specific)

For the 15 most informative features, box plots and KDE curves compare TP, TN, FP, and FN groups:

| Figure | Feature | Key finding |
|--------|---------|-------------|
| `model_evaluation/KDE_TranmembranP.png` | Max TM tendency | TP and FP share a high-propensity peak — TM helices mimic SP cores |
| `model_evaluation/Boxplot_TRTE.png` | Max TM tendency | FP median ≈ TP median; confirms TM confusion |
| `model_evaluation/Boxplot_Chou_Fasman.png` | Mean helix propensity | TP has higher helical tendency; FN are less helical |
| `model_evaluation/Chou_Fasman_BProp.png` | Max β-sheet propensity | FP (TM proteins) show elevated β-sheet — structural signature of TM β-barrels |
| `model_evaluation/Boxplot_ResidueFlex.png` | Max flexibility | FN sequences tend to be more flexible — disordered N-terminus |
| `model_evaluation/Boxplot_PuntaMaritan.png` | Max membrane propensity | FP sequences have elevated membrane propensity vs TP |
| `model_evaluation/boxplot_Bulkiness.png` | Mean bulkiness | TP sequences are bulkier (large hydrophobic residues L, I, V, F) |
| `model_evaluation/Boxplot_ArgosMax.png` | Max Argos hydrophobicity | One of the strongest discriminators: TP >> TN; FN overlaps the low end of TP |

### Sequence logos

Sequences are aligned around the cleavage site (positions −13 to +4) and visualised as information-content logos:

| Figure | Group | What to look for |
|--------|-------|-----------------|
| `model_evaluation/figure12A-1.png` | True Positives | Strong AXA motif at −3/−1; hydrophobic core (L/A/V) from −10 to −4 |
| `model_evaluation/figure12A-2.png` | False Negatives | Weaker version of the TP motif — atypical SPs with less conserved core |
| `model_evaluation/figure12A-3.png` | True Negatives | No structured motif — broad, unspecific residue distribution |
| `model_evaluation/figure12A-4.png` | False Positives | Hydrophobic enrichment at −10 to −4 mimicking the SP core — these are TM helices |

---

## 🔍 Error Analysis — SP-NN (`model_evaluation_DL/`)

The same analyses are repeated for the deep learning model. Comparing SP-NN figures with the SVM equivalents reveals which error patterns are model-independent (a property of the data) vs model-specific.

### Key comparisons

| Analysis | What SP-NN does better |
|----------|----------------------|
| TM helix FP rate | Lower — the LSTM's sequential context helps distinguish some TM helices from true SPs |
| SP length distribution | Narrower FN spread — CNN kernel tuned to the 17-aa canonical SP length handles atypical lengths better |
| Hydrophobicity separation | Cleaner — DL model more robust to low-hydrophobicity SPs |
| Sequence logos | Same canonical AXA motif is learned (TP logo), but FN logo shows slightly stronger residual motif than SVM |

### SP-NN figure map

| Figure | Description |
|--------|-------------|
| `model_evaluation_DL/HD_prediction_1.png` | TM helix prevalence by prediction group |
| `model_evaluation_DL/HD_prediction_2.png` | FP vs TN+FP helix domain comparison |
| `model_evaluation_DL/FPR_Kingdom.png` | Kingdom-level FPR for SP-NN |
| `model_evaluation_DL/Pieplot_species.png` | Taxonomic composition of misclassifications |
| `model_evaluation_DL/Signalength_distribution.png` | SP length distribution — TP vs FN |
| `model_evaluation_DL/Hydrophobicity_SP.png` | Hydrophobicity histogram — TP vs FN |
| `model_evaluation_DL/Boxplot_Hydrophobicity_SP.png` | Hydrophobicity boxplot |
| `model_evaluation_DL/Mean_hydrophobicity.png` | Mean hydrophobicity per group |
| `model_evaluation_DL/AA_frequencies.png` | Full-sequence AA frequency — TP vs FN |
| `model_evaluation_DL/AA_frequencies_SP.png` | SP-region AA frequency — TP vs FN |
| `model_evaluation_DL/AA_frequency_Distribution.png` | 4-group AA frequency comparison |
| `model_evaluation_DL/figure12A-1.png` to `_4.png` | Sequence logos for TP, FN, TN, FP |

---

## 📌 Conclusion

| | Von Heijne | SVM | SP-NN |
|:--|:---:|:---:|:---:|
| **MCC** | 0.688 | 0.808 | **0.902** |
| Feature engineering | Manual (PSWM) | Manual (15 descriptors) | None |
| Training required | ✗ | ✓ | ✓ |
| Sequence representation | 15-residue cleavage window | Full protein (15 features) | N-terminus (90 aa, raw) |
| Architecture | Rule-based | RBF kernel machine | CNN + LSTM + MLP |

The error analyses across all three models consistently point to two root causes of misclassification that are **model-independent** (present in all three):

1. **False Positives:** Transmembrane domain proteins whose N-terminal hydrophobic helix mimics a signal peptide hydrophobic core
2. **False Negatives:** Atypically polar, short, or compositionally unusual signal peptides that deviate from the canonical AXA / hydrophobic-core template

Both are genuine biological challenges, not model defects. Addressing them would require either richer training data (more diverse atypical SPs) or explicit transmembrane-helix discriminators in the feature set.

---

## 🔁 How to Run

Run in this order:

```
1. hyperparameter_tuning.ipynb
   └── Outputs: SignalPeptideSVM.pkl, benchmark_features.npz
       ↓
2. benchmark_test.ipynb
   └── Outputs: SignalPeptideLSTM.pt, model_evaluation_DL/ figures
       ↓
3. model_evaluation_and_plots.ipynb
   └── Outputs: model_evaluation/ figures
```

> ⚠️ `benchmark_test.ipynb` and `model_evaluation_and_plots.ipynb` both require `../2.Data_Preparation/train_bench.tsv`.  
> ⚠️ `model_evaluation_and_plots.ipynb` requires `SignalPeptideSVM.pkl` and `benchmark_features.npz` (produced by `hyperparameter_tuning.ipynb`).

---

## 📦 Required Libraries

```python
torch                # SP-NN inference and TorchScript
ray[tune]            # Hyperparameter search
scikit-learn         # SVM, metrics, preprocessing
scikit-optimize      # BayesSearchCV
biopython            # Physicochemical scales
matplotlib, seaborn  # All plots
logomaker            # Sequence logos
numpy, pandas        # Numerical and tabular data
importnb             # Import .ipynb as modules
```
