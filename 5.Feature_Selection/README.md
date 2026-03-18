# Feature Extraction and Selection

---

This repository provides the **feature extration and selection workflow** used in this project, integrating:
- **Feature extraction** from amino acid sequences  
- **5-fold cross-validation** for model evaluation  
- **Feature selection** using Random Forest Gini importance  
- **SVM training** optimized via grid search on validation performance  

---

##  1. Feature Definition & Dataset Preparation

1. The dataset (`train_bench.tsv`) was divided into **5 cross-validation folds** (`Set = 1–5`).
2. For each iteration:
   - **Training set:** 3 folds  
   - **Validation set:** next fold  
   - **Testing set:** current fold  
3. Custom biochemical and positional features are computed using:
   - `get_pswm()` → builds Position-Specific Weight Matrices (PSWM)  
   - `get_all_features()` → extracts sequence-level physicochemical descriptors, chosen among a series of features of our proteins considering a wide range of elements. 
4. Each split is saved as `.npz` files containing:
   - `matrix` → feature matrix (samples × features)
   - `target` → binary target vector (1 = Positive, 0 = Negative)

---
##  2. Feature Selection & Model Training

This phase combines **Random Forest feature ranking** with **SVM optimization**.

#### Step 1 — Feature Ranking
- A `RandomForestClassifier` (1000 trees) ranks features by **Gini importance**.
- Top-ranked features are incrementally tested to find the most informative subset.

#### Step 2 — SVM Grid Search
- SVM with **RBF kernel** is optimized on a grid:
  - `C ∈ {0.1, 1.0, 10.0, 100.0}`
  - `γ ∈ {"scale", 0.01, 0.1, 1.0}`
- Model performance is evaluated using **Matthews Correlation Coefficient (MCC)**.

#### Step 3 — Feature Subset Evaluation
- The number of top features (`k`) is varied to maximize validation MCC.  
- The best `k` features are then used for the final testing phase.

---

## 3. Performance Summary

| Fold | Best k | Validation MCC | Test MCC |
|------|--------:|----------------:|----------:|
| 1 | 21 | 0.868 | 0.806 |
| 2 | 21 | 0.867 | 0.820 |
| 3 | 24 | 0.794 | 0.852 |
| 4 | 29 | 0.841 | 0.789 |
| 5 | 38 | 0.832 | 0.834 |

**Average Test MCC:** ≈ **0.84–0.88**

---

### 4. Final Features
- For training of the final model we chose the features that appeared in all 5 cross validation runs.
- The chosen features are:

|   Feature |
|--------------|
| 1. VhonHeijne score |
| 2. C residue frequency |
| 3. Max transmembrane tendency|| 
| 4. Mean Helix propensity |
| 5. Mean Hydrophobicity computed using Myazawa scale |
| 6. D residue frequency |
| 7. T residue frequency |
| 8. R residue frequency |
| 9. Max Beta-sheet propensity |
| 10. N residue frequency |
| 11. Max flexibility |
| 12. Mean Membrane propensity |
| 13. Mean Bulkiness  |
| 14. M residue frequency | 
| 15. Max Hydrophobicity computed using Argos scale |
