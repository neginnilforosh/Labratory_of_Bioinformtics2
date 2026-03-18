# The Von Heijne method for SP detection

This repository implements the construction of a **Position-Specific Weight Matrix (PSWM)** using experimentally validated protein sequences. These are compared against background amino acid frequencies derived from SwissProt. A pseudocount of +1 is applied in order to avoid zero probabilities.

The implementation is based on the following parameters:
- **Amino acid window**: [-13, +2] relative to the cleavage site  
- **Background distribution**: SwissProt amino acid frequencies  
- **PSWM**: constructed from positive examples (proteins with confirmed signal peptides)  
- **Objective**: identification of likely cleavage sites through comparison between observed amino acid frequencies and the background model  

---

## PSWM implementation

The notebook [create_pswm.ipynb](./create_pswm.ipynb) is used to construct the Position-Specific Weight Matrix (PSWM) from a training dataset of proteins with annotated signal peptides provided in *train_bench.tsv*.

The input DataFrame is required to include:
- Protein sequences  
- Cleavage site index (SPEnd)  

---

## Evaluation of the Von Heijne method

The script [validation_and_testing_vonheijne](./validation_and_testing_vonheijne) contains the implementation for evaluating the performance of the Von Heijne algorithm in predicting signal peptide cleavage sites.

The evaluation procedure is structured in three main stages:

### 1. Scoring using the PSWM
Each sequence is scanned using a sliding window of approximately 14–15 amino acids. For each window, a score is computed as the sum of log-odds values derived from the PSWM. The window with the highest score is selected as the predicted cleavage site.

### 2. Threshold optimization (validation set)
Precision–Recall curves are computed on the validation dataset. The optimal threshold is determined by maximizing the F1-score.

### 3. Testing and performance metrics
Predictions on the test dataset are evaluated using the following metrics:
- **Matthews Correlation Coefficient (MCC)**  
- **Accuracy (ACC)**  
- **Precision (PPV)**  
- **Recall (Sensitivity, SEN)**  

---

## 5-Fold cross-validation

A more robust evaluation of the Von Heijne method is performed through 5-fold cross-validation on the *train_bench.tsv* dataset using [vonheijne](./vonheijne.ipynb).

For each fold:
- **Data splitting**: one subset is used for testing, one for validation, and the remaining three for training  
- **PSWM computation**: a new matrix is constructed for each iteration using the training subsets  
- **Evaluation**: performance metrics (MCC, PPV, ACC, SEN) are calculated on the test subset  
- **Visualization**: Precision–Recall curves are generated, highlighting the optimal threshold selected from the validation subset  

---

## Results

Following cross-validation, the mean and standard deviation of each metric are reported:

| Metric       | Mean ± Std     |
|--------------|---------------|
| MCC          | 0.663 ± 0.016 |
| Precision    | 0.934 ± 0.003 |
| Accuracy     | 0.691 ± 0.014 |
| Sensitivity  | 0.711 ± 0.029 |
