# **The VonHeijne method for SP detection** 

In this repository we built a **Position-Specific Weight Matrix (PSWM)** using our experimentally validated protein sequences and compared them against background amino acid frequencies (SwissProt). Pseudocount +1 was used to avoid zero probabilities. 
For the method implementation, the following parameters were considered: 
- **Aminoacids Window**: [-13, +2] relative to the cleavage site
- **Background distribution**: SwissProt amino acid frequencies
- **PSWM**: built from positive examples (proteins with confirmed signal peptides)
- **Aim**: identify likely cleavage sites by comparing observed amino acid frequencies with the   background model
---

## **PSWM implementation** 

The script [create_pswm.ipynb](./create_pswm.ipynb) was produced to build a Position-Specific Weight Matrix (PSWM) from a training dataset of proteins with annotated signal peptides in the *train_bench.tsv* file. 
The input DataFrame used must contain sequences and a SPEnd index for the cleavage site.

---

## **Evaluation of the Von Heijne Method**

The script [validation_and_testing_vonheijne](./validation_and_testing_vonheijne) provides code to evaluate the performance of the Von Heijne algorithm for signal peptide cleavage site prediction.
The evaluation is split into three main steps:

### 1. **Scoring with the PSWM**
Each sequence is scanned with a sliding window (14–15 aa length). For each window, a score is computed by summing log-odds values from the Position-Specific Weight Matrix (PSWM).The highest-scoring window is retained as the candidate cleavage site.

### 2.  **Threshold optimization (Validation set)**

Precision–Recall curves are computed on the validation set. The best threshold is chosen by maximizing the F1-score.

### 3.  **Testing and metrics**
On the test set, predictions are evaluated using:
- **Matthews Correlation Coefficient (MCC)**
- **Accuracy (ACC)**
- **Precision (PPV)**
- **Recall (SEN)**

---

## 5-Fold Cross-Validation

A more robust evaluation of the Von Heijne method is performed using a 5-fold cross-validation on the train_bench.tsv dataset using the [vonheijne](./vonheijne.ipynb):
- **Splitting**: In each fold, one set is used for testing, one for validation, and the remaining three sets for training.
- **PSWM computation**: The matrix is built for each iteration using the selected training sets.
- **Metrics**: For each fold, the following metrics are computed on the test set: MCC, PPV, ACC, SEN
- **Visualization**: Precision–Recall curves are plotted for each fold, highlighting the best threshold determined from the validation set.
--- 

## Results 
After cross-validation method, a mean value of each metric was calculated and it is shown below. 

 | Metric      | Mean ± Std  |
|------------|-------------|
| MCC        | 0.663 ± 0.016 |
| Precision  | 0.934 ± 0.003 |
| Accuracy   | 0.691 ± 0.014 |
| Sensitivity| 0.711 ± 0.029 |



