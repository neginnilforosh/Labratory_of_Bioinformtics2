# The Von Heijne Method for Signal Peptide Detection

## Overview
This repository implements the Von Heijne method for predicting signal peptide (SP) cleavage sites using a Position-Specific Weight Matrix (PSWM). The approach is based on statistical enrichment of amino acids around experimentally validated cleavage sites.

The pipeline includes:
1. PSWM construction from annotated sequences  
2. Scoring of candidate cleavage sites  
3. Threshold optimization  
4. Performance evaluation and cross-validation  

---

## Methodological Framework

### 1. Dataset and Preprocessing
The model is trained using a dataset (train_bench.tsv) containing:
- Protein sequences  
- Annotated cleavage site index (SPEnd)  

For each sequence, a fixed window of amino acids is extracted:
- Positions: [-13, +2]  
- Window length: 15 residues  

---

### 2. PSWM Construction

The notebook `create_pswm.ipynb` builds the PSWM.

Steps:
1. Count amino acid frequencies at each position  
2. Apply pseudocounts (+1) to avoid zero probabilities  

   f(i,a) = (count(i,a) + 1) / (N + 20)

3. Normalize using background (SwissProt frequencies)  

4. Compute log-odds matrix:

   PSWM(i,a) = log( f(i,a) / p_background(a) )

---

### 3. Cleavage Site Scoring

Implemented in `validation_and_testing_vonheijne.ipynb`.

Procedure:
- Slide a 15-aa window across the sequence  
- Compute score for each window:

  Score = sum from i = -13 to +2 of PSWM(i, a_i)

- Select highest-scoring window as predicted cleavage site  

---

### 4. Threshold Optimization

- Compute Precision–Recall curve on validation set  
- Select optimal threshold using F1-score:

  F1 = 2 * (Precision * Recall) / (Precision + Recall)

---

### 5. Testing and Metrics

Evaluate predictions using:
- MCC  
- Accuracy  
- Precision  
- Recall (Sensitivity)  

---

### 6. 5-Fold Cross-Validation

- Split dataset into 5 folds  
- For each iteration:
  - 3 folds: training  
  - 1 fold: validation  
  - 1 fold: testing  

- Repeat full pipeline and average results  

---

## Results

| Metric       | Mean ± Std     |
|--------------|---------------|
| MCC          | 0.663 ± 0.016 |
| Precision    | 0.934 ± 0.003 |
| Accuracy     | 0.691 ± 0.014 |
| Sensitivity  | 0.711 ± 0.029 |

---

## Usage

1. Run `create_pswm.ipynb`  
2. Run `validation_and_testing_vonheijne.ipynb`  
3. Run `vonheijne.ipynb`  

---

## Key Idea

The method identifies cleavage sites by maximizing a log-odds score derived from positional amino acid enrichment relative to a background model.
