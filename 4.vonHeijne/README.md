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
The model is trained using a dataset (*train_bench.tsv*) containing:
- Protein sequences  
- Annotated cleavage site index (`SPEnd`)  

For each sequence, a fixed window of amino acids is extracted relative to the cleavage site:
- Positions: **[-13, +2]**
- This results in a window of 15 residues centered around the cleavage region  

This window captures the biological signal known to characterize signal peptides.

---

### 2. PSWM Construction

The notebook `create_pswm.ipynb` builds the Position-Specific Weight Matrix.

#### Step-by-step:
1. **Count frequencies**  
   For each position in the window, count occurrences of each amino acid across all training sequences.

2. **Apply pseudocounts**  
   A value of +1 is added to each count to avoid zero probabilities:
   
   f_{i,a} = \frac{count_{i,a} + 1}{N + 20}
   

3. **Background normalization**  
   Frequencies are compared against SwissProt background amino acid frequencies.

4. **Log-odds transformation**  
   The PSWM is computed as:
   \[
   PSWM_{i,a} = \log \left( \frac{f_{i,a}}{p_a^{background}} \right)
   \]

This matrix encodes positional enrichment or depletion of amino acids.

---

### 3. Cleavage Site Scoring

Implemented in `validation_and_testing_vonheijne.ipynb`.

#### Procedure:
1. A sliding window (~15 amino acids) is moved along each protein sequence  
2. For each window:
   - Extract amino acids at each position  
   - Sum corresponding PSWM log-odds scores:
     \[
     Score = \sum_{i=-13}^{+2} PSWM_{i, a_i}
     \]
3. The window with the **maximum score** is selected as the predicted cleavage site  

This corresponds to a maximum likelihood estimate under the PSWM model.

---

### 4. Threshold Optimization (Validation Phase)

Raw scores must be converted into binary predictions.

#### Steps:
1. Apply scoring to validation set  
2. Generate Precision–Recall (PR) curve  
3. Compute F1-score for different thresholds:
   \[
   F1 = 2 \cdot \frac{Precision \cdot Recall}{Precision + Recall}
   \]
4. Select threshold that maximizes F1-score  

This ensures a balance between false positives and false negatives.

---

### 5. Testing and Performance Evaluation

Using the optimal threshold, predictions are evaluated on the test set.

Metrics used:

- **Matthews Correlation Coefficient (MCC)**  
  Robust metric for imbalanced datasets  

- **Accuracy (ACC)**  
  Overall correctness  

- **Precision (PPV)**  
  Fraction of true positives among predicted positives  

- **Recall (Sensitivity, SEN)**  
  Fraction of true positives detected  

---

### 6. 5-Fold Cross-Validation

Implemented in `vonheijne.ipynb`.

#### Workflow:
1. Split dataset into 5 folds  
2. For each iteration:
   - 3 folds → training  
   - 1 fold → validation  
   - 1 fold → testing  

3. Repeat full pipeline:
   - Build PSWM  
   - Optimize threshold  
   - Evaluate on test set  

4. Aggregate results across folds:
   - Mean  
   - Standard deviation  

This provides a robust estimate of model generalization.

---

## Results

| Metric       | Mean ± Std     |
|--------------|---------------|
| MCC          | 0.663 ± 0.016 |
| Precision    | 0.934 ± 0.003 |
| Accuracy     | 0.691 ± 0.014 |
| Sensitivity  | 0.711 ± 0.029 |



## Key Idea

The Von Heijne method operates as a **probabilistic scoring model**, where:
- Biological signal = positional amino acid bias  
- Detection = maximizing log-odds score relative to background  

It is simple, interpretable, and biologically grounded.

