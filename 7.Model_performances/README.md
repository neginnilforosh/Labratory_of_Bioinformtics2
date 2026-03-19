
# Performance: VonHeijne , SVM, Deep Learning 

Our final models Included:

- **SP-NN (Signal Peptide Neural Network)**
Deep learning model combining CNN, LSTM, and MLP layers for sequence classification.

- **Von Heijne Method**
Classical scoring-based approach used as a simple interpretable baseline.

- **Support Vector Machine (SVM)**
A traditional machine learning model trained on engineered sequence features.

---

## **SVM Optimization**

SVM optimization was developed  via Bayesian search and provides a performance analysis module to evaluate model quality and interpret false positives (FP) and false negatives (FN). 

### **Bayesian search optimization**
The script **hyperparameter_tuning.ipynb** contains the implementation of the Bayesian search using an RBF SVM using BayesSearchCV. It reports best hyperparameters and MCC scores. Bayesian optimization uses past evaluations to predict where the best hyperparameters probably are, wasting  less time testing bad configurations. 

### **Results on Validation Set**
The table below shows the results of the hyperparameter tuining procedure. The MCC obtained suggest a better performance of the SVM method in respect to the classical approach of the VonHeije method. 

  
| Hyperparameter | Value | Description |
|----------------|-------|-------------|
| `svm__C` | **0.6758** | Regularization parameter controlling the trade-off between margin size and classification error |
| `svm__gamma` | **0.0350** | RBF kernel coefficient defining the influence of single training examples |
| `svm__kernel` | **rbf** | Kernel type used by the SVM (Radial Basis Function) |
| **Best MCC (validation)** | **0.8569** | Matthews Correlation Coefficient obtained during cross-validation |

---

## **Performance Analysis**

The script **model_evaluation.ipynb** contains the evaluation procedure that took in account different features to differentiate our results on the SVM model. The script **benchmark_test.ipynb** contains the same procedure applied to the deep learning network. 
During the model evaluation phase, several in-depth analyses were conducted to better understand the model’s behavior and error patterns, focusing on False Negatives (FN) and False Positives (FP).
These analyses aim to identify biological or compositional reasons behind the misclassifications and assess the robustness of the feature set.

---


## **False Negative (FN) Analysis**

- **Taxonomy distribution:**
Pie charts comparing the taxonomy of FN, FP, and total predictions.
→ In particular, we expect to observe many FP among plant sequences due to the presence of transit peptides.

- **Amino acid frequency distribution:**
Histogram comparing the amino acid composition of FN versus True Positives (TP).

- **Signal peptide composition:**
Comparison between the amino acid composition of signal peptides in TP versus FN sequences.

- **Signal peptide length distribution:** 
Length comparison of signal peptides between TP and FN groups.

- **Hydrophobicity distribution:**
Comparison of average hydrophobicity between TP and FN sequences.

- **Feature-wise distributions:**
In general, distributions were generated for all the most informative features to detect systematic differences between FN and TP samples.

---

## **False Positive (FP) Analysis**

- **FPR on transmembrane proteins:**
Calculation of the false positive rate (FPR) specifically for transmembrane domain proteins.

- **FPR with Von Heijne features:**
Same FPR calculation applied to models including the Von Heijne feature set, to evaluate its impact on transmembrane-related misclassifications.

--- 
## **Results** 
 | Metod     |   MCC  |
|------------|--------|
| VonHeijne  | 0.688  |
| SVM        | 0.808  |
| DL         | 0.902  |

The results show a clear performance hierarchy among the tested methods:

- **Von Heijne (MCC = 0.688)**
As expected from a rule-based, biologically inspired scoring approach, performance is moderate. The method captures general signal-peptide tendencies but lacks the flexibility to model complex sequence patterns.

- **SVM (MCC = 0.808)**
The SVM provides a strong improvement over Von Heijne, showing that relatively simple machine-learning methods can effectively capture discriminative sequence features once they are properly encoded.

- **Deep Learning (MCC = 0.902)**
The deep model achieves the best performance by a significant margin. This indicates that the CNN+LSTM architecture is able to learn both local motifs and long-range dependencies within the sequence, outperforming both handcrafted rules and classical ML.

We demonstrated that richer models capture progressively more nuanced sequence information. The MCC of 0.902 suggests that the deep learning approach is highly reliable for this classification task.

