# DATA PREPARATION 

---
## 1. Building Training and Benchmarking Datasets: Splitting

Once the preliminary data are collected, the datasets need to be divided into two distinct subsets:

- **Training set**  
  Used to train the models, optimize hyperparameters, and perform cross-validation.
- **Benchmarking set**  
  Also known as the *holdout dataset*, it is reserved for testing the generalization performance of the models.

 ### Benchmarking Set: Motivation
Cross-validation alone is not sufficient to guarantee an unbiased estimate of generalization performance:

- Hyperparameter tuning through cross-validation and grid search may still introduce **overfitting**.  
- A holdout dataset provides a stronger guarantee of the *never-seen-before* condition.  
- The model tested on the benchmarking set is the one intended for **production** use.  
- During cross-validation, models are trained on slightly different subsets of the training dataset, which may bias results.

---

##  2. Redundancy Reduction

Before splitting the data, it is essential to create a **non-redundant dataset**.  
This involves:
- Controlling **sequence identity** and **alignment coverage**  
- Performing redundancy reduction using clustering tools with **MMSeqs2**  
- Selecting **one representative sequence per cluster** 
Once redundancy has been addressed, the split can safely be performed randomly.

### 3. MMseqs2

For the clusterisation procedure the **MMseqs2** software (version 14.7e284) was used. MMseqs2 is an open-source software use to cluster large databases of sequences. input file were: negative_dataset.fasta and positive_dataset.fasta. For each run two output files were created. The parameters for the clustering procedure were : 30% sequence identity, and 40% coverage. Each of the condition must be satisfied. For the negative set: 
- [neg_cluster.tsv](./neg_cluster/neg_cluster.tsv): .TSV containing two colums ( ID of each sequence in th input file, ID of the representative sequence indetifying the cluster). 
- [neg_rep_seq.fasta](./neg_cluster/neg_rep_seq.tsv): .fasta file containing all the representative sequences, one for each found cluster.
- [neg_all_seq.fasta](./neg_cluster/neg_all_seq.tsv): .fasta file with all the sequences used. 
  The same output files ( [positive_cluster.tsv](./pos_cluster/pos_cluster.tsv), [positive_rep_seq.fasta](./pos_cluster/pos_rep_seq.tsv), [neg_all_seq.fasta](./pos_cluster/pos_all_seq.tsv ) was obtained using as input *positive_dataset.fasta*. 

### 4. Results
results filtered from a proper [script](./output_recap.ipynb) are visualised below: 

 | Non-redundant positives | Non-redundant negatives | N-r negatives with helix transmembrane |
|--------------------------|-------------------------|----------------------------------------|
|          1093            |          8934           |             636                        |

---

## 5.  Data Splitting Strategy

To ensure proper training and unbiased evaluation, the dataset is divided as follows:

- **Training vs. Benchmarking sets**
  - Randomly assign **80%** of both positive and negative sequences to the **training set**.
  - Assign the remaining **20%** of both positive and negative sequences to the **benchmarking (holdout) set**.

- **Cross-validation on training data**
  - Build **5-fold cross-validation subsets** from the training set.
  - Each split preserves the overall **positive/negative ratio**.
  - Store information about the cross-validation subset each protein belongs to, so results remain reproducible.

  Our approach was to generate a .tsv file using the script [prepare_dataset](./prepare_datasets.ipynb),containing representative entries obtained after clustarisation. .tsv contains the following columns

| Column |
|--------------------------|
|        1. EntryID              | 
|        2. OrganismName              | 
|        3. Kingdom              |
|        4. Sequence length              |
|        5. HelixDomain ( True/False for the negative entries, NaN for the positive ones)              |
|        6. Class ( Negative/Positive)              |
|        7. SPstart (defined for positive entries)              |
|        8. SPend (defined for positives entries)              | 
|        9. Set (from 1-5 for the entries of the training sets, Benchmark for the ones in the benchmark set)             |
|        10. Sequence              | 

  
 ## 6. Results 

Resulting numbers are visualised on the table below


|       Dataset          |        Negatives        |     Positives      |  
|-------------------------|-------------------------|-------------------|
|  Training Sets (total)  |          7147           |       874         |  
| Benchmark Set           |          1787           |       219         |
| Total                   |          8934           |       1093        |

---




