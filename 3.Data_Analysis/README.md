# Data Analysis of Signal Peptides

## Introduction

Signal peptides (SPs) are short N-terminal regions that direct newly synthesized proteins toward the secretory pathway. After targeting, these peptides are usually removed by signal peptidases. Because of their biological role, they represent an important signal for distinguishing secreted proteins from proteins that remain in other cellular compartments.

In this project, the data analysis step was carried out to examine the properties of proteins with and without signal peptides before applying predictive methods. The aim was to verify that the dataset was biologically meaningful, to compare the positive and negative classes, and to identify the sequence patterns that could later support approaches such as the **von Heijne method** and **SVM-based classification**.

---

## Workflow Overview

The analysis was performed on the processed dataset generated after the data collection and preprocessing steps. In particular, the input file used in this stage was:

- `../2.Data_Preparation/train_bench.tsv`

This file contains proteins divided into:

- five subsets for training and validation
- one independent benchmark set

The data analysis focused on the following aspects:

- **protein sequence length distribution**
- **signal peptide length distribution**
- **amino acid composition**
- **taxonomic classification**
- **cleavage-site motif visualization**

The analyses and figures were generated in:

- `plots.ipynb`

---

### 1. Protein Length Distribution

**Goal:** Compare the distribution of full protein lengths between proteins with signal peptides and proteins without signal peptides.

**Method:** Sequence lengths were extracted from the processed dataset and compared between positive and negative classes in both the **training** and **benchmark** sets. Since a few very long proteins can dominate the scale, filtered distribution plots were used to make the main patterns easier to interpret.

#### Training set
![Protein length - Training](./1,2.Sequence_lengths_comparison/Distribution_plot_training_distribution.png)

#### Benchmark set
![Protein length - Benchmark](./1,2.Sequence_lengths_comparison/Distribution_plot_benchmark.png)

**Interpretation:**  
The distributions show that the two classes do not have identical length profiles. Negative proteins cover a wider range of sequence lengths and include more long proteins, whereas positive proteins tend to be somewhat more compact. The benchmark set follows the same general trend observed in the training data, indicating that the split remains consistent and suitable for model evaluation.

---

### 2. Signal Peptide Length Distribution

**Goal:** Examine the variability of signal peptide lengths in the positive dataset.

**Method:** For proteins in the positive class, the signal peptide end position (`SPEnd`) was used to estimate signal peptide length. These lengths were then compared between the training and benchmark sets.

![Signal peptide length - Training vs Benchmark](./1,2.Sequence_lengths_comparison/sp_length_train_vs_bench.png)

**Interpretation:**  
Most signal peptides fall within the expected biological range, with the majority concentrated around 20–25 amino acids. This agrees with the known structure of signal peptides and supports the reliability of the annotations used in the dataset. The similarity between training and benchmark sets also suggests that the positive class is consistently represented across the split.

---

### 3. Amino Acid Composition

**Goal:** Identify compositional biases in signal peptide regions compared with a general protein background.

**Method:** Amino acid frequencies were computed from signal peptide sequences and compared with Swiss-Prot reference frequencies. This analysis was used to determine whether the dataset reproduced the expected residue preferences of signal peptides.

![Amino acid composition compared with Swiss-Prot](./3.AA_Comparison/frequencies_sp_train_vs_swissprot.png)

**Interpretation:**  
The composition profile shows a clear enrichment in hydrophobic residues, which is one of the most characteristic properties of signal peptides. This result is biologically expected, since signal peptides must interact with membrane-associated translocation machinery. The observed enrichment confirms that residue composition is likely to be informative for downstream prediction methods.

---

### 4. Taxonomic Classification

**Goal:** Assess whether the dataset remains biologically representative after preprocessing and splitting.

**Method:** Proteins were grouped according to their kingdom labels and their distribution was visualized for the training and benchmark sets.

![Kingdom distribution - Training](./4.Taxonomy_classification/kingdom_barplot_train.png)

**Interpretation:**  
The dataset is largely composed of eukaryotic groups, especially **Metazoa**, with additional representation from **Fungi**, **Viridiplantae**, and other eukaryotic categories. This matches the original data collection strategy. The taxonomic composition observed after preprocessing indicates that the dataset remains biologically coherent and that the benchmark set broadly reflects the same structure as the training data.

---

### 5. Cleavage Site Sequence Logos

**Goal:** Visualize residue conservation around signal peptide cleavage sites.

**Method:** Sequence windows surrounding the annotated cleavage sites were aligned and represented as sequence logos. These logos summarize position-specific residue preferences in the region most directly involved in signal peptide processing.

![Cleavage site sequence logo](./5.SequenceLogo/sequence_logo_training.png)

**Interpretation:**  
The sequence logo shows conserved residue preferences around the cleavage site, especially the presence of small and neutral residues close to the cleavage position. This agrees with the known biological rules governing signal peptide cleavage and provides further evidence that the positive dataset is coherent and suitable for later modeling.

---

## Conclusion

This analysis provided a descriptive overview of the dataset used for signal peptide prediction. The results show that the processed data are biologically plausible and consistent with the known properties of signal peptides. In particular:

- the sequence length distributions remain comparable between training and benchmark sets
- signal peptide lengths fall within the expected biological range
- amino acid composition reflects the strong hydrophobic character of signal peptides
- taxonomic composition remains coherent after preprocessing
- cleavage-site motifs show recognizable conserved patterns

Together, these results confirm that the dataset is suitable for the next stages of the project, including feature extraction, the application of von Heijne-based rules, and machine-learning classification.

---

## Files used in this step

**Main notebook**
- `plots.ipynb`

**Input dataset**
- `../2.Data_Preparation/train_bench.tsv`

**Output figure folders**
- `1,2.Sequence_lengths_comparison/`
- `3.AA_Comparison/`
- `4.Taxonomy_classification/`
- `5.SequenceLogo/`
