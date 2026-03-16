# Data Analysis of Signal Peptide Datasets

## Introduction

This section summarizes the exploratory analysis performed on the curated dataset used for signal peptide prediction. The analysis is based on the files contained in this archive and was designed to assess dataset quality, compare positive and negative classes, and highlight biologically meaningful patterns before downstream modeling.

The dataset used in this step comes from `../2.Data_Preparation/train_bench.tsv` and contains proteins divided into:

- five training subsets (`Set = 1, 2, 3, 4, 5`)
- one independent benchmark set (`Set = Benchmark`)

Each entry includes the following key fields:

- `Class` (`Positive` or `Negative`)
- `SequenceLength`
- `Kingdom`
- `SPStart`
- `SPEnd`
- `Sequence`

The positive class contains proteins annotated with a signal peptide, whereas the negative class contains proteins without a signal peptide.

---

## Workflow Overview

The exploratory analysis in the archive focuses on four main aspects:

1. **Sequence length analysis**
2. **Signal peptide length analysis**
3. **Amino acid composition analysis**
4. **Taxonomic composition analysis**
5. **Cleavage-site motif visualization using sequence logos**

The main analysis code is contained in:

- `plots.ipynb`

The generated figures are stored in:

- `1,2.Sequence_lengths_comparison/`
- `3.AA_Comparison/`
- `4.Taxonomy_classification/`
- `5.SequenceLogo/`

---

## Dataset Summary

The analyzed dataset contains both training and benchmark proteins after preprocessing and redundancy reduction.

### Number of proteins by class

- **Positive proteins:** 1093
- **Negative proteins:** 8934

### Number of proteins by set

- **Training subsets 1–5:** 8021 proteins in total
- **Benchmark set:** 2006 proteins

### Number of proteins by class within benchmark split

- **Benchmark positive:** 219
- **Benchmark negative:** 1787

### Kingdom distribution

The taxonomic annotations in this dataset are restricted to eukaryotic groups and are mainly distributed across:

- **Metazoa**
- **Fungi**
- **Viridiplantae**
- **Other**

This is consistent with the original data collection strategy, which focused on eukaryotic proteins.

---

## 1. Protein Sequence Length Analysis

### Goal

To compare the distribution of full protein lengths between positive and negative proteins and to evaluate how these distributions behave across training and benchmark sets.

### Files

- `1,2.Sequence_lengths_comparison/sequence_length_training.png`
- `1,2.Sequence_lengths_comparison/sequence_length_benchmark.png`
- `1,2.Sequence_lengths_comparison/Distribution_plot_training_distribution.png`
- `1,2.Sequence_lengths_comparison/Distribution_plot_benchmark.png`
- `1,2.Sequence_lengths_comparison/Distribution_plot_set_division.png`
- `1,2.Sequence_lengths_comparison/Distribution_plot_set_division_2.png`
- `1,2.Sequence_lengths_comparison/Distribution_plot_by_set_facet.png`

### Description

The notebook compares protein lengths in the positive and negative classes for both the training and benchmark sets. To improve readability, some visualizations focus on proteins below a chosen length threshold, since a small number of very long proteins would otherwise dominate the scale.

The resulting plots include:

- histograms
- density-style distribution plots
- comparisons across training subsets and benchmark set

### Interpretation

The sequence length distributions show that the two classes are not identical in size profile. In this dataset, negative proteins tend to include more long sequences and a broader range of lengths, while positive proteins are generally more compact. This is consistent with the summary statistics of the dataset, where the negative class has a higher average sequence length than the positive class.

These plots are useful because they confirm that sequence length is a potentially informative variable and that the benchmark set follows a similar overall structure to the training data.

---

## 2. Signal Peptide Length Analysis

### Goal

To examine the distribution of signal peptide lengths in positive proteins and compare these lengths across training subsets and the benchmark set.

### Files

- `1,2.Sequence_lengths_comparison/sp_length_train_vs_bench.png`
- `1,2.Sequence_lengths_comparison/Distribution_plot_SPE_distribution.png`
- `1,2.Sequence_lengths_comparison/Distribution_plot_distribution.png`

### Description

For positive proteins, the signal peptide length is represented by the `SPEnd` value. The notebook extracts this information and visualizes its distribution in the training and benchmark sets.

The plots include comparisons between:

- training vs benchmark
- individual subsets within the training partition
- all sets combined

### Interpretation

The signal peptide lengths are concentrated within the expected biological range. In this dataset, the median signal peptide length is approximately 22 residues, with most values falling around 19–25 residues. This agrees with the known biological structure of signal peptides, which typically consist of a short N-terminal targeting sequence with a hydrophobic core and a cleavage region.

The consistency between training and benchmark distributions suggests that the benchmark set is appropriate for evaluating generalization.

---

## 3. Amino Acid Composition Analysis

### Goal

To compare the amino acid composition of signal peptide regions in the dataset against reference amino acid frequencies from Swiss-Prot and to evaluate residue-category patterns across different subsets.

### Files

- `3.AA_Comparison/frequencies_sp_train_vs_swissprot.png`
- `1,2.Sequence_lengths_comparison/frequencies_sp_train_vs_swissprot.png`
- `1,2.Sequence_lengths_comparison/frequencies_sp_bench_vs_swissprot_aacategories.png`
- `1,2.Sequence_lengths_comparison/AminoAcid_Frequencies_Set.png`
- `1,2.Sequence_lengths_comparison/Frequencies_propierties_Set.png`
- `1,2.Sequence_lengths_comparison/Distribution_Set.png`

### Description

The notebook computes amino acid frequencies from signal peptide sequences and compares them with reference Swiss-Prot frequencies. It also groups amino acids into broader physicochemical categories to assess residue-property composition across the different subsets.

The generated plots examine:

- residue-level frequency enrichment and depletion
- category-level composition
- differences across training subsets and benchmark data

### Interpretation

The signal peptide regions show the expected enrichment in hydrophobic residues, which is a hallmark of signal sequences. This feature reflects the biological role of signal peptides in targeting proteins to the secretory pathway and interacting with membrane-associated translocation machinery.

The comparison with Swiss-Prot background frequencies helps highlight that signal peptides are compositionally distinct from average protein regions. This supports the use of residue composition as an informative feature for both classical and machine-learning-based prediction approaches.

---

## 4. Taxonomic Classification Analysis

### Goal

To evaluate the taxonomic composition of the dataset and verify whether the training and benchmark sets remain comparable from a biological perspective.

### Files

- `4.Taxonomy_classification/kingdom_barplot_train.png`
- `4.Taxonomy_classification/kingdom_barplot_bench.png`
- `4.Taxonomy_classification/kingdom_pie_train.png`
- `4.Taxonomy_classification/kingdom_pie_bench.png`
- `4.Taxonomy_classification/species_barplot_train.png`
- `4.Taxonomy_classification/species_barplot_bench.png`
- `4.Taxonomy_classification/species_pie_train.png`
- `4.Taxonomy_classification/species_pie_bench.png`

### Description

The dataset was grouped by taxonomic labels and visualized using bar plots and pie charts. These figures summarize the distribution of proteins across kingdoms and the most represented organisms.

### Interpretation

The plots confirm that the dataset is dominated by eukaryotic groups, particularly **Metazoa**, followed by **Fungi** and **Viridiplantae**. The benchmark set broadly reflects the same taxonomic structure observed in training.

This is important because it indicates that the split between training and benchmark was not accompanied by a major biological shift in taxonomic composition. As a result, performance differences observed later are less likely to arise from hidden sampling bias.

---

## 5. Cleavage-Site Sequence Logo Analysis

### Goal

To visualize residue conservation around signal peptide cleavage sites and verify whether the dataset reflects known biological cleavage patterns.

### Files

- `5.SequenceLogo/sequence_logo_training.png`
- `5.SequenceLogo/sequence_logo_benchmark.png`
- `5.SequenceLogo/sequence_logo_benchmark_colored.png`
- `5.SequenceLogo/mirrored_logo_fixed.png`

### Description

The sequence logo analysis was generated from aligned windows surrounding the cleavage site of positive proteins. Separate logos were created for the training and benchmark sets.

These logos summarize position-specific residue preferences around the cleavage region and provide an intuitive visual representation of sequence conservation.

### Interpretation

The sequence logos show conserved preferences for small, neutral residues close to the cleavage site, which is in agreement with the classical biological rules described for signal peptide processing. This provides additional evidence that the positive dataset is biologically coherent and suitable for downstream modeling.

The cleavage-site motifs extracted from this dataset are particularly relevant for later implementation of the von Heijne approach, since that method depends directly on position-specific residue patterns around the cleavage position.

---

## Concluding Remarks

The exploratory analysis contained in this archive confirms that the dataset is internally consistent and biologically meaningful. The sequence length distributions, signal peptide length distributions, amino acid composition patterns, taxonomic structure, and cleavage-site sequence logos all align with the expected properties of signal peptide-containing proteins.

Overall, this step provides a strong descriptive foundation for the following stages of the project, including:

- construction of position-specific scoring approaches
- feature extraction
- machine learning classification
- final model evaluation

---

## Files Used in This Step

### Main notebook
- `plots.ipynb`

### Input dataset
- `2.Data_Preparation/train_bench.tsv`

### Output figure folders
- `1,2.Sequence_lengths_comparison/`
- `3.AA_Comparison/`
- `4.Taxonomy_classification/`
- `5.SequenceLogo/`
