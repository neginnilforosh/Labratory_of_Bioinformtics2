# Data Analysis of Signal Peptide Datasets

This section describes the exploratory analysis carried out on the curated dataset before feature extraction and model building. The aim was to check that the training and benchmark sets were biologically coherent, comparable to one another, and suitable for downstream prediction tasks.

All analyses in this step were generated from the processed dataset in `../2.Data_Preparation/train_bench.tsv` and the figures produced by `plots.ipynb`.

---

## What was analyzed

The analysis focused on five questions:

1. Do proteins with and without signal peptides differ in overall sequence length?
2. Are signal peptide lengths consistent across training and benchmark sets?
3. Does the amino acid composition of signal peptides match known biological expectations?
4. Is the dataset taxonomically representative and balanced across the split?
5. Do cleavage-site motifs show the expected conserved patterns?

The figures are grouped in the same order below.

---

## 1. Protein sequence length

The first analysis compares the length of full protein sequences in the training and benchmark sets. Since a small number of very long proteins can dominate the scale, filtered views were also generated to make the main distribution easier to read.

### Training set
![Protein length distribution in training set](./1,2.Sequence_lengths_comparison/sequence_length_training.png)

### Benchmark set
![Protein length distribution in benchmark set](./1,2.Sequence_lengths_comparison/sequence_length_benchmark.png)

### Distribution-focused views
![Training length distribution](./1,2.Sequence_lengths_comparison/Distribution_plot_training_distribution.png)

![Benchmark length distribution](./1,2.Sequence_lengths_comparison/Distribution_plot_benchmark.png)

![Protein length distribution across sets](./1,2.Sequence_lengths_comparison/Distribution_plot_set_division.png)

![Protein length distribution across subsets](./1,2.Sequence_lengths_comparison/Distribution_plot_set_division_2.png)

### Interpretation

These plots show that the positive and negative classes do not have identical length profiles. Negative proteins span a broader range and include more long sequences, while positive proteins tend to be somewhat more compact. The benchmark set follows the same general pattern seen in training, which suggests that the split is structurally consistent and suitable for later evaluation.

---

## 2. Signal peptide length

For positive proteins, the end position of the annotated signal peptide (`SPEnd`) was used to estimate signal peptide length. This analysis was included to verify that the positive class reflects the expected biological range for signal peptides.

### Training vs benchmark comparison
![Signal peptide length in training and benchmark sets](./1,2.Sequence_lengths_comparison/sp_length_train_vs_bench.png)

### Distribution across subsets
![Signal peptide length distribution across sets](./1,2.Sequence_lengths_comparison/Distribution_plot_SPE_distribution.png)

![Signal peptide length overview](./1,2.Sequence_lengths_comparison/Distribution_plot_distribution.png)

### Interpretation

Most signal peptides fall within the expected range, with the bulk of the distribution centered around roughly 20 to 25 residues. The similarity between training and benchmark distributions is important, because it indicates that the evaluation set was not built from a substantially different positive population.

---

## 3. Amino acid composition

The next step examined amino acid composition in signal peptide regions and compared it with reference Swiss-Prot frequencies. Additional plots summarize the same information at the level of broader residue categories and physicochemical groups.

### Comparison with Swiss-Prot background
![Training SP frequencies vs Swiss-Prot](./3.AA_Comparison/frequencies_sp_train_vs_swissprot.png)

![Benchmark SP frequencies vs Swiss-Prot categories](./1,2.Sequence_lengths_comparison/frequencies_sp_bench_vs_swissprot_aacategories.png)

### Composition across subsets
![Amino acid frequencies by set](./1,2.Sequence_lengths_comparison/AminoAcid_Frequencies_Set.png)

![Residue property frequencies by set](./1,2.Sequence_lengths_comparison/Frequencies_propierties_Set.png)

![Residue category distribution by set](./1,2.Sequence_lengths_comparison/Distribution_Set.png)

### Interpretation

The composition analysis highlights the expected enrichment of hydrophobic residues in signal peptide regions. This is one of the defining biological properties of signal peptides and is consistent with their role in membrane targeting and translocation. Compared with a Swiss-Prot background, the difference is clear enough to justify using residue composition and physicochemical features in later prediction models.

---

## 4. Taxonomic composition

To check whether the dataset remained biologically representative after preprocessing and splitting, proteins were grouped by kingdom and by the most represented organisms.

### Kingdom-level distribution
![Kingdom distribution in training set](./4.Taxonomy_classification/kingdom_barplot_train.png)

![Kingdom distribution in benchmark set](./4.Taxonomy_classification/kingdom_barplot_bench.png)

![Kingdom proportions in training set](./4.Taxonomy_classification/kingdom_pie_train.png)

![Kingdom proportions in benchmark set](./4.Taxonomy_classification/kingdom_pie_bench.png)

### Species-level distribution
![Species distribution in training set](./4.Taxonomy_classification/species_barplot_train.png)

![Species distribution in benchmark set](./4.Taxonomy_classification/species_barplot_bench.png)

![Species proportions in training set](./4.Taxonomy_classification/species_pie_train.png)

![Species proportions in benchmark set](./4.Taxonomy_classification/species_pie_bench.png)

### Interpretation

The dataset is dominated by eukaryotic groups, especially Metazoa, with additional representation from Fungi, Viridiplantae, and other eukaryotic categories. The benchmark set broadly mirrors the training set, which reduces the risk that downstream performance differences are caused by an unintended taxonomic shift rather than by the predictive difficulty of the task itself.

---

## 5. Cleavage-site sequence logos

Finally, cleavage-site sequence logos were generated from aligned windows around the annotated cleavage positions. These plots summarize residue preferences in the region most directly related to signal peptide processing.

### Training set logo
![Training cleavage-site sequence logo](./5.SequenceLogo/sequence_logo_training.png)

### Benchmark set logo
![Benchmark cleavage-site sequence logo](./5.SequenceLogo/sequence_logo_benchmark.png)

### Alternative visualizations
![Colored benchmark logo](./5.SequenceLogo/sequence_logo_benchmark_colored.png)

![Mirrored sequence logo](./5.SequenceLogo/mirrored_logo_fixed.png)

### Interpretation

The sequence logos show clear positional preferences around the cleavage site, including enrichment of small, neutral residues near the cleavage position. This agrees with the known biological rules governing signal peptide cleavage and supports the quality of the positive annotations. It also provides a direct bridge to the next stage of the project, where position-specific patterns are used in the von Heijne approach.

---

## Final remarks

Taken together, these analyses show that the dataset is coherent, biologically plausible, and well structured for the following steps of the project. The main trends are consistent with what is expected for signal peptides:

- positive proteins contain realistic signal peptide lengths
- signal peptide regions are compositionally distinct from general protein background
- cleavage-site motifs show recognizable conservation
- training and benchmark sets remain broadly comparable after preprocessing

This step therefore serves as a descriptive validation of the dataset and as a foundation for feature extraction, classical scoring approaches, and machine-learning models.

---

## Files used in this step

**Main notebook**
- `plots.ipynb`

**Input dataset**
- `../2.Data_Preparation/train_bench.tsv`

**Figure folders**
- `1,2.Sequence_lengths_comparison/`
- `3.AA_Comparison/`
- `4.Taxonomy_classification/`
- `5.SequenceLogo/`
