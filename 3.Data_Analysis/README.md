# 📊 Step 3 — Data Analysis

> **What this folder does:** Explores the dataset visually before any modelling begins — checking protein lengths, signal peptide properties, amino acid composition, taxonomic diversity, and sequence patterns around cleavage sites. This step confirms the data is biologically sensible and ready for downstream models.

---

## 🗂️ Files in this Folder

| File | What it does |
|------|-------------|
| `plots.ipynb` | All analyses and figure generation |

### Output figure folders (created when notebook is run)

| Folder | Contents |
|--------|----------|
| `Sequence_lengths_comparison/` | Protein length and SP length distribution plots |
| `3.AA_Comparison/` | Amino acid frequency comparisons |
| `4.Taxonomy_classification/` | Pie charts and bar plots of kingdom/species composition |
| `5.SequenceLogo/` | Sequence logos around the cleavage site |

---

## 📥 Input

All analyses read from one file:

```
../2.Data_Preparation/train_bench.tsv
```

---

## 🔬 Analysis 1 — Protein Sequence Length Distribution

**Question:** Are positive (SP) and negative (non-SP) proteins the same length, or systematically different?

**How it works:** The `SequenceLength` column is extracted and plotted separately for the training set and the benchmark set. To avoid a few very long proteins distorting the scale, the plots focus on sequences ≤2,500 amino acids (outliers are shown separately).

**Figures generated:**

| Figure | Description |
|--------|-------------|
| `Sequence_lengths_comparison/Distribution_plot_training_distribution.png` | Density histogram of protein lengths — Training set, Positive vs Negative |
| `Sequence_lengths_comparison/Distribution_plot_benchmark.png` | Same comparison for the Benchmark set |
| `Sequence_lengths_comparison/Violin_plot_set_division.png` | Violin plot showing length distribution across all 5 folds + benchmark |

**What to look for:** Negative proteins tend to have a wider length range and more very long sequences. Positive proteins (secreted) tend to be more compact. The training and benchmark sets should follow the same pattern — confirming the split is representative.

---

## 🔬 Analysis 2 — Signal Peptide Length Distribution

**Question:** How long are the signal peptides in the dataset? Are the training and benchmark sets consistent?

**How it works:** Only positive-class proteins (those with a signal peptide) are used. The `SPEnd` column gives the cleavage site position, which serves as a proxy for SP length. KDE curves and histogram plots are generated for training vs. benchmark.

**Figures generated:**

| Figure | Description |
|--------|-------------|
| `Sequence_lengths_comparison/sp_length_train_vs_bench.png` | SP length probability distribution — Training vs Benchmark |
| `Sequence_lengths_comparison/Distribution_plot_SPE_distribution.png` | KDE of SP lengths across all 6 sets (Folds 1–5 + Benchmark) |

**What to look for:** Signal peptides should cluster between 15 and 30 amino acids. Both training and benchmark distributions should overlap closely, confirming that the 80/20 split does not bias the SP length profile.

---

## 🔬 Analysis 3 — Amino Acid Composition

**Question:** Do signal peptides have a distinctive amino acid composition compared with typical proteins?

**How it works:**
1. For each positive protein, the subsequence from `SPStart` to `SPEnd` is extracted (the actual signal peptide)
2. The frequency of all 20 amino acids is computed across all signal peptides
3. These frequencies are compared to SwissProt background frequencies (the typical composition of all reviewed proteins)

**SwissProt background frequencies used:**

| AA | Freq | AA | Freq | AA | Freq | AA | Freq |
|----|------|----|------|----|------|----|------|
| A (Ala) | 8.25% | G (Gly) | 7.07% | L (Leu) | 9.64% | S (Ser) | 6.65% |
| R (Arg) | 5.52% | H (His) | 2.27% | K (Lys) | 5.79% | T (Thr) | 5.36% |
| N (Asn) | 4.06% | I (Ile) | 5.90% | M (Met) | 2.41% | W (Trp) | 1.10% |
| D (Asp) | 5.46% | — | — | F (Phe) | 3.86% | Y (Tyr) | 2.92% |
| C (Cys) | 1.38% | — | — | P (Pro) | 4.74% | V (Val) | 6.85% |

Amino acids are grouped by chemical class on the x-axis (colour-coded):

| Colour | Class | Residues |
|--------|-------|---------|
| 🔵 Blue | Nonpolar aliphatic | G, A, V, L, I, M, P |
| 🟠 Orange | Aromatic | F, Y, W |
| 🟢 Green | Polar uncharged | S, T, C, N, Q |
| 🔴 Red | Charged | D, E, K, R |

**Figures generated:**

| Figure | Description |
|--------|-------------|
| `3.AA_Comparison/frequencies_sp_train_vs_swissprot.png` | SP amino acid frequencies (training) vs SwissProt background |
| `Sequence_lengths_comparison/frequencies_sp_bench_vs_swissprot_aacategories.png` | Same comparison grouped by chemical class |
| `Sequence_lengths_comparison/AminoAcid_Frequencies_Set.png` | AA frequencies across all 5 training folds |

**What to look for:** Signal peptides are expected to be strongly enriched in hydrophobic residues (especially L, A, V, I) compared to the SwissProt background. This hydrophobic enrichment is what allows the signal peptide to interact with the membrane translocon machinery. Seeing this clearly confirms the biological relevance of the dataset.

---

## 🔬 Analysis 4 — Taxonomic Classification

**Question:** Which organisms are represented in the dataset? Is the taxonomic composition consistent between training and benchmark sets?

**How it works:** Proteins are grouped by `Kingdom` and `OrganismName`. Because there are hundreds of species, only the **6 most frequent organisms** are kept as individual labels; all others are grouped as "Other". Pie charts and bar plots are generated for both training and benchmark subsets.

**Figures generated:**

| Figure | Description |
|--------|-------------|
| `4.Taxonomy_classification/kingdom_barplot_train.png` | Bar plot of kingdom distribution — Training set |
| `4.Taxonomy_classification/kingdom_barplot_bench.png` | Bar plot of kingdom distribution — Benchmark set |
| `4.Taxonomy_classification/kingdom_pie_train.png` | Pie chart of kingdom fractions — Training set |
| `4.Taxonomy_classification/species_barplot_train.png` | Bar plot of top 6 species + Other — Training |
| `4.Taxonomy_classification/species_pie_bench.png` | Pie chart of species — Benchmark |

**What to look for:** The dataset should be dominated by **Metazoa** (animals), with smaller contributions from **Fungi**, **Viridiplantae** (plants), and other eukaryotic groups. The training and benchmark sets should show similar proportions — if one set has dramatically different species composition, the benchmark would not be a fair test.

---

## 🔬 Analysis 5 — Cleavage Site Sequence Logos

**Question:** Are there conserved sequence patterns around signal peptide cleavage sites?

**How it works:**
1. For each positive protein, the subsequence from **position SPEnd−13 to SPEnd+2** is extracted — a 15-residue window centred just before the cleavage site
2. All such subsequences are aligned and fed to **logomaker** to build a sequence logo
3. The logo is computed using **information content** (measured in bits): positions where one amino acid is strongly preferred show tall stacks; positions with uniform distribution show flat, short stacks

```python
# Extract 15-residue window around cleavage site
for index, row in train_pos.iterrows():
    sequence  = row["Sequence"]
    cleavage  = int(row["SPEnd"])
    window    = sequence[cleavage - 13 : cleavage + 2]  # -13 to +2
    train_seqs.append(window)

# Build information-content logo
logo_matrix = lm.alignment_to_matrix(sequences=train_seqs,
                                      to_type='information',
                                      characters_to_ignore='.-X')
```

**Figures generated:**

| Figure | Description |
|--------|-------------|
| `5.SequenceLogo/sequence_logo_training.png` | Sequence logo from training set sequences |
| `5.SequenceLogo/sequence_logo_benchmark.png` | Sequence logo from benchmark sequences |
| `5.SequenceLogo/sequence_logo_benchmark_colored.png` | Coloured version of the benchmark logo |
| `5.SequenceLogo/mirrored_logo_fixed.png` | Mirror-style logo for extended visualisation |

**What to look for:** The canonical signal peptide cleavage motif is the **AXA rule** — small, neutral amino acids (especially Ala) at positions −3 and −1 relative to the cleavage site. The logo should show:
- High information content at positions −3 and −1 (Ala preferred)
- Hydrophobic residues (L, A, V, I) dominating positions −10 to −4 (the hydrophobic H-region of the SP)
- A drop in conservation immediately after the cleavage site


