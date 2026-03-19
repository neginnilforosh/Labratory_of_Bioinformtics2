# 🧹 Step 2 — Data Preparation

> **What this folder does:** Takes the raw downloaded proteins, removes near-duplicate sequences, then splits everything into a training set (80%) and a held-out benchmark set (20%) ready for model training.

---

## 🗂️ Files in this Folder

| File | What it does |
|------|-------------|
| `recover_ids_clusters.ipynb` | Reads MMseqs2 clustering output and extracts one representative sequence per cluster |
| `prepare_datasets.ipynb` | Randomly splits data into train/benchmark sets, adds sequences, and builds the final `train_bench.tsv` |
| `output_recap.ipynb` | Prints a summary of the final dataset counts |

---

## ❓ Why Is Data Preparation Necessary?

Raw protein databases contain many near-identical sequences (e.g. the same protein from closely related species, or slightly different isoforms). If these near-duplicates appear in both training and testing sets, the model appears to perform better than it really does — it has essentially "seen" the test data. This is called **data leakage**.

The preparation pipeline solves this with three steps:

```
Raw proteins
    │
    ▼
1. Redundancy reduction (MMseqs2 clustering)
    │
    ▼
2. 80/20 train/benchmark split
    │
    ▼
3. 5-fold cross-validation subdivision
    │
    ▼
train_bench.tsv  ← single file used by all downstream models
```

---

## 🧬 Step 1 — Redundancy Reduction with MMseqs2

**MMseqs2** (version 14.7e284) is a fast and memory-efficient tool for clustering large sets of protein sequences.

### Clustering parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Sequence identity | **30%** | Two sequences are in the same cluster only if they share ≥30% identical residues |
| Coverage | **40%** | The alignment must cover ≥40% of both sequences |

Both conditions must be satisfied simultaneously. Using a 30% identity threshold means that even distant homologues are grouped together — preventing the model from simply memorising sequence similarity.

### What happens during clustering

1. MMseqs2 groups all positive sequences into clusters, then all negative sequences into clusters separately
2. One **representative sequence** is chosen per cluster (typically the longest or most central one)
3. All other sequences in the cluster are discarded

### Input/Output files

| Input | Output |
|-------|--------|
| `positive_dataset.fasta` | `pos_cluster/pos_cluster.tsv` — which sequence belongs to which cluster |
| `negative_dataset.fasta` | `pos_cluster/pos_rep_seq.fasta` — one representative sequence per cluster |

The notebook `recover_ids_clusters.ipynb` reads these clustering outputs and creates clean TSV files with only the representative entries:
- `neg_cluster/negative_representative_data.tsv`
- `pos_cluster/positive_representative_data.tsv`

### Numbers after redundancy reduction

| Class | Before clustering | After clustering |
|-------|:-----------------:|:----------------:|
| Positive | 2,932 | **1,093** |
| Negative | 20,615 | **8,934** |
| Negative with N-term helix | 1,384 | **636** |

> Removing ~63% of positive sequences and ~57% of negative sequences shows how much redundancy existed in the raw UniProt downloads.

---

## ✂️ Step 2 — Train / Benchmark Split

The `prepare_datasets.ipynb` notebook performs a **random 80/20 split** on the de-redundant sequences:

```python
# Shuffle randomly first
pos_random = pos_id.sample(frac=1).reset_index(drop=True)
neg_random = neg_id.sample(frac=1).reset_index(drop=True)

# 80% → training,  20% → benchmark
n_pos = int(len(pos_random) * 0.8)  # 874 positives for training
n_neg = int(len(neg_random) * 0.8)  # 7,147 negatives for training

training_pos   = pos_random[:n_pos]
benchmark_pos  = pos_random[n_pos:]

training_neg   = neg_random[:n_neg]
benchmark_neg  = neg_random[n_neg:]
```

A sanity check confirms there is **zero overlap** between training and benchmark sets.

### Why keep a separate benchmark set?

Cross-validation gives an estimate of training performance, but it is still computed on data the model has seen. The benchmark set is:
- **Never seen during training**
- **Never seen during hyperparameter tuning**
- Used only once at the very end to report final model performance

This provides the strongest guarantee of unbiased evaluation.

---

## 📂 Step 3 — 5-Fold Cross-Validation Subdivision

The training set is further subdivided into **5 equal folds** for cross-validation. Each fold maintains the original positive/negative class ratio.

```
Training set (80%)
    │
    ├── Fold 1 ─── Set = "1"
    ├── Fold 2 ─── Set = "2"
    ├── Fold 3 ─── Set = "3"
    ├── Fold 4 ─── Set = "4"
    └── Fold 5 ─── Set = "5"
```

In each cross-validation iteration:
- **3 folds** → training
- **1 fold** → validation (threshold or hyperparameter selection)
- **1 fold** → testing

This rotation ensures every protein is used for testing exactly once.

---

## 📄 Final Output: `train_bench.tsv`

All data — training folds and benchmark — is merged into a single file. This is the **master file used by every subsequent step** in the project.

### Column descriptions

| Column | Type | Description |
|--------|------|-------------|
| `EntryID` | string | UniProt accession code (e.g. `P12345`) |
| `OrganismName` | string | Scientific name of the source organism |
| `Kingdom` | string | Eukaryotic group (Metazoa, Fungi, Viridiplantae, etc.) |
| `SequenceLength` | integer | Total length of the protein in amino acids |
| `HelixDomain` | boolean / NaN | `True` if a transmembrane helix starts within the first 90 residues; `NaN` for positive entries |
| `Class` | string | `Positive` (has SP) or `Negative` (no SP) |
| `SPStart` | integer / NaN | Position where the signal peptide starts (1-indexed); `NaN` for negatives |
| `SPEnd` | integer / NaN | Position where the signal peptide ends / cleavage site; `NaN` for negatives |
| `Set` | string | `"1"` to `"5"` for training folds, `"Benchmark"` for the holdout set |
| `Sequence` | string | Full amino acid sequence |

### Final dataset counts

| Subset | Negatives | Positives | Total |
|--------|:---------:|:---------:|:-----:|
| Training (Folds 1–5) | 7,147 | 874 | 8,021 |
| Benchmark | 1,787 | 219 | 2,006 |
| **Total** | **8,934** | **1,093** | **10,027** |

---

## 🔁 How to Run

Run the notebooks in this order:

1. **`recover_ids_clusters.ipynb`** — must be run after MMseqs2 has already been executed on the FASTA files (MMseqs2 is a command-line tool run outside Jupyter)
2. **`prepare_datasets.ipynb`** — creates `train_bench.tsv`
3. **`output_recap.ipynb`** — verifies the final counts

> ⚠️ Do **not** re-run the randomisation cell (`pos_id.sample(frac=1)`) after the initial run. The split must remain fixed for all downstream experiments to be comparable. The cell is marked with a `# DO NOT RERUN THIS CELL` comment.
