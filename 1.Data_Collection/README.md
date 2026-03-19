# 📥 Step 1 — Data Collection

> **What this folder does:** Downloads protein sequences from the UniProt database and organises them into two groups — proteins **with** a signal peptide (positive set) and proteins **without** one (negative set).

---

## 🗂️ Files in this Folder

| File | What it does |
|------|-------------|
| `get_dataset_pos.ipynb` | Downloads and filters all positive (signal peptide) proteins from UniProt |
| `get_dataset_neg.ipynb` | Downloads and filters all negative (no signal peptide) proteins from UniProt |
| `output_recap.ipynb` | Prints a short summary of how many proteins were collected |

---

## 🧬 What is a Signal Peptide?

A **signal peptide** is a short stretch of amino acids (typically 15–30 residues) at the very beginning (N-terminus) of a protein. Its job is to act as a molecular address label, directing the newly made protein to the secretory pathway so it can be exported outside the cell or embedded in a membrane. After the protein reaches its destination, the signal peptide is cut off by an enzyme called a signal peptidase.

Proteins **without** a signal peptide stay inside the cell in compartments like the nucleus, cytoplasm, or mitochondria.

---

## 🔍 Where Does the Data Come From?

All data come from **UniProt** (Release 2025_03) — the world's most comprehensive and well-annotated protein database. Only manually reviewed entries from **UniProtKB/Swiss-Prot** are used, ensuring high data quality.

Two separate searches were designed using the UniProt Advanced Search interface, then automated via the UniProt REST API.

---

## ✅ Positive Dataset (Proteins WITH a Signal Peptide)

### Search criteria

| Filter | Meaning |
|--------|---------|
| `existence:1` | Protein exists at the protein level (strongest evidence) |
| `length:[40 TO *]` | At least 40 amino acids long |
| `reviewed:true` | Manually reviewed (Swiss-Prot quality) |
| `taxonomy_id:2759` | Belongs to Eukaryota |
| `fragment:false` | Complete sequence, not a fragment |
| `ft_signal_exp:` | Has an **experimentally confirmed** signal peptide annotation |

### API endpoint used

```
https://rest.uniprot.org/uniprotkb/search?format=json&query=
  (existence:1) AND (length:[40 TO *]) AND (reviewed:true) AND
  (taxonomy_id:2759) AND (fragment:false) AND (ft_signal_exp:*)
&size=500
```

Results come back in batches of 500 entries. The notebook handles pagination automatically — it follows the `Link` header in each response to fetch the next batch until all entries are retrieved.

### Additional filtering applied in code

After downloading, two extra checks are applied:
- **Minimum SP length:** entries where the signal peptide is **shorter than 14 residues** are removed (too short to be a reliable signal peptide annotation)
- **Cleavage site must be defined:** entries with no reported cleavage position are removed

### Output files

| File | Contents |
|------|----------|
| `positive_dataset.tsv` | EntryID, OrganismName, Kingdom, SequenceLength, SPStart, SPEnd |
| `positive_dataset.fasta` | Amino acid sequences in FASTA format |

---

## ❌ Negative Dataset (Proteins WITHOUT a Signal Peptide)

### Search criteria

| Filter | Meaning |
|--------|---------|
| `existence:1` | Protein exists at the protein level |
| `length:[40 TO *]` | At least 40 amino acids long |
| `reviewed:true` | Manually reviewed (Swiss-Prot quality) |
| `taxonomy_id:2759` | Belongs to Eukaryota |
| `fragment:false` | Complete sequence |
| `NOT (ft_signal:)` | **Explicitly has no signal peptide** |
| `cc_scl_term_exp:SL-xxxx` | Experimentally localised to one of the intracellular compartments listed below |

### Why are subcellular localisation filters important?

Simply excluding proteins annotated with a signal peptide would not be enough — it might include proteins with no localisation information at all. To ensure the negative set contains proteins that are **genuinely non-secretory**, only proteins experimentally found in specific intracellular compartments are selected:

| Term | Compartment |
|------|-------------|
| SL-0091 | Cytoplasm |
| SL-0191 | Nucleus |
| SL-0173 | Mitochondrion |
| SL-0209 | Plastid |
| SL-0204 | Peroxisome |
| SL-0039 | Cellular membrane |

### Additional check applied in code

For each negative protein, the code checks whether a **transmembrane helix begins within the first 90 residues**. If it does, the `HelixDomain` column is set to `True`. This flag is important later — transmembrane proteins can be confused with signal-peptide proteins because their N-terminal hydrophobic helix looks similar to a signal peptide hydrophobic core.

### Output files

| File | Contents |
|------|----------|
| `negative_dataset.tsv` | EntryID, OrganismName, Kingdom, SequenceLength, HelixDomain |
| `negative_dataset.fasta` | Amino acid sequences in FASTA format |

---

## 🔄 How the API Pagination Works

Both notebooks use the same approach to download all results:

```python
# Retry strategy — handles temporary server errors automatically
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session.mount("https://", HTTPAdapter(max_retries=retries))

# Pagination loop — follows the "next" link until no more pages remain
def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        yield response, total
        batch_url = get_next_link(response.headers)
```

Each response returns up to 500 entries. The `get_next_link()` function reads the HTTP `Link:` header to find the URL for the next batch.

---

## 📊 Results After Collection and Filtering

| Dataset | Count |
|---------|------:|
| Positive proteins (with SP) | **2,932** |
| Negative proteins (without SP) | **20,615** |
| Negative proteins with N-terminal transmembrane helix | **1,384** |

> The large imbalance between positive (2,932) and negative (20,615) sets is intentional and reflects biology — the majority of proteins are intracellular. This imbalance will be addressed in later steps through class weighting and careful evaluation metrics (MCC rather than accuracy).

---

## 🔁 How to Run

1. Open and run `get_dataset_pos.ipynb` — this creates `positive_dataset.tsv` and `positive_dataset.fasta`
2. Open and run `get_dataset_neg.ipynb` — this creates `negative_dataset.tsv` and `negative_dataset.fasta`
3. Open and run `output_recap.ipynb` to verify the counts

> ⚠️ An active internet connection is required. The download may take several minutes depending on server load.
