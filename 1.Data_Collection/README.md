# Data Collection

The first step of the project is **data collection**: for any machine learning procedure, it is necessary to acquire high-quality data to ensure correct model training.  
We retrieved relevant protein sequences from the [UniProt](https://www.uniprot.org/) database, filtering them into two different datasets: a **positive** one containing proteins with the signal peptide, and a **negative** set of proteins without it.  

Our approach was divided into two consecutive steps:
1. Retrieval of proteins using the **web interface**.  
2. Retrieval via the **UniProt API**, starting from the advanced search queries.  

---

## 1. Web Interface Approach

We started using the *Advanced Search* interface in UniProt (Release 2025_03) to filter and retrieve datasets.

### Positive Set

**QUERY:** 
```
(existence:1) AND (length:[40 TO ]) AND (reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) AND (ft_signal_exp:)
```
- *existence:1* → protein existence at protein level (experimental evidence)  
- *length:[40 TO ]* → sequence length ≥ 40 aa  
- *reviewed:true* → only proteins reviewed in UniProtKB/Swiss-Prot  
- *taxonomy_id:2759* → proteins from the taxon Eukaryota  
- *fragment:false* → exclude fragments (only complete proteins)  
- *ft_signal_exp:* → proteins with feature signal peptide (experimental evidence)  

**Result:** 

Curated eukaryotic proteins, ≥40 aa, experimentally confirmed, non-fragment, with an experimentally validated signal peptide. For uniprot the API was retrived. 
The API URL using the search endpoint for positive set was retrived. This endpoint is lighter and returns chunks of 500 at a time and requires pagination: 

https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28existence%3A1%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28reviewed%3Atrue%29+AND+%28taxonomy_id%3A2759%29+AND+%28fragment%3Afalse%29+AND+%28ft_signal_exp%3A*%29%29&size=500

### Negative set

**QUERY:**
 ```
 (existence:1) AND (length:[40 TO ]) AND (reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) NOT (ft_signal:) AND 
 ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))
```
- *existence:1* → protein existence at protein level (experimental evidence)  
- *length:[40 TO ]* → sequence length ≥ 40 aa  
- *reviewed:true* → Swiss-Prot (manually curated)  
- *taxonomy_id:2759* → Eukaryota  
- *fragment:false* → exclude fragments  
- *NOT (ft_signal:)* → no signal peptide  
- *cc_scl_term_exp:*  
  - SL-0091 → Cytoplasm  
  - SL-0191 → Nucleus  
  - SL-0173 → Mitochondrion  
  - SL-0209 → Plastid  
  - SL-0204 → Peroxisome  
  - SL-0039 → Cellular membrane
  
**Result:**  

Curated eukaryotic proteins, ≥40 aa, experimentally confirmed, non-fragment, without signal peptide, localized experimentally to one of the listed compartments.
The API URL using the search endpoint for negative set was downloaded. This endpoint is lighter and returns chunks of 500 at a time and requires pagination.

https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment%3Afalse%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28taxonomy_id%3A2759%29+NOT+%28ft_signal%3A*%29+AND+%28%28cc_scl_term_exp%3ASL-0091%29+OR+%28cc_scl_term_exp%3ASL-0191%29+OR+%28cc_scl_term_exp%3ASL-0173%29+OR+%28cc_scl_term_exp%3ASL-0209%29+OR+%28cc_scl_term_exp%3ASL-0204%29+OR+%28cc_scl_term_exp%3ASL-0039%29%29+AND+%28reviewed%3Atrue%29+AND+%28existence%3A1%29%29&size=500

---

## 2. API approach

Two Python scripts were created for our customised search:

[get_dataset_neg.ipynb](./get_dataset_neg.ipynb)
 → for the positive set

[get_dataset_pos.ipynb](./get_dataset_pos.ipynb)
 → for the negative set

These scripts perform the API calls to UniProt, generate output files in both .tsv and .fasta formats, and apply additional filtering:

- **Positive set**: filtering out proteins with a signal peptide shorter than 14 residues and without a cleavage site.

- **Negative set**: checking for proteins with a transmembrane helix starting within the first 90 residues.

---
  
## 3. Data collection output

After the API calls, our data were saved in:
- **.tsv** format: containing metadata for each protein.
  1. **Positive set:** UniProt accession, organism name, Eukaryotic kingdom, protein length, signal peptide cleavage site position.
  2. **Negative set:** UniProt accession, organism name, Eukaryotic kingdom, protein length, presence of a transmembrane helix within the first 90 residues.

- **.fasta** format: containing the protein sequences.
The number of proteins retrived after the search is showed in the table below. 


| Positive Set | Negative Set | Negative with HD | 
|--------------|--------------|------------------|
|  2932        |    20615     |      1384        |
































