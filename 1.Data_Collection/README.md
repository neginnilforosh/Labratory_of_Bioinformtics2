Data Collection

Data collection represents the first stage of the project and is a critical prerequisite for any machine learning workflow. The quality, reliability, and consistency of the input data directly influence the performance and robustness of the predictive models. For this reason, protein sequences were retrieved from the UniProt￼ database and organized into two distinct datasets: a positive dataset, composed of proteins containing a signal peptide, and a negative dataset, composed of proteins lacking a signal peptide.

The data collection strategy was carried out in two sequential phases:
	1.	Identification of relevant proteins through the UniProt web interface
	2.	Programmatic retrieval through the UniProt API, based on the advanced search queries previously defined

⸻

1. Retrieval through the UniProt Web Interface

As an initial step, the Advanced Search function provided by UniProt (Release 2025_03) was used to define the selection criteria for both datasets.

Positive Dataset

Query

(existence:1) AND (length:[40 TO ]) AND (reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) AND (ft_signal_exp:)

The search parameters were selected as follows:
	•	existence:1: proteins with evidence at the protein level
	•	length:[40 TO ]: proteins with sequence length greater than or equal to 40 amino acids
	•	reviewed:true: manually reviewed entries from UniProtKB/Swiss-Prot
	•	taxonomy_id:2759: proteins belonging to the taxonomic group Eukaryota
	•	fragment:false: exclusion of fragmentary or incomplete sequences
	•	ft_signal_exp:: proteins annotated with an experimentally validated signal peptide

Result

This query returned curated eukaryotic proteins of at least 40 amino acids in length, supported by experimental evidence, non-fragmentary, and annotated with an experimentally confirmed signal peptide.
Starting from this search, the corresponding UniProt API endpoint for the positive dataset was obtained. The search endpoint is lightweight, returns results in batches of 500 entries, and therefore requires pagination.

https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28existence%3A1%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28reviewed%3Atrue%29+AND+%28taxonomy_id%3A2759%29+AND+%28fragment%3Afalse%29+AND+%28ft_signal_exp%3A*%29%29&size=500

Negative Dataset

Query

(existence:1) AND (length:[40 TO ]) AND (reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) NOT (ft_signal:) AND 
((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))

The search parameters were selected as follows:
	•	existence:1: proteins with evidence at the protein level
	•	length:[40 TO ]: proteins with sequence length greater than or equal to 40 amino acids
	•	reviewed:true: manually reviewed entries from UniProtKB/Swiss-Prot
	•	taxonomy_id:2759: proteins belonging to Eukaryota
	•	fragment:false: exclusion of fragmentary or incomplete sequences
	•	NOT (ft_signal:): exclusion of proteins annotated with a signal peptide
	•	cc_scl_term_exp: experimentally validated subcellular localization terms, including:
	•	SL-0091: Cytoplasm
	•	SL-0191: Nucleus
	•	SL-0173: Mitochondrion
	•	SL-0209: Plastid
	•	SL-0204: Peroxisome
	•	SL-0039: Cellular membrane

Result

This query returned curated eukaryotic proteins of at least 40 amino acids in length, experimentally supported, non-fragmentary, lacking a signal peptide, and experimentally localized to one of the selected cellular compartments.
The corresponding UniProt API endpoint for the negative dataset was then derived from the web query. As in the previous case, the endpoint returns data in batches of 500 entries and requires pagination.

https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment%3Afalse%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28taxonomy_id%3A2759%29+NOT+%28ft_signal%3A*%29+AND+%28%28cc_scl_term_exp%3ASL-0091%29+OR+%28cc_scl_term_exp%3ASL-0191%29+OR+%28cc_scl_term_exp%3ASL-0173%29+OR+%28cc_scl_term_exp%3ASL-0209%29+OR+%28cc_scl_term_exp%3ASL-0204%29+OR+%28cc_scl_term_exp%3ASL-0039%29%29+AND+%28reviewed%3Atrue%29+AND+%28existence%3A1%29%29&size=500


⸻

2. Retrieval through the UniProt API

To automate the retrieval process and refine the datasets, two Jupyter notebooks were developed:
	•	get_dataset_neg.ipynb￼
Used for the positive dataset
	•	get_dataset_pos.ipynb￼
Used for the negative dataset

These notebooks perform the API requests to UniProt, export the results in both .tsv and .fasta formats, and apply additional filtering criteria.

The additional filtering steps were as follows:
	•	Positive dataset: proteins were excluded if the signal peptide was shorter than 14 residues or if no cleavage site was reported
	•	Negative dataset: proteins were checked for the presence of a transmembrane helix beginning within the first 90 residues

⸻

3. Data Collection Outputs

Following the API-based retrieval, the datasets were stored in two output formats:

Tab-Separated Files (.tsv)

These files contain the metadata associated with each protein entry.
	•	Positive dataset
	•	UniProt accession
	•	Organism name
	•	Eukaryotic kingdom
	•	Protein length
	•	Signal peptide cleavage-site position
	•	Negative dataset
	•	UniProt accession
	•	Organism name
	•	Eukaryotic kingdom
	•	Protein length
	•	Presence of a transmembrane helix within the first 90 residues

FASTA Files (.fasta)

These files contain the amino acid sequences of the retrieved proteins.

⸻

4. Summary of Retrieved Data

The number of proteins obtained after the retrieval and filtering steps is summarized below.

Positive Set	Negative Set	Negative with HD
2932	20615	1384


⸻

5. Concluding Remarks

This data collection procedure ensured the construction of two high-quality and biologically relevant datasets suitable for downstream analysis and model development. By combining manual query design through the UniProt web interface with automated retrieval via the UniProt API, it was possible to obtain curated protein sets while maintaining reproducibility and control over the selection criteria.
