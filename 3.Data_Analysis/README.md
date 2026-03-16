# Protein Sequence Length and Signal Peptide Analysis

To explore and visualize several properties of protein sequences in both **training** and **benchmarking** datasets the script *plots.ipynb* has been produced.
The script analyses focus on sequence lengths, signal peptide (SP) properties, amino acid composition, taxonomic distribution, and SP cleavage sites.

---

## Aim
This analysis is helpful to: 
- Detect dataset **biases** (length, amino acid composition, taxonomy).
- Assess whether sequence properties are **informative features** for classification tasks.
- Provide **biological insights** into signal peptide characteristics.
- Evaluate differences between **training** and **benchmarking** datasets.

---

### Requirements
To run the notebook, the following libraries are needed. 
  - `pandas`
  - `numpy`
  - `matplotlib`
  - `seaborn`
  - `logomaker` (for sequence logos)

Install them with:
```bash
pip install pandas numpy matplotlib seaborn logomaker
```
---
## Analysis

The following elements were the subject of our analysis

### **1. Protein Length Distribution**
- **Description**: Compares the lengths of protein sequences between **positive** and **negative** classes.  The goal was to identify potential differences in sequence length distributions across classes in both training and benchmarking datasets.

### **2. Signal Peptide (SP) Length Distribution**
- **Description**: Examines the lengths of signal peptides (SPs), to investigate whether SP lengths differ between classes or datasets.

### **3. Amino Acid Composition Analysis**
- **Description**: Compares the amino acid (AA) composition of SP sequences against a **reference background distribution** (SwissProt database, see [SwissProt statistics](https://web.expasy.org/docs/relnotes/relstat.html)), to highlight amino acid biases characteristic of SPs compared to general proteins.

### **4. Taxonomic Classification**
- **Description**: Analyzes the taxonomic origin of proteins at both **kingdom** and **species** levels. to explore how sequences are distributed across taxa, checking for possible dataset biases

### **5. Sequence Logos of SP Cleavage Sites**
- **Description**: Visualizes sequence motifs at SP cleavage sites. to identify conserved motifs and patterns around cleavage sites.
---


