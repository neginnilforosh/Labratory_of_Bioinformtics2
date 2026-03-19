# 🧠 Step 6 — Deep Learning (SP-NN)

> **What this folder does:** Builds, trains, and evaluates a hybrid **CNN-LSTM neural network** that predicts signal peptides directly from raw amino acid sequences — no manual feature engineering required.

---

## 🗂️ Files in this Folder

| File | What it does |
|------|-------------|
| `lstm.ipynb` | Defines the model architecture, encoding utilities, training infrastructure, and runs automated hyperparameter optimisation via Ray Tune |
| `benchmark_test.ipynb` | Self-contained notebook: trains the final model with the best config, exports it to TorchScript, evaluates it on the benchmark set, and generates all error analysis figures |

> **Note:** `benchmark_test.ipynb` does **not** import anything from `lstm.ipynb`. All required classes and functions are redefined internally, so it can be run independently.

### Output figures (saved to `model_evaluation_DL/`)

| Figure | Description |
|--------|-------------|
| `Confusion_Matrix.png` | Heatmap showing TP, TN, FP, FN counts |
| `HD_prediction_1.png` | Transmembrane helix prevalence across TP/TN/FP/FN groups |
| `HD_prediction_2.png` | TM helix prevalence: FP only vs all Negatives |
| `SP_length_distribution.png` | SP length histogram — True Positives vs False Negatives |
| `Hydrophobicity_SP.png` | Hydrophobicity distribution (SP region) — TP vs FN |
| `Boxplot_Hydrophobicity_SP.png` | Boxplot version of the hydrophobicity comparison |
| `Mean_hydrophobicity.png` | Mean hydrophobicity bar chart per prediction group |
| `AA_frequencies_TP_vs_FN.png` | Full-sequence amino acid frequency — TP vs FN |
| `AA_frequencies_SP.png` | SP-region amino acid frequency — TP vs FN |
| `AA_frequency_Distribution.png` | 4-group AA frequency: FP / FN / all Positives / all Negatives |
| `FPR_Kingdom.png` | False positive rate pie chart per kingdom |
| `Pieplot_species.png` | Taxonomic composition of all Positives / FN / FP |
| `SeqLogo_IC_1.png` to `_4.png` | Sequence logos (information content) for TP, FN, TN, FP |
| `SeqLogo_Prob_1.png` to `_4.png` | Sequence logos (residue probability) for TP, FN, TN, FP |

---

## 🔤 Step 1 — Sequence Encoding

Before a protein sequence can be fed to a neural network, it must be converted to numbers. Two preprocessing steps are applied.

### 1a. Length standardisation (target: 90 residues)

Signal peptides are always at the N-terminus (the very beginning of the protein). The first 90 amino acids contain all the biologically relevant information, so every sequence is standardised to exactly 90 residues:

| Case | What happens |
|------|-------------|
| Sequence shorter than 90 | Pad with `X` characters appended to the end |
| Sequence longer than 90 | Truncate — keep only the first 90 residues |

The padding character `X` is treated as "unknown" and encoded as an all-zero vector — it carries no information and does not confuse the model.

### 1b. One-hot encoding (21 dimensions per position)

Each amino acid is converted to a **21-dimensional binary vector**:

```
Alphabet: A C D E F G H I K L M N P Q R S T V W Y X
Index:    0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

Example: Leucine (L) → [0,0,0,0,0,0,0,0,0, 1 ,0,0,0,0,0,0,0,0,0,0,0]
Example: Padding (X) → [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 1 ]
```

A full sequence of 90 residues becomes a **90 × 21 matrix** — this is the direct input to the CNN.

---

## 🏗️ Step 2 — Model Architecture (SP-NN)

SP-NN is a three-block hybrid architecture. Each block has a specific biological motivation:

```
Input  [Batch, 90, 21]   ← 90 positions, 21 channels per position
  │
  ▼
┌───────────────────────────────────────────────────────────────┐
│  BLOCK 1: Conv1D (Local Motif Detection)                      │
│  kernel=17, 64 feature maps, padding='same'                   │
│  [Batch, 90, 21]  →  [Batch, 90, 64]                         │
└───────────────────────────────────────────────────────────────┘
  │
  ▼
┌───────────────────────────────────────────────────────────────┐
│  BLOCK 2: LSTM × 2 (Sequential Context)                       │
│  hidden size=128, batch_first=True                            │
│  Takes last hidden state: [Batch, 90, 64]  →  [Batch, 128]   │
│  + BatchNorm1D(128)                                           │
└───────────────────────────────────────────────────────────────┘
  │
  ▼
┌───────────────────────────────────────────────────────────────┐
│  BLOCK 3: MLP Head (Classification)                           │
│  Linear(128→256) → ReLU → Dropout(0.498)                     │
│  Linear(256→128) → ReLU → Dropout(0.498)                     │
│  Linear(128→ 64) → ReLU → Dropout(0.498)  ← bottleneck       │
│  Linear( 64→1024)→ ReLU → Dropout(0.498)  ← expand           │
│  Linear(1024→  1) → Sigmoid                                   │
└───────────────────────────────────────────────────────────────┘
  │
  ▼
P(SP) ∈ [0, 1]  →  if > 0.5: Positive (has SP) else: Negative
```

### Why each block is designed this way

**Block 1 — Conv1D with kernel size 17:**
Signal peptides have a hydrophobic core that is typically 7–15 residues long. A convolution kernel of size 17 spans this entire region in one operation, allowing the CNN to detect the hydrophobic stretch as a pattern regardless of its exact starting position. 64 parallel feature maps capture diverse local sequence patterns simultaneously. `padding='same'` keeps the sequence length at 90 after convolution.

**Block 2 — LSTM × 2:**
After the CNN extracts local features, the LSTM processes the resulting 90-position feature sequence from left to right (N-terminus to position 90), building up a representation of the full sequential context. Only the **last hidden state** `out[:, -1, :]` is kept — after reading all 90 positions, this single vector summarises the entire sequence. Stacking 2 LSTM layers allows the model to build more abstract representations.

**Block 3 — Hourglass MLP:**
The `[256 → 128 → 64 → 1024]` pattern compresses features to a 64-dimensional bottleneck and then expands to 1024 before the final classification layer. This "compress then expand" structure was identified as optimal by the hyperparameter search — projecting compressed features into a higher-dimensional space improves class separation.

---

## 🛠️ Step 3 — Training Infrastructure

### Memory-efficient dataset loading (`SignalDataset`)

The one-hot encoded dataset is large. To avoid loading everything into RAM at once, `SignalDataset` stores only references to NumPy arrays. Tensors are created **only when a batch is requested**:

```python
class SignalDataset(Dataset):
    def __init__(self, X, y):
        self.X = X  # numpy array reference — not yet a Tensor
        self.y = y

    def __getitem__(self, idx):
        # Convert to Tensor only at access time
        x = torch.from_numpy(self.X[idx].copy()).float()
        y = torch.tensor(self.y[idx]).float().view(1)
        return x, y
```

### Training loop (`train_val`)

Three key mechanisms ensure stable, generalising training:

**1. Gradient clipping:**
LSTMs are prone to "exploding gradients" — if gradients grow too large, weight updates become destabilising. Gradient clipping caps the gradient vector's magnitude before each weight update:
```python
torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
```

**2. Early stopping:**
After every epoch, the validation MCC is computed. If it fails to improve for 20 consecutive epochs, training stops early:
```python
if val_MCC > best_val_MCC:
    best_val_MCC = val_MCC
    best_state_dict = model.state_dict()   # save the best weights
    epochs_without_improvement = 0
else:
    epochs_without_improvement += 1
    if epochs_without_improvement >= patience:
        break
```

**3. Best-model tracking:**
The function returns the `state_dict` (weights) from the **best epoch**, not the last epoch. This ensures the returned model generalises well even if later epochs start to overfit.

---

## ⚙️ Step 4 — Hyperparameter Optimisation (`lstm.ipynb`)

Manual hyperparameter selection for a model this complex would be impractical. **Ray Tune** automates this with 15 random trials, each evaluated with 5-fold cross-validation.

### Why all code is inside the trainable function

Ray Tune spawns isolated worker processes. If `SP_NN`, `SignalDataset`, etc. are defined outside the trainable function, the worker cannot import them (there is no `lstm.py` module). All class and function definitions are therefore placed **inside** `_tune_trainable` as closures, so Ray serialises them together.

### Memory sharing

The one-hot encoded dataset (~100+ MB) is stored once in Ray's shared memory:

```python
data_ref = ray.put(all_data)   # store once
del all_data                   # free local RAM
gc.collect()

# Inside each worker:
all_data = ray.get(config["data_ref"])  # cheap reference, not a copy
```

### Search space

| Parameter | What it controls | Range explored |
|-----------|-----------------|----------------|
| `num_layers` | Depth of the MLP head | {2, 3, 4, 5} |
| `hidden_sizes` | Width of each MLP layer | Sampled from {1024,512,256,128,64,32} |
| MLP topology | Funnel (descending) vs Flexible (shuffled) | 50/50 random |
| `dropout` | Dropout probability | [0.1, 0.5] uniform |
| `lr` | Learning rate | [1e-4, 1e-3] log-uniform |
| `batch_size` | Mini-batch size | {10, 20} |
| `num_lstm_layers` | Stacked LSTM layers | {1, 2} |
| `lstm_hidden_size` | LSTM hidden state dimension | {64, 128} |

Ray minimises `-mean_MCC` (which is equivalent to maximising mean_MCC across 5 folds).

---

## 📊 Step 5 — Best Configuration Found

| Parameter | Value |
|-----------|-------|
| LSTM layers | 2 |
| LSTM hidden size | 128 |
| MLP hidden layers | `[256, 128, 64, 1024]` |
| Dropout | 0.498 |
| Learning rate | 2.86 × 10⁻⁴ |
| Batch size | 20 |
| Optimiser | Adam |
| Loss function | BCELoss |

### Benchmark performance

| Metric | Value |
|--------|-------|
| **MCC** | **0.9019** |
| Decision threshold | 0.5 |

---

## 💾 Step 6 — TorchScript Model Export

The trained model is saved using **JIT tracing**, producing a standalone `.pt` file:

```python
model.eval()
dummy_input  = torch.randn(1, 90, 21).to(device)
traced_model = torch.jit.trace(model, dummy_input)
torch.jit.save(traced_model, "SignalPeptideLSTM.pt")
```

**Why not just save the weights (`state_dict`)?**

| `state_dict` | TorchScript `.pt` |
|---|---|
| Requires the original Python class definition to load | Self-contained — no source code needed |
| Python only | Works in Python, C++, Java, mobile |
| Dropout enabled at load time | `model.eval()` before tracing permanently disables Dropout |

**Loading the saved model:**
```python
model = torch.jit.load("SignalPeptideLSTM.pt")
model.eval()
```

---

## 🔎 Step 7 — Error Analysis (`benchmark_test.ipynb`)

After evaluation, every benchmark entry is labelled:

| Label | Meaning |
|-------|---------|
| **TP** | True Positive — signal peptide correctly identified |
| **TN** | True Negative — non-signal correctly rejected |
| **FP** | False Positive — non-signal wrongly predicted as signal peptide |
| **FN** | False Negative — real signal peptide missed by the model |

The analysis investigates *why* each error type occurs:

- **FP analysis:** Are false positives enriched for transmembrane domain proteins? Are certain taxonomic kingdoms over-represented?
- **FN analysis:** Do missed signal peptides have shorter or more polar sequences than correctly detected ones? Are their hydrophobic cores weaker?
- **Sequence logos:** Align sequences around the cleavage site to visualise what TP, FN, TN, FP sequences look like structurally.

