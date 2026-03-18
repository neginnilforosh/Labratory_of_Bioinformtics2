# Deep Learning — Hybrid CNN-LSTM Signal Peptide Classifier

This repository implements a hybrid **CNN-LSTM neural network** for binary signal peptide prediction from protein sequences. The pipeline covers data preprocessing, model definition, hyperparameter optimisation, benchmark evaluation, and detailed error analysis.

The code is organised into two self-contained Jupyter notebooks:

| Notebook | Role |
|---|---|
| `lstm.ipynb` | Defines all encoding utilities, the model architecture, the training/inference functions, and runs automated hyperparameter optimisation via Ray Tune |
| `benchmark_test.ipynb` | Fully self-contained — redefines all shared classes internally, trains the final model using the best configuration, exports it to TorchScript, and performs complete benchmark evaluation and error analysis |

> **Note:** `benchmark_test.ipynb` does **not** import from `lstm.ipynb`. All required classes and functions are defined directly inside it, making it runnable independently.

---

## Input Data

| File | Description |
|---|---|
| `train_bench.tsv` | Tab-separated file containing protein sequences, binary class labels (`Positive` / `Negative`), fold assignments (`Set` column: `1`–`5` and `Benchmark`), and metadata columns (`Kingdom`, `OrganismName`, `HelixDomain`, `SPStart`, `SPEnd`) |

---

## Pipeline Overview

### 1. Data Preprocessing

**File:** both notebooks — `create_one_hot_sets()` / `create_one_hot_sets_with_benchmark()`

The raw protein sequences undergo two standardisation steps before being fed to the network:

**Sequence length standardisation (fixed length: 90)**
- Sequences longer than 90 residues are **truncated**, keeping the N-terminal region which contains the signal peptide
- Sequences shorter than 90 residues are **padded** by appending `'X'` characters until length 90 is reached
- 90 residues was chosen because signal peptides are always located at the N-terminus and are typically 15–30 residues long; 90 captures all biologically relevant positional information

**One-hot encoding — `one_hot_encoding(seq)`**
- Each amino acid is converted to a 21-dimensional binary vector
- Alphabet order: `A C D E F G H I K L M N P Q R S T V W Y X`
- The 21st channel (`X`) encodes both the padding token and any unknown residues
- A sequence of length 90 becomes a **90 × 21 matrix** — the direct input to the CNN
- Unknown characters not in the alphabet map to an all-zero vector (safe fallback)

**Padding helper — `_increase_lenseq(seq, target=90)`**
- Appends `'X' * (90 - len(seq))` to reach the target length
- `X` encodes as an all-zero vector, so padding is informationally inert

**Dataset splits used in `benchmark_test.ipynb`:**

| Index | Fold | Role |
|---|---|---|
| 0–3 | 1–4 | Training (concatenated) |
| 4 | 5 | Validation — used only for early stopping |
| 5 | Benchmark | Final evaluation — never seen during training or HPO |

---

### 2. Dataset Handling — `SignalDataset`

`SignalDataset` is a custom PyTorch `Dataset` implementing a **lazy loading** strategy:

- Stores raw NumPy arrays (`self.X`, `self.y`) without converting them to tensors at initialisation
- Converts a sample to a PyTorch tensor **only when the DataLoader requests it** inside `__getitem__`
- This keeps peak RAM usage low for large one-hot encoded datasets
- `.copy()` is applied before tensor conversion to produce a writable array, preventing a PyTorch `UserWarning` that occurs when data originates from Ray's read-only shared memory
- Labels are reshaped with `.view(1)` to match the model's output shape `[B, 1]`, which is required by `BCELoss`

---

### 3. Model Architecture — `SP_NN`

The model is a three-block hybrid architecture:

```
Input  [Batch, 90, 21]
  │
  ▼
┌─────────────────────────────────────────┐
│  Conv1D Block                           │
│  kernel=17, 64 feature maps             │  ← local motif detection
│  padding='same'  →  [B, 90, 64]        │
└─────────────────────────────────────────┘
  │
  ▼
┌─────────────────────────────────────────┐
│  LSTM Block                             │
│  stacked layers, hidden=128             │  ← sequential context
│  → last hidden state  [B, 128]         │
└─────────────────────────────────────────┘
  │
  ▼
  BatchNorm1d
  │
  ▼
┌─────────────────────────────────────────┐
│  MLP Head  (dynamic depth)              │
│  Linear → ReLU → Dropout  × n layers   │  ← classification
│  → Linear → Sigmoid → probability      │
└─────────────────────────────────────────┘
```

**Block 1 — CNN (local feature extraction)**
- `nn.Conv1d` with `kernel_size=17` and `padding='same'`
- Kernel size 17 spans the typical length of a signal peptide hydrophobic core (7–15 residues), allowing the convolution to detect it in a single operation
- 64 feature maps capture diverse local sequence patterns simultaneously
- `padding='same'` preserves the sequence length (90 → 90) after convolution

**Block 2 — LSTM (sequential processing)**
- Reads the 64-channel CNN feature map as a sequence of 90 time steps
- Supports stacked layers (`num_lstm_layers`); inter-layer dropout is applied only when `num_lstm_layers > 1` (PyTorch requirement)
- Only the **last hidden state** `out[:, -1, :]` is used — this summarises the entire 90-residue sequence into a single vector after reading all positions

**Block 3 — MLP classification head (dynamic)**
- Architecture is driven by the `hidden_sizes` list — any depth and width can be tested without code changes
- Each hidden layer: `Linear → ReLU → Dropout`
- Final layer: `Linear → Sigmoid` producing a probability in (0, 1)
- `BatchNorm1d` is applied to the LSTM output before the MLP to stabilise the input distribution

**Constructor parameters:**

| Parameter | Type | Description |
|---|---|---|
| `input_size` | int | One-hot dimension — always 21 |
| `hidden_sizes` | list | Width of each MLP hidden layer, e.g. `[256, 128, 64, 1024]` |
| `lstm_hidden_size` | int | LSTM hidden-state dimension |
| `num_lstm_layers` | int | Number of stacked LSTM layers |
| `output_size` | int | 1 for binary classification |
| `dropout_p` | float | Dropout probability for MLP and inter-LSTM layers |

---

### 4. Training Loop — `train_val`

```python
train_val(model, train_loader, val_loader, optimizer, criterion,
          epochs, patience, scorer, init_best_score, output_transform, verbose)
```

The training loop incorporates three stability and generalisation mechanisms:

**Gradient clipping**
- `torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)` is applied after every backward pass
- LSTMs are prone to exploding gradients over long sequences; capping the gradient vector norm to 1.0 prevents weight updates from destabilising training

**Early stopping**
- Validation MCC is computed after every epoch
- If it fails to improve for `patience` consecutive epochs, training halts early
- This prevents overfitting and saves compute time on CPU

**Best-model tracking**
- The `state_dict` of the highest-scoring epoch is saved in RAM throughout training
- The function returns this checkpoint — not the weights from the final epoch
- This ensures the returned model is always the best-generalising version, even if later epochs overfit

**Evaluation metric — MCC (Matthews Correlation Coefficient)**
- Preferred over accuracy because the dataset is class-imbalanced
- MCC accounts for all four cells of the confusion matrix (TP, TN, FP, FN) and remains informative when one class dominates

**Returns:** `state_dict` of the best-performing epoch

---

### 5. Inference — `test`

```python
score, all_preds = test(model, test_loader, scorer, output_transform)
```

- Sets `model.eval()` to disable Dropout and freeze BatchNorm running statistics, ensuring deterministic results
- Uses `torch.no_grad()` to skip gradient computation — reduces memory usage and speeds up inference
- Default threshold: **0.5** (probabilities > 0.5 → class 1)
- **Returns two values:** the MCC scalar and the full list of binary predictions for every sample — both are needed for subsequent error analysis

---

### 6. Hyperparameter Optimisation — `lstm.ipynb`

Automated HPO is performed using **Ray Tune** with random search over 15 trials.

**Why all definitions live inside `_tune_trainable`**

Ray spawns isolated worker processes with no access to the notebook kernel. If `SP_NN`, `SignalDataset` etc. are referenced from outside `_tune_trainable`, Ray attempts `import lstm` in each worker — which fails with `ModuleNotFoundError` since no `lstm.py` exists. The fix is to redefine all required classes and functions **inside** `_tune_trainable` as closures, so Ray serialises them together with the function.

**Memory optimisation — Ray shared object store**
- The encoded dataset (~100+ MB) is stored once in Ray's shared memory via `ray.put(all_data)`
- Workers receive a lightweight reference (`data_ref`) and retrieve the data with `ray.get(config["data_ref"])`
- This prevents per-worker dataset copies that would cause out-of-memory errors
- The local `all_data` variable is deleted and `gc.collect()` is called immediately after `ray.put()` to free RAM

**Search space:**

| Parameter | Distribution |
|---|---|
| `num_layers` | choice [2, 3, 4, 5] |
| `hidden_sizes` | sampled via `_generate_hidden_sizes` |
| `dropout` | uniform [0.1, 0.5] |
| `lr` | log-uniform [1e-4, 1e-3] |
| `batch_size` | choice [10, 20] |
| `num_lstm_layers` | choice [1, 2] |
| `lstm_hidden_size` | choice [64, 128] |

**MLP topology sampling — `_generate_hidden_sizes`**

Two design philosophies are explored at random (50/50):
- **Funnel:** layer sizes in descending order (e.g. `512 → 256 → 64`) — compresses features progressively, stable
- **Flexible:** layer sizes shuffled randomly (e.g. `128 → 512 → 64`) — allows intermediate expansion, more expressive

**Reporting:** Ray minimises `loss` by default, so `-mean_mcc` is reported as `loss` to maximise MCC.

---

### 7. Benchmark Evaluation & Error Analysis — `benchmark_test.ipynb`

After loading the JIT-traced model, every benchmark entry is annotated with its prediction outcome:

| Label | Meaning |
|---|---|
| **TP** | True Positive — signal peptide correctly identified |
| **TN** | True Negative — non-signal correctly classified |
| **FP** | False Positive — non-signal wrongly predicted as signal |
| **FN** | False Negative — real signal peptide missed by the model |

**Analyses performed and figures saved to `model_evaluation/`:**

| Figure | Description |
|---|---|
| `Confusion_Matrix.png` | Heatmap of TP/TN/FP/FN counts |
| `HD_prediction_1.png` | HelixDomain protein rate per prediction category |
| `HD_prediction_2.png` | HelixDomain rate: FP-only vs all Negatives |
| `SP_length_distribution.png` | Signal peptide length histogram — TP vs FN |
| `Hydrophobicity_SP.png` | Miyazawa hydrophobicity histogram — TP vs FN |
| `Boxplot_Hydrophobicity_SP.png` | Hydrophobicity boxplot — TP vs FN |
| `Mean_hydrophobicity.png` | Mean hydrophobicity bar chart per category |
| `AA_frequencies_TP_vs_FN.png` | Full-sequence AA frequency — TP vs FN |
| `AA_frequencies_SP.png` | SP-region AA frequency — TP vs FN |
| `AA_frequency_Distribution.png` | AA frequency — FP / FN / all Positives / all Negatives |
| `SeqLogo_IC_1–4.png` | Sequence logos (information content) per category |
| `SeqLogo_Prob_1–4.png` | Sequence logos (probability) per category |
| `FPR_Kingdom.png` | False positive rate pie chart per taxonomic kingdom |
| `Pieplot_species.png` | Taxonomic composition: all Positives / FN / FP |

---

### 8. Model Export — TorchScript JIT Tracing

The trained model is exported using `torch.jit.trace`:

```python
model.eval()
dummy_input  = torch.randn(1, 90, 21).to(device)
traced_model = torch.jit.trace(model, dummy_input)
torch.jit.save(traced_model, "SignalPeptideLSTM.pt")
```

**Why JIT tracing over `state_dict`:**
- The `.pt` file contains **both architecture and weights** — no source code is needed to load it
- Compatible with non-Python environments (C++, Java, mobile)
- `model.eval()` before tracing permanently disables Dropout in the saved model, ensuring deterministic inference

**Loading the saved model:**
```python
model = torch.jit.load("SignalPeptideLSTM.pt")
model.eval()
```

---

## Final Results

### Best Hyperparameter Configuration

Identified via Ray Tune random search (15 trials), maximising average 5-fold MCC:

| Parameter | Value |
|---|---|
| LSTM layers | 2 |
| LSTM hidden size | 128 |
| MLP hidden layers | `[256, 128, 64, 1024]` |
| Dropout | 0.498 |
| Learning rate | 2.86 × 10⁻⁴ |
| Batch size | 20 |

**MLP topology note:** The selected `[256 → 128 → 64 → 1024]` architecture compresses features to a 64-dimensional bottleneck before expanding to 1024. This hourglass pattern suggests that projecting compressed representations into a higher-dimensional space improves class separation at the final linear layer.

### Benchmark Performance

| Metric | Value |
|---|---|
| **MCC** | **0.9019** |
| Loss function | BCELoss |
| Optimiser | Adam |
| Threshold | 0.5 |


