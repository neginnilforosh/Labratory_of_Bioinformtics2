# Deep Learning : LSTM model

This repository implements a hybrid **CNN-LSTM neural network** for our signal peptide prediction problem.
The code is organized into two main scripts:

1. **`benchmark_test.ipyb`** – Handles **data preprocessing**, fold splitting, sequence padding/truncation, and one-hot encoding of protein sequences.
   This script prepares all datasets required for training, validation, and benchmarking.

3. **`lstm.ipyb`** – Contains the **core deep learning model** (`SP_NN`), dataset handling (`SignalDataset`),
   training and validation loop (`train_val`), testing (`test`), and model export to TorchScript. This script performs the model learning, evaluation, and saving for deployment.

---

## Pipeline Overview

### 1️. Data Preprocessing
The first step is the **cleaning, standardization, and encoding** of sequences. Our input file is the `train_bench.tsv` (protein sequences, labels, fold IDs)
- **Sequence length standardization**:
  - Sequences shorter than 90 residues are padded with `'X'`  
  - Sequences longer than 90 residues are truncated (N-terminal preserved)
- **One-hot encoding**:
  - Converts sequences to 90×21 matrices (20 amino acids + padding)
- Splits data into folds:
  - Training folds: 1–4  
  - Validation fold: 5  
  - Benchmark fold: 6

**Key functions**:

- `create_one_hot_sets(dataset)`: main preprocessing and encoding pipeline  
- `_increase_lenseq(seq)`: pads sequences to length 90  
- `one_hot_encoding(seq)`: converts amino acid sequences to one-hot matrices  

---

### 2️. Dataset Handling 

For Dataset Handling the `SignalDataset` class:
  - Stores NumPy arrays of features and labels
  - Converts samples to PyTorch tensors **only when requested** (lazy loading)
  - Compatible with `DataLoader` for batching

---

### 3️. Model Definition 
Our model is a Hybrid CNN_LSTM model.

**Architecture**:

1. **CNN Block**  
   - Input: 21 channels (one-hot)  
   - 1D convolution, 64 feature maps, kernel size 17  

2. **LSTM Block**  
   - Processes CNN features sequentially  
   - Supports stacked layers  
   - Output summarized to last time step  

3. **MLP Head**  
   - Dynamic fully connected layers (configurable via `hidden_sizes`)  
   - ReLU + Dropout for regularization  
   - Sigmoid for binary output

---

### 4️. Training Loop

For the Training of our model, we used the function `train_val(model, train_loader, val_loader, ...)`
  - Trains model using BCELoss and Adam optimizer  
  - **Early stopping** based on MCC  
  - **Gradient clipping** to stabilize LSTM training  
  - Returns **best model weights** (state dict)

The validation performed every epoch to monitor MCC: training is stopped if no improvement for `patience` epochs

---

### 5️. Testing 

We created the function `test(model, test_loader, ...)`  
- Evaluates the model on the **benchmark fold**  
- Outputs:
  - MCC (Matthews Correlation Coefficient)  
  - Binary predictions (`0/1`)  
The function uses a threshold of 0.5 for classifying probabilities

---
 

## Result Final And Model Export

We use d**TorchScript JIT tracing** to create a standalone `.pt` model that Includes both **architecture and weights**, can be loaded in Python or other environments without the original code
The final optimized configuration (from 5-fold cross-validation) includes:

- LSTM: 2 layers, 128 hidden units each  
- MLP hidden layers: `[256, 128, 64, 1024]`  
- Dropout: 0.5  
- Learning rate: 2.8e-4  
- Batch size: 20  

This architecture compresses and expands features before classification, improving separation of classes, with a final MCC on benchmark set: 0.9018711109052984




