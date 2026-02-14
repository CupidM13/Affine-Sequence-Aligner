# Affine Sequence Aligner

A comprehensive bioinformatics tool written in Python that performs pairwise nucleotide sequence alignment using **Affine Gap Penalties**. 

This tool calculates both **Global** (Gotoh/Needleman-Wunsch) and **Local** (Smith-Waterman) alignments simultaneously, providing a detailed report on the biological relationship between two sequences.

## 🧬 Key Features

* **Affine Gap Scoring:** Implements biologically accurate gap penalties distinguishing between gap opening ($\sigma$) and gap extension ($\epsilon$).
* **Dual Mode:** Automatically computes and reports both Global and Local alignments in a single run.
* **Contextual Visualization:**
    * Standard alignment blocks (80 characters width).
    * "Global Context" view for local alignments (shows where the match occurs within the full sequence).
* **Robust Input Handling:** Supports standard FASTA files (including multiline sequences).
* **Automated Reporting:** Saves formatted results to a timestamped file in the `/outputs` directory.

## ⚙️ Algorithms & Logic

Unlike standard linear gap penalties ($Score = k \times gap$), this tool uses the **Affine Gap Penalty Model**, which is more biologically realistic as insertions/deletions often occur in contiguous chunks.

### Formula
$$Penalty(k) = \sigma + \epsilon \cdot (k-1)$$

Where:
* $\sigma$ (**Sigma**): Gap Opening Penalty (High cost to start a gap).
* $\epsilon$ (**Epsilon**): Gap Extension Penalty (Lower cost to extend a gap).
* $k$: Length of the gap.

### Matrices
The tool utilizes **Gotoh's Algorithm** with three Dynamic Programming matrices:
1.  **Middle ($M$):** Best score ending in a Match or Mismatch (Diagonal).
2.  **Lower ($L$):** Best score ending in a gap in Sequence 2 (Vertical).
3.  **Upper ($U$):** Best score ending in a gap in Sequence 1 (Horizontal).

## 📂 Project Structure

```text
Affine-Sequence-Aligner/
├── data/               # Place your input FASTA files here
│   ├── seq1.fasta
│   └── seq2.fasta
├── outputs/            # Results are automatically saved here
├── src/                # Source code
│   └── affine_aligner.py
├── README.md           # This documentation
└── requirements.txt    # Dependencies
```

## 🚀 Installation

1. **Clone the repository**:

2. **Install Dependencies:**
This tool requires `numpy` for efficient matrix operations.
```bash
pip install numpy
```



## 💻 Usage

**Prerequisite:** Ensure your FASTA files are saved somewhere inside your project folder (e.g., inside the `data/` folder).

Run the tool from your project root using the command line:

```bash
python src/affine_aligner.py <file1_path> <file2_path>

```
Example:
```bash
python src/affine_aligner.py data/seq1.fasta data/seq2.fasta

```

### ⚠️ CRITICAL: SCORING PARAMETERS ⚠️

When prompted for scoring parameters, you **MUST enter negative numbers** for penalties. The algorithm treats positive numbers as rewards.

* **Match:** Positive (e.g., `1.0` or `5.0`)
* **Mismatch:** **NEGATIVE** (e.g., `-1.0` or `-4.0`)
* **Sigma (Gap Open):** **NEGATIVE** (e.g., `-5.0` or `-10.0`)
* **Epsilon (Gap Extend):** **NEGATIVE** (e.g., `-1.0` or `-0.5`)

*> Failure to enter negative values will result in incorrect alignments.*

## 📄 Output Format

Results are saved to the `outputs/` folder with a filename format:
`{seq1}_{seq2}_{DDHHMMSS}.txt`

The report includes:

1. **Global Alignment:** The best end-to-end alignment.
2. **Local Alignment:** The best matching substring.
3. **Local Context View:** The local alignment shown relative to the full original sequences (useful for seeing where a domain lies within a gene).

## 🧪 Example

**Input:**

* Match: `1`
* Mismatch: `-1`
* Sigma: `-5`
* Epsilon: `-1`

**Console Output:**

```text
Calculating Global Alignment...
Calculating Local Alignment...

============================================================
SUCCESS!
Global Score: -3.0
Local Score:  12.0
Report saved: outputs/seqA_seqB_141230.txt
============================================================

```

## 📝 Author

* **Cupid Moolchandani**
* Bioinformatics Coding Challenge Submission

```

```
