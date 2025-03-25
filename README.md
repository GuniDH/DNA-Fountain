# DNA Fountain Encoding & Decoding by Guni

## Overview

This project implements **DNA Fountain Encoding and Decoding**, a method for storing digital data in DNA sequences. The system **transforms binary data into DNA oligomers** (short sequences of nucleotides) and enables reconstruction using a randomized encoding scheme.
This is an improved version of the normal DNA fountain because it has memory extension with special bases M,K.
DNA Fountain is an erasure-correcting code designed specifically for DNA storage systems. It addresses challenges like sequence dropout and synthesis/sequencing errors by incorporating:

1. A Luby Transform code-like approach with a robust soliton distribution
2. Efficient binary-to-DNA encoding and DNA-to-binary decoding
3. XOR-based payload generation to improve error resilience

## How It Works

1. **Encoding Process**:
   - Binary input is segmented into fixed-size chunks
   - For each "tip" (oligomer):
     - A random subset of segments is selected based on the degree distribution
     - The segments are XORed together to create a payload
     - A unique barcode, seed, and degree value are added to the payload
     - The binary sequence is converted to DNA bases
   - Multiple copies of each tip are created to ensure redundancy

2. **Decoding Process**:
   - DNA sequences are converted back to binary
   - Metadata (barcode, seed, degree) is extracted from each tip
   - Tips with degree 1 (containing exactly one segment) are identified
   - As segments are recovered, they're removed from other tips through XOR operations
   - This process continues until all original segments are recovered
   - The original binary input is reassembled

## Implementation Details

The implementation includes:

- A `dna_fountain` class that handles both encoding and decoding operations
- Customizable parameters for barcode size, seed values, segment sizes, etc.
- DNA base mapping system (using A, C, G, T and special markers M, K)
- Robust Soliton distribution for optimizing tip (oligomer) degrees

DNA-based storage offers:
- **High-density data storage**: Can store vast amounts of data in a tiny physical space.
- **Long-term stability**: DNA is more durable than traditional storage media.
- **Energy efficiency**: Requires minimal power compared to data centers.

## Features

- **DNA Fountain Encoding**
  - Converts binary data into short **DNA oligomers** using **randomized tip selection**.
  - Uses **XOR-based redundancy** to enhance data recovery.
- **Lossless Data Recovery**
  - Ensures **reconstruction of original binary data** from DNA sequences.
- **Randomized Encoding for Security**
  - Uses **seed-based segment selection** for randomized encoding.

## Challenges

- **Expensive DNA synthesis & sequencing.**
- **Slow read/write speeds compared to traditional storage.**

## Usage

Run the script with:

```python
python main.py
```

The script includes a sample binary input string for demonstration, but can be modified to work with any binary input.

## Author

Guni


