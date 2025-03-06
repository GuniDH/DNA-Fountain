# DNA Fountain Encoding & Decoding by Guni

## Overview

This project implements **DNA Fountain Encoding and Decoding**, a method for storing digital data in DNA sequences. The system **transforms binary data into DNA oligomers** (short sequences of nucleotides) and enables reconstruction using a randomized encoding scheme.

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

## How DNA Encoding Works

1. **Binary Data Preparation**
   - The input data is split into **4-bit segments**.
2. **Tip Generation (Oligomers)**
   - Each segment is randomly assigned to one or more DNA tips.
   - XOR-based redundancy ensures error correction.
3. **Base Mapping**
   - Binary pairs (00, 01, 10, 11) are mapped to **nucleotide bases (A, T, C, G)**.
4. **Final DNA Encoding**
   - The tips are combined to form DNA sequences ready for storage.

## How DNA Decoding Works

1. **Read the stored DNA tips.**
2. **Reverse the base mapping** to retrieve binary values.
3. **Use the original seed-based selection** to reconstruct missing data.
4. **Recover the original binary sequence**.

## Challenges

- **Expensive DNA synthesis & sequencing.**
- **Slow read/write speeds compared to traditional storage.**

---
### Author
**Guni**  


