✅ README Template: Program 2 — Sparse Matrix Data Structures & MV Products
📌 Overview

This project implements multiple sparse matrix data structures and uses them to compute matrix–vector products efficiently.
The assignment covers COO, CSR, Modified CSR, Ellpack–Itpack, Symmetric CSR, and Compressed Diagonal formats.

🧠 Key Concepts

Sparse matrix compression

Row pointers (IA), column indices (JA), and compressed arrays (AA)

Efficient matrix–vector multiplication

Handling unordered COO input

Exploiting matrix symmetry and bandedness

📦 Implemented Data Structures & Algorithms
✔️ 1. COO → CSR → Matrix–Vector Product

Standard CSR format:

AA = nonzero values

JA = column indices

IA = row pointer

✔️ 2. Modified CSR

Stores:

Diagonal elements separately

Off-diagonal values and pointers

Reduces branching during multiplication

✔️ 3. COO → Ellpack–Itpack

Fixed number of columns per row

Useful for vectorized computation

✔️ 4. Unordered COO → Relaxed CSR

Handles input where row order is arbitrary

Adds “elbow room” for insertions

✔️ 5. Symmetric CSR

Stores only lower triangular part

Automatically adds symmetric contributions during MV product

✔️ 6. Compressed Diagonal Storage (CDS)

Extracts the main diagonal and offset diagonals

For banded matrices with bandwidth 
𝑘
k

🧪 Experimental Design

Random sparse matrices tested across many dimensions

Accuracy verified by comparing with dense multiplication

Performance improvement measured qualitatively

📈 Results Summary

All formats produce correct products

Symmetric CSR and CDS significantly reduce storage

Ellpack performs consistently for evenly distributed sparsity
