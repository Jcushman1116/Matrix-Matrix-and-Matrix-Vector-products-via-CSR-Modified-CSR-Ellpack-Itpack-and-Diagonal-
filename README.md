# Sparse Matrix–Vector Multiplication Using Structured Storage Formats

## Overview
This project implements and evaluates multiple algorithms for **sparse matrix–vector multiplication** using different compressed storage formats. The goal is to reduce unnecessary computation and memory usage by exploiting sparsity and structure in matrices.

All algorithms are implemented in **MATLAB** and tested against MATLAB’s built-in routines to evaluate **accuracy**, **runtime performance**, and **scalability** as the matrix dimension `n` increases.

---

## Tasks and Algorithms Implemented

### Task 1: In-Row-Order COO Matrix–Vector Products

#### (a) COO → CSR Conversion and Multiplication
Converts a sparse matrix from **Coordinate (COO)** format to **Compressed Sparse Row (CSR)** format and computes:

    z = A v

- Builds row pointer array `IA`
- Uses dot products per row
- Efficient for row-based access
- Time complexity: **O(nnz)**

---

#### (b) COO → Modified CSR Conversion and Multiplication
Stores diagonal and off-diagonal elements separately to reduce storage overhead.

- Diagonal elements stored explicitly
- Off-diagonal elements stored row-wise
- Optimized storage for sparse matrices
- Demonstrates space–time tradeoffs

---

#### (c) COO → Ellpack–Itpack Conversion and Multiplication
Converts COO format into **Ellpack–Itpack** representation.

- Uses fixed-width row storage
- Pads rows with zeros as needed
- Efficient for vectorized computation
- Time complexity depends on maximum nonzeros per row

---

### Task 2: Unordered COO → Relaxed CSR
Handles COO matrices where row indices are **not sorted**.

- Sorts row indices
- Builds CSR structure with **elbow room** for insertion
- Computes matrix–vector product using relaxed CSR
- Demonstrates dynamic sparse matrix storage

---

### Task 3: Symmetric CSR Storage
Implements CSR storage for **symmetric matrices** by storing **only the lower triangular part**.

- Reduces memory requirements
- Exploits symmetry during multiplication
- Automatically adds mirrored contributions
- Accurate for large matrix sizes

---

### Task 4: Compressed Diagonal Storage
Handles matrices with **three nonzero diagonals** (main, sub, super).

- Automatically detects diagonal offset `k`
- Stores diagonals in a compact array
- Uses offset indexing for multiplication

**Note:**  
- Works correctly for `k = 1`
- Produces incorrect results for `k > 1` due to an indexing issue

---

## Experimental Design

- Matrix and vector entries generated synthetically
- Dimensions tested up to `n = 1000`
- Performance metrics:
  - **Runtime**
  - **Relative error** compared to MATLAB built-in routines

### Relative Error Definition

    ||x_comp − x_true|| / ||x_true||

---

## Results and Observations

- All Task 1 algorithms produce correct results with low relative error
- MATLAB’s built-in routines outperform custom implementations in runtime
- Error increases with problem size due to floating-point accumulation
- Symmetric CSR achieves correct results with reduced storage
- Compressed diagonal method fails for offsets greater than 1

---

## Files Included

### Task 1
- `COO_to_CSR_Mv_mult.m`
- `COO_to_modified_CSR_Mv_mult.m`
- `COO_to_ELL_Mv_mult.m`
- `task1_driver.m`

### Task 2
- `COO_to_relaxed_CSR_Mv_mult.m`
- `task2_driver.m`

### Task 3
- `Task3.m`
- `task3_driver.m`

### Task 4
- `task4.m`
- `task4_driver.m`

Each task is implemented in a separate MATLAB file with an accompanying testing script.

---

## Key Takeaways

- Sparse storage formats significantly reduce unnecessary computation
- Algorithmic efficiency depends heavily on matrix structure
- Built-in MATLAB routines remain highly optimized
- Correct indexing is critical in compressed storage algorithms
