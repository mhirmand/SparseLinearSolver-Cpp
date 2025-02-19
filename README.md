# C++ SparseLinearSolver

This C++ project demonstrates the solution of sparse linear systems using Gaussian Elimination. The implementation uses the **Compressed Row Storage (CRS)** format for efficient storage and manipulation of sparse matrices. This code is designed for **demonstration purposes** and may not be optimal for solving large systems.

---

## Table of Contents
1. [Purpose](#purpose)
2. [Usage](#usage)
   - [Example](#example)
3. [Implementation Details](#implementation-details)
   - [SparseMatrix Class](#sparsematrix-class)
   - [Solver Class](#solver-class)
4. [Time and Memory Complexity](#time-and-memory-complexity)
5. [Limitations](#limitations)
6. [License](#license)

---

## Purpose

The purpose of this code is to demonstrate how to solve sparse linear systems of equations using Gaussian Elimination. The implementation is designed to handle various edge cases, such as:
- Singular matrices
- Zero matrices
- Rectangular matrices

The code serves best as a learning tool for understanding sparse matrix storage and linear system solvers.

---

## Usage

To use the solver, include the `SparseLinearSolver.h` header in your project and link against the `SparseLinearSolver.cpp` implementation. The main functionality is encapsulated in the `SparseSolver` namespace, with the `SparseMatrix` class representing the sparse matrix and the `Solver` class providing the solving capabilities.

### Example

```cpp
#include "SparseLinearSolver.h"

int main() {
  int nError = 0;

  const int n = 3;   // 3x3 system
  const int nnz = 7; // Number of non-zero elements

  SparseSolver::SparseMatrix A(n, n, nnz);
  double* b = new double[n] {0.0};
  b[0] = 8.0;
  b[1] = -11.0;
  b[2] = -3.0;

	// define the matrix entries 
  int rows[] = { 0, 0, 0, 1, 1, 1, 2 };
  int cols[] = { 0, 1, 2, 0, 1, 2, 2 };
  double vals[] = { 2.0, 1.0, -1.0, -3.0, -1.0, 2.0, 2.0 };

	// Add non-zero values (we could do this element by element instead of giving a list of values)
  A.addValue(rows, cols, vals, nnz);

	// Finalize the matrix
	A.finalize();

  try {
    // Solve the system
    SparseSolver::Solver::solve(A, b);
	
	std::cout << "Solution: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;

    // Check the solution
    assert(std::fabs(b[0] - 1.5) < epsilon);
    assert(std::fabs(b[1] - 3.5) < epsilon);
    assert(std::fabs(b[2] + 1.5) < epsilon);
    std::cout << "Test Simple System: Passed" << std::endl;
  }
  catch (const std::runtime_error& e) {
    std::cout << "Test Simple Matrix: Failed (" << e.what() << ")" << std::endl;
    nError = 1;
  }

	// clean memory
	delete[] b;
	b = nullptr;

  return nError;    

    return 0;
}

## Implementation Details

### SparseMatrix Class

The `SparseMatrix` class represents a sparse matrix in **Compressed Row Storage (CRS)** format. It provides the following functionality:

- **Constructor**: Initializes the matrix with the given number of rows, columns, and non-zero elements.
- **addValue**: Adds non-zero values to the matrix.
- **finalize**: Converts the matrix to CRS format for efficient storage and computation.

### Key Members:
- `values`: Array of non-zero values.
- `rowIndex`: Array of row pointers.
- `colIndex`: Array of column indices.
- `numRows`, `numCols`, `numNonZero`: Dimensions and number of non-zero elements.

---

# Solver Class

The `Solver` class provides the functionality to solve the linear system \( Ax = b \) using Gaussian Elimination. It consists of two main steps:

1. **Forward Elimination**: Transforms the matrix into an upper triangular form.
2. **Backward Substitution**: Solves for the unknowns using back substitution.

### Key Methods:
- `solve`: Solves the system \( Ax = b \).
- `forwardElimination`: Performs forward elimination to reduce the matrix to upper triangular form.
- `backwardSubstitution`: Solves the system using back substitution.

---

# Time and Memory Complexity

## Time Complexity
- **Gaussian Elimination**: The algorithm has a time complexity of \( O(n^3) \) for dense matrices. For sparse matrices, the complexity depends on the sparsity pattern but is generally lower.
- **Forward Elimination**: \( O(n^3) \) in the worst case.
- **Backward Substitution**: \( O(n^2) \).

## Memory Complexity
- The **Compressed Row Storage (CRS)** format ensures efficient memory usage:
  - \( O(n + nnz) \) space complexity, where \( n \) is the number of rows and \( nnz \) is the number of non-zero elements.
