#include "SparseLinearSolver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace SparseSolver {

  // constructor for SparseMatrix class
  SparseMatrix::SparseMatrix(size_t rows, size_t cols, size_t nnz)
    : numRows(rows), numCols(cols), numNonZero(nnz), finalized(false), currentNNZ(0) {
    values = new double[numNonZero];
    rowIndex = new int[numRows + 1] {0};
    colIndex = new int[numNonZero];
  }

  // deconstructor for SparseMatrix class
  SparseMatrix::~SparseMatrix() {
    delete[] values;
    delete[] rowIndex;
    delete[] colIndex;
  }

  // method to add value to the matrix
  void SparseMatrix::addValue(size_t row, size_t col, double value) {

    // check if the matrix is finalized
    if (finalized) {
      throw std::logic_error("Cannot add values after finalizing the matrix.");
    }

    // check if too many non-zero elements are added
    if (currentNNZ >= numNonZero) {
      throw std::overflow_error("Exceeded the allocated non-zero elements.");
    }

    // add the value to the matrix
    values[currentNNZ] = value;
    colIndex[currentNNZ] = col;
    rowIndex[row + 1]++;
    currentNNZ++;
  }

  // method to finalize the matrix
  // it transitions the matrix from its user - friendly input format to 
  // the compact Compressed Row Storage(CRS) format.
  void SparseMatrix::finalize() {
    if (finalized) return;

    for (size_t i = 1; i <= numRows; ++i) {
      rowIndex[i] += rowIndex[i - 1];
    }
    finalized = true;
  }

  // solve the linear system Ax = b
  void Solver::solve(const SparseMatrix& A, const double* b, double* x) {
    // Step 1: Create a mutable copy of `b`
    double* bcopy = new double[A.numRows];
    memcpy(bcopy, b, sizeof(*(A.values)) * A.numRows);

    size_t n = A.numRows;

    // Forward elimination
    for (size_t k = 0; k < n; ++k) {
      // Find diagonal element in row k
      double diag = 0.0;
      size_t diagIndex = 0;
      for (int idx = A.rowIndex[k]; idx < A.rowIndex[k + 1]; ++idx) {
        if (A.colIndex[idx] == k) {
          diag = A.values[idx];
          diagIndex = idx;
          break;
        }
      }
      if (std::fabs(diag) < 1e-10) {
        throw std::runtime_error("Matrix is singular or nearly singular.");
      }

      // Normalize the k-th row
      for (int idx = A.rowIndex[k]; idx < A.rowIndex[k + 1]; ++idx) {
        A.values[idx] /= diag;
      }
      bcopy[k] /= diag;

      // Eliminate below
      for (size_t i = k + 1; i < n; ++i) {
        double factor = 0.0;
        for (int idx = A.rowIndex[i]; idx < A.rowIndex[i + 1]; ++idx) {
          if (A.colIndex[idx] == k) {
            factor = A.values[idx];
            break;
          }
        }

        if (std::fabs(factor) > 1e-10) {
          for (int idx = A.rowIndex[i]; idx < A.rowIndex[i + 1]; ++idx) {
            int col = A.colIndex[idx];
            A.values[idx] -= factor * A.values[A.rowIndex[k] + (col - k)];
          }
          bcopy[i] -= factor * bcopy[k];
        }
      }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
      x[i] = bcopy[i];
      for (int idx = A.rowIndex[i]; idx < A.rowIndex[i + 1]; ++idx) {
        int col = A.colIndex[idx];
        if (col > i) {
          x[i] -= A.values[idx] * x[col];
        }
      }
    }

    // Step 4: Clean up
    delete[] bcopy;
  }

  // backward substitution
  void Solver::forwardSubstitution(const SparseMatrix& A, double*& b) {

    return;
  }


  // forward substitution
  void Solver::backwardSubstitution(const SparseMatrix& A, const double* b, double* &x) {


    return;
	}
}
