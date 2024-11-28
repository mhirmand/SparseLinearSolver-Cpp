#include "SparseLinearSolver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace SparseSolver {

  // constructor for SparseMatrix class
  SparseMatrix::SparseMatrix(int rows, int cols, int nnz)
    : numRows(rows), numCols(cols), numNonZero(nnz), finalized(false), currentNNZ(0) {
    values = nullptr;
    colIndex = nullptr;
    rowIndex = nullptr;
    if (numNonZero > 0)
    {
      values = new double[numNonZero] {0.0};
      colIndex = new int[numNonZero];
    }
    if (numRows > 0)
      rowIndex = new int[numRows + 1] {0};
    
  }

  // deconstructor for SparseMatrix class
  SparseMatrix::~SparseMatrix() {
    delete[] values;
    values = nullptr;

    delete[] rowIndex;
    rowIndex = nullptr;

    delete[] colIndex;
    colIndex = nullptr;
  }

  // method to add value to the matrix
  void SparseMatrix::addValue(int* rows, int* cols, double* vals, int nVals) {

    // check if the matrix is finalized
    if (finalized) {
      throw std::logic_error("Cannot add values after finalizing the matrix.");
    }

    // check if too many non-zero elements are added
    if (currentNNZ >= numNonZero) {
      throw std::overflow_error("Exceeded the allocated non-zero elements.");
    }

    // add the value to the matrix
    int r = 0;
    int c = 0;
    double val = 0.0;
    for (int i = 0; i < nVals; i++) {
      r = rows[i];
      c = cols[i];
      val = vals[i];
      if (r >= numRows || c >= numCols) {
        throw std::out_of_range("Row or column index out of range.");
      }
      values[currentNNZ] = val;
      colIndex[currentNNZ] = c;
      rowIndex[r + 1]++;
      currentNNZ++;
    }
  }

  // method to finalize the matrix
  // it transitions the matrix from its user - friendly input format to 
  // the compact Compressed Row Storage(CRS) format.
  void SparseMatrix::finalize() {
    if (finalized) return;

    for (int i = 1; i <= numRows; ++i) {
      rowIndex[i] += rowIndex[i - 1];
    }
    finalized = true;
  }

  // solve the linear system Ax = b
  void Solver::solve(SparseMatrix& A, double *& b) {
    
    // Step 1: Perform forward elimination
    forwardElimination(A, b);

    // Step 2: Perform backward substitution
    backwardSubstitution(A, b);

  }

  // backward substitution
  void Solver::forwardElimination(const SparseMatrix& A, double*& b) {

    int n = A.numRows;

    // Forward elimination
    for (int k = 0; k < n; ++k) {
      // Find diagonal element in row k
      double diag = 0.0;
      int diagIndex = 0;
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
      b[k] /= diag;

      // Eliminate below
      for (int i = k + 1; i < n; ++i) {
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
          b[i] -= factor * b[k];
        }
      }
    }
    return;
  }

  // forward substitution
  void Solver::backwardSubstitution(const SparseMatrix& A, double*& x) {

    int n = A.numRows;

    // Backward substitution
    for (int i = n - 1; i >= 0; --i) {
      for (int idx = A.rowIndex[i]; idx < A.rowIndex[i + 1]; ++idx) {
        int col = A.colIndex[idx];
        if (col > i) {
          x[i] -= A.values[idx] * x[col];
        }
      }
    }
    return;
  }
}
