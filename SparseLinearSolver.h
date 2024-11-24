#pragma once

#ifndef SPARSE_LINEAR_SOLVER_H
#define SPARSE_LINEAR_SOLVER_H

#include <stddef.h>

namespace SparseSolver {

  class SparseMatrix {
  public:
    double* values;    // Non-zero values
    int* rowIndex;     // Row pointers
    int* colIndex;     // Column indices
    size_t numRows;    // Number of rows
    size_t numCols;    // Number of columns
    size_t numNonZero; // Number of non-zero values

    SparseMatrix(size_t rows, size_t cols, size_t nnz);
    ~SparseMatrix();

    void addValue(size_t row, size_t col, double value);
    void finalize();

  private:

    bool finalized;
    size_t currentNNZ;
  };

  class Solver {
  public:
    // Solver() : tolZero(1e-10) {};
    static void solve(const SparseMatrix& A, const double* b, double* x);
  private:
    // double tolZero;
    // void setTol(const double tol) { tolZero = tol; };
    static void forwardSubstitution(const SparseMatrix& A, double*& b);
    static void backwardSubstitution(const SparseMatrix& A, const double* b, double*& x);
  };

}

#endif // SPARSE_LINEAR_SOLVER_H
