#include <iostream>
#include "SparseLinearSolver.h"

using namespace SparseSolver;

int main() {
  const size_t n = 3; // 3x3 system
  const size_t nnz = 7; // Number of non-zero elements

  SparseMatrix A(n, n, nnz);
  double b[] = { 8, -11, -3 };
  double x[n] = { 0 };

  size_t rows[] = { 0, 0, 0, 1, 1, 1, 2 };
  size_t cols[] = { 0, 1, 2, 0, 1, 2, 2 };
  double vals[] = { 2.0, 1.0, -1.0, -3.0, -1.0, 2.0, 2.0 };

  // Add non-zero values
  A.addValue(rows, cols, vals, 7);

  // Finalize the matrix
  A.finalize();

  tolZero = 1e-8;

  try {

    Solver::solve(A, b);

    std::cout << "Solution:" << std::endl;
    for (size_t i = 0; i < n; ++i) {
      std::cout << x[i] << " ";
    }
    std::cout << std::endl;
  }
  catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
