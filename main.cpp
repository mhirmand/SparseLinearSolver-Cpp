#include <iostream>
#include "SparseLinearSolver.h"

int main() {
  const size_t n = 3; // 3x3 system
  const size_t nnz = 7; // Number of non-zero elements

  SparseSolver::SparseMatrix A(n, n, nnz);
  double b[] = { 8, -11, -3 };
  double x[n] = { 0 };

  // Add non-zero values
  A.addValue(0, 0, 2.0);
  A.addValue(0, 1, 1.0);
  A.addValue(0, 2, -1.0);
  A.addValue(1, 0, -3.0);
  A.addValue(1, 1, -1.0);
  A.addValue(1, 2, 2.0);
  // A.addValue(2, 1, 1.0);
  A.addValue(2, 2, 2.0);

  // Finalize the matrix
  A.finalize();

  try {
    SparseSolver::Solver::solve(A, b, x);

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
