#include <iostream>
#include <cassert>
#include <cmath>
#include "SparseLinearSolver.h"

#define EPSILON 1e-10

using namespace SparseSolver;

int testSimpleSystem() {

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

    // Check the solution
    assert(std::fabs(b[0] - 1.5) < EPSILON);
    assert(std::fabs(b[1] - 3.5) < EPSILON);
    assert(std::fabs(b[2] + 1.5) < EPSILON);
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
}

int testSingularMatrix() {
  int nError = 0;

  const int n = 2;
  const int nnz = 2;

  SparseSolver::SparseMatrix A(n, n, nnz);
  double* b = new double[2] {0.0};
  b[0] = 1.0;
  b[1] = 2.0;
  int rows[2] = { 0, 0 };
  int cols[2] = { 0, 1 };
  double vals[2] = { 1.0, 2.0 };

  // Define a singular matrix (second row is all zero)
  A.addValue(rows, cols, vals, nnz);

  // FInalize the matrix
  A.finalize();

  try {
    SparseSolver::Solver::solve(A, b);
    assert(false); // Should not reach here
    std::cout << "Test Singular System: Failed" << std::endl;
    nError = 1;
  }
  catch (const std::runtime_error& e) {
    std::cout << "Test Singular Matrix: Passed (" << e.what() << ")" << std::endl;
  }

  delete[] b; 
  b = nullptr;


  return nError;
}

int testZeroMatrix() {

  int nError = 0;
  const int n = 3;
  const int nnz = 0;

  SparseSolver::SparseMatrix A(n, n, nnz);

  double* b = new double[3] {0.0};
  b[0] = 1.0;
  b[1] = 2.0;
  b[2] = 3.0;

  A.finalize();

  try {
    SparseSolver::Solver::solve(A, b);
    assert(false); // Should not reach here
    std::cout << "Test Zero System: Failed" << std::endl;
    nError = 1;
  }
  catch (const std::runtime_error& e) {
    std::cout << "Test Zero Matrix: Passed (" << e.what() << ")" << std::endl;
  }

  // clean memory
  delete[] b;
  b = nullptr;

  return nError;
}

int testRectangularMatrix() {

  int nError = 0;

  const int nrows = 3;
  const int ncols = 2; // Rectangular matrix (non-square)
  const int nnz = 4;

  SparseSolver::SparseMatrix A(nrows, ncols, nnz);
  double* b = new double[3] {0.0};
  b[0] = 1.0;
  b[1] = 2.0;
  b[2] = 3.0;

  int rows[4] = { 0, 1, 2, 2 };
  int cols[4] = { 0, 1, 0, 1 };
  double vals[4] = { 1.0, 1.0, 1.0, 1.0 };

  // Add values to the matrix
  A.addValue(rows, cols, vals, nnz);

  A.finalize();

  try {
    SparseSolver::Solver::solve(A, b);
    assert(false); // Should not reach here
    std::cout << "Test Rectangular System: Failed" << std::endl;
    nError = 1;
  }
  catch (const std::runtime_error& e) {
    std::cout << "Test Rectangular Matrix: Passed (" << e.what() << ")" << std::endl;
  }

  delete[] b;
  b = nullptr;

  return nError;
}

int testDiagonalSparseSystem() {

  int nError = 0;

  const int n = 5;
  const int nnz = 5;

  SparseSolver::SparseMatrix A(n, n, nnz);
  double* b = new double[5] {0.0};
  b[0] = 1.0;
  b[1] = 2.0;
  b[2] = 3.0;
  b[3] = 4.0;
  b[4] = 5.0;

  int rows[5] = { 0, 1, 2, 3, 4 };
  int cols[5] = { 0, 1, 2, 3, 4 };
  double vals[5] = { 10.0, 10.0, 10.0, 10.0, 10.0 };

  // Define a sparse diagonal matrix
  A.addValue(rows, cols, vals, nnz);

  A.finalize();


  try
  {
  // Solve the system
  SparseSolver::Solver::solve(A, b);

  // Check the solution
  assert(std::fabs(b[0] - 0.1) < EPSILON);
  assert(std::fabs(b[1] - 0.2) < EPSILON);
  assert(std::fabs(b[2] - 0.3) < EPSILON);
  assert(std::fabs(b[3] - 0.4) < EPSILON);
  assert(std::fabs(b[4] - 0.5) < EPSILON);

  std::cout << "Test Diagonal Sparse System: Passed" << std::endl;
  }
  catch (const std::runtime_error& e)
  {
    std::cout << "Test Diagonal System: Failed(" << e.what() << ")" << std::endl;
    nError = 1;
  }

  delete[] b;
  b = nullptr;

  return nError;
}

std::pair<int, int> testAllCases() {
  std::pair<int, int> results;

  int failed = 0, ntests = 0;
  failed += testSimpleSystem();
  ntests++;

  //
  failed += testSingularMatrix();
  ntests++;

  //
  failed += testZeroMatrix();
  ntests++;

  //
  failed += testRectangularMatrix();
  ntests++;

  //
  failed += testDiagonalSparseSystem();
  ntests++;

  //
	results.first = ntests;
	results.second = failed;
	return results;
}

int main() {

	
	auto results = testAllCases();
  std::cout << (results.first - results.second) << "/" << results.first << " tests passed!" << std::endl;
  return 0;
}
