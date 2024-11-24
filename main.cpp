#include <iostream>
#include "SparseLinearSolver.h"

// using namespace SparseSolver;

int main() {
	const size_t n = 3; // 3x3 system
	const size_t nnz = 7; // Number of non-zero elements

	SparseSolver::SparseMatrix A(n, n, nnz);
	double* b = new double[n];
	b[0] = 8.0;
	b[1] = -11.0;
	b[2] = -3.0;

	// double x[n] = { 0 };
	int rows[] = { 0, 0, 0, 1, 1, 1, 2 };
	int cols[] = { 0, 1, 2, 0, 1, 2, 2 };
	double vals[] = { 2.0, 1.0, -1.0, -3.0, -1.0, 2.0, 2.0 };

	// Add non-zero values
	A.addValue(rows, cols, vals, 7);

	// Finalize the matrix
	A.finalize();

	try {
		SparseSolver::Solver::solve(A, b);

		std::cout << "Solution:" << std::endl;
		for (size_t i = 0; i < n; ++i) {
			std::cout << b[i] << " ";
		}
		std::cout << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}

	return 0;
}
