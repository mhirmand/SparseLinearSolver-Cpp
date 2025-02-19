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
		int numRows;       // Number of rows
		int numCols;       // Number of columns
		int numNonZero;    // Number of non-zero values

		SparseMatrix(int rows, int cols, int nnz);
		~SparseMatrix();

		// Add values to the matrix at the specified row and column indices.
		void addValue(int* row, int* col, double* values, int nVals = 1); 

		// Finalize the matrix after adding all the values.
		void finalize(); 

		bool isFinalized() const { return finalized; }

	private:

		bool finalized;
		int currentNNZ;
	};

	class Solver {
	public:
		// Solve the linear system Ax = b using Gaussian elimination.
		//		- The solution is stored in the b array. 
		//		- A is modified in place.
		static void solve(SparseMatrix& A, double *& b, const double tolerance = 1e-10);
	private:
		// Perform forward elimination to convert the matrix to upper triangular form.
		static void forwardElimination(const SparseMatrix& A, double*& b, const double tolerance = 1e-10);

		// Perform backward substitution to solve the system.
		static void backwardSubstitution(const SparseMatrix& A, double*& b);

	};
}

#endif // SPARSE_LINEAR_SOLVER_H