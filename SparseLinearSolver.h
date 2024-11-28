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
		int numRows;    // Number of rows
		int numCols;    // Number of columns
		int numNonZero; // Number of non-zero values

		SparseMatrix(int rows, int cols, int nnz);
		~SparseMatrix();

		void addValue(int* row, int* col, double* value, int nVals = 1);
		void finalize();

	private:

		bool finalized;
		int currentNNZ;
	};

	class Solver {
	public:
		static void solve(SparseMatrix& A, double *& b);
	private:
		static void forwardElimination(const SparseMatrix& A, double*& b);
		static void backwardSubstitution(const SparseMatrix& A, double*& b);
	};

	

}

#endif // SPARSE_LINEAR_SOLVER_H