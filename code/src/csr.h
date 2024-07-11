#ifndef _CSR_
#define _CSR_
#include <string>
#include <vector>
using namespace std;

typedef struct 
{
	int row;
	int col;
	double value;
} CSR;

struct row_equal
{
	explicit row_equal(const int s) : row(s) {}

	bool operator () (const CSR& csr) const
	{
		return csr.row == row;
	}

	int row;
};

double sqrNormRow(int& m, int*& row_idx, int*& cols, double*& values);

double dotProductCSR(int& m, int*& row_idx, int*& cols, double*& values, double*& vector);

void scaleVecLine(int& m, int*& row_idx, int*& cols, double*& values, double& scale, double*& vector);

void scaleVecLinePosProj(int& m, int*& row_idx, int*& cols, double*& values, double& scale, double*& vector);

void scaleVecLineBoxProj(int& m, int*& row_idx, int*& cols, double*& values, double& scale, double*& vector);

void scaleDivideVecLine(int& m, int div, int*& row_idx, int*& cols, double*& values, double& scale, double*& vector);

void scaleVecLinePartial(int& m, int i, int j, int*& row_idx, int*& cols, double*& values, double& scale, double*& vector);

void scaleNewVecLine(int& m, int*& row_idx, int*& cols, double*& values, double& scale, double*& vector, double*& vector_out);

void deleteRowsZeros(int& M, int*& row_idx);

void reorderRows(int M, int N, vector<int> row_order, int* row_idx_in, int* cols_in, double* values_in, int*& row_idx_out, int*& cols_out, double*& values_out);

void importSparseMatrixBIN(int M, int N, int& NNZ, string filename_row_idx, string filename_cols, string filename_values, int*& row_idx, int*& cols, double*& values);

#endif