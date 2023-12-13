#include "aux_func.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LeastSquaresConjugateGradient;
using namespace std;

// https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html

int main(int argc, char *argv[]) {

	if(argc != 5) {
		cout << "Incorrect number of arguments: Corret usage is ";
		cout << "'./bin/cgls_eigen.exe <data_set> <M> <N> <it_max>'" << endl;
		exit(1);
	}

	int M;
	int N;
	int it_max;
	double** A_matrix;
	double* b_vector;
	double* x_vector;

	M = atoi(argv[2]);
	N = atoi(argv[3]);

	it_max = atoi(argv[4]);

	string matrix_type = argv[1];
	string filename_A;
	string filename_b;
	string filename_x;
	if (matrix_type.compare("dense") == 0) {
		filename_A = "../data/dense/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A_matrix, b_vector, x_vector);
	}
	else if (matrix_type.compare("dense_norm") == 0) {
		filename_A = "../data/dense_norm/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense_norm/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense_norm/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A_matrix, b_vector, x_vector);
	}
	else if (matrix_type.compare("dense_coherent") == 0) {
		filename_A = "../data/dense_coherent/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense_coherent/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense_coherent/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A_matrix, b_vector, x_vector);
	}
	else if (matrix_type.compare("dense_ls") == 0) {
		filename_A = "../data/dense_ls/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense_ls/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense_ls/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A_matrix, b_vector, x_vector);
	}
	else {
		cout << "Incorrect number of arguments: Corret usage is ";
		cout << "'./bin/cgls_eigen.exe <data_set> <M> <N> <it_max>'" << endl;
		exit(1);
	}

	MatrixXd m(M,N);
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				m(i,j) = A_matrix[i][j];

	VectorXd b(M);
		for (int i = 0; i < M; i++)
			b(i) = b_vector[i];

	VectorXd x(N);
		for (int i = 0; i < N; i++)
			x(i) = x_vector[i];

	delete[] A_matrix[0];
	delete[] A_matrix;
	delete[] b_vector;
	delete[] x_vector;

	LeastSquaresConjugateGradient<MatrixXd> cgls;

	double time_compute_start = omp_get_wtime();

	cgls.compute(m);

	double time_compute_end = omp_get_wtime();
	double time_compute = time_compute_end-time_compute_start;

	cgls.setMaxIterations(it_max);
	VectorXd x_sol = cgls.solve(b);

	cgls.setTolerance(cgls.error());

	double start = omp_get_wtime();

	x_sol = cgls.solve(b);

	double end = omp_get_wtime();

	auto diff = x_sol - x;
	auto res = m*x_sol-b;

	cout << M << " " << N << " " << diff.squaredNorm() << " " << res.squaredNorm() << " ";
	
	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);
	cout.precision(15);

	cout << end - start + time_compute << endl;
	
	return 0;

}