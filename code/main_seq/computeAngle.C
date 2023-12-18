#include "kacz.h"
#include "aux_func.h"
#include <iostream>
#include <math.h>
#include <omp.h>
#include <climits>
using namespace std;

int main (int argc, char *argv[]) {

	if(argc != 4) {
		cout << "Incorrect number of arguments: Corret usage is ";
		cout << "'./bin/computeAngle.exe <data_set> <M> <N>'" << endl;
		exit(1);
	}

	int M;
	int N;
	double* b;
	double* x;
	double* x_sol;
	double** A;

	M = atoi(argv[2]);
	N = atoi(argv[3]);

	double start = omp_get_wtime();

	string matrix_type = argv[1];
	string filename_A;
	string filename_b;
	string filename_x;
	if (matrix_type.compare("dense") == 0) {
		filename_A = "../data/dense/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A, b, x);
	}
	else if (matrix_type.compare("dense_coherent") == 0) {
		filename_A = "../data/dense_coherent/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense_coherent/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense_coherent/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A, b, x);
	}
	else {
		cout << "Incorrect number of arguments: Corret usage is ";
		cout << "'./bin/computeAngle.exe <method> <data_set> <n_runs> <eps> <M> <N>'" << endl;
		exit(1);
	}

	double angle;
	double max_angle = 0;
	for (int i = 0; i < M-1; i++) {
		for (int j = i+1; j < M; j++) {
			angle = dotProduct(A[i], A[j], N)/sqrt(sqrNorm(A[i], N))/sqrt(sqrNorm(A[j], N));
			angle = acos(angle);
			if (abs(angle) > max_angle)
				max_angle = abs(angle);
		}
	}

	cout << max_angle << endl;

	max_angle = 0;
	for (int i = 0; i < M-1; i++) {
		angle = dotProduct(A[i], A[i+1], N)/sqrt(sqrNorm(A[i], N))/sqrt(sqrNorm(A[i+1], N));
		angle = acos(angle);
		if (abs(angle) > max_angle)
			max_angle = abs(angle);
	}

	cout << max_angle << endl;

	delete[] A[0];
	delete[] A;
	delete[] b;
	delete[] x;
	delete[] x_sol;

	return 0;
}
