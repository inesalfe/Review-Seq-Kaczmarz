#include "kacz.h"
#include "aux_func.h"
#include <iostream>
#include <math.h>
#include <omp.h>
using namespace std;

int main (int argc, char *argv[]) {

	if(argc != 7) {
		cout << "Incorrect number of arguments: Corret usage is ";
		cout << "'./bin/solveDense.exe <method> <data_set> <n_runs> <eps> <M> <N>'" << endl;
		exit(1);
	}

	int M;
	int N;
	double* b;
	double* x;
	double* x_sol;
	double** A;

	int n_runs = atoi(argv[3]);
	double eps = atof(argv[4]);

	M = atoi(argv[5]);
	N = atoi(argv[6]);

	double start = omp_get_wtime();

	string matrix_type = argv[2];
	string filename_A;
	string filename_b;
	string filename_x;
	if (matrix_type.compare("dense") == 0) {
		filename_A = "../data/dense/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A, b, x);
	}
	else if (matrix_type.compare("dense_norm") == 0) {
		filename_A = "../data/dense_norm/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense_norm/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense_norm/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A, b, x);
	}
	else if (matrix_type.compare("dense_coherent") == 0) {
		filename_A = "../data/dense_coherent/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense_coherent/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense_coherent/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A, b, x);
	}
	else if (matrix_type.compare("dense_ls") == 0) {
		filename_A = "../data/dense_ls/A_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_b = "../data/dense_ls/b_" + to_string(M) + "_" + to_string(N) + ".bin";
		filename_x = "../data/dense_ls/x_" + to_string(M) + "_" + to_string(N) + ".bin";
		importDenseSystemBIN(M, N, filename_A, filename_b, filename_x, A, b, x);
	}
	else {
		cout << "Incorrect number of arguments: Corret usage is ";
		cout << "'./bin/solveDense.exe <method> <data_set> <n_runs> <eps> <M> <N>'" << endl;
		exit(1);
	}

	long long it;

	double stop;

	string alg = argv[1];
	if (alg.compare("RK") == 0) {
		x_sol = RK_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = RK_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("NSSRK") == 0) {
		x_sol = NSSRK_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = NSSRK_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("GSSRK") == 0) {
		x_sol = GSSRK_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = GSSRK_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("GRK") == 0) {
		x_sol = GRK_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = GRK_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("REK") == 0) {
		x_sol = REK_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = REK_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("RGS") == 0) {
		x_sol = RGS_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = RGS_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("RK_cyclic") == 0) {
		x_sol = RK_cyclic_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = RK_cyclic_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("RK_rand") == 0) {
		x_sol = RK_rand_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = RK_rand_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("RK_norep_rand") == 0) {
		x_sol = RK_norep_rand_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = RK_norep_rand_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("RK_norep_rand_noshuffle") == 0) {
		x_sol = RK_norep_rand_noshuffle_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = RK_norep_rand_noshuffle_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("RK_quasirand_halton") == 0) {
		x_sol = RK_quasirand_halton_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = RK_quasirand_halton_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else if (alg.compare("RK_quasirand_sobol") == 0) {
		x_sol = RK_quasirand_sobol_errorSC(M, N, A, b, x, eps, it, n_runs);
		// stop = omp_get_wtime();
		// cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
		x_sol = RK_quasirand_sobol_maxit(M, N, A, b, it, n_runs);
		stop = omp_get_wtime();
		cout << sqrNormDiff(x_sol, x, N) << " " << stop - start << endl;
	}
	else {
		cout << "Error: Invalid algorithm." << endl;
		delete[] A[0];
		delete[] A;
		delete[] b;
		delete[] x;
		return 0;
	}

	delete[] A[0];
	delete[] A;
	delete[] b;
	delete[] x;
	delete[] x_sol;

	return 0;
}
