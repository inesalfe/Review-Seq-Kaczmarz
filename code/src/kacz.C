#include "kacz.h"
#include "aux_func.h"
#include "sobol.h"
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <omp.h>
using namespace std;

double* RK_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator;
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			generator.seed(run+1);
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < it; i++)
				line = dist(generator);
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				line = dist(generator);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				line = dist(generator);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* NSSRK_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	int old_line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator;
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		line = 0;
		old_line = 0;
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			generator.seed(run+1);
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < it; i++) {
				old_line = line;
				line = dist(generator);
				while (line == old_line) {
					line = dist(generator);
				}
			}
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				old_line = line;
				line = dist(generator);
				while (line == old_line) {
					line = dist(generator);
				}
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				old_line = line;
				line = dist(generator);
				while (line == old_line) {
					line = dist(generator);
				}
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* GSSRK_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double * aux = new double[(long)M*(long)M];
	double ** G = new double*[(long)M];
	for (long i = 0; i < M; i++) {
		G[i] = &aux[i * M];
		for (long j = 0; j < M; j++)
			G[i][j] = dotProduct(A[i], A[j], N);
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	unordered_set<int> select_set;
	for (int i = 0; i < M; i++)
		select_set.insert(i);

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator;
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			generator.seed(run+1);
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < it; i++) {
				line = dist(generator);
				while(select_set.count(line) == 0)
					line = dist(generator);
			}
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				line = dist(generator);
				while(select_set.count(line) == 0)
					line = dist(generator);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
				select_set.erase(line);
				for (int j = 0; j < M; j++)
					if(G[line][j] != 0 && select_set.count(j) == 0)
						select_set.insert(j);
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				line = dist(generator);
				while(select_set.count(line) == 0)
					line = dist(generator);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
				select_set.erase(line);
				for (int j = 0; j < M; j++)
					if(G[line][j] != 0 && select_set.count(j) == 0)
						select_set.insert(j);
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	delete[] G[0];
	delete[] G;
	return x_sol;
}

double* GRK_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	double sqr_matrixNorm = 0;
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] sqrNorm_line;
			exit(1);
		}
		sqr_matrixNorm += sqrNorm_line[i];
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	double* res_samp = new double[M];
	double normRes;
	double maxVal;
	double eps_k;
	double temp;
	vector<double> probs(M);
	discrete_distribution<> dist;

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator;
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			generator.seed(run+1);
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < it; i++)
				line = dist(generator);
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				normRes = 0;
				maxVal = 0;
				for (int i = 0; i < M; i++) {
					res_samp[i] = b[i] - dotProduct(A[i], x_k, N);
					temp = res_samp[i]*res_samp[i];
					if (temp/sqrNorm_line[i] > maxVal)
						maxVal = temp/sqrNorm_line[i];
					normRes += temp;
				}
				eps_k = 0.5*(maxVal/normRes + 1/sqr_matrixNorm);
				maxVal = 0;
				for (int i = 0; i < M; i++) {
					if (res_samp[i]*res_samp[i] < eps_k*normRes*sqrNorm_line[i])
						res_samp[i] = 0;
					maxVal += res_samp[i]*res_samp[i];
				}
				for (int i = 0; i < M; i++)
					probs[i] = res_samp[i]*res_samp[i]/maxVal;
				dist = discrete_distribution<>(probs.begin(), probs.end());
				line = dist(generator);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				normRes = 0;
				maxVal = 0;
				for (int i = 0; i < M; i++) {
					res_samp[i] = b[i] - dotProduct(A[i], x_k, N);
					temp = res_samp[i]*res_samp[i];
					if (temp/sqrNorm_line[i] > maxVal)
						maxVal = temp/sqrNorm_line[i];
					normRes += temp;
				}
				eps_k = 0.5*(maxVal/normRes + 1/sqr_matrixNorm);
				maxVal = 0;
				for (int i = 0; i < M; i++) {
					if (res_samp[i]*res_samp[i] < eps_k*normRes*sqrNorm_line[i])
						res_samp[i] = 0;
					maxVal += res_samp[i]*res_samp[i];
				}
				for (int i = 0; i < M; i++)
					probs[i] = res_samp[i]*res_samp[i]/maxVal;
				dist = discrete_distribution<>(probs.begin(), probs.end());
				line = dist(generator);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	delete[] res_samp;
	delete[] sqrNorm_line;
	return x_sol;
}

double* RK_cyclic_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				line = it%M;
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				line = it%M;
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	delete[] sqrNorm_line;
	return x_sol;
}

double* RK_rand_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	mt19937_64 gen(0);
	uniform_int_distribution<> dis(0, M-1);

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			gen.seed(run+1);
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < it; i++)
				line = dis(gen);
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				line = dis(gen);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				line = dis(gen);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	delete[] sqrNorm_line;
	return x_sol;
}

double* RK_norep_rand_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	vector<int> samp_line(M);
	for (int i = 0; i < M; i++)
		samp_line[i] = i;
	mt19937 rng(1);

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				line = it%M;
				if (line == 1)
					shuffle(begin(samp_line), end(samp_line), rng);
				line = samp_line[line];
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				line = it%M;
				if (line == 1)
					shuffle(begin(samp_line), end(samp_line), rng);
				line = samp_line[line];
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	delete[] sqrNorm_line;
	return x_sol;
}

double* RK_norep_rand_noshuffle_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	vector<int> samp_line(M);
	for (int i = 0; i < M; i++)
		samp_line[i] = i;
	mt19937 rng(1);

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		shuffle(begin(samp_line), end(samp_line), rng);
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				line = samp_line[it%M];
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				line = samp_line[it%M];
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	delete[] sqrNorm_line;
	return x_sol;
}

double* RK_quasirand_halton_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	double phi = 1.6180339887498948482;
	double alpha = 1.0/phi;
	double seed;
	mt19937_64 gen(0);
	uniform_real_distribution<> dist(0.0, 1.0);
	double integral;

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		gen.seed(run);
		seed = dist(gen);
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				line = (int)(modf(seed+alpha*it, &integral)*M);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				line = (int)(modf(seed+alpha*it, &integral)*M);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	delete[] sqrNorm_line;
	return x_sol;
}

double* RK_quasirand_sobol_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	long long it_min;
	long long it_mid;
	long long it_max;
	double* x_it_min = new double[N];
	double* x_it_max = new double[N];
	double* x_it_mid = new double[N];
	double* res_it_max = new double[M];
	double* res_it_mid = new double[M];

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator;
		for (int i = 0; i < N; i++) {
			x_it_min[i] = 0;
		}
		it_min = 0;
		it_max = max_it;
		start = omp_get_wtime();
		while(it_max > it_min) {
			generator.seed(run+1);
			it_mid = it_min+(it_max-it_min)/2;
			it = it_min;
			for (int i = 0; i < N; i++) {
				x_k[i] = x_it_min[i];
			}
			while(it < it_mid) {
				it++;
				line = (int)(sobol::sample(it, run)*M);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_mid[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_mid[i] = b[i] - dotProduct(A[i], x_it_mid, N);
			while(it < it_max) {
				it++;
				line = (int)(sobol::sample(it, run)*M);
				scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
				for (int i = 0; i < N; i++) {
					x_k[i] += scale * A[line][i];
				}
			}
			for (int i = 0; i < N; i++) {
				x_it_max[i] = x_k[i];
			}
			for (int i = 0; i < M; i++)
				res_it_max[i] = b[i] - dotProduct(A[i], x_it_max, N);
			if(sqrNorm(res_it_max, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_max[i];
				}
				it_min = it_max;
				it_max = 10*it_max;
			}
			else if(sqrNorm(res_it_mid, M) > eps) {
				for (int i = 0; i < N; i++) {
					x_it_min[i] = x_it_mid[i];
				}
				it_min = it_mid;
			}
			else {
				for (int i = 0; i < N; i++) {
					x_it_max[i] = x_it_mid[i];
				}
				it_max = it_mid;
			}
			if (it_max - it_min <= 1) {
				it = it_max;
				for (int i = 0; i < N; i++) {
					x_k[i] = x_it_max[i];
				}
				break;
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	cout << M << " " << N << " " << duration << " ";
	cout << avg_it << " ";

	delete[] x_k;
	delete[] x_it_min;
	delete[] x_it_max;
	delete[] x_it_mid;
	delete[] res_it_max;
	delete[] res_it_mid;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << sqrNorm(res, M) << " ";
	delete[] res;
	delete[] sqrNorm_line;
	return x_sol;
}

double* RK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			line = dist(generator);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* NSSRK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	int old_line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		old_line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			old_line = line;
			line = dist(generator);
			while (line == old_line) {
				line = dist(generator);
			}
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* GSSRK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double** G = new double*[M];
	for (size_t i = 0; i < M; i++) {
		G[i] = new double[M];
		for (int j = 0; j < M; j++)
			G[i][j] = dotProduct(A[i], A[j], N);
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	unordered_set<int> select_set;
	for (int i = 0; i < M; i++)
		select_set.insert(i);

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			line = dist(generator);
			while(select_set.count(line) == 0)
				line = dist(generator);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
			select_set.erase(line);
			for (int j = 0; j < M; j++)
				if(G[line][j] != 0 && select_set.count(j) == 0)
					select_set.insert(j);
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	for (int i = 0; i < M; i++) {
		delete[] G[i];
	}
	delete[] G;
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* GRK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	double sqr_matrixNorm = 0;
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			delete[] sqrNorm_line;
			exit(1);
		}
		sqr_matrixNorm += sqrNorm_line[i];
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double* res_samp = new double[M];
	double normRes;
	double maxVal;
	double eps_k;
	double temp;
	vector<double> probs(M);
	discrete_distribution<> dist;

	double start;
	double stop;
	double duration = 0;		
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			normRes = 0;
			maxVal = 0;
			for (int i = 0; i < M; i++) {
				res_samp[i] = b[i] - dotProduct(A[i], x_k, N);
				temp = res_samp[i]*res_samp[i];
				if (temp/sqrNorm_line[i] > maxVal)
					maxVal = temp/sqrNorm_line[i];
				normRes += temp;
			}
			eps_k = 0.5*(maxVal/normRes + 1/sqr_matrixNorm);
			maxVal = 0;
			for (int i = 0; i < M; i++) {
				if (res_samp[i]*res_samp[i] < eps_k*normRes*sqrNorm_line[i])
					res_samp[i] = 0;
				maxVal += res_samp[i]*res_samp[i];
			}
			for (int i = 0; i < M; i++)
				probs[i] = res_samp[i]*res_samp[i]/maxVal;
			dist = discrete_distribution<>(probs.begin(), probs.end());
			line = dist(generator);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	delete[] res_samp;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* REK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			exit(1);
		}
	}

	discrete_distribution<> dist_lines(sqrNorm_line.begin(), sqrNorm_line.end());

	vector<double> sqrNorm_col(N);
	for (int i = 0; i < N; i++) {
		sqrNorm_col[i] = sqrNormMatrixCol(A, i, M);
	}

	discrete_distribution<> dist_cols(sqrNorm_col.begin(), sqrNorm_col.end());

	double* x_k = new double[N];
	double* z_k = new double[M];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	int col;
	int counter;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		for (int i = 0; i < M; i++) {
			z_k[i] = b[i];
		}
		line = 0;
		col = 0;
		it = 0;
		counter = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			col = dist_cols(generator);
			scale = dotProductCol(A, col, z_k, M)/sqrNorm_col[col];
			for (int i = 0; i < M; i++) {
				z_k[i] -= scale * A[i][col];
			}
			line = dist_lines(generator);
			scale = (b[line]-z_k[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	delete[] z_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* RGS_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			exit(1);
		}
	}

	vector<double> sqrNorm_col(N);
	for (int i = 0; i < N; i++) {
		sqrNorm_col[i] = sqrNormMatrixCol(A, i, M);
	}

	discrete_distribution<> dist_cols(sqrNorm_col.begin(), sqrNorm_col.end());

	double* x_k = new double[N];
	double* r_k = new double[M];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int col;
	int counter;

	double start;
	double stop;
	double duration = 0;
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		for (int i = 0; i < M; i++) {
			r_k[i] = b[i];
		}
		it = 0;
		col = 0;
		counter = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			col = dist_cols(generator);
			scale = dotProductCol(A, col, r_k, M)/sqrNorm_col[col];
			x_k[col] += scale;
			for (int i = 0; i < M; i++) {
				r_k[i] -= scale * A[i][col];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	delete[] r_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* RK_cyclic_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;		
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			line = it%M;
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* RK_rand_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;		
	long long avg_it = 0;

	mt19937_64 gen(0);
	uniform_int_distribution<> dis(0, M-1);

	for(int run = 0; run < n_runs; run++) {
		gen.seed(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			line = dis(gen);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";
	
	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* RK_norep_rand_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;		
	long long avg_it = 0;

	vector<int> samp_line(M);
	for (int i = 0; i < M; i++)
		samp_line[i] = i;
	mt19937 rng(1);

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			line = it%M;
			if (line == 1)
				shuffle(begin(samp_line), end(samp_line), rng);
			line = samp_line[line];
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* RK_norep_rand_noshuffle_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	vector<int> samp_line(M);
	for (int i = 0; i < M; i++)
		samp_line[i] = i;
	mt19937 rng(1);

	double start;
	double stop;
	double duration = 0;		
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		shuffle(begin(samp_line), end(samp_line), rng);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			line = samp_line[it%M];
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;

		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* RK_quasirand_halton_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double phi = 1.6180339887498948482;
	double alpha = 1.0/phi;
	double seed;
	mt19937_64 gen(0);
	uniform_real_distribution<> dist(0.0, 1.0);
	double integral;

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;		
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		gen.seed(run);
		seed = dist(gen);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			line = (int)(modf(seed+alpha*it, &integral)*M);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* RK_quasirand_sobol_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] x;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;

	double start;
	double stop;
	double duration = 0;		
	long long avg_it = 0;

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(1) {
			it++;
			line = (int)(sobol::sample(it, run)*M);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			if (sqrNormDiff(x_k, x, N) < eps)
				break;
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
		avg_it += it;
	}
	avg_it /= n_runs;
	it = avg_it;
	// cout << M << " " << N << " " << duration << " ";
	// cout << avg_it << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	// double* res = new double[M];
	// for (int i = 0; i < M; i++)
	// 	res[i] = b[i] - dotProduct(A[i], x_sol, N);
	// cout << sqrNorm(res, M) << " ";
	// delete[] res;
	return x_sol;
}

double* RK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			line = dist(generator);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* NSSRK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	int old_line;
	long long it;

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		old_line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			old_line = line;
			line = dist(generator);
			while (line == old_line) {
				line = dist(generator);
			}
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* GSSRK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	discrete_distribution<> dist(sqrNorm_line.begin(), sqrNorm_line.end());

	double** G = new double*[M];
	for (size_t i = 0; i < M; i++) {
		G[i] = new double[M];
		for (int j = 0; j < M; j++)
			G[i][j] = dotProduct(A[i], A[j], N);
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	double start;
	double stop;
	double duration = 0;

	unordered_set<int> select_set;
	for (int i = 0; i < M; i++)
		select_set.insert(i);

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			line = dist(generator);
			while(select_set.count(line) == 0)
				line = dist(generator);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
			select_set.erase(line);
			for (int j = 0; j < M; j++)
				if(G[line][j] != 0 && select_set.count(j) == 0)
					select_set.insert(j);
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	for (int i = 0; i < M; i++) {
		delete[] G[i];
	}
	delete[] G;
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* GRK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	double* sqrNorm_line = new double[M];
	double sqr_matrixNorm = 0;
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] sqrNorm_line;
			exit(1);
		}
		sqr_matrixNorm += sqrNorm_line[i];
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	double* res_samp = new double[M];
	double normRes;
	double maxVal;
	double eps_k;
	double temp;
	vector<double> probs(M);
	discrete_distribution<> dist;

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			normRes = 0;
			maxVal = 0;
			for (int i = 0; i < M; i++) {
				res_samp[i] = b[i] - dotProduct(A[i], x_k, N);
				temp = res_samp[i]*res_samp[i];
				if (temp/sqrNorm_line[i] > maxVal)
					maxVal = temp/sqrNorm_line[i];
				normRes += temp;
			}
			eps_k = 0.5*(maxVal/normRes + 1/sqr_matrixNorm);
			maxVal = 0;
			for (int i = 0; i < M; i++) {
				if (res_samp[i]*res_samp[i] < eps_k*normRes*sqrNorm_line[i])
					res_samp[i] = 0;
				maxVal += res_samp[i]*res_samp[i];
			}
			for (int i = 0; i < M; i++)
				probs[i] = res_samp[i]*res_samp[i]/maxVal;
			dist = discrete_distribution<>(probs.begin(), probs.end());
			line = dist(generator);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	delete[] res_samp;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* REK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	discrete_distribution<> dist_lines(sqrNorm_line.begin(), sqrNorm_line.end());

	vector<double> sqrNorm_col(N);
	for (int i = 0; i < N; i++) {
		sqrNorm_col[i] = sqrNormMatrixCol(A, i, M);
	}

	discrete_distribution<> dist_cols(sqrNorm_col.begin(), sqrNorm_col.end());

	double* x_k = new double[N];
	double* z_k = new double[M];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	int col;
	long long it;
	int counter;

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		for (int i = 0; i < M; i++) {
			z_k[i] = b[i];
		}
		line = 0;
		col = 0;
		it = 0;
		counter = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			col = dist_cols(generator);
			scale = dotProductCol(A, col, z_k, M)/sqrNorm_col[col];
			for (int i = 0; i < M; i++) {
				z_k[i] -= scale * A[i][col];
			}
			line = dist_lines(generator);
			scale = (b[line]-z_k[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	delete[] z_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* RGS_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	vector<double> sqrNorm_line(M);
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			exit(1);
		}
	}

	vector<double> sqrNorm_col(N);
	for (int i = 0; i < N; i++) {
		sqrNorm_col[i] = sqrNormMatrixCol(A, i, M);
	}

	discrete_distribution<> dist_cols(sqrNorm_col.begin(), sqrNorm_col.end());

	double* x_k = new double[N];
	double* r_k = new double[M];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	long long it;
	int col;
	int counter;

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		mt19937 generator(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		for (int i = 0; i < M; i++) {
			r_k[i] = b[i];
		}
		it = 0;
		col = 0;
		counter = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			col = dist_cols(generator);
			scale = dotProductCol(A, col, r_k, M)/sqrNorm_col[col];
			x_k[col] += scale;
			for (int i = 0; i < M; i++) {
				r_k[i] -= scale * A[i][col];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	delete[] r_k;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* RK_cyclic_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			line = it%M;
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* RK_rand_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	double start;
	double stop;
	double duration = 0;

	mt19937_64 gen(0);
	uniform_int_distribution<> dis(0, M-1);

	for(int run = 0; run < n_runs; run++) {
		gen.seed(run+1);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			line = dis(gen);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";
	
	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* RK_norep_rand_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	double start;
	double stop;
	double duration = 0;

	vector<int> samp_line(M);
	for (int i = 0; i < M; i++)
		samp_line[i] = i;
	mt19937 rng(1);

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			line = it%M;
			if (line == 1)
				shuffle(begin(samp_line), end(samp_line), rng);
			line = samp_line[line];
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* RK_norep_rand_noshuffle_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	vector<int> samp_line(M);
	for (int i = 0; i < M; i++)
		samp_line[i] = i;
	mt19937 rng(1);

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		shuffle(begin(samp_line), end(samp_line), rng);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			line = samp_line[it%M];
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;

		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* RK_quasirand_halton_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double phi = 1.6180339887498948482;
	double alpha = 1.0/phi;
	double seed;
	mt19937_64 gen(0);
	uniform_real_distribution<> dist(0.0, 1.0);
	double integral;

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		gen.seed(run);
		seed = dist(gen);
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			line = (int)(modf(seed+alpha*it, &integral)*M);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}

double* RK_quasirand_sobol_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs) {

	double* sqrNorm_line = new double[M];
	for (int i = 0; i < M; i++) {
		sqrNorm_line[i] = sqrNorm(A[i], N);
		if (sqrNorm_line[i] == 0) {
			cout << "Invalid input: matrix with zero norm line" << endl;
			delete[] A[0];
			delete[] A;
			delete[] b;
			delete[] sqrNorm_line;
			exit(1);
		}
	}

	double* x_k = new double[N];
	double* x_sol = new double[N];
	for (int i = 0; i < N; i++) {
		x_sol[i] = 0;
	}
	double scale;
	int line;
	long long it;

	double start;
	double stop;
	double duration = 0;

	for(int run = 0; run < n_runs; run++) {
		for (int i = 0; i < N; i++) {
			x_k[i] = 0;
		}
		line = 0;
		it = 0;
		start = omp_get_wtime();
		while(it < max_it) {
			it++;
			line = (int)(sobol::sample(it, run)*M);
			scale = (b[line]-dotProduct(A[line], x_k, N))/sqrNorm_line[line];
			for (int i = 0; i < N; i++) {
				x_k[i] += scale * A[line][i];
			}
		}
		stop = omp_get_wtime();
		duration += stop - start;
		for (int i = 0; i < N; i++) {
			x_sol[i] += x_k[i];
		}
	}
	cout << M << " " << N << " " << duration << " ";

	delete[] x_k;
	delete[] sqrNorm_line;
	for (int i = 0; i < N; i++) {
		x_sol[i] /= n_runs;
	}
	double* res = new double[M];
	for (int i = 0; i < M; i++)
		res[i] = b[i] - dotProduct(A[i], x_sol, N);
	cout << max_it << " " << sqrNorm(res, M) << " ";
	delete[] res;
	return x_sol;
}