#ifndef _KACZ_
#define _KACZ_

double* RK_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* NSSRK_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* GSSRK_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* GRK_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* RK_cyclic_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* RK_rand_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* RK_norep_rand_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* RK_norep_rand_noshuffle_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* RK_quasirand_halton_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* RK_quasirand_sobol_resSC(int M, int N, double**& A, double*& b, double eps, int max_it, long long& it, int n_runs);

double* RK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* NSSRK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* GSSRK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* GRK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* REK_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* RGS_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* RK_cyclic_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* RK_rand_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* RK_norep_rand_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* RK_norep_rand_noshuffle_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* RK_quasirand_halton_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* RK_quasirand_sobol_errorSC(int M, int N, double**& A, double*& b, double*& x, double eps, long long& it, int n_runs);

double* RK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* NSSRK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* GSSRK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* GRK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* REK_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* RGS_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* RK_cyclic_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* RK_rand_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* RK_norep_rand_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* RK_norep_rand_noshuffle_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* RK_quasirand_halton_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

double* RK_quasirand_sobol_maxit(int M, int N, double**& A, double*& b, int max_it, int n_runs);

#endif