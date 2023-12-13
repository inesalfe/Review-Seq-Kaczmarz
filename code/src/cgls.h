#ifndef _CGLS_
#define _CGLS_

#include <vector>

double* cglsSolve(int M, int N, int& it, double**& A, double*& b, double eps);

#endif
