# A Survey of a Class of Iterative Row-Action Methods: The Kaczmarz Method

This repository containts the source code for testing several variations of the Kaczmarz method.

The paper describing the obtained results was submitted to SIAM Review and its preprint is available in ArXiv.

# Requirements

* C++11
* Eigen Library (version 3.4)

### Python Packages (for plots)

* numpy
* matplotlib
* sys
* scipy
* scikit-spatial

# Usage

Importante note: Run all commands inside the code directory.

Compile all files:

```
make
```

Generate Consistent Data Sets:

```
./bin/genConsistDataSets.exe
```

Generate Inconsistent Data Sets:

```
./bin/genLSDataSets.exe
```

Run Kaczmarz and Variations:

```
./bin/solveDense.exe <method> <data_set> <n_runs> <eps> <M> <N>
```

### Parameter options for the algorithms

__<method>__

Options  | Algorithm Name
------------- | -------------
RK  | Randomized Kaczmarz
NSSRK  | Non-Repetitive Selective Set Randomized Kaczmarz
GSSRK  | Gramian Selective Set Randomized Kaczmarz
GRK  | Greedy Randomized Kaczmarz
REK  | Randomized Extended Kaczmarz
RGS  | Randomized Gauss-Seidel
RK_cyclic  | Cyclic Kaczmarz
RK_rand  | Simple Randomized Kaczmarz
RK_norep_rand  | Simple Randomized Kaczmarz Without Replacement (with shuffling)
RK_norep_rand_noshuffle  | Simple Randomized Kaczmarz Without Replacement (without shuffling)
RK_quasirand_halton  | Simple Randomized Kaczmarz with Quasirandom Numbers from Halton Sequence
RK_quasirand_sobol  | Simple Randomized Kaczmarz with Quasirandom Numbers from Sobol Sequence

__<data_set>__

All data sets are dense and entries are always generated from a normal distribution.

Options  | Average | Standard Deviation
------------- | ------------- | -------------
dense  | -5 to 5 | 1 to 20
dense_norm  | 0 | 1
dense_coherent  | 2 | 20
dense_ls  | -5 to 5 | 1 to 20

For "dense_coherent" only 5 elements change in every row.
For "dense_ls" gaussian noise was added to the right hand side of the every equation in the system.

__<n_runs>__

Value used for simulations: 10.

__<eps>__

Value used for simulations: $10^{-8}$.

__<M> and <N>__

M value Options: 2000, 4000, 20000, 40000, 80000, 160000.
N value Options: 50, 100, 200, 500, 750, 1000, 2000, 4000, 10000, 20000.

Overdetermined systems were created such that M/N > 2.
