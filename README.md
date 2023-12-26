# A Survey of a Class of Iterative Row-Action Methods: The Kaczmarz Method

This repository containts the source code for testing several variations of the Kaczmarz method.

The paper describing the obtained results was submitted to SIAM Review and its preprint is available in ArXiv.

# Requirements

* C++11
* Eigen Libreary (version 3.4)

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

<method>


