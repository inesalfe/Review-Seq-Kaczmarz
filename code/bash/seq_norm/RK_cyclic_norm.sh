#!/bin/bash

# bash bash/seq_norm/RK_cyclic_norm.sh

export OMP_NUM_THREADS=1

N=(50 100 200 500 750 1000 2000 4000 10000 20000)
M=(2000 4000 20000 40000 80000 160000)

echo "--- Remove Old Files ---"

rm outputs/seq_norm/RK_cyclic_norm.txt

echo "--- RK_cyclic Dense Norm ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe RK_cyclic dense_norm 10 1E-8 $m $n >> outputs/seq_norm/RK_cyclic_norm.txt
			echo "RK_cyclic dense_norm $m $n"
		fi
	done
done