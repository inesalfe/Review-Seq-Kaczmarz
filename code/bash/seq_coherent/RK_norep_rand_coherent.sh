#!/bin/bash

# bash bash/seq_coherent/RK_norep_rand_coherent.sh

export OMP_NUM_THREADS=1

N=(1000)
M=(4000 20000 40000 80000 160000)

echo "--- Remove Old Files ---"

rm outputs/seq_coherent/RK_norep_rand_coherent.txt

echo "--- RK_norep_rand_coherent Dense ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe RK_norep_rand dense_coherent 10 1E-8 $m $n >> outputs/seq_coherent/RK_norep_rand_coherent.txt
			echo "RK_norep_rand dense_coherent $m $n"
		fi
	done
done