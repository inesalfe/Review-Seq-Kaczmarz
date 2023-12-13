#!/bin/bash

# bash bash/seq/NSSRK_N_1000.sh

export OMP_NUM_THREADS=1

N=(1000)
M=(2000 4000 20000 40000 80000)

echo "--- Remove Old Files ---"

rm outputs/seq/NSSRK_N_1000.txt

echo "--- NSSRK Dense ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe NSSRK dense 10 1E-8 $m $n >> outputs/seq/NSSRK_N_1000.txt
			echo "NSSRK dense $m $n"
		fi
	done
done