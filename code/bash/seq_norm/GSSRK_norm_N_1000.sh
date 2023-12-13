#!/bin/bash

# bash bash/seq_norm/GSSRK_norm_N_1000.sh > outputs/progress/GSSRK_norm_N_1000.txt &

export OMP_NUM_THREADS=1

N=(1000)
M=(2000 4000 20000 40000 80000)

echo "--- Remove Old Files ---"

rm outputs/seq_norm/GSSRK_norm_N_1000.txt

echo "--- GSSRK Dense Norm ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe GSSRK dense_norm 10 1E-8 $m $n >> outputs/seq_norm/GSSRK_norm_N_1000.txt
			echo "GSSRK dense_norm $m $n"
		fi
	done
done