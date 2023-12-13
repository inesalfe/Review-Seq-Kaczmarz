#!/bin/bash

# bash bash/seq_norm/cgls_norm_it_N_1000.sh

echo "--- Remove Old Files ---"

rm outputs/seq_norm/cgls_norm_it_N_1000.txt

N=(1000)
M=(4000 20000 40000 80000 160000)

echo "--- Dense Norm CGLS ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
		then	
			export OMP_NUM_THREADS=1
			str=$(./bin/cgls_eigen_it_error.exe dense_norm $m $n 1E-8)
			IFS=' ' read -r -a arr <<< "$str"
			it=${arr[2]}
			wait
			echo "$m $n $it" >> outputs/seq_norm/cgls_norm_it_N_1000.txt
		fi
	done
done