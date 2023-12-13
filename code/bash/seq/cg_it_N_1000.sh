#!/bin/bash

# bash bash/seq/cg_it_N_1000.sh

echo "--- Remove Old Files ---"

rm outputs/seq/cg_it_N_1000.txt

N=(1000)
M=(4000 20000 40000 80000 160000)

echo "--- Dense CG ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
		then	
			export OMP_NUM_THREADS=1
			str=$(./bin/cg_eigen_it_error.exe dense $m $n 1E-8)
			IFS=' ' read -r -a arr <<< "$str"
			it=${arr[2]}
			wait
			echo "$m $n $it" >> outputs/seq/cg_it_N_1000.txt
		fi
	done
done