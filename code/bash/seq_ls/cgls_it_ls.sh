#!/bin/bash

# bash bash/seq_ls/cgls_it_ls.sh

echo "--- Remove Old Files ---"

rm outputs/seq_ls/cgls_it_ls.txt

# N=(50 100 200 500 750 1000 2000 4000 10000 20000)
# M=(2000 4000 20000 40000 80000 160000)

N=(50 100 200 500 750 1000 2000 4000 10000)
M=(2000 4000 20000 40000 80000)

echo "--- Dense CGLS ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
		then	
			export OMP_NUM_THREADS=1
			str=$(./bin/cgls_eigen_it_error.exe dense_ls $m $n 1E-8)
			IFS=' ' read -r -a arr <<< "$str"
			it=${arr[2]}
			wait
			echo "$m $n $it" >> outputs/seq_ls/cgls_it_ls.txt
		fi
	done
done