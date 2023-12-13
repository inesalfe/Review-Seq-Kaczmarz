#!/bin/bash

# bash bash/seq/RK.sh

export OMP_NUM_THREADS=1

N=(50 100 200 500 750 1000 2000 4000 10000 20000)
M=(2000 4000 20000 40000 80000 160000)

echo "--- Remove Old Files ---"

rm outputs/seq/RK.txt

echo "--- RK Dense ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe RK dense 10 1E-8 $m $n >> outputs/seq/RK.txt
			echo "RK dense $m $n"
		fi
	done
done