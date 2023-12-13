#!/bin/bash

# bash bash/seq/REK_ls.sh

export OMP_NUM_THREADS=1

# N=(50 100 200 500 750 1000 2000 4000 10000 20000)
# M=(2000 4000 20000 40000 80000 160000)

N=(50 100 200 500 750 1000 2000 4000 10000)
M=(2000 4000 20000 40000 80000)

echo "--- Remove Old Files ---"

rm outputs/seq_ls/REK_ls.txt

echo "--- REK Dense Error ---"

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe REK dense_ls 10 1E-8 $m $n >> outputs/seq_ls/REK_ls.txt
			echo "REK dense_ls $m $n"
		fi
	done
done