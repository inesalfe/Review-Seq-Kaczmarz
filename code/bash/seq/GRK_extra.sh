#!/bin/bash

# bash bash/seq/GRK_extra.sh > outputs/progress/GRK_extra.txt &

export OMP_NUM_THREADS=1

echo "--- Remove Old Files ---"

rm outputs/seq/GRK_extra.txt

echo "--- GRK Dense ---"

N=(500)
M=(4000 20000 40000 80000)

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe GRK dense 10 1E-8 $m $n >> outputs/seq/GRK_extra.txt
			echo "GRK dense $m $n"
		fi
	done
done

N=(2000)
M=(20000 40000 80000)

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe GRK dense 10 1E-8 $m $n >> outputs/seq/GRK_extra.txt
			echo "GRK dense $m $n"
		fi
	done
done