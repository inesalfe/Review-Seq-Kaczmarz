#!/bin/bash

# bash bash/seq/GSSRK_extra.sh > outputs/progress/GSSRK_extra.txt &

export OMP_NUM_THREADS=1

echo "--- Remove Old Files ---"

rm outputs/seq/GSSRK_extra.txt

echo "--- GSSRK Dense ---"

N=(500)
M=(4000 20000 40000 80000)

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe GSSRK dense 10 1E-8 $m $n >> outputs/seq/GSSRK_extra.txt
			echo "GSSRK dense $m $n"
		fi
	done
done

N=(2000)
M=(20000 40000 80000)

for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
			then
			./bin/solveDense.exe GSSRK dense 10 1E-8 $m $n >> outputs/seq/GSSRK_extra.txt
			echo "GSSRK dense $m $n"
		fi
	done
done