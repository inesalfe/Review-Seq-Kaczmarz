#!/bin/bash

# bash bash/seq_norm/cgls_norm_N_1000.sh

loop_variable=10
N=(1000)
M=(4000 20000 40000 80000 160000)

echo "--- Remove Old Files ---"

rm outputs/seq_norm/cgls_norm_N_1000.txt

echo "--- Dense Norm CGLS ---"

line=0
for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
		then	
			export OMP_NUM_THREADS=1
			((line++))
			str="$(sed "${line}q;d" outputs/seq_norm/cgls_norm_it_N_1000.txt)"
			IFS=' ' read -r -a arr <<< "$str"
			it=${arr[2]}
			wait
			sum=0.0
			sume=0.0
			sumr=0.0
			for (( i=0; i<loop_variable; i++ )) do
				export OMP_NUM_THREADS=1
				str=$(./bin/cgls_eigen.exe dense_norm $m $n $it)
				IFS=' ' read -r -a arr <<< "$str"
				time=${arr[4]}
				error=${arr[2]}
				res=${arr[3]}
				wait
				sum=$(awk "BEGIN {print $sum+$time; exit}")
				sume=$(awk "BEGIN {print $sume+$error; exit}")
				sumr=$(awk "BEGIN {print $sumr+$res; exit}")
			done
			avge=$(awk "BEGIN {print $sume/$loop_variable; exit}")
			avgr=$(awk "BEGIN {print $sumr/$loop_variable; exit}")
			echo "dense $m $n - it: $it - time: $sum - avg error: $avge - avg res: $avgr"
			echo "$m $n $it $sum $avge $avgr" >> outputs/seq_norm/cgls_norm_N_1000.txt
		fi
	done
done