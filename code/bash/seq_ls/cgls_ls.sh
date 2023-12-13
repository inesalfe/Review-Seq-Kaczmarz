#!/bin/bash

# bash bash/seq_ls/cgls_ls.sh

loop_variable=10
# N=(50 100 200 500 750 1000 2000 4000 10000 20000)
# M=(2000 4000 20000 40000 80000 160000)

N=(50 100 200 500 750 1000 2000 4000 10000)
M=(2000 4000 20000 40000 80000)

echo "--- Remove Old Files ---"

rm outputs/seq_ls/cgls_ls.txt

echo "--- Dense CGLS ---"

line=0
for m in ${M[@]}; do
	for n in ${N[@]}; do
		if [ "$m" -gt "$(( 2*n ))" ]
		then	
			export OMP_NUM_THREADS=1
			((line++))
			str="$(sed "${line}q;d" outputs/seq_ls/cgls_it_ls.txt)"
			IFS=' ' read -r -a arr <<< "$str"
			it=${arr[2]}
			wait
			sum=0.0
			sume=0.0
			sumr=0.0
			for (( i=0; i<loop_variable; i++ )) do
				export OMP_NUM_THREADS=1
				str=$(./bin/cgls_eigen.exe dense_ls $m $n $it)
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
			echo "$m $n $it $sum $avge $avgr" >> outputs/seq_ls/cgls_ls.txt
		fi
	done
done