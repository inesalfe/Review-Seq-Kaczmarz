#!/bin/bash

# bash bash/seq/CG.sh > outputs/progress/CG.txt &

bash bash/seq/cg_it_N_1000.sh > outputs/progress/cg_it_N_1000.txt &
wait
bash bash/seq/cg_N_1000.sh > outputs/progress/cg_N_1000.txt &
wait
bash bash/seq/cgls_it_N_1000.sh > outputs/progress/cgls_it_N_1000.txt &
wait
bash bash/seq/cgls_N_1000.sh > outputs/progress/cgls_N_1000.txt &
wait