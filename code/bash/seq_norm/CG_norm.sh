#!/bin/bash

# bash bash/seq_norm/CG_norm.sh > outputs/progress/CG_norm.txt &

bash bash/seq_norm/cg_norm_it_N_1000.sh > outputs/progress/cg_norm_it_N_1000.txt &
wait
bash bash/seq_norm/cg_norm_N_1000.sh > outputs/progress/cg_norm_N_1000.txt &
wait
bash bash/seq_norm/cgls_norm_it_N_1000.sh > outputs/progress/cgls_norm_it_N_1000.txt &
wait
bash bash/seq_norm/cgls_norm_N_1000.sh > outputs/progress/cgls_norm_N_1000.txt &
wait