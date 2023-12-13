#!/bin/bash

# bash bash/seq_norm/RK_var_norm.sh > outputs/progress/RK_var_norm.txt &

bash bash/seq_norm/RK_norm.sh > outputs/progress/RK_norm.txt &
wait
bash bash/seq_norm/RK_cyclic_norm.sh > outputs/progress/RK_cyclic_norm.txt &
wait
bash bash/seq_norm/RK_rand_norm.sh > outputs/progress/RK_rand_norm.txt &
wait
bash bash/seq_norm/RK_norep_rand_norm.sh > outputs/progress/RK_norep_rand_norm.txt &
wait
bash bash/seq_norm/RK_norep_rand_noshuffle_norm.sh > outputs/progress/RK_norep_rand_noshuffle_norm.txt &
wait
bash bash/seq_norm/RK_quasirand_halton_norm.sh > outputs/progress/RK_quasirand_halton_norm.txt &
wait
bash bash/seq_norm/RK_quasirand_sobol_norm.sh > outputs/progress/RK_quasirand_sobol_norm.txt &
wait
bash bash/seq_norm/NSSRK_norm_N_1000.sh > outputs/progress/NSSRK_norm_N_1000.txt &
wait
bash bash/seq_norm/GSSRK_norm_N_1000.sh > outputs/progress/GSSRK_norm_N_1000.txt &
wait
bash bash/seq_norm/GRK_norm_N_1000.sh > outputs/progress/GRK_norm_N_1000.txt &
wait