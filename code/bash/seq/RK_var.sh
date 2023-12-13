#!/bin/bash

# bash bash/seq/RK_var.sh > outputs/progress/RK_var.txt &

bash bash/seq/RK.sh > outputs/progress/RK.txt &
wait
bash bash/seq/RK_cyclic.sh > outputs/progress/RK_cyclic.txt &
wait
bash bash/seq/RK_rand.sh > outputs/progress/RK_rand.txt &
wait
bash bash/seq/RK_norep_rand.sh > outputs/progress/RK_norep_rand.txt &
wait
bash bash/seq/RK_norep_rand_noshuffle.sh > outputs/progress/RK_norep_rand_noshuffle.txt &
wait
bash bash/seq/RK_quasirand_halton.sh > outputs/progress/RK_quasirand_halton.txt &
wait
bash bash/seq/RK_quasirand_sobol.sh > outputs/progress/RK_quasirand_sobol.txt &
wait
bash bash/seq/NSSRK_N_1000.sh > outputs/progress/NSSRK_N_1000.txt &
wait
bash bash/seq/GSSRK_N_1000.sh > outputs/progress/GSSRK_N_1000.txt &
wait
bash bash/seq/GRK_N_1000.sh > outputs/progress/GRK_N_1000.txt &
wait