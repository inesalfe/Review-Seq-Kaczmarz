#!/bin/bash

# bash bash/seq_coherent/RK_var_coherent.sh > outputs/progress/RK_var_coherent.txt &

bash bash/seq_coherent/RK_coherent.sh > outputs/progress/RK_coherent.txt &
wait
bash bash/seq_coherent/RK_cyclic_coherent.sh > outputs/progress/RK_cyclic_coherent.txt &
wait
bash bash/seq_coherent/RK_rand_coherent.sh > outputs/progress/RK_rand_coherent.txt &
wait
bash bash/seq_coherent/RK_norep_rand_coherent.sh > outputs/progress/RK_norep_rand_coherent.txt &
wait
bash bash/seq_coherent/RK_norep_rand_noshuffle_coherent.sh > outputs/progress/RK_norep_rand_noshuffle_coherent.txt &
wait
bash bash/seq_coherent/RK_quasirand_halton_coherent.sh > outputs/progress/RK_quasirand_halton_coherent.txt &
wait
bash bash/seq_coherent/RK_quasirand_sobol_coherent.sh > outputs/progress/RK_quasirand_sobol_coherent.txt &
wait