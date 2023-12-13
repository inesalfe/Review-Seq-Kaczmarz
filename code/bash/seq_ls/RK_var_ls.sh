#!/bin/bash

# bash bash/seq_ls/RK_var_ls.sh > outputs/progress/RK_var_ls.txt &

bash bash/seq_ls/REK_ls.sh > outputs/progress/REK_ls.txt &
wait
bash bash/seq_ls/RGS_ls.sh > outputs/progress/RGS_ls.txt &
wait
bash bash/seq_ls/cgls_it_ls.sh > outputs/progress/cgls_it_ls.txt &
wait
bash bash/seq_ls/cgls_ls.sh > outputs/progress/cgls_ls.txt &
wait