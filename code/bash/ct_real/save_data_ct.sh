#!/bin/bash

# bash bash/ct_real/save_data_ct.sh

./bin/CK_box_proj_csr_max_it_data.exe ct_real 10 100703 50625 2000000 1000 1
./bin/CK_csr_max_it_data.exe ct_real 10 100703 50625 2000000 1000 1
./bin/RK_box_proj_csr_max_it_data.exe ct_real 10 100703 50625 2000000 1000 1
./bin/RK_csr_max_it_data.exe ct_real 10 100703 50625 2000000 1000 1
./bin/SRKWOR_box_proj_csr_max_it_data.exe ct_real 10 100703 50625 2000000 1000 1
./bin/SRKWOR_csr_max_it_data.exe ct_real 10 100703 50625 2000000 1000 1