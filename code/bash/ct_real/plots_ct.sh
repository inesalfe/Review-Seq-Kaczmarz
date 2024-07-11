#!/bin/bash

# bash bash/ct_real/plots_ct.sh

python3 plots/ct_real/CK.py ct_real 100703 50625 2000000 1
python3 plots/ct_real/CK_box_proj.py ct_real 100703 50625 2000000 1
python3 plots/ct_real/RK.py ct_real 100703 50625 2000000 1
python3 plots/ct_real/RK_box_proj.py ct_real 100703 50625 2000000 1
python3 plots/ct_real/SRKWOR.py ct_real 100703 50625 2000000 1
python3 plots/ct_real/SRKWOR_box_proj.py ct_real 100703 50625 2000000 1

# 1891001 12.2924
# 511001 4.99425
# 745001 9.2523
# 591001 4.9138
# 579001 11.4398
# 680001 4.95699