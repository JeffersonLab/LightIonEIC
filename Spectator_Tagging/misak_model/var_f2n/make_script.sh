# !/bin/bash
#

touch 'MISAK_DATA_1.OUT'
touch 'MISAK_DATA_2.OUT'
touch 'MISAK_DATA_3.OUT'
gfortran -o semi_eic semi_eic.f c_res.f
./semi_eic
