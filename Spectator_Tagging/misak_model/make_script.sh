# !/bin/bash
#

touch 'MISAK_DATA.OUT'
gfortran -o semi_eic semi_eic.f
./semi_eic
