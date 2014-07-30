#!/bin/bash
#

touch 'MISAK_DATA_1.OUT'
touch 'MISAK_DATA_2.OUT'
touch 'MISAK_DATA_3.OUT'
gfortran -o semi_eic semi_eic.f tag_user_res.f
./semi_eic
