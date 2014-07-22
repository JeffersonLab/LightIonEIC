#!/bin/bash

touch EVTP1.OUT
#touch EVTP2.OUT
#touch EVTP3.OUT
#touch EVTP4.OUT
#touch EVTP5.OUT
#touch EVTP6.OUT
#touch EVTP7.OUT

gfortran -c *.f
bash tag_app_evtp.link
./tag_app_evtp
rm *.o
