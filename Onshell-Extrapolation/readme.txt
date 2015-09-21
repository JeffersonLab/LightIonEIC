This directory contains code to generate electron-deuteron scattering cross section data from
Misak's and Christian's models and was also designed to be used with the event generator
for extract F2N from virtual data.   This allowed for testing of different order fits
as well as different t' ranges.

The the directory is divided into the sub-directories christians_model misak_model and data_reading.
The cross section data is generated from modified fortran code located in the christians_model and 
misaks_model directories. There are simple bash scripts titled make_script.sh which can be run with

bash make_script.sh

These will compile (using gfortran) the fortran code in that directory, link files, execute and then cleanup.

The data is outputted to text files which are read by ROOT macros (written in C++) in data_reading. 

Within these 3 main sub-directories there are other subdirectories which contain different versions of the code
to do different things. The subdirectories within christians_model and misaks_model are read by the directory of
the same name in data_reading. 

The var_f2n directories are used to just plot the cross section data and perhaps look at scaling with different f2ns.
Note that in data_reading/var_f2n/ the c_and_m_comparison.cpp macro doesn't read the data from the output files directly 
but rather reads the root files in that directory, so you must call c_model_f2n_reader.cpp and m_model_f2n_reader.cpp to 
generate the root
files before calling c_and_m_comparison.cpp.

The polarized directory contains the newest version of christians code with polarization features. It has not been plotted
or analyzed in data_reading yet because so far the asymmetry measurement was just a single value. 

The var_x directory is a leftover from when I was playing around with comparing the 2 models for different Bjorken X's. 
I left the misak_model version of it because it shows how to change around different parameters in that code, but it is
not actually used. 

The extract_f2n directories are used to play around with the extraction of F2N from cross section data. The data in the output files
is not just cross section but rather cross section / pole factor (part of the extraction process is dividing by pole factor). 



Probably the most useful, interesting and complex file is data_reading/extract_f2n/c_model_extraction_mc.cpp. 
It reads the cross section / POLE data outputted by christians_model, it can use monte carlo methods to 
add a random gaussian error to the data and then do many fits of different orders and t' ranges for different 
iterations of the randomization.It is able to do the 3 different methods of extraction of F2N 
(model independent, model dependent and blind test). 
    Note that the blind test requires there to be 2 different input data files, one with the nominal F2N to generate the
theoretical shape and the other with the blind F2N.     
    Note that there is a weird memory issue when doing the model dependent or blind test modes that I wasn't able to solve. Those
modes use a lot of RAM to run (around 5GB for a 500 iteration run) so doing more iterations is not advised unless you have the RAM
for it.
    Also if you want to rerun the c_model_extraction_mc.cpp macro you must quit ROOT and the start it again. If you run it twice in 
the same session it will give you a seg fault. I believed this has something to do with the Audo_dict files that ROOT
generates and leaves in the directory.

I have tried to comment this file well, but as you can see it was written in a somewhat sloppy manner because of time constraints. 
However it does work well if you know how to use it. 

