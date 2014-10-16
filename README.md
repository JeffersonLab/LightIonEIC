Polarized & unpolarized light ions will play a major role in the physics program of an electron ion collider.   In this git repository, the codes that are being used to study the electron light-ion collision are stored and documented. 

Collider - This code uses simulates electron ion collisions, including the cross angle, as well as a D(e,e'p_recoil)X cross section calculation.  The output of this code is setup to be passed to the full detector simulation.  This code has been cross checked against the fixed target code where in the limit of a zero momentum hardon, both codes give the same answer.

Fixed-Target - The fixed target code simulates the D(e,e'p_recoil) reaction for a fixed target.  This code can be is to test the collider code where in the limit of a zero momentum, intial deuteron, both codes should give the same answer.

Hyde - The Hyde directory contains code to calculate the electron-ion collision in the lab frame and then boost to a frame where the ion has zero momentum.   This "fix target" frame is useful for many codes.   The code can also do the transformation back to the lab frame.

Lorentz Transform - This is a basic utility to make the Lorentz transforms from the lab frame to the target rest frame.   

Spectator Tagging - He we took the physics model of Christian Weiss and studied the precision of on-shell extrapolations given different assumptions for data quality (randon and systematic errors) as well as studying the effect of doing a model independent vs. model depenedent extraction.

-

git clone https://github.com/JeffersonLab/LightIonEIC

