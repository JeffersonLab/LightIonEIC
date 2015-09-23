Light Ion Spectator Tagging Project
====================================
Spectator tagging with polarized & unpolarized light ions will play a major role in the physics program of the future electron ion collider.   In an effect to bettter quantify the impact this physics would have, a collaborative effort of theorists and experimentalists proposed a Jefferson Lab's Laboratory-Directed Research and Development (LDRD) FY14/15 Project "Physics potential of polarized light ions with EIC at JLab." The resources developed (theoretical models, computer codes, event generators) as part of this project are available here for scientific applications and further development.  More details about the physics of semi-inclusive scattering, as well as a collection of results and publications can be found at https://www.jlab.org/theory/tag/.

Event-Generator
---------------
This code uses simulates electron ion collisions, including the cross angle, as well as a D(e,e'p_recoil)X cross section calculation.  The output of this code is setup to be passed to the full detector simulation.  This code has been cross checked against a fixed target code, found in the toolbox directory, where in the limit of a zero momentum hardon, both codes give the same answer.   Detailed documentaion for the event generator can be found at https://github.com/JeffersonLab/LightIonEIC/blob/master/Event-Generator/Documenation/STEGs.pdf.

Onshell-Extrapolation
---------------------
Using a physics model, we  studied the precision of on-shell extrapolations given different assumptions for data quality (randon and systematic errors) as well as studying the effect of doing a model independent vs. model depenedent extraction.  

Theory-Codes
------------
Symbolic links to the theory codes currently included with the event generator.

Toolbox
-------

Smaller codes that were developed during the project for debugging:

* Fixed-Target - The fixed target code simulates the D(e,e'p_recoil) reaction for a fixed target.  This code can be is to test the collider code where in the limit of a zero momentum, intial deuteron, both codes should give the same answer.

* Lorentz-Transform - A basic utility to make the Lorentz transforms from the lab frame to the target rest frame.   

* Spectator-Tagging - This directory contains code started by Charles Hyde to calculate the electron-ion collision in the lab frame and then boost to a frame where the ion has zero momentum.   This "fix target" frame is useful for many codes.   The code can also do the transformation back to the lab frame.



-

git clone https://github.com/JeffersonLab/LightIonEIC

----

