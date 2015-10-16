# molecular-cogs

## UNDER CONSTRUCTION 

This source code (*.py and *nb files) is an implementation of a technique of identifying molecular cogs.
That is, two sepparate parts of a molecule which produce the 

To run all simulations run:
./run_all

To analyze the results, run the analyze.nb Mathematica notebook (you'll need a licence for that, though, but Wolfram CDF Viewer is a free-to-use program for viewing these Mathematica files).
The notebook relies on additional packages (*.m files).


The analysis depends on a slightly modified version of the Open Babel package, accessible at: https://github.com/ponadto/openbabel
The difference lies in two files: /src/forcefield/forcefieldgaff.cpp and /src/forcefield/forcefieldgaff.h where the "GetVariousForces" function was added.
It extracts forces arising from various interaction types: ele, vdw, bond, angl, tors.
These forces are then analyzed



Questions? Contact me: ponadto@gmail.com
