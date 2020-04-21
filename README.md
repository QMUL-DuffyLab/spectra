PLACEHOLDER README
==================


I haven't yet written any of the convenience stuff which will make this code usable by anyone but me.
I also have not made it particularly readable yet so bear with me lol, make it work then make it nice

throughout this repo I've labelled pigments by ligand codes from RCSB PDB; e.g. chlorophyll a is CLA, chlorophyll b is CHL, fucoxanthin is A86, etc. etc.

lineshapes for individual pigments are in the lineshape folder; the only reason to mess with these is if the parameters change.

in the LHCII directory are chlorophyll TrESP files from Renger and Muh's paper; the final eight columns are q(0,0) and q(1,1) (S0 and S1) charges obtained using three different DFT functionals and HF-CIS.
I do not know enough about quantum chemistry to understand what the differences between these are in detail.

for the couplings between pigments do 

`gfortran -g coupling_calc.f08 -llapack -std=f2008 -ffree-form -Wall -fcheck=bounds -o coupling_calc`

to make; I haven't bothered with a makefile yet because it's one f08 file. Also requires gfortran and LAPACK obviously.

running `./scripts/run.py` will run the fortran code on a given Hamiltonian - add `-i FCP` to run it on FCP. I'll add more convenience stuff soon. it also generates input files for the spectra-calculating code which can be run standalone by e.g.

`cd spectra && make cleaner && make && cd .. && ./spectra/exec_spectra in/input_spectra.dat in/lineshapes.1 > out/exponent_output.dat`
