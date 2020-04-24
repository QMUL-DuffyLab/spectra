PLACEHOLDER README
==================


I haven't yet written any of the convenience stuff which will make this code usable by anyone but me.
I also have not made it particularly readable yet so bear with me lol, make it work then make it nice

basic requirements are `gcc, LAPACK, FFTW`

throughout this repo I've labelled pigments by ligand codes from RCSB PDB; e.g. chlorophyll a is CLA, chlorophyll b is CHL, fucoxanthin is A86, etc. etc.

lineshapes for individual pigments are in the lineshape folder; the only reason to mess with these is if the parameters change.

to calculate couplings between pigments and the resulting absorption spectrum just do `./scripts/run.py`; the default is just to run for one frame of LHCII but that can be changed quite easily, hopefully the script is readable enough to figure that out. To make the coupling code do

`gfortran -g coupling_calc.f08 -llapack -std=f2008 -ffree-form -Wall -fcheck=bounds -o coupling_calc`

to make; I haven't bothered with a makefile yet because it's one f08 file. The lineshape and spectra-calculating codes have makefiles so just `cd` to the lower directories and do `make cleaner && make` or whatever.

A note on input file setup: the big folder of LHCII snapshots isn't included here because it's like 100MB or whatever, but `scripts/fix_files.sh` is set up to generate input files in the right way. The fortran code requires four `fmt=%16.8e` columns where the first three are coordinates of an atom and the fourth is its TrEsp charge; the script pulls coordinates from a PDB file and stitches them together with the corresponding TrEsp charges, as long as they're stored in the form 

`ATOM_NUMBER   ATOM_NAME   PARTIAL_CHARGE`

and the atom names match those in the PDB. It does this by generating dicts and checking keys in both files.
