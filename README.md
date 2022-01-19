PLACEHOLDER README
==================

this code takes MD data of some pigments along with corresponding spectroscopic data and calculates couplings, excitonic states, absorption and fluorescence spectra and population dynamics using Redfield theory and VERA.
it's not particularly convenient to use (or well written!) yet, pls bear with me :)

throughout this repo I've labelled pigments by ligand codes from RCSB PDB; e.g. chlorophyll a is CLA, chlorophyll b is CHL, fucoxanthin is A86, etc. etc.

lineshapes for individual pigments are in the lineshape folder; the only reason to mess with these is if the parameters or spectral desnsity ansatzes change.

REQUIREMENTS
------------

basic requirements are `gcc/g++, gfortran, LAPACK, FFTW`.
The fortran code is all strictly `f2008`, not that it's very complicated or requires any specific F2008 behaviour, and I think the rest is strict `C++11`; some of it is just old `C99` code with the malloc calls changed to be C++ compliant.
For the plotting scripts you need python 3.7, numpy and matplotlib set up with latex support - I use latex in the line and axis labels because I like them to look nice.
If you have all those things you might need to add something like
```
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{sfmath}'
```
I think that preamble should be escaped correctly but I have it set in my matplotlibrc so I'm not sure.

COMPILATION
-----------

for making the required programs it should be sufficient (is for me) to do `make` in the root folder; it loops over the three subdirs and runs the makefile in each.
to force recompilation of everything do `make clean && make all`.

RUNNING
-------

just do `./scripts/run.py -f NUM` in the root folder for whichever frame number `NUM` you want to look at.
Running it with `-f 0` will loop over all the frames it finds so you can run it easily in batch mode.

the python script is clever enough to notice if the required lineshape data isn't there and will attempt to create it if it isn't.
there are also two parameters which affect everything; T, the temperature we're simulating at and τ, the length (in femtoseconds) that we've calculated line-broadening functions for.
the defaults are set in the run script as 300K and 2048 (setting 2048 instead of, say, 2000, makes basically no difference. but I get to choose because I'm writing the code.); the script checks whether these match the values in the `lineshape` directory, specifically in `lineshape/in/protocol`.
if they don't, it recalculates the line-broadening functions with whatever T, τ you give it.
it doesn't check for the executable files being there because if they're not it'll just fail; you can run `make clean && make all` yourself.

A note on input file setup: the big folder of LHCII snapshots isn't included here because it's like 100MB or whatever, but `scripts/fix_files.sh` is set up to generate input files in the right way.
The fortran code requires four `fmt=%16.8e` columns where the first three are coordinates of an atom and the fourth is its TrEsp charge; the script pulls coordinates from a PDB file and stitches them together with the corresponding TrEsp charges, as long as they're stored in the form 

`ATOM_NUMBER   ATOM_NAME   PARTIAL_CHARGE`

and the atom names match those in the PDB.
It does this by generating dicts and checking keys in both files.

The script `scripts/plot_aw.py` is run after each frame in the run script and plots A(w) and F(w), with corresponding experimental data ripped from the supplementary material in [Tomas Mancal's arxiv paper](http://arxiv.org/abs/1512.00887) (ultimately from a different paper but I forget where).
It also plots them together on the same axis without doing any normalisation to see what the raw data looks like.
Finally the script `scripts/average_spectra.py` loops over all 1000 frames and sums A(w) and F(w), then plots their averages in the same way - reason for this being that an individual frame of MD data can have weird instantaneous couplings due to the movement of the pigments.

The file `vera.def` contains all the parameters we use to describe the carotenoids - changes to these are very non-trivial and the parameters there are taken from fits to TA data of lutein in pyridine, so unless you very much know what you're doing I wouldn't mess with them.
