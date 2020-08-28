#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
Eventually - should take keyword arguments e.g. S1-Qy, S2-Qy etc. etc. and generate control file.
Also, should generate the diagonal of the matrix using oscillator strengths for each given pigment - can do this by reading strings from control file but I suppose we won't need to if we're generating the control file here.
"""

import os
import numpy as np # this is extremely wasteful to read in one float but i want to make sure everything else works first
import argparse
import subprocess
import pigment_data

parser = argparse.ArgumentParser(description="Generate control and osc. strength files")
parser.add_argument("-i", "--input_dir", default='LHCII',
        help="Relative path to input directory.")
parser.add_argument("-o", "--output_dir", default='out',
        help="Relative path to output directory.")
parser.add_argument("-f", "--frame", default=1,
        help="MD frame to calculate for - pass 0 to loop over all frames")
parser.add_argument("-T", "--temperature", type=float, default=300.0,
        help="The temperature you want to calculate the spectra for")
parser.add_argument("-t", "--tau", type=int, default=2048,
        help='''The length of the lineshape function - this is in 
              picoseconds. Default is 2048 since it's a power of 2,
              which makes the FFTs slightly more optimised later.
              ''')

args = parser.parse_args()

def get_pigments(input_dir):
    pigment_dirs = []
    numbers = []
    # programatically get pigment name/numbers
    for item in os.scandir(input_dir):
        if item.is_dir():
            '''
            filter(str.isdigit, str(item)) makes item into
            a string instead of a dirEntry and then makes a list
            of every digit in the string. Then we join them with
            nulls to make a number, take the final 3 (ligand codes
            for FX and DD have numbers in), make it an int, append to
            the list of numbers, and sort both arrays based on
            number. don't need numbers so just return pigments
            '''
            code = str(''.join(filter(str.isdigit, str(item))))[-3:]
            numbers.append(int(code))
            pigment_dirs.append(item.name)

    numbers, pigment_dirs = zip(*sorted(zip(numbers, pigment_dirs)))
    return pigment_dirs

def construct_input_files(pigment_dirs, direc, snapshot_number, protein,
    recalc_lineshapes):
    # fortran won't create the directory; do it here
    output_path = "{}/{}/{}".format(direc, args.input_dir, snapshot_number)
    os.makedirs(output_path, exist_ok=True)
    # there must be a nicer way of doing this but i can't think of it:
    # different information needs to be printed to the file based on
    # file name. Maybe a pair of dicts and then a comprehension
    input_file     = "{}/pigments.{}".format(output_path, snapshot_number) 
    energy_file    = "{}/ei.txt".format(output_path) 
    lifetimes_file = "{}/lifetimes.txt".format(output_path) 
    lambda_file    = "{}/lambda.txt".format(output_path) 
    gnt_file       = "{}/gnt.txt".format(output_path) 
    dipole_file    = "{}/dipoles.txt".format(output_path) 
    lineshape_file = "{}/lineshapes.{}".format(output_path, snapshot_number) 
    f = open(input_file, "w")
    g = open(energy_file, "w")
    h = open(lifetimes_file, "w")
    j = open(lambda_file, "w")
    k = open(gnt_file, "w")
    l = open(lineshape_file, "w")
    m = open(dipole_file, "w")
    for p in pigment_dirs:
        gt = "lineshape/out/{}_gt.dat".format(p[0:3])
        if not os.path.isfile(gt):
            print("Lineshape data does not exist for ligand {}. Generating now.".format(p[0:3]))
            os.system("cd lineshape && ./test ./in/prot ./in/{}.def".format(p[0:3]))

        reorg = np.loadtxt("lineshape/out/{}_lambda.dat".format(p[0:3]))[0]
        lineshape = "lineshape/in/{}.def".format(p[0:3])
        print("{}/{}/frame{}.csv".format(input_dir, p, snapshot_number), file=f)
        if protein is 'FCP':
            state = "S{}".format(snapshot_number)
            print(pigment_data.pigment_data[p[0:3]][state]["energy"], file=g)
            print(pigment_data.pigment_data[p[0:3]][state]["lifetime"], file=h)
        else:
            if p in pigment_data.site_energies.keys():
                print(pigment_data.site_energies[p], file=g)
                # print(pigment_data.novod_energies[p], file=g)
            else:
                print(pigment_data.pigment_data[p[0:3]]["S1"]["energy"], file=g)

            print(pigment_data.pigment_data[p[0:3]]["S1"]["lifetime"], file=h)
            print(pigment_data.pigment_data[p[0:3]]["D"], file=m)

        print(reorg, file=j)
        print(gt, file=k)
        print(lineshape, file=l)

    f.close()
    g.close()
    h.close()
    j.close()
    return (input_file, output_path)

with open("./lineshape/in/prot") as f:
    lineshape_dict = dict([tuple(line.rstrip().split(" = ")) for line in f.readlines()])

recalc_lineshapes = (abs(float(lineshape_dict["T"]) - args.temperature) > 1E-9) or (int(lineshape_dict["ns"]) != args.tau)
print(recalc_lineshapes)
input_dir  = os.path.join(os.getcwd(), args.input_dir)
output_dir = os.path.join(os.getcwd(), args.output_dir)
pigment_dirs = get_pigments(input_dir)

# this is so ugly lol needs tidying up in future
def run_frame(i, do_plots):
    input_file, output_path = construct_input_files(pigment_dirs, output_dir, i, args.input_dir, recalc_lineshapes) # NB: assumes input_dir is just the name of the protein
    print("Calculating for frame {}.\n\n".format(output_path))
    print("./couplings/coupling_calc {} {} {}".format(input_file, output_path, output_path))
    print("./spectra/exec_spectra {} {}".format("in/input_spectra.dat", "{}/lineshapes.{}".format(output_path, i)))
    os.system("./couplings/coupling_calc {} {} {}".format(input_file, output_path, output_path))
    os.system("./spectra/exec_spectra {} {}".format("in/input_spectra.dat", "{}/lineshapes.{}".format(output_path, i)))
    if do_plots is not 0:
        os.system("python ./scripts/plot_aw.py -f {}".format(i))
        os.system("python ./scripts/plot_chiw.py -d {}".format(output_path))

if int(args.frame) == 0:
    for i in range(1000):
        run_frame(i + 1, 0) # range starts from 0

else:
    run_frame(args.frame, 1)

