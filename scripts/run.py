#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
Eventually - should take keyword arguments e.g. S1-Qy, S2-Qy etc. etc. and generate control file.
Also, should generate the diagonal of the matrix using oscillator strengths for each given pigment - can do this by reading strings from control file but I suppose we won't need to if we're generating the control file here.
"""

import os
import time
import numpy as np # this is extremely wasteful to read in one float but i want to make sure everything else works first
import argparse
import subprocess
import pigment_data

parser = argparse.ArgumentParser(description="Generate control and osc. strength files")
parser.add_argument("-i", "--input_dir", default='structures/LHCII',
        help="Relative path to input directory.")
parser.add_argument("-o", "--output_dir", default='out',
        help="Relative path to output directory.")
parser.add_argument("-f", "--frame", default=1,
        help="MD frame to calculate for - pass 0 to loop over all frames")
parser.add_argument("-ps", "--plot", type=int, default=0,
        help="Plot spectra - default is yes, do -p 0 to disable")
parser.add_argument("-pc", "--plot_chiw", type=int, default=0,
        help="Plot exciton spectra - default is yes, do -p 0 to disable")
parser.add_argument("-pr", "--protocol", default="in/protocol",
        help="Protocol file")
parser.add_argument("-T", "--temperature", type=float, default=300.0,
        help="The temperature you want to calculate the spectra for")
parser.add_argument("-t", "--tau", type=int, default=2048,
        help='''The length of the lineshape function - this is in 
              picoseconds. Default is 2048 since it's a power of 2,
              which makes the FFTs slightly more optimised later.
              ''')
parser.add_argument("-c", "--chl_ansatz", type=int, default=0,
        help='''The ansatz to use for chlorophyll a and b: 0 = OBO,
        1 = RENGER, 2 = BIG
        ''')

args = parser.parse_args()

def get_pigments(input_dir):
    pigment_dirs = []
    numbers = []
    # uncomment to do carotenoids/chl c for FCP as well
    # pigment_names = ['CLA', 'CHL', 'KC' , 'NEX', 'LUT', 'XAT', 'A86']
    pigment_names = ['CLA', 'CHL', 'LUT']
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
            if (any(p in str(item) for p in pigment_names)):
                code = str(''.join(filter(str.isdigit, str(item))))[-3:]
                numbers.append(int(code))
                pigment_dirs.append(item.name)

    numbers, pigment_dirs = zip(*sorted(zip(numbers, pigment_dirs)))
    return pigment_dirs

def construct_input_files(pigment_dirs, direc, snapshot_number, protein,
    recalc_lineshapes):
    # fortran won't create the directory; do it here
    output_path = "{}/{}/{}".format(direc, protein, snapshot_number)
    os.makedirs(output_path, exist_ok=True)
    # there must be a nicer way of doing this but i can't think of it:
    # different information needs to be printed to the file based on
    # file name. Maybe a pair of dicts and then a comprehension
    input_file     = "{}/pigments.{}".format(output_path, snapshot_number) 
    energy_file    = "{}/en.txt".format(output_path) 
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
            os.system("./lineshape/test {} lineshape/in/{}.def".format(args.protocol, p[0:3]))

        reorg = np.loadtxt("lineshape/out/{}_lambda.dat".format(p[0:3]))[0]
        lineshape = "lineshape/in/{}.def".format(p[0:3])
        print("{}/{}/frame{}.csv".format(input_dir, p, snapshot_number), file=f)
        if protein is 'FCP':
            state = "S{}".format(snapshot_number)
            print(pigment_data.pigment_data[p[0:3]][state]["energy"], file=g)
            print(pigment_data.pigment_data[p[0:3]][state]["lifetime"], file=h)
        else:
            if p in pigment_data.site_energies.keys():
                # print(pigment_data.site_energies[p], file=g)
                print(pigment_data.novod_energies[p], file=g)
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
    k.close()
    l.close()
    m.close()
    if recalc_lineshapes:
        recalc_lineshapes = False

    return (input_file, output_path)

def run_frame(i, plot_spectra, plot_excitons):
    protein = args.input_dir.split('structures/')[-1]
    input_file, output_path = construct_input_files(pigment_dirs, output_dir, i, protein, recalc_lineshapes) # NB: assumes input_dir is just the name of the protein
    print("Calculating for frame {}.\n\n".format(output_path))
    print("./couplings/coupling_calc {} {} {} {}".format(input_file, output_path, args.temperature, args.tau))
    print("./spectra/exec_spectra {} {} {}".format("in/input_spectra.dat", args.protocol, "{}/lineshapes.{}".format(output_path, i)))
    os.system("./couplings/coupling_calc {} {} {} {}".format(input_file, output_path, args.temperature, args.tau))
    ti = time.time_ns()
    os.system("./spectra/exec_spectra {} {} {}".format("in/input_spectra.dat", args.protocol, "{}/lineshapes.{}".format(output_path, i)))
    tf = time.time_ns()
    print("Time taken (ms): {:6.3f}".format(float(tf - ti) / 1000000))
    if plot_spectra != 0:
        os.system("python ./scripts/plot_aw.py -d {} -f {}".format(output_path, i))

    if plot_excitons != 0:
        os.system("python ./scripts/plot_chiw.py -d {} -n {}".format(output_path, len(pigment_dirs)))

input_dir  = os.path.join(os.getcwd(), args.input_dir)
output_dir = os.path.join(os.getcwd(), args.output_dir)
pigment_dirs = get_pigments(input_dir)

with open(args.protocol) as f:
    lineshape_dict = dict([tuple(line.rstrip().split(" = ")) for line in f.readlines()])

recalc_lineshapes = (abs(float(lineshape_dict["T"]) - args.temperature) > 1E-9) or (int(lineshape_dict["ns"]) != args.tau) or (int(lineshape_dict["chl_ansatz"]) != args.chl_ansatz)
if recalc_lineshapes:
    print("Temperature/tau parameters given to script don't match those in lineshape folder. Recalculating lineshapes.")
    pigments = np.unique([p[0:3] for p in pigment_dirs])
    # make sure we don't get caught in a loop
    with open(args.protocol, 'w') as n:
        n.write("T = {}\n".format(args.temperature))
        n.write("ns = {}\n".format(args.tau))
        n.write("chl_ansatz = {}".format(args.chl_ansatz))
        recalc_lineshapes = False

    for p in pigments:
        print("Recalculating lineshape for pigment {}\n".format(p))
        os.system("./lineshape/test {} lineshape/in/{}.def".format(args.protocol, p))

if int(args.frame) == 0:
    t0 = time.time_ns()
    # get the numbers - all pigment dirs have the same numbers by construction
    numbers = [fn[5:-4] for fn in os.listdir("{}/{}".format(input_dir, pigment_dirs[0]))]
    for i in range(len(numbers)):
        run_frame(numbers[i], 0, 0) # range starts from 0

    t1 = time.time_ns()
    print("Total time taken (s): {:6.3f}".format(float(t1 - t0) / 1000000000))
else:
    run_frame(args.frame, args.plot, args.plot_chiw)

