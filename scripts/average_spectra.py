#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
For LHCII data only atm - average the A(w) and F(w) data over every snapshot,
write them out, plot them.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from plot_aw import plot_aw_fw

parser = argparse.ArgumentParser(description="Generate control and osc. strength files")
parser.add_argument("-i", "--input_dir", default='out/LHCII',
        help="Relative path to input directory.")
parser.add_argument("-p", "--protein", default='LHCII',
        help="the protein we're looking at ")
parser.add_argument("-r", "--recalc", type=int, default=1,
        help="Recalculate the average - 1 if yes, 0 if no")
args = parser.parse_args()

initial_data = np.loadtxt("{}/1/aw.dat".format(args.input_dir))
aws = np.zeros_like(initial_data)
fws = np.zeros_like(initial_data)
jij = np.zeros_like(np.loadtxt("{}/1/J_ij.out".format(args.input_dir)))

num_frames = 1000 # NB: this should probably be calculated
aw_max = np.zeros(num_frames)
fw_max = np.zeros(num_frames)
taus = np.zeros(num_frames)
avg_tau = 0.0

if args.recalc is 1:
    print("Summing A(w) and F(w) per frame")
    for i in range(num_frames):
        if (i % 100) is 0:
            print(".", end='', flush=True)

        direc = "{}/{}".format(args.input_dir, i + 1)
        aw_temp = np.loadtxt("{}/aw.dat".format(direc))
        fw_temp = np.loadtxt("{}/fw.dat".format(direc))
        taus[i] = np.loadtxt("{}/tau.dat".format(direc))
        aw_max[i] = aw_temp[np.argmax(aw_temp[:, 1])][0]
        fw_max[i] = fw_temp[np.argmax(fw_temp[:, 1])][0]
        jij = jij + np.loadtxt("{}/J_ij.out".format(direc))
        aws = aws + aw_temp
        fws = fws + fw_temp
        avg_tau = avg_tau + taus[i]

    print("\nDone.")
    aws = aws / float(num_frames)
    fws = fws / float(num_frames)
    jij = jij / float(num_frames)
    avg_tau = avg_tau / float(num_frames)
else:
    aws = np.loadtxt("{}/aw_average.dat".format(args.input_dir))
    fws = np.loadtxt("{}/fw_average.dat".format(args.input_dir))
    jij = np.loadtxt("{}/jij_average.dat".format(args.input_dir))

np.savetxt("{}/aw_max.dat".format(args.input_dir), aw_max)
np.savetxt("{}/fw_max.dat".format(args.input_dir), fw_max)
np.savetxt("{}/aw_average.dat".format(args.input_dir), aws)
np.savetxt("{}/fw_average.dat".format(args.input_dir), fws)
np.savetxt("{}/jij_average.dat".format(args.input_dir), jij)
np.savetxt("{}/tau_average.dat".format(args.input_dir), np.array([avg_tau, np.std(taus)]))
print("Standard deviation of A(w) max = {} (cm^[-1])".format(np.std(aw_max)))
print("Standard deviation of F(w) max = {} (cm^[-1])".format(np.std(fw_max)))
print("Standard deviation of <τ> = {}".format(np.std(taus)))

# experimental data: this filename construction's ugly
if (args.protein is 'LHCII'):
    aw_exp = np.loadtxt("out/{}/aw_exp.dat".format(args.protein), skiprows=1)
    fw_exp = np.loadtxt("out/{}/fw_exp.dat".format(args.protein), skiprows=1)
else:
    aw_exp = np.zeros_like(aws)
    fw_exp = np.zeros_like(aws)

draw_maximums = (True, True, True)
plot_aw_fw(aws, fws, aw_exp, fw_exp, draw_maximums, args.input_dir)
