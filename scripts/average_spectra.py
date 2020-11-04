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

if args.recalc is 1:
    print("Summing A(w) and F(w) per frame")
    for i in range(1000):
        if (i % 100) is 0:
            print(".", end='', flush=True)

        direc = "{}/{}".format(args.input_dir, i + 1)
        aws = aws + np.loadtxt("{}/aw.dat".format(direc))
        fws = fws + np.loadtxt("{}/fw.dat".format(direc))
        jij = jij + np.loadtxt("{}/J_ij.out".format(direc))

    print("\nDone.")
    aws = aws / 1000.
    fws = fws / 1000.
    jij = jij / 1000.
else:
    aws = np.loadtxt("{}/aw_average.dat".format(args.input_dir))
    fws = np.loadtxt("{}/fw_average.dat".format(args.input_dir))
    jij = np.loadtxt("{}/jij_average.dat".format(args.input_dir))

np.savetxt("{}/aw_average.dat".format(args.input_dir), aws)
np.savetxt("{}/fw_average.dat".format(args.input_dir), fws)
np.savetxt("{}/jij_average.dat".format(args.input_dir), jij)

# delete 0 wavenumber rows to prevent divide-by-zero warning
aws = np.delete(aws, 0, 0)
fws = np.delete(fws, 0, 0)
# convert to wavelength (nm)
aws[:, 0] = 10000000/aws[:, 0]
fws[:, 0] = 10000000/fws[:, 0]

# experimental data: this filename construction's ugly
if (args.protein is 'LHCII'):
    aw_exp = np.loadtxt("out/{}/aw_exp.dat".format(args.protein), skiprows=1)
    fw_exp = np.loadtxt("out/{}/fw_exp.dat".format(args.protein), skiprows=1)

fig, ax = plt.subplots()
ax.set_xlim([600, 700])
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(aws[:, 0], aws[:, 1]/np.max(aws[:, 1]),
         label=r'$ A(\omega) $')
if (args.protein is 'LHCII'):
    plt.plot(aw_exp[:, 0], aw_exp[:, 1]/np.max(aw_exp[:, 1]),
             label=r'Experiment')

plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("{}/aw_average.pdf".format(args.input_dir))

fig, ax = plt.subplots()
ax.set_xlim([660, 780])
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(fws[:, 0], fws[:, 1]/np.max(fws[:, 1]),
         label=r'$ F(\omega) $')
if (args.protein is 'LHCII'):
    plt.plot(fw_exp[:, 0], fw_exp[:, 1]/np.max(fw_exp[:, 1]),
             label=r'Experiment')

plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("{}/fw_average.pdf".format(args.input_dir))

fig, ax = plt.subplots()
ax.set_xlim([600, 780])
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(aws[:, 0], aws[:, 1], label=r'$ A(\omega) $')
plt.plot(fws[:, 0], fws[:, 1], label=r'$ F(\omega) $')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("{}/aw_fw_average.pdf".format(args.input_dir))
