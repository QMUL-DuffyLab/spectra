#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
For LHCII data only atm - average the A(w) and F(w) data over every snapshot,
write them out, plot them.
"""

import os
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from plot_aw import plot_aw_fw

parser = argparse.ArgumentParser(description="Generate control and osc. strength files")
parser.add_argument("-i", "--input_dir", default='out/LHCII',
        help="Relative path to input directory.")
parser.add_argument("-p", "--protein", default='LHCII',
        help="the protein we're looking at ")
parser.add_argument("-c", "--compress", type=int, default=1,
        help="Compress each frame of data - 1 if yes, 0 if no")
parser.add_argument("-r", "--recalc", type=int, default=1,
        help="Recalculate the average - 1 if yes, 0 if no")
args = parser.parse_args()

# scandir to get the directories only
dirs = [f for f in os.scandir(args.input_dir) if f.is_dir()]
# the frames all have numbers for directory names - ignore any letters
numbers = [f.name for f in dirs if not any(c.isalpha() for c in f.name)]
# sort the numbers
numbers = sorted(numbers, key=lambda x: int(x))

initial_data = np.loadtxt("{}/{}/aw.dat".format(args.input_dir, numbers[0]))
aws = np.zeros_like(initial_data)
fws = np.zeros_like(initial_data)
jij = np.zeros_like(np.loadtxt("{}/{}/J_ij.out".format(args.input_dir, numbers[0])))

theta_shape = (len(numbers), ) + np.shape(jij)

gs_pops = np.zeros((len(numbers), 3))
aw_max = np.zeros(len(numbers))
fw_max = np.zeros(len(numbers))
taus = np.zeros((len(numbers), 2))
thetas = np.zeros(theta_shape)
jijs = np.zeros(theta_shape)
j620612 = np.zeros(len(numbers))
kappas = np.zeros(theta_shape)
final_pop = np.zeros((len(numbers), 2))
avg_tau = 0.0
curdir = os.getcwd()
os.chdir(args.input_dir)

if args.recalc is 1:
    print("Summing A(w) and F(w) per frame")
    for i, number in enumerate(numbers):
        if (i % 100) is 0:
            print(".", end='', flush=True)

        direc = "{}/{}".format(os.getcwd(), number)
        aw_temp = np.loadtxt("{}/aw.dat".format(direc))
        fw_temp = np.loadtxt("{}/fw.dat".format(direc))
        taus[i, 0] = int(number)
        if os.path.isfile("{}/lifetime.dat".format(direc)):
            taus[i, 1] = np.loadtxt("{}/lifetime.dat".format(direc))
        elif os.path.isfile("{}/tau.dat".format(direc)):
            taus[i, 1] = np.loadtxt("{}/tau.dat".format(direc))
        else:
            taus[i, 1] = 0.0

        aw_max[i] = aw_temp[np.argmax(aw_temp[:, 1])][0]
        fw_max[i] = fw_temp[np.argmax(fw_temp[:, 1])][0]
        jijs[i]   = np.loadtxt("{}/J_ij.out".format(direc))
        j620612[i]= jijs[i][14][11]
        thetas[i] = np.loadtxt("{}/theta.dat".format(direc))
        kappas[i] = np.loadtxt("{}/kappa.dat".format(direc))
        aws = aws + aw_temp
        fws = fws + fw_temp
        avg_tau = avg_tau + taus[i, 1]
        final_pop[i, 0] = int(number)
        final_pop[i, 1] = np.loadtxt("{}/final_pop.dat".format(direc))
        pop_at_tau = np.loadtxt("{}/pop_at_tau.dat".format(direc))
        if (len(pop_at_tau) == 49):
            gs_pops[i, 0] = pop_at_tau[0]
            gs_pops[i, 1] = pop_at_tau[15]
            gs_pops[i, 2] = pop_at_tau[32]
        else:
            gs_pops[i, 0] = pop_at_tau[0]
            gs_pops[i, 1] = pop_at_tau[15]

        if (args.compress == 1):
            subprocess.run(["zip", "-rm", "{}.zip".format(direc), "{}".format(number)], check=True)

    print("\nDone.")
    aws = aws / float(len(numbers))
    fws = fws / float(len(numbers))
    jij = jij / float(len(numbers))
    avg_tau = avg_tau / float(len(numbers))
else:
    aws = np.loadtxt("{}/aw_average.dat".format(args.input_dir))
    fws = np.loadtxt("{}/fw_average.dat".format(args.input_dir))
    jij = np.loadtxt("{}/jij_average.dat".format(args.input_dir))
    tau = np.loadtxt("{}/tau_average.dat".format(args.input_dir))


os.chdir(curdir)
np.savetxt("{}/aw_max.dat".format(args.input_dir), aw_max)
np.savetxt("{}/fw_max.dat".format(args.input_dir), fw_max)
np.savetxt("{}/aw_average.dat".format(args.input_dir), aws)
np.savetxt("{}/fw_average.dat".format(args.input_dir), fws)
np.savetxt("{}/tau_average.dat".format(args.input_dir), np.array([avg_tau, np.std(taus[:, 1])]))
np.savetxt("{}/tau.dat".format(args.input_dir), taus)
np.savetxt("{}/gs_pops_at_tau.dat".format(args.input_dir), gs_pops)
np.savetxt("{}/final_pop.dat".format(args.input_dir), final_pop)
np.savetxt("{}/j620612.dat".format(args.input_dir), j620612)
print("<τ> = {}".format(avg_tau))
print("Standard deviation of A(w) max = {} (cm^[-1])".format(np.std(aw_max)))
print("Standard deviation of F(w) max = {} (cm^[-1])".format(np.std(fw_max)))
print("Standard deviation of <τ> = {}".format(np.std(taus[:, 1])))

jij_avg   = np.mean(jijs, axis=0)
jij_std   = np.std(jijs, axis=0)
theta_avg = np.mean(thetas, axis=0)
theta_std = np.std(thetas, axis=0)
kappa_avg = np.mean(kappas, axis=0)
kappa_std = np.std(kappas, axis=0)
np.savetxt("{}/jij_average.dat".format(args.input_dir), jij_avg)
np.savetxt("{}/jij_std.dat".format(args.input_dir), jij_std)
np.savetxt("{}/theta_average.dat".format(args.input_dir), theta_avg)
np.savetxt("{}/theta_std.dat".format(args.input_dir), theta_std)
np.savetxt("{}/kappa_average.dat".format(args.input_dir), kappa_avg)
np.savetxt("{}/kappa_std.dat".format(args.input_dir), kappa_std)

# experimental data: this filename construction's ugly
if (args.protein is 'LHCII'):
    aw_exp = np.loadtxt("out/{}/aw_exp.dat".format(args.protein), skiprows=1)
    fw_exp = np.loadtxt("out/{}/fw_exp.dat".format(args.protein), skiprows=1)
else:
    aw_exp = np.zeros_like(aws)
    fw_exp = np.zeros_like(aws)

draw_maximums = (True, True, True)
plot_aw_fw(aws, fws, aw_exp, fw_exp, draw_maximums, args.input_dir)

fig, ax = plt.subplots()
plt.title(r'Lifetime: avg = $ {} $, $ \sigma = {} $'.format(avg_tau, np.std(taus[:, 1])))
plt.xlabel("Frame")
plt.ylabel(r'$ \left< \tau \right>$')
plt.plot(taus[:, 0], taus[:, 1])
plt.savefig("{}/tau.pdf".format(args.input_dir))
