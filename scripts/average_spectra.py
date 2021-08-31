#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
"""

import os
import argparse
import subprocess
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from plot import plot_all

parser = argparse.ArgumentParser(description="Generate control and osc. strength files")
parser.add_argument("-i", "--input_dir", default='out/LHCII',
        help="Relative path to input directory.")
parser.add_argument("-p", "--protein", default='LHCII',
        help="the protein we're looking at ")
parser.add_argument("-c", "--compress", type=int, default=0,
        help="Compress each frame of data - 1 if yes, 0 if no")
parser.add_argument("-d", "--delete", type=int, default=0,
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
eig = np.zeros_like(jij)
lut620len = np.zeros((len(numbers), ) + np.shape(np.loadtxt("{}/{}/lut620len.out".format(args.input_dir, numbers[0]), skiprows=1)))

theta_shape = (len(numbers), ) + np.shape(jij)
exc_620_coupling = np.zeros((len(numbers), np.shape(jij)[0]))
eigvals = np.zeros_like(exc_620_coupling)

gs_pops = np.zeros((len(numbers), 3))
aw_max = np.zeros(len(numbers))
fw_max = np.zeros(len(numbers))
taus = np.zeros((len(numbers), 2))
thetas = np.zeros(theta_shape)
jijs = np.zeros(theta_shape)
eigs = np.zeros(theta_shape)
rijsq = np.zeros(theta_shape)
musq620 = np.zeros(len(numbers))
j620cl = np.zeros((len(numbers), 3))
j620ex = np.zeros((len(numbers), 6))
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

        aw_max[i]       = aw_temp[np.argmax(aw_temp[:, 1])][0]
        fw_max[i]       = fw_temp[np.argmax(fw_temp[:, 1])][0]
        jijs[i]         = np.loadtxt("{}/J_ij.out".format(direc))
        eigs[i]         = np.loadtxt("{}/eigvecs.out".format(direc))
        j620cl[i][0]    = jijs[i][14][9]
        j620cl[i][1]    = jijs[i][14][10]
        j620cl[i][2]    = jijs[i][14][11]
        musq620[i]      = np.sum(np.square(np.loadtxt("{}/mu_site.out".format(direc))[14, :]))
        com             = np.loadtxt("{}/c_o_m.dat".format(direc))
        for j in range(np.shape(jij)[0]):
            for k in range(np.shape(jij)[0]):
                rijsq[i][j][k] = np.sum(np.square(com[j] - com[k]))

        # we want the coupling between lut 1 and the cluster
        # in the site basis, then the top two? excitons
        # by their participation in the cluster, the couplings
        # to those, and their energies
        cluster_part = [eigs[i][j][9]**2  +
                        eigs[i][j][10]**2 + 
                        eigs[i][j][11]**2
                        for j in range(np.shape(jij)[0])]
        exc_620_coupling[i] = np.matmul(jijs[i][14, :], eigs[i])
        # -1:-3:-1 gives us the most cluster-y one and then
        # the second most cluster-y, bc of the way argsort works
        cluster_eigvecs = np.argsort(cluster_part)[-1:-3:-1]
        # print(exc_620_coupling)
        # print(eig[:, cluster_eigvecs[0]])
        eigvals[i] = np.loadtxt("{}/eigvals.out".format(direc))
        j620ex[i][0] = eigvals[i][cluster_eigvecs[0]]
        j620ex[i][1] = cluster_part[cluster_eigvecs[0]]
        j620ex[i][2] = exc_620_coupling[i][cluster_eigvecs[0]]
        j620ex[i][3] = eigvals[i][cluster_eigvecs[1]]
        j620ex[i][4] = cluster_part[cluster_eigvecs[1]]
        j620ex[i][5] = exc_620_coupling[i][cluster_eigvecs[1]]

        # we want to average the square of the eigenvector components
        # since the squares represent the actual participation coefficients
        # likewise the couplings - only squared couplings appear in the rates
        eigs[i] = np.square(eigs[i])
        exc_620_coupling[i] = np.square(exc_620_coupling[i])
        lut620len[i] = np.loadtxt("{}/lut620len.out".format(direc),skiprows=1)

        thetas[i] = np.loadtxt("{}/theta.dat".format(direc))
        kappas[i] = np.loadtxt("{}/kappa.dat".format(direc))
        aws = aws + aw_temp
        fws = fws + fw_temp
        avg_tau = avg_tau + taus[i, 1]
        final_pop[i, 0] = int(number)
        final_pop[i, 1] = np.loadtxt("{}/final_pop.dat".format(direc))
        pop_at_tau = np.loadtxt("{}/gs_pops_at_tau.dat".format(direc))
        if (len(pop_at_tau) == 3):
            gs_pops[i, 0] = pop_at_tau[0]
            gs_pops[i, 1] = pop_at_tau[1]
            gs_pops[i, 2] = pop_at_tau[2]
        else:
            gs_pops[i, 0] = pop_at_tau[0]
            gs_pops[i, 1] = pop_at_tau[1]

        if (args.compress == 1):
            subprocess.run(["zip", "-rm", "{}.zip".format(direc), "{}".format(number)], check=True)

        if (args.delete == 1):
            subprocess.run(["rm", "-rf", "{}".format(number)], check=True)

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
np.savetxt("{}/j620cl.dat".format(args.input_dir), j620cl)
np.savetxt("{}/j620ex.dat".format(args.input_dir), j620ex)
np.savetxt("{}/musq620.dat".format(args.input_dir), musq620)
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
eig_avg   = np.mean(eigs, axis=0)
eig_std   = np.std(eigs, axis=0)
eigvals_avg   = np.mean(eigvals, axis=0)
eigvals_std   = np.std(eigvals, axis=0)
exc_620_avg = np.mean(exc_620_coupling, axis=0)
exc_620_std = np.std(exc_620_coupling, axis=0)
rmsd_avg    = np.sqrt(np.mean(rijsq, axis=0))
rmsd_std    = np.std(rijsq, axis=0)
lut620len_avg   = np.mean(lut620len, axis=0)
lut620len_std   = np.std(lut620len, axis=0)
np.savetxt("{}/jij_average.dat".format(args.input_dir), jij_avg)
np.savetxt("{}/jij_std.dat".format(args.input_dir), jij_std)
np.savetxt("{}/theta_average.dat".format(args.input_dir), theta_avg)
np.savetxt("{}/theta_std.dat".format(args.input_dir), theta_std)
np.savetxt("{}/kappa_average.dat".format(args.input_dir), kappa_avg)
np.savetxt("{}/kappa_std.dat".format(args.input_dir), kappa_std)
np.savetxt("{}/eig_average.dat".format(args.input_dir), eig_avg)
np.savetxt("{}/eig_std.dat".format(args.input_dir), eig_std)
np.savetxt("{}/eigvals_average.dat".format(args.input_dir), eigvals_avg)
np.savetxt("{}/eigvals_std.dat".format(args.input_dir), eigvals_std)
np.savetxt("{}/exc_620_average.dat".format(args.input_dir), exc_620_avg)
np.savetxt("{}/exc_620_std.dat".format(args.input_dir), exc_620_std)
np.savetxt("{}/rmsd_average.dat".format(args.input_dir), rmsd_avg)
np.savetxt("{}/rmsd_std.dat".format(args.input_dir), rmsd_std)
np.savetxt("{}/lut620len_average.dat".format(args.input_dir), lut620len_avg)
np.savetxt("{}/lut620len_std.dat".format(args.input_dir), lut620len_std)

plot_all(args.input_dir)
