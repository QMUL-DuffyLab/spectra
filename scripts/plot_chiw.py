#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", nargs="+",
                    help="List of data files")
parser.add_argument("-o", "--output_file", help="Output filename")
parser.add_argument("-d", "--dir", help="Root directory of data")
parser.add_argument("-c", "--cmap", default='viridis', help="Colour map to use")

args = parser.parse_args()

print("Plotting χ_i(w)")

N = 14 # make this a parameter later i guess lol
filenames = ["{}/chi_i_{:02d}.dat".format(args.dir, i + 1) for i in range(N)]
mu = np.loadtxt("{}/mu_exciton.out".format(args.dir))
mu_sq = np.array([mu[i, 0]**2 + mu[i, 1]**2 + mu[i, 2]**2 for i in range(N)]).flatten()
eigvals = np.loadtxt("{}/eigvals.out".format(args.dir))
wn = (np.loadtxt(filenames[0])[:, 0].flatten())
# wl = 10000000/wn
chi = np.array([np.column_stack(np.loadtxt(fn)[:, 1]).flatten() for fn in filenames])
chi = np.column_stack((wn, chi.T))

fc = cm.viridis(np.linspace(0, 1, N))
fig, ax = plt.subplots(figsize=(10,6))
# ax.set_xlim([580, 750])
ax.set_xlim([14000, 17000])
plt.xlabel(r'Wavenumber $ (\text{cm}^{-1}) $')
plt.ylabel(r'Intensity (abu)')
for i in range(N):
    label = r'$ \chi_{' + "{}".format(i + 1) + r'}(\omega) \; |\mu|^2 = ' + "{:6.3f}".format(mu_sq[i]) + r' \; \epsilon = ' + "{:5d}".format(int(eigvals[i])) + r' $'
    plt.plot(chi[:, 0], chi[:, i + 1] * mu_sq[i], color=fc[i], label=label)
    ax.fill(chi[:, 0], chi[:, i + 1] * mu_sq[i], color=fc[i], alpha=0.25)

plt.legend(fontsize=10)
plt.savefig("{}/chi_w.pdf".format(args.dir))
plt.close()

# now for chi_bar which makes up F(w)
filenames = ["{}/chi_bar_i_{:02d}.dat".format(args.dir, i + 1) for i in range(N)]
chi = np.array([np.column_stack(np.loadtxt(fn)[:, 1]).flatten() for fn in filenames])
chi = np.column_stack((wn, chi.T))
fig, ax = plt.subplots(figsize=(10,6))
# ax.set_xlim([580, 750])
ax.set_xlim([14000, 17000])
plt.xlabel(r'Wavenumber $ (\text{cm}^{-1}) $')
plt.ylabel(r'Intensity (abu)')
for i in range(N):
    label = r'$ \bar{\chi}_{' + "{}".format(i + 1) + r'}(\omega) \; |\mu|^2 = ' + "{:6.3f}".format(mu_sq[i]) + r' \; \epsilon = ' + "{:5d}".format(int(eigvals[i])) + r' $'
    plt.plot(chi[:, 0], chi[:, i + 1] * mu_sq[i], color=fc[i], label=label)
    ax.fill(chi[:, 0], chi[:, i + 1] * mu_sq[i], color=fc[i], alpha=0.25)

plt.legend(fontsize=10)
plt.savefig("{}/chi_bar_w.pdf".format(args.dir))
plt.close()
