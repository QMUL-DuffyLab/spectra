#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", nargs="+",
                    help="List of data files")
parser.add_argument("-n", type=int, help="args.number of pigments")
parser.add_argument("-o", "--output_file", help="Output filename")
parser.add_argument("-d", "--dir", help="Root directory of data")
parser.add_argument("-c", "--cmap", default='viridis', help="Colour map to use")

args = parser.parse_args()

print("Plotting χ_i(w)")

filenames = ["{}/chi_i_{:02d}.dat".format(args.dir, i + 1) for i in range(args.n)]
mu = np.loadtxt("{}/mu_exciton.out".format(args.dir))
mu_sq = np.array([mu[i, 0]**2 + mu[i, 1]**2 + mu[i, 2]**2 for i in range(args.n)]).flatten()
eigvals = np.loadtxt("{}/eigvals.out".format(args.dir))
wn = (np.loadtxt(filenames[0])[:, 0].flatten())
# wl = 10000000/wn
chi = np.array([np.column_stack(np.loadtxt(fn)[:, 1]).flatten() for fn in filenames])
chi = np.column_stack((wn, chi.T))

fc = cm.viridis(np.linspace(0, 1, args.n))
fig, ax = plt.subplots(figsize=(10,6))
# ax.set_xlim([580, 750])
ax.set_xlim([14000, 17000])
plt.xlabel(r'Wavenumber $ (\text{cm}^{-1}) $')
plt.ylabel(r'Intensity (abu)')
for i in range(args.n):
    label = r'$ \chi_{' + "{}".format(i + 1) + r'}(\omega) \; |\mu|^2 = ' + "{:6.3f}".format(mu_sq[i]) + r' \; \epsilon = ' + "{:5d}".format(int(eigvals[i])) + r' $'
    plt.plot(chi[:, 0], chi[:, i + 1] * mu_sq[i], color=fc[i], label=label)
    ax.fill(chi[:, 0], chi[:, i + 1] * mu_sq[i], color=fc[i], alpha=0.25)

plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("{}/chi_w.pdf".format(args.dir))
plt.close()

# now for chi_bar which makes up F(w)
filenames = ["{}/chi_bar_i_{:02d}.dat".format(args.dir, i + 1) for i in range(args.n)]
chi = np.array([np.column_stack(np.loadtxt(fn)[:, 1]).flatten() for fn in filenames])
chi = np.column_stack((wn, chi.T))
fig, ax = plt.subplots(figsize=(10,6))
# ax.set_xlim([580, 750])
ax.set_xlim([14000, 17000])
plt.xlabel(r'Wavenumber $ (\text{cm}^{-1}) $')
plt.ylabel(r'Intensity (abu)')
for i in range(args.n):
    label = r'$ \bar{\chi}_{' + "{}".format(i + 1) + r'}(\omega) \; |\mu|^2 = ' + "{:6.3f}".format(mu_sq[i]) + r' \; \epsilon = ' + "{:5d}".format(int(eigvals[i])) + r' $'
    plt.plot(chi[:, 0], chi[:, i + 1] * mu_sq[i], color=fc[i], label=label)
    ax.fill(chi[:, 0], chi[:, i + 1] * mu_sq[i], color=fc[i], alpha=0.25)

plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("{}/chi_bar_w.pdf".format(args.dir))
plt.close()

fig, ax = plt.subplots(figsize=(10,6))
ax.set_ylim([0, 0.25])
pop_full = np.loadtxt("{}/populations.dat".format(args.dir))
pop = pop_full
dec = np.argmin(np.abs(np.sum(pop[:, 1:], axis=1) - (1./np.e)))
print(pop[dec])
print(np.abs(np.sum(pop[dec, 1:]) - (1./np.e)))
print("Fluorescence lifetime is {:8.3f} ps.".format(pop[dec, 0]))
plt.xlabel(r'Time (ps) ')
plt.ylabel(r'Population (rel. units)')
plt.title('Populations vs t - fluorescence lifetime = {:8.3f} ps'.format(pop[dec, 0]))
for i in range(args.n):
    plt.plot(pop[:, 0], pop[:, i + 1], color=fc[i])
    # ax.fill(pop[:, 0], pop[:, i + 1], color=fc[i], alpha=0.25)

plt.tight_layout()
plt.savefig("{}/populations.pdf".format(args.dir))
plt.close()
