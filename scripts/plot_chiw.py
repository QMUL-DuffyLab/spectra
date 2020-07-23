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
parser.add_argument("-c", "--cmap", default='viridis', help="Colour map to use")
parser.add_argument("-x", "--x_label",
                    help="String for x-axis")
parser.add_argument("-y", "--y_label",
                    help="String for y-axis")
parser.add_argument("-z", "--z_label",
                    help="String for z-axis")

args = parser.parse_args()

N = 14 # make this a parameter later i guess lol
root_dir = "out/LHCII"
frame = 1
filenames = ["{}/{}/chi_i_{:02d}.dat".format(root_dir, frame, i + 1) for i in range(N)]
mu = np.loadtxt("{}/{}/mu_exciton.out".format(root_dir, frame))
mu_sq = np.array([mu[i, 0]**2 + mu[i, 1]**2 + mu[i, 2]**2 for i in range(N)]).flatten()
eigvals = np.loadtxt("{}/{}/eigvals.out".format(root_dir, frame))
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

plt.legend()
plt.savefig("{}/{}/chi_w.pdf".format(root_dir, frame))
plt.close()

# now for chi_bar which makes up F(w)
filenames = ["{}/{}/chi_bar_i_{:02d}.dat".format(root_dir, frame, i + 1) for i in range(N)]
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

plt.legend()
plt.savefig("{}/{}/chi_bar_w.pdf".format(root_dir, frame))
plt.close()
