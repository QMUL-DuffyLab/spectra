#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

parser = argparse.ArgumentParser()
parser.add_argument("-n", type=int, default=15, help="args.number of pigments")
# parser.add_argument("-o", "--output_file", help="Output filename")
parser.add_argument("-d", "--dir", help="Root directory of data")
parser.add_argument("-c", "--cmap", default='viridis', help="Colour map to use")

args = parser.parse_args()

fc = cm.viridis(np.linspace(0, 1, args.n))
fig, ax = plt.subplots(figsize=(10,6))
# ax.set_xlim([580, 750])
ax.set_xlim([14000, 17000])

fig, ax = plt.subplots(3, sharex=True, figsize=(12,12))
ax[0].set_ylim([0, 1.0])
ax[1].set_ylim([0, 0.25])
ax[2].set_ylim([0, 0.55])
ax[0].set_xlim([0, 1000])
ax[1].set_xlim([0, 1000])
ax[2].set_xlim([0, 1000])
pop = np.loadtxt("{}/populations.dat".format(args.dir))
gs_pop = np.loadtxt("{}/gs_populations.dat".format(args.dir))
plt.xlabel(r'Time (ps) ')
ax[0].set_ylabel(r'Population (rel. units)')
ax[1].set_ylabel(r'Population (rel. units)')
ax[2].set_ylabel(r'Population (rel. units)')
ax[0].set_title(r'Ground states')
ax[1].set_title(r'Chlorophylls')
ax[2].set_title(r'Carotenoids')

# redfield_gs = np.zeros(len(pop[:, 0]))
# gs_620 = np.zeros(len(pop[:, 0]))
# gs_621 = np.zeros(len(pop[:, 0]))
# for i in range(len(pop[:, 0])):
#     redfield_gs[i] = redfield_gs[i] + pop[i, 1]
#     for j in range(16):
#         gs_620[i] = gs_620[i] + pop[i, j + 16]
#         gs_621[i] = gs_621[i] + pop[i, j + 64]

# ax[0].plot(pop[:, 0], redfield_gs, color='C0', label='Chls')
# ax[0].plot(pop[:, 0], gs_620, color='C1', label='620')
# ax[0].plot(pop[:, 0], gs_621, color='C2', label='621')
# ax[0].plot(pop[:, 0], 1 - (redfield_gs + gs_620 + gs_621), color='C3', label='Total excited state pop')
ax[0].plot(gs_pop[:, 0], 1 - gs_pop[:, 1], color='C4', label='Total excited state population')
ax[0].axhline(y=1/np.e, ls='--', color='k', label=r'$ 1 / e $')

for i in range(15):
    ax[1].plot(pop[:, 0], pop[:, i + 1], color=fc[i])
    # ax[0].fill(pop[:, 0], pop[:, i + 1], color=fc[i], alpha=0.25)

ax[2].plot(pop[:, 0], pop[:, 16], label=r'620 $S0_{00}$')
ax[2].plot(pop[:, 0], pop[:, 64], label=r'621 $S0_{00}$')
ax[2].plot(pop[:, 0], pop[:, 32], label=r'620 $S1_{00}$')
ax[2].plot(pop[:, 0], pop[:, 80], label=r'621 $S1_{00}$')
ax[2].plot(pop[:, 0], pop[:, 33], label=r'620 $S1_{10}$')
ax[2].plot(pop[:, 0], pop[:, 81], label=r'621 $S1_{10}$')

# print(pop[np.argmin(abs(chl_sum - 1/np.e)), 0])
lifetime = pop[np.argmin(abs((gs_pop[:, 1]) - (1 - 1/np.e))), 0]
with open("{}/lifetime.dat".format(args.dir), 'w') as f:
    f.write("{:16.8e}".format(lifetime))

fig.suptitle("Populations - lifetime = {}".format(lifetime))
ax[0].legend(fontsize=16)
# ax[1].legend(fontsize=16)
ax[2].legend(fontsize=16)
# plt.tight_layout()
plt.savefig("{}/populations.pdf".format(args.dir))
plt.close()
