#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
For LHCII data only atm - average the A(w) and F(w) data over every snapshot,
write them out, plot them.
"""

import os
import numpy as np
import argparse
import subprocess
import pigment_data

root_dir = "out/LHCII"
initial_data = np.loadtxt("{}/1/aw.dat".format(root_dir))
aws = np.zeros_like(initial_data)
fws = np.zeros_like(initial_data)

for i in range(1000):
    direc = "{}/{}".format(root_dir, i + 1)
    aws = aws + np.loadtxt("{}/aw.dat".format(direc))
    fws = fws + np.loadtxt("{}/fw.dat".format(direc))

aws = aws / 1000.
fws = fws / 1000.

np.savetxt("{}/aw_average.dat".format(root_dir), aws)
np.savetxt("{}/fw_average.dat".format(root_dir), fws)

# delete 0 wavenumber rows to prevent divide-by-zero warning
aws = np.delete(aws, 0, 0)
fws = np.delete(fws, 0, 0)
# convert to wavelength (nm)
aws[:, 0] = 1000000/aws[:, 0]
fws[:, 0] = 1000000/fws[:, 0]
# plot
aw_exp = np.loadtxt("{}/aw_exp.dat", skiprows=1)
fw_exp = np.loadtxt("{}/fw_exp.dat", skiprows=1)

fig, ax = plt.subplots()
ax.set_xlim([580, 700])
plt.grid()
plt.legend()
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(aw[:, 0], aw[:, 1], label=r'$ A(\omega) $')
plt.plot(aw_exp[:, 0], aw_exp[:, 1], label=r'Experiment')
plt.savefig("{}/aw_average.pdf".format(root_dir))

fig, ax = plt.subplots()
ax.set_xlim([640, 780])
plt.grid()
plt.legend()
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(fw[:, 0], fw[:, 1], label=r'$ F(\omega) $')
plt.plot(fw_exp[:, 0], fw_exp[:, 1], label=r'Experiment')
plt.savefig("{}/fw_average.pdf".format(root_dir))

fig, ax = plt.subplots()
ax.set_xlim([580, 780])
plt.grid()
plt.legend()
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(aw[:, 0], aw[:, 1], label=r'$ A(\omega) $')
plt.plot(fw[:, 0], fw[:, 1], label=r'$ F(\omega) $')
plt.savefig("{}/aw_fw_average.pdf".format(root_dir))
