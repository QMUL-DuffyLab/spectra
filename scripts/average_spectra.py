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

root_dir = "out/LHCII"
initial_data = np.loadtxt("{}/1/aw.dat".format(root_dir))
aws = np.zeros_like(initial_data)
fws = np.zeros_like(initial_data)

resum = 1
if resum is not 0:
    print("Summing A(w) and F(w) per frame")
    for i in range(1000):
        if (i % 100) is 0:
            print(".", end='', flush=True)

        direc = "{}/{}".format(root_dir, i + 1)
        aws = aws + np.loadtxt("{}/aw.dat".format(direc))
        fws = fws + np.loadtxt("{}/fw.dat".format(direc))

    print("Done.")
    aws = aws / 1000.
    fws = fws / 1000.
else:
    aws = np.loadtxt("{}/aw_average.dat".format(root_dir))
    fws = np.loadtxt("{}/fw_average.dat".format(root_dir))

np.savetxt("{}/aw_average.dat".format(root_dir), aws)
np.savetxt("{}/fw_average.dat".format(root_dir), fws)

# delete 0 wavenumber rows to prevent divide-by-zero warning
aws = np.delete(aws, 0, 0)
fws = np.delete(fws, 0, 0)
# convert to wavelength (nm)
aws[:, 0] = 10000000/aws[:, 0]
fws[:, 0] = 10000000/fws[:, 0]

# plot
aw_exp = np.loadtxt("{}/aw_exp.dat".format(root_dir), skiprows=1)
fw_exp = np.loadtxt("{}/fw_exp.dat".format(root_dir), skiprows=1)

fig, ax = plt.subplots()
ax.set_xlim([580, 700])
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(aws[:, 0], aws[:, 1]/np.max(aws[:, 1]),
         label=r'$ A(\omega) $')
plt.plot(aw_exp[:, 0], aw_exp[:, 1]/np.max(aw_exp[:, 1]),
         label=r'Experiment')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("{}/aw_average.pdf".format(root_dir))

fig, ax = plt.subplots()
ax.set_xlim([640, 780])
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(fws[:, 0], fws[:, 1]/np.max(fws[:, 1]),
         label=r'$ F(\omega) $')
plt.plot(fw_exp[:, 0], fw_exp[:, 1]/np.max(fw_exp[:, 1]),
         label=r'Experiment')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("{}/fw_average.pdf".format(root_dir))

fig, ax = plt.subplots()
ax.set_xlim([580, 780])
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.plot(aws[:, 0], aws[:, 1], label=r'$ A(\omega) $')
plt.plot(fws[:, 0], fws[:, 1], label=r'$ F(\omega) $')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("{}/aw_fw_average.pdf".format(root_dir))
