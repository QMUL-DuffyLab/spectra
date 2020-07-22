#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
from numpy import loadtxt, max

parser = argparse.ArgumentParser(description="Plot A(w)")
parser.add_argument("-f", "--frame", default=1,
                    help="MD frame to calculate for - pass 0 to \
                    loop over all frames")

args = parser.parse_args()

exp_data = loadtxt('out/LHCII/mancal_2016_qy.dat')
aw_data = loadtxt("out/LHCII/{}/aw.dat".format(args.frame))

plt.plot(10000000/aw_data[:, 0], aw_data[:, 1]/max(aw_data[:, 1]),
         label=r'$ A(\omega) $')
plt.plot(exp_data[:, 0], exp_data[:, 1]/max(exp_data[:, 1]),
         label=r'Mancal $ Q_y $')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'$ A(\omega) $ (abu)')

ax = plt.gca()
ax.set_xlim([580, 700])
plt.grid()
plt.legend()
plt.savefig("out/LHCII/{}/aw_new.pdf".format(args.frame))
plt.close()

fw_data = loadtxt("out/LHCII/{}/fw.dat".format(args.frame))
exp_data = loadtxt('out/LHCII/fw_lhcii_mancal.csv')

plt.plot(10000000/fw_data[:, 0], fw_data[:, 1]/max(fw_data[:, 1]),
         label=r'$ F(\omega) $')
plt.plot(exp_data[:, 0], exp_data[:, 1]/max(exp_data[:, 1]),
         label='Mancal')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'$ F(\omega) $ (abu)')

ax = plt.gca()
ax.set_xlim([640, 780])
plt.grid()
plt.legend()
plt.savefig("out/LHCII/{}/fw_new.pdf".format(args.frame))

fig, ax = plt.subplots()
ax.set_xlim([580, 780])
plt.plot(10000000/aw_data[:, 0], aw_data[:, 1],
         label=r'$ A(\omega) $')
plt.plot(10000000/fw_data[:, 0], fw_data[:, 1],
         label=r'$ F(\omega) $')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.grid()
plt.legend()
plt.savefig("out/LHCII/{}/aw_fw.pdf".format(args.frame))
