#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Plot A(w)")
parser.add_argument("-f", "--frame", default=1,
                    help="MD frame to calculate for - pass 0 to \
                    loop over all frames")

args = parser.parse_args()

print("Plotting A(w) and F(w)")

aw_data = np.loadtxt("out/LHCII/{}/aw.dat".format(args.frame))
fw_data = np.loadtxt("out/LHCII/{}/fw.dat".format(args.frame))

aw_exp = np.loadtxt('out/LHCII/aw_exp.dat', skiprows=1)
fw_exp = np.loadtxt('out/LHCII/fw_exp.dat', skiprows=1)

'''
the next two lines switch from wavenumbers to a wavelength
in nanometres. NB: this is a hacky way of doing it!
I deliberately delete the first row, since it contains the
point at zero wavenumber. This is because it'll throw a
divide-by-zero RuntimeWarning otherwise. this is fine for our
purposes here, because we're only showing from ~550-750nm anyway
'''
aw_data = np.delete(aw_data, 0, 0)
fw_data = np.delete(fw_data, 0, 0)
aw_data[:, 0] = 10000000/aw_data[:, 0]
fw_data[:, 0] = 10000000/fw_data[:, 0]
print("Max of A(w):     {}".format(aw_data[np.argmax(aw_data[:, 1])]))
print("Max of A(w) exp: {}".format(aw_exp[np.argmax(aw_exp[:, 1])]))
plt.plot(aw_data[:, 0], aw_data[:, 1]/max(aw_data[:, 1]),
         label=r'$ A(\omega) $')
plt.plot(aw_exp[:, 0], aw_exp[:, 1]/max(aw_exp[:, 1]),
         label=r'Mancal $ Q_y $')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'$ A(\omega) $ (abu)')

ax = plt.gca()
ax.set_xlim([550, 750])
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("out/LHCII/{}/aw.pdf".format(args.frame))
plt.close()


plt.plot(fw_data[:, 0], fw_data[:, 1]/max(fw_data[:, 1]),
         label=r'$ F(\omega) $')
plt.plot(fw_exp[:, 0], fw_exp[:, 1]/max(fw_exp[:, 1]),
         label='Mancal')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'$ F(\omega) $ (abu)')
print("Max of F(w):     {}".format(fw_data[np.argmax(fw_data[:, 1])]))
print("Max of F(w) exp: {}".format(fw_exp[np.argmax(fw_exp[:, 1])]))

ax = plt.gca()
ax.set_xlim([600, 850])
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("out/LHCII/{}/fw.pdf".format(args.frame))

fig, ax = plt.subplots()
ax.set_xlim([550, 850])
plt.plot(aw_data[:, 0], aw_data[:, 1],
         label=r'$ A(\omega) $')
plt.plot(fw_data[:, 0], fw_data[:, 1],
         label=r'$ F(\omega) $')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Intensity (abu)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("out/LHCII/{}/aw_fw.pdf".format(args.frame))
