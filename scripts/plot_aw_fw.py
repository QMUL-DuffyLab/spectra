import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--ligand", help="ligand code")
parser.add_argument("-a", "--aw_file", help="A(w) input file")
parser.add_argument("-f", "--fw_file", help="F(w) input file")
parser.add_argument("-oa", "--output_aw_file", help="Output A(w) file")
parser.add_argument("-of", "--output_fw_file", help="Output F(w) file")
args = parser.parse_args()

if (args.ligand is not None):
    args.aw_file = "out/{}_Aw_file.dat".format(args.ligand)
    args.fw_file = "out/{}_Fw_file.dat".format(args.ligand)
    args.output_aw_file = "out/{}_Aw.pdf".format(args.ligand)
    args.output_fw_file = "out/{}_Fw.pdf".format(args.ligand)
else:
    args.ligand = args.aw_file[4:7]

aw = np.loadtxt(args.aw_file)
aw = aw[aw[:, 0].argsort()]
fw = np.loadtxt(args.fw_file)
fw = fw[fw[:, 0].argsort()]

exp_file = "in/{0}_exp.dat".format(args.ligand)
aw_exp = np.loadtxt(exp_file)

fig, ax = plt.subplots()
plt.grid()
plt.xlabel('$ \omega $')
plt.ylabel(r'$ A(\omega) $')
# plt.xlim((-2000.0, 2000.0))
plt.plot(aw[:, 0] + 2600, (aw[:, 1] - np.min(aw[:, 1])), label="My result", lw=2.5)
plt.plot(aw_exp[:, 0], 1.7*aw_exp[:, 1], label="Exp result", lw=1.5)
plt.legend()
plt.savefig(args.output_aw_file)
plt.close()

fig, ax = plt.subplots()
plt.grid()
plt.xlabel('$ \omega $')
plt.ylabel(r'$ F(\omega) $')
# plt.xlim((-1000.0, 1000.0))
plt.plot(fw[:, 0], fw[:, 1], label="My result")
plt.legend()
plt.savefig(args.output_fw_file)
plt.close()
