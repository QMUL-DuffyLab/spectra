import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--aw_file", help="A(w) input file")
parser.add_argument("-f", "--fw_file", help="F(w) input file")
parser.add_argument("-oa", "--output_aw_file", help="Output A(w) file")
parser.add_argument("-of", "--output_fw_file", help="Output F(w) file")
args = parser.parse_args()

aw = np.loadtxt(args.aw_file)
aw = aw[aw[:, 0].argsort()]
fw = np.loadtxt(args.fw_file)
fw = fw[fw[:, 0].argsort()]

# testing for comparison
aw_chris = np.loadtxt("/Users/cgray/Downloads/Duffy_Fitting_Code/A_spec_theor_CAR_2MODE.txt")

fig, ax = plt.subplots()
plt.grid()
plt.xlabel('$ \omega $')
plt.ylabel(r'$ A(\omega) $')
# plt.xlim((-1000.0, 1000.0))
plt.plot(aw[:, 0], aw[:, 1], label="My result")
plt.plot(aw_chris[:, 0], aw_chris[:, 1], label="Chris result")
plt.legend()
plt.savefig(args.output_aw_file)
plt.close()

fig, ax = plt.subplots()
plt.grid()
plt.xlabel('$ \omega $')
plt.ylabel(r'$ F(\omega) $')
plt.xlim((-1000.0, 1000.0))
plt.plot(fw[:, 0], fw[:, 1], label="My result")
plt.legend()
plt.savefig(args.output_fw_file)
plt.close()
