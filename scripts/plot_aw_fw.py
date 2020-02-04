import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

def plot(filename, *args, **kwargs):
    print(kwargs)
    fig, ax = plt.subplots()
    plt.grid()
    if 'xlabel' in kwargs:
        plt.xlabel(kwargs['xlabel'])

    if 'ylabel' in kwargs:
        plt.ylabel(kwargs['ylabel'])

    for i in range(0, len(args) - 1, 2):
        if 'label' in kwargs:
            plt.plot(args[i], args[i + 1], label=kwargs['label'][i//2])
        else:
            plt.plot(args[i], args[i + 1])

    fig.tight_layout()
    if 'label' in kwargs:
        plt.legend()

    plt.savefig(filename, bbox_inches='tight')
    plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--ligand", help="ligand code")
parser.add_argument("-aw", "--aw_file", help="A(w) input file")
parser.add_argument("-fw", "--fw_file", help="F(w) input file")
parser.add_argument("-at", "--at_file", help="A(w) input file")
parser.add_argument("-ft", "--ft_file", help="F(w) input file")
parser.add_argument("-ag", "--aw_graph", help="Output A(w) file")
parser.add_argument("-fg", "--fw_graph", help="Output F(w) file")
args = parser.parse_args()

if (args.ligand is not None):
    args.aw_file = "out/{}_Aw.dat".format(args.ligand)
    args.fw_file = "out/{}_Fw.dat".format(args.ligand)
    args.at_file = "out/{}_At.dat".format(args.ligand)
    args.ft_file = "out/{}_Ft.dat".format(args.ligand)
    args.aw_graph = "out/{}_Aw.pdf".format(args.ligand)
    args.fw_graph = "out/{}_Fw.pdf".format(args.ligand)
else:
    args.ligand = args.aw_file[4:7]

aw = np.loadtxt(args.aw_file)
aw = aw[aw[:, 0].argsort()]
fw = np.loadtxt(args.fw_file)
fw = fw[fw[:, 0].argsort()]

exp_file = "in/{0}_exp.dat".format(args.ligand)
aw_exp = np.loadtxt(exp_file)

plot(args.output_aw_file, aw[:, 0], (aw[:, 1] - np.min(aw[:, 1])),
        aw_exp[:, 0], (aw_exp[:, 1] - np.min(aw_exp[:, 1])),
        xlabel=r'$ \omega $', ylabel=r'$ A(\omega) $',
        label=[r'My result', r'Exp. result'])
plot(args.output_fw_file, fw[:, 0], (fw[:, 1] - np.min(fw[:, 1])),
        xlabel=r'$ \omega $', ylabel=r'$ F(\omega) $',
        label=[r'My result', r'Exp. result'])
