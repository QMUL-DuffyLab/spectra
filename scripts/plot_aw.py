#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np

def plot_aw_fw(aw, fw, aw_exp, fw_exp, draw_maximums, path):
    '''
    Take simulated absorption and fluorescence spectra (aw and fw)
    along with experimental equivalents if desired; output plots of
    A(w), F(w) and the pair together. draw_maximums should be a 
    list or tuple of logicals which are used to determine whether to
    draw in vertical lines at the maximum in the A(w)/F(w)/combined plots.
    There's definitely a neater/more pythonic way of doing this but ¯\_(ツ)_/¯
    '''
    '''
    the first four lines lines switch from wavenumbers to a wavelength
    in nanometres. NB: this is a hacky way of doing it!
    I deliberately delete the first row, since it contains the
    point at zero wavenumber. This is because it'll throw a
    divide-by-zero RuntimeWarning otherwise. this is fine for our
    purposes here, because we're only showing from ~550-750nm anyway
    '''
    aw = np.delete(aw, 0, 0)
    fw = np.delete(fw, 0, 0)
    aw[:, 0] = 10000000/aw[:, 0]
    fw[:, 0] = 10000000/fw[:, 0]
    # Plot A(\omega)
    aw_max = aw[np.argmax(aw[:, 1])]
    print("Max of A(w):     {}".format(aw_max))
    print("Max of A(w) exp: {}".format(aw_exp[np.argmax(aw_exp[:, 1])]))

    plt.plot(aw[:, 0], aw[:, 1]/max(aw[:, 1]), label=r'$ A(\omega) $')
    if (draw_maximums[0] is True):
        plt.axvline(x=aw_max[0], linestyle='dashed', color='k',
                    label=r'max = $ {:8.3f} $'.format(aw_max[0]))

    # {a/f}w_exp are set to zero below if there's no experimental data
    if (np.any(aw_exp)):
        plt.plot(aw_exp[:, 0], aw_exp[:, 1]/max(aw_exp[:, 1]),
                 label=r'Mancal $ Q_y $')

    plt.xlabel(r'Wavelength (nm)')
    plt.ylabel(r'$ A(\omega) $ (abu)')
    ax = plt.gca()
    ax.set_xlim([550, 750])
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("{}/aw.pdf".format(path))
    plt.close()

    # Plot F(\omega)
    fw_max = fw[np.argmax(fw[:, 1])]
    print("Max of F(w):     {}".format(fw_max))
    print("Max of F(w) exp: {}".format(fw_exp[np.argmax(fw_exp[:, 1])]))
    plt.plot(fw[:, 0], fw[:, 1]/max(fw[:, 1]), label=r'$ F(\omega) $')
    if (draw_maximums[1] is True):
        plt.axvline(x=fw_max[0], linestyle='dashed', color='k', 
                    label=r'max = $ {:8.3f} $'.format(fw_max[0]))

    if (np.any(fw_exp)):
        plt.plot(fw_exp[:, 0], fw_exp[:, 1]/max(fw_exp[:, 1]),
                 label='Mancal')

    plt.xlabel(r'Wavelength (nm)')
    plt.ylabel(r'$ F(\omega) $ (abu)')
    ax = plt.gca()
    ax.set_xlim([550, 850])
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("{}/fw.pdf".format(path))

    # combined
    fig, ax = plt.subplots()
    ax.set_xlim([550, 850])
    plt.plot(aw[:, 0], aw[:, 1], label=r'$ A(\omega) $')
    plt.plot(fw[:, 0], fw[:, 1], label=r'$ F(\omega) $')
    if (draw_maximums[2] is True):
        plt.axvline(x=aw_max[0], linestyle='dashed', color='k',
                    label=r' $ A(\omega) $ max = $ {:8.3f} $'.format(aw_max[0]))
        plt.axvline(x=fw_max[0], linestyle='dashed', color='k',
                    label=r' $ F(\omega) $ max = $ {:8.3f} $'.format(fw_max[0]))

    plt.xlabel(r'Wavelength (nm)')
    plt.ylabel(r'Intensity (abu)')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("{}/aw_fw.pdf".format(path))

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description="Plot A(w)")
    parser.add_argument("-d", "--dir", default='out/LHCII',
                        help="directory of A(w) and F(w) data to plot")
    parser.add_argument("-f", "--frame", default=1,
                        help="MD frame to calculate for - pass 0 to \
                        loop over all frames")

    args = parser.parse_args()

    print("Plotting A(w) and F(w)")

    aw_data = np.loadtxt("{}/aw.dat".format(args.dir))
    fw_data = np.loadtxt("{}/fw.dat".format(args.dir))

    if 'LHCII' in args.dir:
        aw_exp = np.loadtxt('out/LHCII/aw_exp.dat', skiprows=1)
        fw_exp = np.loadtxt('out/LHCII/fw_exp.dat', skiprows=1)
    else:
        aw_exp = np.zeros_like(aw_data)
        fw_exp = np.zeros_like(aw_data)

    draw_maximums = (True, True, True)
    plot_aw_fw(aw_data, fw_data, aw_exp, fw_exp, draw_maximums, args.dir)
