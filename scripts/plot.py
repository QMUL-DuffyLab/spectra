#!/usr/bin/env python3
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plot_aw_fw(path):
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
    draw_maximums = (False)
    aw = np.loadtxt("{}/aw_average.dat".format(path))
    fw = np.loadtxt("{}/fw_average.dat".format(path))
    aw = np.delete(aw, 0, 0)
    fw = np.delete(fw, 0, 0)
    aw[:, 0] = 10000000/aw[:, 0]
    fw[:, 0] = 10000000/fw[:, 0]
    # Plot A(\omega)
    aw_max = aw[np.argmax(aw[:, 1])]
    print("Max of A(w):     {}".format(aw_max))

    plt.plot(aw[:, 0], aw[:, 1]/max(aw[:, 1]), label=r'$ A(\omega) $')
    if (draw_maximums is True):
        plt.axvline(x=aw_max[0], linestyle='dashed', color='k',
                    label=r'max = $ {:8.3f} $'.format(aw_max[0]))

    plt.xlabel(r'Wavelength (nm)')
    plt.ylabel(r'$ A(\omega) $ (abu)')
    ax = plt.gca()
    ax.set_xlim([600, 700])
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("{}/aw.pdf".format(path))
    plt.close()

    # Plot F(\omega)
    fw_max = fw[np.argmax(fw[:, 1])]
    print("Max of F(w):     {}".format(fw_max))
    plt.plot(fw[:, 0], fw[:, 1]/max(fw[:, 1]), label=r'$ F(\omega) $')
    if (draw_maximums is True):
        plt.axvline(x=fw_max[0], linestyle='dashed', color='k', 
                    label=r'max = $ {:8.3f} $'.format(fw_max[0]))

    plt.xlabel(r'Wavelength (nm)')
    plt.ylabel(r'$ F(\omega) $ (abu)')
    ax = plt.gca()
    ax.set_xlim([650, 750])
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("{}/fw.pdf".format(path))

    # combined
    fig, ax = plt.subplots()
    ax.set_xlim([600, 750])
    plt.plot(aw[:, 0], aw[:, 1], label=r'$ A(\omega) $')
    plt.plot(fw[:, 0], fw[:, 1], label=r'$ F(\omega) $')
    if (draw_maximums is True):
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
    plt.close()

def plot_heatmap(data, out, *args, kwargs={}, plot_kwargs={}, cbar_kwargs={}):
    '''
    Easily plot a heatmap with controllable axis labels, tick labels, colour
    bar details and various other things. Note that you have to construct the
    kwargs dictionaries yourself and change/delete keys as needed. This is
    because normally you'd just have one kwargs with all the keyword arguments
    in it, but I wanted to be able to set plot parameters and colour bar
    parameters separately, and you can't pull out chunks of the normal kwargs
    dictionary to pass them to the plot or colour bar as **kwargs.
    '''
    fig, ax = plt.subplots(figsize=(12, 12))
    if 'title' in kwargs:
        plt.title(kwargs['title'])
        
    if 'cmap' in kwargs:
        cm = kwargs['cmap']
    else:
        cm = "RdBu_r"
            
    X, Y = np.meshgrid(np.arange(np.shape(data)[0] + 1), np.arange(np.shape(data)[1] + 1))
    # - 0.5 for X and Y to centre the data on the ticks
    im = ax.pcolormesh(X - 0.5, Y - 0.5, data, cmap=cm, ec="#BBBBBB", **plot_kwargs)
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'])
    
    if 'ticklabels' in kwargs:
        kwargs['x_ticklabels'] = kwargs['ticklabels']
        kwargs['y_ticklabels'] = kwargs['ticklabels']
    
    if 'y_ticklabels' in kwargs:
        ax.set_yticks(np.arange(len(kwargs['y_ticklabels'])))
        ax.set_yticklabels(kwargs['y_ticklabels'])
        
    if 'x_ticklabels' in kwargs:
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.set_xticks(np.arange(len(kwargs['x_ticklabels'])))
        ax.set_xticklabels(kwargs['x_ticklabels'])
        plt.setp(ax.get_xticklabels(), rotation=90, ha="left", va="center", rotation_mode="anchor")
        
    cax = plt.axes([0.125, 0.05, 0.775, 0.055])
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal', **cbar_kwargs)
            
    if 'cbar_label' in kwargs:
        cbar.set_label(kwargs['cbar_label']) 
        
    if 'cbar_ticklabels' in kwargs:
        cbar.ax.set_xticklabels(kwargs['cbar_ticklabels'])
        
    fig.savefig(out)
    plt.close()

def plot_all(path):
    plot_aw_fw(path)

    figsize=(10,8)
    sns.set(font_scale=2) 
    pigments = [r'Chl \textit{b} 601', r'Chl \textit{a} 602', r'Chl \textit{a} 603', r'Chl \textit{a} 604', r'Chl \textit{b} 605', r'Chl \textit{b} 606', r'Chl \textit{b} 607', r'Chl \textit{b} 608', r'Chl \textit{b} 609', r'Chl \textit{a} 610', r'Chl \textit{a} 611', r'Chl \textit{a} 612', r'Chl \textit{a} 613', r'Chl \textit{a} 614', r'Lut 620', r'Lut 621']
    valfmt = "{:3.1f}"
    valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    taus = np.loadtxt("{}/{}".format(path, "tau.dat"))
    tau_avg = np.loadtxt("{}/tau_average.dat".format(path))
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(r'$ \left< \tau \right> = {:8.3f} $ ps, $ \sigma = {:8.3f} $ ps'.format(tau_avg[0], tau_avg[1]))
    plt.xlabel("Frame")
    plt.ylabel(r'Lifetime $ \tau $ (ps)')
    plt.plot(taus[:, 0], taus[:, 1])
    plt.tight_layout()
    plt.savefig("{}/tau.pdf".format(path))
    plt.close()

    ''' lifetime histogram'''
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(taus[:, 1], edgecolor='k', bins=50, histtype='step')
    ax.set_xlabel(r' Lifetime $ \tau $')
    ax.set_ylabel(r'Counts')
    ax.set_title(r'$ \left< \tau \right> = {:8.3f} $ ps, $ \sigma = {:8.3f} $ ps'.format(tau_avg[0], tau_avg[1]))
    plt.tight_layout()
    plt.savefig("{}/tau_hist.pdf".format(path))
    plt.close()

    '''620-612 coupling'''
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(np.loadtxt("{}/j620cl.dat".format(path))[:, 2])
    ax.set_ylabel(r'620-612 coupling $ (cm^{-1}) $')
    ax.set_xlabel(r'Frame')
    plt.tight_layout()
    fig.savefig("{}/j620612.pdf".format(path))
    plt.close()

    '''620-612 coupling vs. lifetime'''
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(np.loadtxt("{}/tau.dat".format(path))[:, 1], np.loadtxt("{}/j620cl.dat".format(path))[:, 2])
    ax.set_ylabel(r'$ J_{620-612} (\text{cm}^{-1}) $')
    ax.set_xlabel(r'$ \tau $ (ps)')
    plt.tight_layout()
    fig.savefig("{}/j620612_tau.pdf".format(path))
    plt.close()

    '''620 squared dipole moment'''
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(np.loadtxt("{}/musq620.dat".format(path)))
    ax.set_xlabel("Frame")
    ax.set_ylabel(r'$ \mu^{2}_{\text{LUT620}} $')
    plt.tight_layout()
    fig.savefig("{}/musq620.pdf".format(path))
    plt.close()

    ''' 
    big plot with participation / coupling data
    this plot is an absolute fucking nightmare!!!
    '''
    avg_part = np.loadtxt("{}/eig_average.dat".format(path))
    exc_620_coupling = np.loadtxt("{}/exc_620_average.dat".format(path))
    eigvals = np.loadtxt("{}/eigvals_average.dat".format(path))
    eigvals_std = np.loadtxt("{}/eigvals_std.dat".format(path))

    edgecolour='#BBBBBB'
    #cm = sns.color_palette("_r", as_cmap=True)
    cm = sns.light_palette("firebrick", as_cmap=True)
    #cm = sns.dark_palette("seagreen", as_cmap=True)

    avg_part = avg_part[:14, :14]
    fig, ax = plt.subplots(figsize=(12, 12))
    sidey=np.arange(15)
    X, Y = np.meshgrid(sidey, sidey)
    # - 0.5 for X and Y to centre the data at the exciton numbers - this would start them at 0 though, so shift the tick labels by 1
    im = ax.pcolormesh(X - 0.5, Y - 0.5, np.transpose(avg_part), cmap=cm, ec=edgecolour)
    ax.set_yticks(np.arange(0, np.shape(avg_part)[0], dtype=int))
    ax.set_yticklabels("{:5.0f}".format(eigval) for eigval in eigvals[:14])
    #ax.set_yticks(np.arange(0, np.shape(avg_part)[0], dtype=int))
    #ax.set_yticklabels([str(i + 1) for i in np.arange(0, np.shape(avg_part)[0], dtype=int)])
    ax.set_xticks(np.arange(len(pigments[:14])))
    ax.set_xlabel("Pigment")
    ax.set_ylabel(r'$ \left< E_{\,\text{\huge{exciton}}} \right> (cm^{-1}) $')
    ax.set_xticklabels(pigments[:14])
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left", va="center", rotation_mode="anchor")
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    # for some reason "location" does not work as a kw arg here. the matplotlib docs say it should. fucking useless shit
    cax = plt.axes([0.125, 0.03, 0.675, 0.055])
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    cbar.set_label("Exciton participation")
    sideax = plt.axes([0.82, 0.1, 0.09, 0.8])
    #sideax.plot(exc_620_coupling[:14])
    sidex=np.linspace(-0.5, 0.5, num=2, endpoint=True)
    sidey=np.arange(15)
    X, Y = np.meshgrid(sidex, sidey)
    sideax.pcolormesh(X, Y, np.transpose(exc_620_coupling[np.newaxis, :14]), cmap=cm, ec=edgecolour)
    sideax.set_xticks([])
    sideax.set_yticks([])
    #sideax.yaxis.set_ticks_position('right')
    sideax.set_title(r'$ \left< J^{2}_{\text{\Large{exc}},\, \text{\Large{Lut 620}}} \right>  $', pad=12.0, fontsize=28)
    for i in range(len(exc_620_coupling[:14])):
        sideax.text(0, i + 0.5, "{:3.1f}".format(exc_620_coupling[i]), ha="center", va="center", color="k")
        
    # plt.tight_layout()
    plt.savefig("{}/exciton_heatmap.pdf".format(path), bbox_inches='tight')
    plt.close()

        
    kwargs = {
        'cmap' : sns.color_palette("flare", as_cmap=True),
        'ticklabels' : pigments,
        'xlabel' : r'Pigment',
        'ylabel' : r'Pigment',
        'cbar_label' : r'RMSD $ (\text{\AA}) $',
    }

    rmsd = np.loadtxt("{}/rmsd_average.dat".format(path))
    rmsd = np.ma.masked_where(rmsd < 0.1, rmsd)
    plot_heatmap(rmsd, "{}/rmsd.pdf".format(path), kwargs=kwargs)

    theta = np.loadtxt("{}/theta_average.dat".format(path))

    plot_kwargs = {
        'vmin' : 0.,
        'vmax' : np.pi,
    }

    cbar_kwargs = {
        'ticks' : [0., np.pi/4., np.pi/2., 3. * np.pi / 4., np.pi],
    }

    kwargs['cmap'] = "RdBu_r"
    kwargs['cbar_ticklabels'] = [r'$ 0 $', r'$ \frac{\pi}{4} $', r'$ \frac{\pi}{2} $', r'$ \frac{3\pi}{4} $', r'$ \pi $']
    kwargs['cbar_label'] = r'$ \theta $'
    plot_heatmap(theta, "{}/theta.pdf".format(path), kwargs=kwargs, plot_kwargs=plot_kwargs, cbar_kwargs=cbar_kwargs)

    jij = np.loadtxt("{}/jij_average.dat".format(path))
    kwargs['cmap'] = "RdBu_r"
    kwargs['cbar_label'] = r'$ J_{ij} (\text{cm}^{-1}) $'
    plot_kwargs['vmin'] = -50.
    plot_kwargs['vmax'] = +50.
    del kwargs['cbar_ticklabels']
    plot_heatmap(np.ma.masked_where(jij > 1000.0, jij), "{}/jij.pdf".format(path), kwargs=kwargs, plot_kwargs=plot_kwargs)

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description="Plot A(w)")
    parser.add_argument("-d", "--dir", default='out/LHCII',
                        help="directory of A(w) and F(w) data to plot")
    parser.add_argument("-f", "--frame", default=1,
                        help="MD frame to calculate for - pass 0 to \
                        loop over all frames")

    args = parser.parse_args()

    plot_all(args.dir)

