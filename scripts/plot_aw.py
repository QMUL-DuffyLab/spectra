#!/usr/bin/env python3
from numpy import loadtxt,max;
import matplotlib.pyplot as plt;

data = loadtxt('out/LHCII/1/aw.dat')
exp_data = loadtxt('out/LHCII/mancal_2016_qy.dat')
plt.plot(10000000/data[:,0], data[:,1]/max(data[:,1]), label='LHCII chls')
plt.plot(exp_data[:,0], exp_data[:,1]/max(exp_data[:,1]), label='Mancal $ Q_y $')
plt.xlabel(r'Wavelength (nm)'); plt.ylabel(r'$ A(\omega) $ (abu)')
ax = plt.gca(); ax.set_xlim([580,700]); 
plt.grid(); 
plt.savefig('out/LHCII/1/aw_new.pdf');
