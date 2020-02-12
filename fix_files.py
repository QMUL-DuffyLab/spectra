#!/usr/bin/env python3

import os
import glob
import numpy as np

filelist = glob.glob('LHCII' + '/**/*.csv', recursive=True)

for item in filelist:
    arr = np.loadtxt(item)
    print(item)
    np.savetxt(item, arr, fmt='%016.8e')
