
import os
import numpy as np

from scipy.signal import convolve

from astropy.io import fits
from astropy.io import ascii

import matplotlib.pyplot as plt


cwd = os.getcwd()
spectrum_list = []

for file in sorted(os.listdir(cwd)):
    if '.dat'in file:
        spectrum_list.append(file)

#for i in range(len(spectrum_list)):
for i in range(1):

    espectrum_name = spectrum_list[i]
    espectrum_path = f'{cwd}/{espectrum_name}'


    data = ascii.read(espectrum_path)

    wl_values = data['col1']
    int_values = data['col2']

    plt.figure()
    plt.plot(wl_values,int_values)
    plt.show()