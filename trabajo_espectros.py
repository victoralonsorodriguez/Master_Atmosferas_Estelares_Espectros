
# Imoprting packages
import os

from astropy.io import fits
from astropy.io import ascii

from scipy.signal import find_peaks

import matplotlib.pyplot as plt

import pdb

'''#-----CODE-----#'''

# Obtaining the files to work with
cwd = os.getcwd()
spectrum_list = []

for file in sorted(os.listdir(cwd)):
    if '.dat'in file:
        spectrum_list.append(file)


# Analyzing each spectrum separately
#for i in range(len(spectrum_list)):
for i in range(1):

    # Selecting the spectum
    espectrum_name = spectrum_list[i]
    espectrum_path = f'{cwd}/{espectrum_name}'

    # Obtaining the data of each spectum 
    data = ascii.read(espectrum_path)

    # separating the wavelenght values from 
    # the intensisty values
    wl_values = data['col1']
    int_values = data['col2']


    #-----PEAKS-----#
    # Selecting the upper peaks of the spectum
    # This give us the position of the peak inside
    # the data array
    peaks_upper, _ = find_peaks(int_values)

    # Selecting the lower peaks of the spectum
    peaks_lower, _ = find_peaks(-int_values)

    # Selecting the peaks whose have a small diference
    # between them
    peaks_diff = int_values[peaks_upper] - int_values[peaks_lower]
    
    # Creating a mask for those values
    peaks_mask = peaks_diff<0.01

    peaks_upper_good = peaks_upper[peaks_mask]
    peaks_lower_good = peaks_lower[peaks_mask]

    # Plotting the spectrum
    plt.figure()

    # Plotting the raw spectum
    plt.plot(wl_values,int_values)

    # Plotting the peaks
    plt.plot(wl_values[peaks_upper_good], int_values[peaks_upper_good], "x")
    plt.plot(wl_values[peaks_lower_good], int_values[peaks_lower_good], "x")

    # Showing the plot
    plt.show()