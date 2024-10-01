
# Imoprting packages
import os

import numpy as np

from astropy.io import fits
from astropy.io import ascii

import scipy as scp
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
for i in range(len(spectrum_list)):
#for i in range(1):

    # Selecting the spectum
    espectrum_name = spectrum_list[i]
    espectrum_path = f'{cwd}/{espectrum_name}'

    # Obtaining the data of each spectum 
    data = ascii.read(espectrum_path)

    # separating the wavelenght values from 
    # the intensisty values
    wl_values = np.array(data['col1'])
    int_values = np.array(data['col2'])


    #-----PEAKS-----#
    # Selecting the upper peaks of the spectum
    # This give us the position of the peak inside
    # the data array
    peaks_upper, _ = find_peaks(int_values)

    wl_values_peaks_upper = wl_values[peaks_upper]
    int_values_peaks_upper = int_values[peaks_upper]

    # Selecting the lower peaks of the spectum
    peaks_lower, _ = find_peaks(-int_values)

    # Selecting the peaks whose have a small diference
    # between them
    peaks_diff = int_values[peaks_upper] - int_values[peaks_lower]
    
    # Creating a mask for those peaks wich their 
    # difference is less than dif_value
    diff_value = 0.025
    peaks_mask = peaks_diff < diff_value

    # selecting those good peaks positions
    peaks_upper_good = peaks_upper[peaks_mask]
    peaks_lower_good = peaks_lower[peaks_mask]

    # Wavelenght values for thos peaks
    wl_values_good = wl_values[peaks_upper_good]

    # Wavelenght values for those peaks
    wl_values_upper_good = wl_values[peaks_upper_good]
    wl_values_lower_good = wl_values[peaks_lower_good]

    # Intensity values for those peaks
    int_values_upper_good = int_values[peaks_upper_good]
    int_values_lower_good = int_values[peaks_lower_good]

    # Mean value between the upper and lower peaks
    peaks_good_mean = np.mean( np.array([int_values_upper_good,int_values_lower_good]), axis=0 )

    # Apllying a mean filter to reduce the noise 
    # between peaks
    peaks_good_median = scp.signal.medfilt(peaks_good_mean, kernel_size = 5)

    #-----FITTING-----#

    # Now that points to fit are selected we can fit a 
    # polynomial function to those peaks

    fitting_points_wl = wl_values_upper_good
    
    fitting_points_int = int_values_upper_good



    # Fitting the continuous 
    fitting_function = np.polyfit(fitting_points_wl,fitting_points_int,4)

    # Obtaining the fitting function
    fitting_polynomial = np.poly1d(fitting_function)

    # Creating the fitting data
    fitting_intensity = fitting_polynomial(wl_values)

    # Normalizing the spectrum
    normalize_spectrum = int_values / fitting_intensity

    #-----PLOTTING-----#

    # Plotting the raw spectrum
    # along with the fitting function
    plt.figure(1)

    # Plotting the raw spectum
    plt.plot(wl_values,int_values)

    # Plotting the point selected to fit
    plt.plot(fitting_points_wl, fitting_points_int, "x")

    # Plotting the fitting function
    plt.plot(wl_values,fitting_intensity)



    # Plotting the normalized spectrum
    plt.figure(2)

    # Plotting the raw spectum
    plt.plot(wl_values,normalize_spectrum)
    plt.ylim(bottom=0,top=1.5)

    # Showing the plots
    plt.show()

