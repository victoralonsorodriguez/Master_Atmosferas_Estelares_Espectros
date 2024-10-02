# Imoprting packages
import os

import numpy as np

from scipy.signal import find_peaks
from scipy.signal import savgol_filter

import matplotlib.pyplot as plt

from astropy.stats import sigma_clip

'''#-----CODE-----#'''

def plot_spec(wl_values, it_values, continuum=None, lines=None):
    """
    Receives the intensity and wavelength data and makes a plot.

    Parameters
    ----------
    wl_values : np.array
        Wavelength values.
        
    it_values : np.array
        Intensity values.
    
    continuum : 
        Fit to the continuum.
        
    lines : 
        Set of spectral lines to be plotted on top of the spectrum.

    Returns
    -------
    None.

    """
    
def norm_spec(wl_values, it_values):
    """
    Receives intensity and wavelength data of an spectrum and fits a function to the 
    continuum.
    
    Parameters
    ----------
    wl_values : np.array
        Wavelength values.
        
    it_values : np.array
        Intensity values.

    Returns
    -------
    Function fitted to the continuum, normalized values of intensity.

    """

def main():
    """
    Main function. Starts by looking for available spectra in the working directory.

    Returns
    -------
    None.

    """
    # Obtaining the files to work with
    cwd = os.getcwd()
    spectrum_list = []

    for file in sorted(os.listdir(cwd)):
        if '.dat'in file:
            spectrum_list.append(file)
    
    for i, spec in enumerate(spectrum_list):
        print(f"{i}: {spec}")
    
    spectrum_path = f'{cwd}/{spectrum_list[2]}'
    
    # Loading data:
    data = np.loadtxt(spectrum_path)
    wl_values = data[:,0] # A
    it_values = data[:,1] # arbitrary units
    nu_values = 3e+8/(wl_values*10**-10) # Hz
    
    resolution = abs(wl_values[1]-wl_values[0])
    
    # Raw spectrum
    plt.figure(dpi=700)
    plt.plot(wl_values, it_values, color="red", linewidth=0.2)
    
    # A typical line width is defined. It can be 10 GHz for example. This will be used
    # for getting the maxima.
    typical_width_freq = 10e+9
    typical_width_wl = (typical_width_freq*3e+8/((max(nu_values))**2))*10**10
    # When looking for the maxima, a window of a certain size will be used. This size 
    # is chosen considering the typical width of a line in the spectrum:
    window_size = 600*typical_width_wl/resolution

    #Maxima are selected:
    peaks, _ = find_peaks(it_values, distance=window_size)
    
    iterations = 2 # This will be an argument of a function
    
    for i in range(iterations):
        
        # Wavelength and intensity values for the peaks:
        wl_values_peaks = wl_values[peaks]
        it_values_peaks = it_values[peaks]
    
        # A polinomium is fitted:
        fitted_pol = np.polyfit(wl_values_peaks,it_values_peaks,5)
        fitted_intensity = np.polyval(fitted_pol,wl_values)
        
        if i<iterations-1:
            # Calculate the difference between maxima and fitted continuum
            difference = it_values_peaks - fitted_intensity[peaks]
            
            # Asymetric sigma clipping to filter the peaks, removing the ones
            # far from the fitted polynomium
            clipped_differences = sigma_clip(difference, sigma_upper=7, sigma_lower=3)  # Change sigma value as needed
            # Keep only the maxima that pass the sigma clipping
            peaks = peaks[~clipped_differences.mask]
            wl_values_peaks = wl_values[peaks]
            it_values_peaks = it_values[peaks]

    plt.plot(wl_values,fitted_intensity, "green")
    #plt.plot(wl_values_peaks, it_values_peaks,"x")

   

