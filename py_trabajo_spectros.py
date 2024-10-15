# Imoprting packages
import os
import shutil
import sys

import numpy as np

from scipy.signal import find_peaks
#from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline

import matplotlib.pyplot as plt

from astropy.stats import sigma_clip

from py_parametros_dict import fiting_parameters


'''#-----CODE-----#'''

def plot_spec(wl_values, it_values, path, color="blue", continuum=None, lines=None, points=None):
    """
    Receives the intensity and wavelength data and makes a plot.

    Parameters
    ----------
    wl_values : np.array
        Wavelength values.
        
    it_values : np.array
        Intensity values.
        
    path : str
        Path of the figure that will be produced and saved.
    
    color : str
        Color of the plot.
    
    continuum : 
        Fit to the continuum.
        
    lines : 
        Set of spectral lines to be plotted on top of the spectrum.
    
    points : tuple
        Set of points used in the fit.

    Returns
    -------
    None.

    """
    # Enable LaTeX rendering
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')  # Use a serif font for LaTeX rendering
    # Define the LaTeX preamble with siunitx
    plt.rcParams['text.latex.preamble'] = r'''
                \usepackage{siunitx}
                \sisetup{
                  detect-family,
                  separate-uncertainty=true,
                  output-decimal-marker={.},
                  exponent-product=\cdot,
                  inter-unit-product=\cdot,
                }
                \DeclareSIUnit{\cts}{cts}
                '''
    # Figure is created:
    plt.figure(dpi=700)
    plt.plot(wl_values, it_values, color=color, linewidth=0.2)
    plt.fill_between(wl_values, it_values, 0, color=color, alpha=0.2, edgecolor=None)
    
    if continuum is not None:
        plt.plot(wl_values, continuum, "green", linewidth=0.7)
        
    # Save image as SVG
    plt.savefig(f'{path}.svg', format='svg')
    plt.savefig(f'{path}.pdf', format='pdf') 
    
    
def norm_spec(wl_values, it_values, filter_iterations=3, s=0.01):
    """
    Receives intensity and wavelength data of an spectrum and fits a function to the 
    continuum.
    
    Parameters
    ----------
    wl_values : np.array
        Wavelength values.
        
    it_values : np.array
        Intensity values.
    
    filter_iterations : int
        Number of iterations for the maxima filtering process.
    
    s : float
        Positive smoothing factor used in UnivariateSpline.

    Returns
    -------
    Intensity of the continuum, normalized intensity and the indexes of the maxima used
    to fit the continuum.

    """
    # Resolution of the spectrum:
    resolution = abs(wl_values[1]-wl_values[0])
    
    # Frequency values:
    nu_values = 3e+8/(wl_values*10**-10) # Hz
    
    # When looking for the maxima, a window of a certain size will be used. This size 
    # is chosen considering the typical width of a line in the spectrum:
    typical_width_freq = 10e+9
    typical_width_wl = (typical_width_freq*3e+8/((max(nu_values))**2))*10**10
    # Window size:
    window_size = 600*typical_width_wl/resolution
    
    #Peaks are selected:
    peaks, _ = find_peaks(it_values, distance=window_size)
    
    for i in range(filter_iterations):
        
        # Wavelength and intensity values for the peaks:
        wl_values_peaks = wl_values[peaks]
        it_values_peaks = it_values[peaks]
    
        # A polinomium is fitted:
        fitted_pol = np.polyfit(wl_values_peaks,it_values_peaks,5)
        # Intensity values of the provisional continuum:
        fitted_intensity = np.polyval(fitted_pol,wl_values)
        
        if i<filter_iterations-1:
            # Calculate the difference between maxima and fitted continuum:
            difference = it_values_peaks - fitted_intensity[peaks]
            
            # Asymetric sigma clipping to filter the peaks, removing the ones
            # far from the fitted polynomium
            clipped_differences = sigma_clip(difference, sigma_upper=7, sigma_lower=3)  # Change sigma value as needed
            # Keep only the maxima that pass the sigma clipping
            peaks = peaks[~clipped_differences.mask]
            wl_values_peaks = wl_values[peaks]
            it_values_peaks = it_values[peaks]
    
    # Definitive cubic spline fit:
    cs = UnivariateSpline(wl_values_peaks, it_values_peaks, k=3, s=s)
    # Intensity of the continuum
    continuum_it = cs(wl_values)
    # Normalized intensities:
    normalized_it = it_values/continuum_it
    
    return continuum_it, normalized_it, peaks

#def main(run_version,spectrum_list,i=None,spec_name=None,spec_smooth_dict=None,spec_iter_dict=None):
def main():
    """
    Main function. Starts by looking for available spectra in the working directory.

    Returns
    -------
    None.

    """

    # Spectrum to analyse
    if run_version == '-m':
        # User selection
        spec_num = int(input("### Spectrum to be analysed (number): "))
        
    else:
        # Automatic mode
        spec_num = spect_pos

    spectrum_name = spectrum_list[spec_num]
    spectrum_path = f'{cwd}/{spectrum_name}'
 
    
    print(f"Loading spectrum {spectrum_name}\n")
    
    # A directory where figures will be saved is created:
    dirname = spectrum_name.split('.')[0]
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)
        
    # Loading data:
    data = np.loadtxt(spectrum_path)
    wl_values = data[:,0] # Angstrons
    it_values = data[:,1] # arbitrary units
    
    #########################################################################
    
    # Raw spectrum is plotted:
    plot_spec(wl_values, it_values, os.path.join(dirname, "01_raw_spectrum"))
    
    #########################################################################
    
    # Normalization:
    if run_version == '-m':
        # User selection of fitting parameters
        smoothing_parameter = float(input("### Indicate smoothing parameter for spline fit: "))
        print("")
        iterations = int(input("### Indicate number of filtering iterations: "))
    
    else:
        # Automatic selection  of fitting parameters
        smoothing_parameter = spec_smooth_dict[f'{spectrum_name}']
        iterations = spec_iter_dict[f'{spectrum_name}']
    
    it_continuum, it_normalized, used_maxima = norm_spec(wl_values, it_values, 
                                                         s=smoothing_parameter, 
                                                         filter_iterations=iterations)
    
    # Spectrum with continuum is produced:
    plot_spec(wl_values, it_values, os.path.join(dirname, "02_continuum"), continuum=it_continuum)
    
    # Normalized spectrum is produced:
    plot_spec(wl_values, it_normalized, os.path.join(dirname, "03_normalized_spectrum"))


if __name__ == "__main__":

    # Obtaining the files to work with
    cwd = os.getcwd()
    spectrum_list = []

    for file in sorted(os.listdir(cwd)):
        if '.dat' in file:
            spectrum_list.append(file)
    print("Available spectra:\n")
    for i, spec in enumerate(spectrum_list):
        print(f"{i}: {spec}")
    print("\n")

    # Selecting the automatic mode or the manual mode by arguments
    if len(sys.argv) == 1:

        # By default the running version is the automatic
        sys.argv.append('-a')

    run_version = sys.argv[1]


    if run_version == '-a':

        # Loading fitting parameters dictionaries
        spec_smooth_dict, spec_iter_dict = fiting_parameters()

        # Analyzing each spectrum
        for spect_pos, spec_name in enumerate(spectrum_list):
  
            main()

    else:

        main()


