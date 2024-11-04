# Imoprting packages
import os
import shutil
import sys

import pdb

import numpy as np
import pandas as pd

from scipy.signal import find_peaks
#from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline

import matplotlib.pyplot as plt

from itertools import cycle

from astropy.stats import sigma_clip

from py_dict_fitting_params import fiting_parameters
from py_dict_fitting_params import tipos_estrella
from py_dict_wavelenght import wl_lines


'''#-----CODE-----#'''

def plot_spec(wl_values, it_values, path, color="black", continuum=None, lines=None, points=None, peaks=[],lines_dict=None, hline=[]):
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
    plt.figure(figsize=(8,3))
    plt.plot(wl_values, it_values, color=color, linewidth=0.3)
    #plt.fill_between(wl_values, it_values, 0, color=color, alpha=0.2, edgecolor=None)
    
    if continuum is not None:
        plt.plot(wl_values, continuum, "green", linewidth=0.7)

    #wl_lines_keys = wl_lines_dict.keys()

    

    if lines_dict is not None:
        
        lines_dict_plot = lines_dict
    
    else:
        lines_dict_plot = wl_lines_dict
        
        
    
    lines_dict_plot_keys = lines_dict_plot.keys()
        
    colours = ['indianred','red','orange','gold','greenyellow','lime','green','cyan',
               'dodgerblue','navy','darkviolet','purple','fuchsia','hotpink','crimson', 'red']

    plt.savefig(f'{path}.pdf', dpi=300, format='pdf') 

    if 'normalized' in path:
        for key_pos,key in enumerate(lines_dict_plot_keys):
            
            plt.vlines(x=lines_dict_plot[key], ymin=min(it_values), ymax=max(it_values),color=colours[key_pos],label=key,linewidth=0.7)
            
            for wavelength in lines_dict_plot[key]:
                plt.annotate(key, xy=(wavelength, 1), xycoords=("data", "axes fraction"),va='bottom', ha='center', color=colours[key_pos])
                plt.axvline(x=wavelength, color=colours[key_pos],linewidth=0.7)

        if len(hline) != 0:
            for line in hline:
                plt.axhline(y = line, color = 'r', linestyle = '-',linewidth=0.5) 

        if len(peaks) != 0:
            wl_values_peaks = wl_values[peaks]
            it_values_peaks = it_values[peaks]

            plt.scatter(wl_values_peaks,it_values_peaks,s=10,marker='.')

        plt.legend(loc='lower right')    
        plt.savefig(f'{path}.pdf', dpi=300, format='pdf') 
        
    plt.show()

    # Save image as SVG
    #plt.savefig(f'{path}.svg', format='svg')
    
    #plt.close()

def spec_tower(norm_specs, path, color="black", names=None, annotation_wl=None,lines_dict=None, factor=3):
    """
    Receives several normalized spectra and plots them on top of each other

    Parameters
    ----------
    norm_specs: list
        List of tuples, each one corresponding to a spectrum (wl, intensity)
        
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
    plt.figure(figsize=(8,3))
    for i,spec in enumerate(norm_specs):
        plt.plot(spec[0], spec[1]+i/factor, color='blue' if i==0 else color, linewidth=0.6 if i==0 else 0.3)
        
        if names is not None:
            plt.annotate(names[i], xy=(annotation_wl , 1+i/factor+0.1+1/(11*factor)), 
                 fontsize=12, fontweight='bold', color='blue' if i==0 else color,
                 ha='left', va='top', backgroundcolor="white", bbox=None)
    plt.xlim(min(norm_specs[1][0]), max(norm_specs[1][0]))
    #plt.ylim(-1,6)

    

    if lines_dict is not None:
        
        lines_dict_plot = lines_dict
    
    else:
        lines_dict_plot = wl_lines_dict
        
        
    
    lines_dict_plot_keys = lines_dict_plot.keys()
        
    colours = ['indianred','red','orange','yellow','lime','green','cyan',
               'dodgerblue','navy','darkviolet','purple','fuchsia','hotpink','crimson', 'red']

    plt.savefig(f'{path}.pdf', dpi=300, format='pdf') 

    if 'normalized' in path or 'referencia' in path:
        for key_pos,key in enumerate(lines_dict_plot_keys):
            
            for wavelength in lines_dict_plot[key]:
                plt.annotate(key, xy=(wavelength, 1), xycoords=("data", "axes fraction"),va='bottom', ha='center', color=colours[key_pos])
                plt.axvline(x=wavelength, color=colours[key_pos],linewidth=0.7)


        plt.legend(loc='lower right')    
        plt.savefig(f'{path}.pdf', dpi=300, format='pdf') 
        
    plt.show()

    # Save image as SVG
    #plt.savefig(f'{path}.svg', format='svg')
    
    #plt.close()

def finding_peaks(wl_values, it_values, window_size_factor=450):

    '''
    Receives intensity and wavelength data of an spectrum and searches
    for the peaks

    '''

    # Resolution of the spectrum:
    resolution = abs(wl_values[1]-wl_values[0])
    
    # Frequency values:
    nu_values = 3e+8/(wl_values*10**-10) # Hz
    
    # When looking for the maxima, a window of a certain size will be used. This size 
    # is chosen considering the typical width of a line in the spectrum:
    typical_width_freq = 10e+9
    typical_width_wl = (typical_width_freq*3e+8/((max(nu_values))**2))*10**10
    # Window size:
    window_size = window_size_factor*typical_width_wl/resolution
    
    #Peaks are selected:
    peaks, _ = find_peaks(it_values, distance=window_size)

    return peaks

    
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
    # Finding the peaks
    peaks = finding_peaks(wl_values, it_values)

    # Fitting the continuum
    for i in range(filter_iterations):
        
        # Wavelength and intensity values for the peaks:
        wl_values_peaks = wl_values[peaks]
        it_values_peaks = it_values[peaks]
    
        # A polinomium is fitted:
        fitted_pol = np.polyfit(wl_values_peaks,it_values_peaks,4)
        # Intensity values of the provisional continuum:
        fitted_intensity = np.polyval(fitted_pol,wl_values)
        
        if i<filter_iterations-1:
            # Calculate the difference between maxima and fitted continuum:
            difference = it_values_peaks - fitted_intensity[peaks]
            
            # Asymetric sigma clipping to filter the peaks, removing the ones
            # far from the fitted polynomium
            clipped_differences = sigma_clip(difference, sigma_upper=3, sigma_lower=3)  # Change sigma value as needed
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



def radialv_correct(wl_values, it_values, rv_threshold):
    """
    Receives a spectrum and uses the spectral H and He lines to calculate the
    radial velocity and generate a corrected spectrum.

    Parameters
    ----------
    wl_values : np.array
        Wavelength values.
        
    it_values : np.array
        Intensity values.

    rv_threshold: float
        Intensity threshold for selecting lines.
    
    Returns
    -------
    Wavelength values of the corrected spectrum.

    """
    # radial velocities computed with different lines will be saved here:
    radial_velocities = []
    # we take HI lines only:
    lines = wl_lines()['H I']+wl_lines()['He I']+wl_lines()['He II']
    
    for rest_wavelength in [wl for wl in lines if wl > 4500]:
        # Range around the rest-frame line where a peak in the observed spectrum will be searched:
        window_size = 2e-3*rest_wavelength
        # List of index in the range.
        wl_range = (wl_values >= rest_wavelength - window_size) & (wl_values <= rest_wavelength + window_size)
        # If the lines are strong enough, they are used to calculate the radial velocity:
        if np.max((abs(it_values-1))[wl_range]) > rv_threshold:
            #print(rest_wavelength)
            observed_wavelength = wl_values[wl_range][np.argmax((abs(it_values-1))[wl_range])]
            v_over_c = (observed_wavelength - rest_wavelength) / rest_wavelength
            radial_velocities.append(v_over_c)
    
    # We filter the list of radial velocites to remove outliers
    clipped_velocities = sigma_clip(radial_velocities, sigma=1.5, cenfunc='median') 
    #print(clipped_velocities)
    # Average of the good values for the velocites:
    radial_velocity = np.mean(clipped_velocities)
    print(f"radial velocity: {radial_velocity}")
    
    # The spectrum is now corrected:
    corrected_wl_values = wl_values/(1+radial_velocity)
    
    return corrected_wl_values
    


def matching_lines(wl_values, it_values,lines_dict,spectrum_type,path):

    '''
    Function to match the theorical lines with the spectrum
    '''

    # 1. Seleccionar los minimos
    # 2. Comprobar que peaks corresponden con lineas
    # 3. Refinar el match cteniendo en cuenta criterios de altura de las lineas
    # 4. Seleccionar las keys de las lineas

    # Finding the peaks
    peak_min = finding_peaks(wl_values, -it_values,window_size_factor=20)
    peaks_max = finding_peaks(wl_values, it_values,window_size_factor=20)

    peaks_all = np.concatenate((peak_min, peaks_max))

    # Selecting peaks in a range around linea wavelength
    wl_deviation = 1.5 # nm

    if spectrum_type == 'frias':

        wl_deviation = 0.75


    peaks_range = []
    wl_lines_keys = lines_dict.keys()

    for key_pos,key in enumerate(wl_lines_keys):
        for line_wl in lines_dict[key]:
            for peak_num in peaks_all:
                if line_wl-wl_deviation<wl_values[peak_num] and wl_values[peak_num]<line_wl+wl_deviation:

                    peaks_range.append(peak_num)

    # Obtaining a mask from standard deviation
    clipped_int = sigma_clip(it_values[peaks_range], sigma=0.5, cenfunc='median',masked=True) 
    peaks_sigma_mask = np.ma.getmask(clipped_int)

    peaks_sigma = (np.array(peaks_range))[peaks_sigma_mask]

    # Selecting peaks below the mean intensity value
    peaks_below = []
    mtp = 0.015 # mean tolerance parameter

    for peak_num in peaks_all:      

        if spectrum_type=='frias':
            if it_values[peak_num] < (np.mean(it_values)-(np.mean(it_values))*5*mtp):
                peaks_below.append(peak_num)

        elif it_values[peak_num] < (np.mean(it_values)-(np.mean(it_values))*mtp) or it_values[peak_num] > (np.mean(it_values)+(np.mean(it_values))*mtp):
            peaks_below.append(peak_num)

    new_lines_dict = {}

    # creating the Pandas dataframe to store the line information
    df = pd.DataFrame(columns = ['Element','Wavelength','Intensity']) 

    # Selecting the keys and values to plot
    for key_pos,key in enumerate(wl_lines_keys):
        for line_wl in lines_dict[key]:
            for peak_num in peaks_below:
                if line_wl-wl_deviation<wl_values[peak_num] and wl_values[peak_num]<line_wl+wl_deviation:
                    
                    if key not in new_lines_dict:
                        new_lines_dict[key] = [line_wl]

                        external_file(key,f'{line_wl:.3f}',f'{it_values[peak_num]:.3f}',df,path)

                    else:
                        if line_wl not in new_lines_dict[key]:
                            new_lines_dict[key].append(line_wl)

                            external_file(key,f'{line_wl:.3f}',f'{it_values[peak_num]:.3f}',df,path)

                    
    # Peaks selected that matched the lines
    peaks = peaks_all
    peaks = peaks_range

    lines_dict = new_lines_dict

    return peaks, lines_dict


def external_file(element,wl,it,df,path):

    csv_file = open(f'{path}','w+')

    df.loc[len(df)]=[element,wl,it]

    # Sort the dataframe by element
    df.sort_values(by=['Element'],ascending=True, inplace=True, ignore_index=True)
    df_string = df.to_csv(header=True, index=False, sep=',')
    csv_file.write(df_string)
    csv_file.close()


#------------RUNNING THE CODE------------#
    

#def main(run_version,spectrum_list,i=None,spec_name=None,spec_smooth_dict=None,spec_iter_dict=None):
def main():
    """
    Main function. Starts by looking for available spectra in the working directory.

    Returns
    -------
    None.

    """

    # Spectrum to analyse
    if run_version == '-m' or run_version == '-e':
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
    
    elif run_version == '-a' or  run_version == '-e':

        # Automatic selection  of fitting parameters
        if spectrum_name in calientes:
            spectrum_type = 'calientes'
        if spectrum_name in frias:
            spectrum_type = 'frias'
            
        smoothing_parameter = spec_smooth_dict[f'{spectrum_type}']
        iterations = spec_iter_dict[f'{spectrum_type}']
        rv_threshold = rv_thresholds_dict[f'{spectrum_type}']
    
    it_continuum, it_normalized, used_maxima = norm_spec(wl_values, it_values, 
                                                         s=smoothing_parameter, 
                                                         filter_iterations=iterations)
    
    # Spectrum with continuum is produced:
    plot_spec(wl_values, it_values, os.path.join(dirname, "02_continuum"), continuum=it_continuum)
    
    # Normalized spectrum is produced:
    plot_spec(wl_values, it_normalized, os.path.join(dirname, "03_normalized_spectrum"))

    # Radial velocity is corrected:
    corrected_wl_values = radialv_correct(wl_values, it_normalized, rv_threshold)
    plot_spec(corrected_wl_values, it_normalized, os.path.join(dirname, "04_normalized_spectrum_RF"))

    # Matching the lines
    path = os.path.join(dirname, f"csv_line_information.csv")
    peaks_matched,lines_matched_dict = matching_lines(corrected_wl_values, it_normalized,wl_lines_dict,spectrum_type,path)

    mtp = 0.015
    plot_spec(corrected_wl_values, it_normalized ,os.path.join(dirname, "05_normalized_lines_peaks"),peaks=peaks_matched)
    
    path = os.path.join(dirname, "06_normalized_lines_matched")
    plot_spec(corrected_wl_values, it_normalized ,path, lines_dict=lines_matched_dict)

    
    ######################################################################################################
    # Para producir el plot con las estrellas de referencia para el tipo espectral. 4471 HeI / 4541 HeII #
    ######################################################################################################
    if spectrum_name == "EstrellaProblema1.dat": # La estrella caliente tipo O
        # Estrellas de referencia para el tipo espectral:
        referencias_o_spectral = [x for x in os.listdir("Referencias_O_espectral") if ".dat" in x]
        sorted_ref = sorted(referencias_o_spectral , key=lambda x: int(x.split("_")[0]))
        espectro_problema = (corrected_wl_values[(corrected_wl_values >= 4460) & (corrected_wl_values <= 4560)],
                             it_normalized[(corrected_wl_values >= 4460) & (corrected_wl_values <= 4560)])
        norm_specs =[espectro_problema]
        names = ["Estrella problema 1"]
        for ref in sorted_ref:
            # Loading data:
            data = np.loadtxt(os.path.join("Referencias_O_espectral", ref))
            ref_wl_values = data[:,0] # Angstrons
            ref_it_values = data[:,1] # arbitrary units
            # Radial velocity is corrected:
            corrected_ref_wl_values = radialv_correct(ref_wl_values, ref_it_values, rv_threshold)
            
            # Valores de intensidad en el rango 4460-4560, para ver las lineas 4471 HeI / 4541 HeII:
            xx = corrected_ref_wl_values[(corrected_ref_wl_values >= 4460) & (corrected_ref_wl_values <= 4560)]
            yy = ref_it_values[(corrected_ref_wl_values >= 4460) & (corrected_ref_wl_values <= 4560)]
            norm_specs.append((xx,yy))
            #print(ref)
            names.append(ref[2:-4])
        spec_tower(norm_specs, os.path.join(dirname, "07_referencia_espectral"), names=names, annotation_wl=4490, lines_dict=lines_matched_dict)
        
    ######################################################################################################
    # Para producir el plot con las estrellas de referencia para el tipo de luminosidad Si IV 4116/He I 4121
    ######################################################################################################
    if spectrum_name == "EstrellaProblema1.dat": # La estrella caliente tipo O
        wl_min=4300
        wl_max=4510
        # Estrellas de referencia para el tipo espectral:
        referencias_o_lum = [x for x in os.listdir("Referencias_O_luminosidad") if ".dat" in x]
        sorted_ref = sorted(referencias_o_lum , key=lambda x: int(x.split("_")[0]))
        espectro_problema = (corrected_wl_values[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)],
                             it_normalized[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)])
        norm_specs =[espectro_problema]
        names = ["Estrella problema 1"]
        for ref in sorted_ref:
            # Loading data:
            data = np.loadtxt(os.path.join("Referencias_O_luminosidad", ref))
            ref_wl_values = data[:,0] # Angstrons
            ref_it_values = data[:,1] # arbitrary units
            # Radial velocity is corrected:
            corrected_ref_wl_values = radialv_correct(ref_wl_values, ref_it_values, rv_threshold)
            
            # Valores de intensidad en el rango 4460-4560, para ver las lineas 4471 HeI / 4541 HeII:
            xx = corrected_ref_wl_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            yy = ref_it_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            norm_specs.append((xx,yy))
            #print(ref)
            names.append(ref[2:-4])
        spec_tower(norm_specs, os.path.join(dirname, "08_referencia_luminosidad"), names=names, annotation_wl=4350, lines_dict=lines_matched_dict)
    
    ######################################################################################################
    # Para producir el plot con las estrellas de referencia para el tipo de luminosidad H I
    ######################################################################################################
    if spectrum_name == "EstrellaProblema1.dat": # La estrella caliente tipo O
        wl_min=5650
        wl_max=5750
        # Estrellas de referencia para el tipo espectral:
        referencias_o_lum = [x for x in os.listdir("Referencias_O_luminosidad") if ".dat" in x]
        sorted_ref = sorted(referencias_o_lum , key=lambda x: int(x.split("_")[0]))
        espectro_problema = (corrected_wl_values[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)],
                             it_normalized[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)])
        norm_specs =[espectro_problema]
        names = ["Estrella problema 1"]
        for ref in sorted_ref:
            # Loading data:
            data = np.loadtxt(os.path.join("Referencias_O_luminosidad", ref))
            ref_wl_values = data[:,0] # Angstrons
            ref_it_values = data[:,1] # arbitrary units
            # Radial velocity is corrected:
            corrected_ref_wl_values = radialv_correct(ref_wl_values, ref_it_values, rv_threshold)
            
            # Valores de intensidad en el rango 4460-4560, para ver las lineas 4471 HeI / 4541 HeII:
            xx = corrected_ref_wl_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            yy = ref_it_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            norm_specs.append((xx,yy))
            #print(ref)
            names.append(ref[2:-4])
        spec_tower(norm_specs, os.path.join(dirname, "09_referencia_luminosidad"), names=names, annotation_wl=5670, lines_dict=lines_matched_dict)
    
        
    ######################################################################################################
    # Para producir el plot con las estrellas de referencia para el tipo de luminosidad H I
    ######################################################################################################
    if spectrum_name == "EstrellaProblema1.dat": # La estrella caliente tipo O
        wl_min=4600
        wl_max=4700
        # Estrellas de referencia para el tipo espectral:
        referencias_o_lum = [x for x in os.listdir("Referencias_O_luminosidad") if ".dat" in x]
        sorted_ref = sorted(referencias_o_lum , key=lambda x: int(x.split("_")[0]))
        espectro_problema = (corrected_wl_values[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)],
                             it_normalized[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)])
        norm_specs =[espectro_problema]
        names = ["Estrella problema 1"]
        for ref in sorted_ref:
            # Loading data:
            data = np.loadtxt(os.path.join("Referencias_O_luminosidad", ref))
            ref_wl_values = data[:,0] # Angstrons
            ref_it_values = data[:,1] # arbitrary units
            # Radial velocity is corrected:
            corrected_ref_wl_values = radialv_correct(ref_wl_values, ref_it_values, rv_threshold)
            
            # Valores de intensidad en el rango 4460-4560, para ver las lineas 4471 HeI / 4541 HeII:
            xx = corrected_ref_wl_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            yy = ref_it_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            norm_specs.append((xx,yy))
            #print(ref)
            names.append(ref[2:-4])
        spec_tower(norm_specs, os.path.join(dirname, "10_referencia_luminosidad"), names=names, annotation_wl=4650, lines_dict=lines_matched_dict)
        
    ######################################################################################################
    # Para producir el plot con las estrellas de referencia para el tipo de luminosidad H I
    ######################################################################################################
    if spectrum_name == "EstrellaProblema1.dat": # La estrella caliente tipo O
        wl_min=4080
        wl_max=4130
        # Estrellas de referencia para el tipo espectral:
        referencias_o_lum = [x for x in os.listdir("Referencias_O_luminosidad") if ".dat" in x]
        sorted_ref = sorted(referencias_o_lum , key=lambda x: int(x.split("_")[0]))
        espectro_problema = (corrected_wl_values[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)],
                             it_normalized[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)])
        norm_specs =[espectro_problema]
        names = ["Estrella problema 1"]
        for ref in sorted_ref:
            # Loading data:
            data = np.loadtxt(os.path.join("Referencias_O_luminosidad", ref))
            ref_wl_values = data[:,0] # Angstrons
            ref_it_values = data[:,1] # arbitrary units
            # Radial velocity is corrected:
            corrected_ref_wl_values = radialv_correct(ref_wl_values, ref_it_values, rv_threshold)
            
            # Valores de intensidad en el rango 4460-4560, para ver las lineas 4471 HeI / 4541 HeII:
            xx = corrected_ref_wl_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            yy = ref_it_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            norm_specs.append((xx,yy))
            #print(ref)
            names.append(ref[2:-4])
        spec_tower(norm_specs, os.path.join(dirname, "11_referencia_luminosidad"), names=names, annotation_wl=4081, lines_dict=lines_matched_dict)

        
    ######################################################################################################
    # Para producir el plot con las estrellas de referencia para el tipo espectral
    ######################################################################################################
    if spectrum_name == "EstrellaProblema2.dat": # La estrella caliente tipo O
        wl_min=4040
        wl_max=4120
        # Estrellas de referencia para el tipo espectral:
        referencias_o_lum = [x for x in os.listdir("Referencias_G_espectral") if ".dat" in x]
        sorted_ref = sorted(referencias_o_lum , key=lambda x: int(x.split("_")[0]))
        espectro_problema = (corrected_wl_values[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)],
                             it_normalized[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)])
        norm_specs =[espectro_problema]
        names = ["Estrella problema 2"]
        for ref in sorted_ref:
            # Loading data:
            data = np.loadtxt(os.path.join("Referencias_G_espectral", ref))
            ref_wl_values = data[:,0] # Angstrons
            ref_it_values = data[:,1] # arbitrary units
            # Radial velocity is corrected:
            corrected_ref_wl_values = radialv_correct(ref_wl_values, ref_it_values, rv_threshold)
            
            # Valores de intensidad en el rango 4460-4560, para ver las lineas 4471 HeI / 4541 HeII:
            xx = corrected_ref_wl_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            yy = ref_it_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            norm_specs.append((xx,yy))
            #print(ref)
            names.append(ref[2:-4])
        spec_tower(norm_specs, os.path.join(dirname, "07_referencia_espectral"), names=names, annotation_wl=4050, lines_dict=lines_matched_dict, factor=1)
        
    ######################################################################################################
    # Para producir el plot con las estrellas de referencia para el tipo espectral
    ######################################################################################################
    if spectrum_name == "EstrellaProblema2.dat": # La estrella caliente tipo O
        wl_min=3900
        wl_max=4100
        # Estrellas de referencia para el tipo espectral:
        referencias_o_lum = [x for x in os.listdir("Referencias_G_espectral") if ".dat" in x]
        sorted_ref = sorted(referencias_o_lum , key=lambda x: int(x.split("_")[0]))
        espectro_problema = (corrected_wl_values[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)],
                             it_normalized[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)])
        norm_specs =[espectro_problema]
        names = ["Estrella problema 2"]
        for ref in sorted_ref:
            # Loading data:
            data = np.loadtxt(os.path.join("Referencias_G_espectral", ref))
            ref_wl_values = data[:,0] # Angstrons
            ref_it_values = data[:,1] # arbitrary units
            # Radial velocity is corrected:
            corrected_ref_wl_values = radialv_correct(ref_wl_values, ref_it_values, rv_threshold)
            
            # Valores de intensidad en el rango 4460-4560, para ver las lineas 4471 HeI / 4541 HeII:
            xx = corrected_ref_wl_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            yy = ref_it_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            norm_specs.append((xx,yy))
            #print(ref)
            names.append(ref[2:-4])
        spec_tower(norm_specs, os.path.join(dirname, "08_referencia_espectral"), names=names, annotation_wl=3950, lines_dict=lines_matched_dict, factor=1)

    ######################################################################################################
    # Para producir el plot con las estrellas de referencia para el tipo espectral. banda G
    ######################################################################################################
    if spectrum_name == "EstrellaProblema2.dat": # La estrella caliente tipo O
        wl_min=3900
        wl_max=4500
        # Estrellas de referencia para el tipo espectral:
        referencias_o_lum = [x for x in os.listdir("Referencias_G_espectral") if ".dat" in x]
        sorted_ref = sorted(referencias_o_lum , key=lambda x: int(x.split("_")[0]))
        espectro_problema = (corrected_wl_values[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)],
                             it_normalized[(corrected_wl_values >= wl_min) & (corrected_wl_values <= wl_max)])
        norm_specs =[espectro_problema]
        names = ["Estrella problema 2"]
        for ref in sorted_ref:
            # Loading data:
            data = np.loadtxt(os.path.join("Referencias_G_espectral", ref))
            ref_wl_values = data[:,0] # Angstrons
            ref_it_values = data[:,1] # arbitrary units
            # Radial velocity is corrected:
            corrected_ref_wl_values = radialv_correct(ref_wl_values, ref_it_values, rv_threshold)
            
            # Valores de intensidad en el rango 4460-4560, para ver las lineas 4471 HeI / 4541 HeII:
            xx = corrected_ref_wl_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            yy = ref_it_values[(corrected_ref_wl_values >= wl_min) & (corrected_ref_wl_values <= wl_max)]
            norm_specs.append((xx,yy))
            #print(ref)
            names.append(ref[2:-4])
            
        lines = {key: lines_matched_dict[key] for key in ['Ca II', 'H I'] if key in lines_matched_dict}
        spec_tower(norm_specs, os.path.join(dirname, "09_referencia_espectral"), names=names, annotation_wl=4140, lines_dict=lines, factor=1)
       

    

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

    wl_lines_dict = wl_lines()

    if run_version == '-a':

        # Loading fitting parameters dictionaries
        spec_smooth_dict, spec_iter_dict, rv_thresholds_dict = fiting_parameters()
        calientes, frias = tipos_estrella()

        # Analyzing each spectrum
        for spect_pos, spec_name in enumerate(spectrum_list):
  
            main()

    else:

        calientes, frias = tipos_estrella()
        spec_smooth_dict, spec_iter_dict, rv_thresholds_dict = fiting_parameters()
        main()

