import numpy as np
import astropy.io.fits as pf
import os
import shutil

cwd = os.getcwd()
spec_list = os.listdir(cwd)
print(spec_list)


for spectrum in [x for x in spec_list if ".fits" in x]:

    
    ### We read the header, where we can find info about the target (name, coordinates...) and the observations (date, exposure time, sky conditions...) ###
    headerc = pf.getheader(spectrum)
    starname = headerc['OBJECT'] 
    
    ### We read the data ###
    specc = pf.getdata(spectrum)
    wave = specc['wave']
    flux = specc['flux_norm'] #If you want it NORMALIZED
    #flux = specc['flux_tac'] #If you want it NOT normalized
    
    ### We save wavelength and flux in a 2-column .dat file with the same original name ###
    np.savetxt(spectrum[:-5]+'.dat', np.array([wave, flux]).T, fmt='%f', delimiter='\t')
    os.remove(spectrum)
    print(f"{spectrum} saved")

