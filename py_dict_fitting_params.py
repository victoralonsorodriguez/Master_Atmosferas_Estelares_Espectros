
def fiting_parameters():

    spectrum_smoothing = {'EstrellaProblema1.dat':0.012,
                        'EstrellaProblema2.dat':0.010
                        }
    

    spectrum_iterations = {'EstrellaProblema1.dat':5,
                        'EstrellaProblema2.dat':2
                        }
    
    radial_velocity_threshold = {'EstrellaProblema1.dat':0.8,
                        'EstrellaProblema2.dat':0.4
                        }
    return spectrum_smoothing,spectrum_iterations,radial_velocity_threshold
