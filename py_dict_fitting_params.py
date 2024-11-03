
def fiting_parameters():

    spectrum_smoothing = {'calientes':0.012,
                          'frias':0.010
                        }
    

    spectrum_iterations = {'calientes':5,
                           'frias':2
                        }
    
    radial_velocity_threshold = {'calientes':0.1,
                                 'frias':0.6
                                 }
    return spectrum_smoothing,spectrum_iterations,radial_velocity_threshold


def tipos_estrella():
    
    calientes = ['EstrellaProblema1.dat',
                ]
    
    frias = ['EstrellaProblema2.dat',
            ]

    return calientes, frias
