
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
                 '1_HD 42088_O6V((f))z.dat',
                 '2_HD 199579_O6.5V((f))z.dat',
                 '3_HD 46573_O7V((f))z.dat',
                 '4_HD 164492_O7.5Vz.dat',
                 '5_HD 41161_O8Vn.dat',
                 '6_HD 214680_O9V.dat'
                ]
    
    frias = ['EstrellaProblema2.dat',
            ]

    return calientes, frias
