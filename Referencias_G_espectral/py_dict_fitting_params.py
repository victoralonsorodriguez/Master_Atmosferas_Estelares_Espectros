
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
                 '6_HD 214680_O9V.dat',
                 '1_HD 192639_O7.5Iabf.dat',
                 '2_HD 34656_O7.5II(f).dat',
                 '3_HD 186980_O7.5III((f)).dat',
                 '4_HD 164492_O7.5Vz.dat'
                ]
    
    frias = ['EstrellaProblema2.dat',
             '1_HD141004_G0V.dat',
             '2_HD209750_G2I.dat',
             '3_HD50806_G5V.dat',
             '4_HD101501_G8V.dat',
             '5_HD48329_G8I.dat',
             '6_HD3651_K0V.dat',
             '7_HD200905_K4.5I.dat',
             '8_HD201091_K5V.dat'
            ]

    return calientes, frias
