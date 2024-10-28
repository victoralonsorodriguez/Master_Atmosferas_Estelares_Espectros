
def wl_lines():

    # Obtained from https://linelist.pa.uky.edu/newpage/
    # range between 400 and 6500 A

    wv_by_elements = {
        'Ca I':[4226.1,4260.636],
        'Ca II':[3934.777],
        'Fe I':[4045.03974,4046.95504,4064.8481,4072.66979,4144.58312,
                4172.07716,4215.53935,4271.07578,4383.99000,5250.8353,5259.0459],
        'H I':[4102.892,4341.684,4862.683,6564.610],
        'He I':[4010.3898875,4025.1170140,4027.3221728,4121.9733341,
               4144.9276331,4170.1466968,4389.3811248,4472.3340652,
               4714.457832,4923.3050507,5017.0769260,5049.146012,5877.2271632],
        'He II':[4026.739,4101.198,4201.015,4339.891,4542.864,4687.02,4860.677,
                5413.03,5753.673],
        'Mg I':[5168.7605,5174.1251,5185.0479],
        'Mg II':[4482.383],
        'N II':[3995.373],
        'N III':[4634.24,4640.7,4642.59],
        'O II': [4416.1384],
        'Si II':[4129.218],
        'Si III':[4553.898],
        #'SiIII':[4011.281,4031.889,4066.525,4072.623,4081.944,4102.87,4103.581,
        #         4112.33,4116.648,4122.812,4145.183,4553.898]
        'Si IV':[4090.016,4117.265]
        #'Sr II':[]

    }

    return wv_by_elements