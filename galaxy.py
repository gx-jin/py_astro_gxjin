# -*- coding: utf-8 -*-
"""
    Copyright (C) 2019-2023, Gaoxiang Jin

    E-mail: gx-jin@outlook.com

    Updated versions of the codes are available from github pages:
    https://github.com/gx-jin/py_astro_gxjin

    This software is provided as is without any warranty whatsoever.
    Permission to use, for non-commercial purposes is granted.
    Permission to modify for personal or internal use is granted,
    provided this copyright and disclaimer are included unchanged
    at the beginning of the file. All other rights are reserved.
    In particular, redistribution of the code is not allowed.
"""

from uncertainties import unumpy
# import sys
import numpy as np

###############################################################################
# Constant

c = 2.99792458e5  # speed of light, km/s

################################################################################


def sfr_sfms(mass, z=0.0, imf='Chabrier'):
    
    # TODO: description, Speagle et al. 2014
    
    if imf == 'Chabrier':
        convert = 0. + np.log10(0.63 / 0.67)
    elif imf == 'Kroupa':
        convert = 0.
    elif imf == 'Salpeter':
        convert = 0. + np.log10(1.0 / 0.67)
    else:
        raise RuntimeError("Wrong IMF name! Use 'Chabrier', 'Kroupa', 'Salpeter'.")

    # else:
    #     raise RuntimeError("'fha', 'ld', and 'error' must have same shapes!.")
    
    sfr = unumpy.log10(ha_array * 4.0 * np.pi * np.power(ld, 2)) + convert
    
    return unumpy.nominal_values(sfr), unumpy.std_devs(sfr)