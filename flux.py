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


def sfr_from_fha(fha, ld, error=None, imf='Chabrier'):
    
    # TODO: description
    
    ld = np.float64(ld)
    
    if imf == 'Chabrier':
        convert = - 41.27 + np.log10(0.63 / 0.67)
    elif imf == 'Kroupa':
        convert = - 41.27
    elif imf == 'Salpeter':
        convert = - 41.27 + np.log10(1.0 / 0.67)
    else:
        raise RuntimeError("Wrong IMF name! Use 'Chabrier', 'Kroupa', 'Salpeter'.")
    
    if error is None:
        ha_array = unumpy.uarray(fha, 0)
    else:  # (fha.shape == error.shape) & (fha.shape == ld.shape):
        ha_array = unumpy.uarray(fha, error)
    # else:
    #     raise RuntimeError("'fha', 'ld', and 'error' must have same shapes!.")
    
    sfr = unumpy.log10(ha_array * 4.0 * np.pi * np.power(ld, 2)) + convert
    
    return unumpy.nominal_values(sfr), unumpy.std_devs(sfr)


def ha_agncorr_lines(lines, error=None):
    
    # TODO: description and raise errors
      
    if (len(lines) == 6) & (len(error) == 6):
        fhb = unumpy.uarray(lines[0], error[0])
        fo3 = unumpy.uarray(lines[1], error[1])
        fha = unumpy.uarray(lines[2], error[2])
        fn2 = unumpy.uarray(lines[3], error[3])
        fs2 = unumpy.uarray(lines[4], error[4]) + unumpy.uarray(lines[5], error[5])
        o3hb = np.log10(unumpy.nominal_values(fo3 / fhb))
        n2ha = np.log10(unumpy.nominal_values(fn2 / fha))
        s2ha = np.log10(unumpy.nominal_values(fs2 / fha))
        
        p1 = 0.63 * n2ha + 0.51 * s2ha + 0.59 * o3hb  # pca component from Ji+2019
        fsfr = 1.0 - np.clip(0.14 * p1**2 + 0.96 * p1 + 0.47, 0.0, 0.9)
        
        cha_unlim = (fha / fhb / 2.86)**unumpy.uarray(2.62404, 0.6304)  # Dominguez+2013 ApJ Eq.6
        # Av range: 0-3
        cha_val = np.clip(unumpy.nominal_values(cha_unlim), 1.0, 10**1.3)
        cha = unumpy.uarray(cha_val, unumpy.std_devs(cha_unlim))
        
        fha_all = fha * cha
    
    else:
        raise RuntimeError("'lines' and 'error' must have shapes of 6!.")    
    # if error is None:
    #     ha_array = unumpy.uarray(fha, 0)
    # elif ~((fha.shape == 6) & (error.shape == 6)):
    #     raise RuntimeError("'fha' and 'error' must include 6 lines!.")
    #     ha_array = unumpy.uarray(fha, error)
    # elif (fha[0].shape == ld.shape) & (error[0].shape == ld.shape):
    #     raise RuntimeError("'fha', 'ld', and 'error' must have same shapes!.")

    return unumpy.nominal_values(fha_all) * fsfr, unumpy.std_devs(fha_all), fsfr
    
    
# TODO def flux_dustcorr_hahb(fold, fha, fhb, fold_err=0.0, fha_err=0.0, fhb_err=0.0,): 

# TODO def flux_units(fold, unitsold, unitsnew, fold_err=0.0, ): 

