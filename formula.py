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

# from uncertainties import unumpy
# import sys
import numpy as np

###############################################################################


def bpt_n2_ka03(x):
    y = 0.61 / (x - 0.05) + 1.3
    y[x >= 0.05] = np.nan
    return y
    
    
def bpt_n2_ke01(x):
    y = 0.61 / (x - 0.47) + 1.19
    y[x >= 0.47] = np.nan
    return y


def bpt_s2_ke01(x):
    y = 0.72 / (x - 0.32) + 1.3
    y[x >= 0.32] = np.nan
    return y


def l144_sfr_jin(x):
    y = 1.161 * x + 22.011
    y[x >= 1.6] = np.nan
    y[x <= -1.4] = np.nan
    return y
