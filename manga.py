# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 01:24:50 2023
@author: Gaoxiang Jin

    Copyright (C) 2019-2023, Gaoxiang Jin

    E-mail: gx-jin@outlook.com

    Updated versions of the codes are available from github pages:
    https://

    This software is provided as is without any warranty whatsoever.
    Permission to use, for non-commercial purposes is granted.
    Permission to modify for personal or internal use is granted,
    provided this copyright and disclaimer are included unchanged
    at the beginning of the file. All other rights are reserved.
    In particular, redistribution of the code is not allowed.

"""

import os
import requests
import shutil

###############################################################################
# Constant

c = 2.99792458e5  # speed of light, km/s

################################################################################


def download_logcube(plate, ifudesign, daptype='SPX', save_dir='./'):
    """Download the MaNGA DAP LOGCUBE file to given path, version SDSS DR17,
        SSP MaSTAR.

    Args:
        plate (int): MaNGA PLATE
        ifudesign (int): MaNGA IFUDESIGN 
        daptype (str, optional): 'SPX', 'VOR10', or 'HYB10'. Defaults to 'SPX'.
        save_dir (str, optional): Directory for saving. Defaults to current path.
    """
     
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; \
                rv:80.0) Gecko/20100101 Firefox/80.0'}
            
    if not os.path.exists(save_dir):
        raise RuntimeError("No such directory: " + str(save_dir))
    else:
        save_loc = save_dir + '/manga-' + str(plate) + '-' + str(ifudesign) \
                   + '-LOGCUBE-' + daptype + '-MILESHC-MASTARSSP.fits.gz'
        if os.path.exists(save_loc):
            return
        logcube_url = 'https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/SPX-MILESHC-MASTARSSP/' \
            + str(plate) + '/' + str(ifudesign) + '/manga-' + str(plate) + '-' + str(ifudesign) \
            + '-LOGCUBE-' + daptype + '-MILESHC-MASTARSSP.fits.gz'
        with requests.get(logcube_url, headers=headers, stream=True) as r:
            if r.status_code == requests.codes.ok:
                with open(save_loc, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)
        return

