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

import os
import requests
import shutil
from astropy.io import fits
import numpy as np

###############################################################################
# Constant

c = 2.99792458e5  # speed of light, km/s

################################################################################


def download_logcube(plate, ifudesign, 
                     daptype='SPX', save_dir='./', quiet=True, ):
    """Download the MaNGA DAP LOGCUBE file to given path, version SDSS DR17,
        SSP MaSTAR.

    Args:
        plate (int): MaNGA PLATE
        ifudesign (int): MaNGA IFUDESIGN 
        daptype (str, optional): 'SPX', 'VOR10', or 'HYB10'. Default = 'SPX'.
        save_dir (str, optional): Directory for saving. Default = current path.
        quiet (bool, optional): Print result or not. Default is not.
    """
     
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; \
                rv:80.0) Gecko/20100101 Firefox/80.0'}
            
    if not os.path.exists(save_dir):
        raise RuntimeError("No such directory: " + str(save_dir))
    else:
        save_loc = f'{save_dir}/manga-{plate}-{ifudesign}-LOGCUBE-{daptype}-MILESHC-MASTARSSP.fits.gz'
        if os.path.exists(save_loc):
            if ~quiet:
                print('File existed')
        logcube_url = f'https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/{daptype}-MILESHC-MASTARSSP/{plate}/{ifudesign}/manga-{plate}-{ifudesign}-LOGCUBE-{daptype}-MILESHC-MASTARSSP.fits.gz'
        with requests.get(logcube_url, headers=headers, stream=True) as r:
            if r.status_code == requests.codes.ok:
                with open(save_loc, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)
        if ~quiet:
            print(f"File downloaded: {save_loc}")


def stack_logcube(mode, region,
                  plate=None, ifudesign=None, cubefile=None, mapsfile=None,
                  rmin=0., rmax=-1., pamin=0., pamax=360.):
    """Stack DAP LOGCUBE spectra based on given morphology mode.

    Args:
        mode (str): 'online' or 'local'.
        region (str): Type of stacked region. Supported: 'all', 'cen', 'rering', 'pasector'.
        plate (int, optional): Plate number, only for 'online' mode. Defaults to None.
        ifudesign (int, optional): IFU number, only for 'online' mode. Defaults to None.
        cubefile (str, optional): File path for LOGCUBE file, only for 'local' mode. Defaults to None.
        mapsfile (str, optional): File path for MAPS file, only for 'local' mode. Defaults to None.
        rmin (float, optional): For 'rering' or 'pasector' mode, the inner radius of the ring or sector region. Defaults to 0.0.
        rmax (float, optional): For 'rering' or 'pasector' mode, the outer radius of the ring region. Will be changed to resolution if remax<fwhm. Defaults to -1.0.
        pamin (float, optional): For 'pasector' mode, the start PA in clockwise direction. Defaults to 0.0.
        pamax (float, optional): For 'pasector' mode, the end PA in clockwise direction. Defaults to 360.0.

    Returns:
        turple with 4: observed wavelength, stacked spectra, 1 sigma error, number of spectra used for stacking along the wavelength
    """    
    
    if ~(np.sum(region == np.array(['all', 'cen', 'rering', 'pasector'])) == 1):
        raise RuntimeError("Stack region should be 'tot', 'cen', 'rering', or 'pasector'!")
    
    if mode == 'online':
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; \
                   rv:80.0) Gecko/20100101 Firefox/80.0'}
        locurl = f'https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/SPX-MILESHC-MASTARSSP/{plate}/{ifudesign}/manga-{plate}-{ifudesign}'          
        with requests.get(locurl + '-LOGCUBE-SPX-MILESHC-MASTARSSP.fits.gz',
                          headers=headers) as cubefile:
            if cubefile.status_code == requests.codes.ok:
                hducube = fits.open(cubefile.content)
            else:
                raise RuntimeError(f'No LOGCUBE file for {plate}-{ifudesign}')
        with requests.get(locurl + '-MAPS-SPX-MILESHC-MASTARSSP.fits.gz',
                          headers=headers) as mapsfile:
            if cubefile.status_code == requests.codes.ok:
                hdumaps = fits.open(mapsfile.content)
            else:
                raise RuntimeError('No MAPS file for ' + str(plate) + '-' + str(ifudesign)) 
    elif mode == 'local':
        try:    
            hducube = fits.open(cubefile)
            hdumaps = fits.open(mapsfile)
        except Exception:
            raise RuntimeError("Wrong filename!")
    else:
        raise RuntimeError("'mode' should be 'online' or 'local'!")
    try:  
        wave = hducube['WAVE'].data
        flux = hducube['FLUX'].data
        ivar = hducube['IVAR'].data
        mask = hducube['MASK'].data
        remap = hdumaps['SPX_ELLCOO'].data[1, :, :]
        # kpcmap = hdumaps['SPX_ELLCOO'].data[2, :, :] / 0.7
        pamap = hdumaps['SPX_ELLCOO'].data[3, :, :]
        rfwhm = hdumaps[0].header['RFWHM']
    except Exception:
        raise RuntimeError("Not valid LOGCUBE or MAPS file!")

    try:
        bta = 1.0 - hdumaps[0].header['ECOOELL']
    except Exception:
        bta = 1.0
        raise RuntimeWarning("No valid 'Re' measurement!")

    rmax = rmax if rmax > 0 else np.nanmax(remap)
    arcmap = hdumaps['SPX_ELLCOO'].data[0, :, :]
    roundmap = np.sqrt((np.cos(pamap / 180.0 * np.pi) * arcmap)**2 + 
                       (np.sin(pamap / 180.0 * np.pi) * arcmap * bta)**2)

    if region == 'all':
        rcut = remap >= 0
    elif region == 'cen':
        rcut = roundmap <= rfwhm
    elif (region == 'rering'):
        if rmax < rfwhm:
            raise RuntimeWarning("'rmax' is smaller than the resolution!")
        rcut = (remap <= rmax) & (remap >= rmin)
    elif region == 'pasector':
        if (pamin >= 0.) & (pamax <= 360.) & (pamax > pamin):
            if rmax < rfwhm:
                raise RuntimeWarning("'rmax' is smaller than the resolution! Return empty spectra.")            
            rcut = (remap <= rmax) & (remap >= rmin) & (pamap >= pamin) & (pamap <= pamax) & (roundmap >= rfwhm)
        else:
            raise RuntimeWarning("Invalid 'pamin' or 'pamax'! Define them in clockwise direction.")
    
    nstack = np.sum(rcut)
    if nstack == 0:
        raise RuntimeWarning("No spectra found in defined region! Return empty spectra.")
        return (wave, np.zeros_like(wave) - np.nan, np.zeros_like(wave) - np.nan, np.zeros_like(wave, dtype=int))
    else:
        flux[mask > 0] = np.nan
        ivar[mask > 0] = np.nan
        goodspec = np.nansum(~np.isnan(flux[:, rcut]), axis=1)
        fluxstack = np.nansum(flux[:, rcut], axis=1)
        errstack = np.sqrt(np.nansum(1. / ivar[:, rcut], axis=1)) 
        return wave, fluxstack, errstack, goodspec
    
    hducube.close()
    hdumaps.close()
        
        
# TODO ppxf_miles_fitting_lines(spectra, wave, z,)
