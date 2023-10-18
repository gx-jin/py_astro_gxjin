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
import matplotlib.pyplot as plt
import urllib
from PIL import Image

###############################################################################
# Constant

c = 2.99792458e5  # speed of light, km/s

################################################################################


# def download_file(url):
#     headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; \
#                 rv:80.0) Gecko/20100101 Firefox/80.0'}
#     with requests.get(url, headers=headers, stream=True) as r:
#         return r
#     r.close()
                

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
          
    if not os.path.exists(save_dir):
        raise RuntimeError("No such directory: " + str(save_dir))
    else:
        save_loc = f'{save_dir}/manga-{plate}-{ifudesign}-LOGCUBE-{daptype}-MILESHC-MASTARSSP.fits.gz'
        if os.path.exists(save_loc):
            if ~quiet:
                print('File existed')
        else:
            logcube_url = f'https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/{daptype}-MILESHC-MASTARSSP/{plate}/{ifudesign}/manga-{plate}-{ifudesign}-LOGCUBE-{daptype}-MILESHC-MASTARSSP.fits.gz'
            headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; \
                rv:80.0) Gecko/20100101 Firefox/80.0'}
            with requests.get(logcube_url, headers=headers, stream=True) as r:
                if r.status_code == requests.codes.ok:
                    with open(save_loc, 'wb') as f:
                        shutil.copyfileobj(r.raw, f)
            if not quiet:
                print(f"File downloaded: {save_loc}")


def stack_logcube(mode, region,
                  plate=None, ifudesign=None, cubefile=None, mapsfile=None,
                  rmin=0., rmax=-1., pamin=0., pamax=360., recustom=-1.0,):
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
        recustom (float, optional): For 'rering' mode, the 1 unit radius in arcsec. Must be positive number.
        
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
        hducube = fits.open(cubefile)
        hdumaps = fits.open(mapsfile)
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
        rfwhm = hdumaps[0].header['RFWHM'] / 2.0
    except Exception:
        raise RuntimeError("Not valid LOGCUBE or MAPS file!")

    try:
        bta = 1.0 - hdumaps[0].header['ECOOELL']
        # redap = hdumaps[0].header['REFF']
        padap = hdumaps[0].header['ECOOPA']
    except Exception:
        bta = 1.0
        # redap = 1.0
        padap = 0.0
        print("Warning: No valid 'Re' measurement!")

    arcmap = hdumaps['SPX_ELLCOO'].data[0, :, :]
    mapmaxarc, mapmaxre = np.nanmax(arcmap), np.nanmax(remap)
    rmax = rmax if rmax > 0 else np.max([mapmaxarc, mapmaxre])
    roundmap = np.sqrt((np.cos(pamap / 180.0 * np.pi) * arcmap)**2 + 
                       (np.sin(pamap / 180.0 * np.pi) * arcmap * bta)**2)

    if region == 'all':
        rcut = remap >= 0
    elif region == 'cen':
        rcut = roundmap <= 1.5  # rfwhm
    elif (region == 'rering'):
        if (rmax * mapmaxarc / mapmaxre) < rfwhm:
            print("Warning: 'rmax' is smaller than the resolution!")
        if recustom > 0:
            remap = arcmap / recustom
        rcut = (remap <= rmax) & (remap >= rmin)
    elif region == 'pasector':
        if (pamax > pamin):
            if (rmax * mapmaxarc / mapmaxre) < rfwhm:
                print("Warning: 'rmax' is smaller than the resolution! Return empty spectra.")
            pamap = np.mod(pamap + padap, 360.0) - pamin        
            rcut = (remap <= rmax) & (remap >= rmin) & (roundmap >= rfwhm) & ((pamap >= 0) & (pamap <= (pamax - pamin)) | (pamap >= 180.0) & (pamap <= (pamax - pamin + 180.0)))
        else:
            raise RuntimeError("Invalid 'pamin' or 'pamax'! Define them in clockwise direction.")
    
    nstack = np.sum(rcut)
    if nstack == 0:
        print("Warning: No spectra found in defined region! Return empty spectra.")
        return (wave, np.zeros_like(wave) - np.nan, np.zeros_like(wave) - np.nan, np.zeros_like(wave, dtype=int))
    else:
        flux[mask > 0] = np.nan  # mask 5578 skyline
        ivar[mask > 0] = np.nan
        goodspec = np.sum(~np.isnan(flux[:, rcut]), axis=1)
        goodspec[(wave >= 5570) & (wave <= 5586)] = 0
        fluxstack = np.nanmean(flux[:, rcut], axis=1) * goodspec
        errstack = np.sqrt(np.nanmean(1. / ivar[:, rcut], axis=1) * goodspec) 
        fluxstack[(wave >= 5570) & (wave <= 5586)] = np.nan
        errstack[(wave >= 5570) & (wave <= 5586)] = np.nan
        return wave, fluxstack, errstack, goodspec
    
    hducube.close()
    hdumaps.close()  


def plot_gri_image(plate, ifudesign):
    url = f'https://data.sdss.org/sas/dr17/manga/spectro/redux/v3_1_1/{plate}/images/{ifudesign}.png'
    image = Image.open(urllib.request.urlopen(url))
    plt.figure(figsize=(3, 3), layout='compressed', frameon=False)
    plt.imshow(image)
    plt.axis('off')
    plt.show()
        
         
def plot_manga_edge(ax, ifu, color='k', ls='--', lw=2, ):
    asx1 = np.linspace(-17, -8.5, 10)
    asx2 = np.linspace(-8.5, 8.5, 10)
    asx3 = np.linspace(8.5, 17, 100)
    asy1 = (asx1) * (3**0.5) + 29.44486
    asy2 = (-3**0.5) * (asx1) - 29.44486
    asy3 = (asx2) * 0 + 14.72243
    asy4 = (asx2) * 0 - 14.72243
    asy5 = (asx3) * (-3**0.5) + 29.44486
    asy6 = (asx3) * (3**0.5) - 29.44486

    bsx1 = np.linspace(-14.5, -7.25, 10)
    bsx2 = np.linspace(-7.25, 7.25, 10)
    bsx3 = np.linspace(7.25, 14.5, 10)
    bsy1 = (bsx1) * (3**0.5) + 25.1147367
    bsy2 = (-3**0.5) * (bsx1) - 25.1147367
    bsy3 = (bsx2) * 0 + 12.557368
    bsy4 = (bsx2) * 0 - 12.557368
    bsy5 = (bsx3) * (-3**0.5) + 25.1147367
    bsy6 = (bsx3) * (3**0.5) - 25.1147367

    csx1 = np.linspace(-12, -6, 100)
    csx2 = np.linspace(-6, 6, 100)
    csx3 = np.linspace(6, 12, 100)
    csy1 = (csx1) * (3**0.5) + 20.7846
    csy2 = (-3**0.5) * (csx1) - 20.7846
    csy3 = (csx2) * 0 + 10.3923
    csy4 = (csx2) * 0 - 10.3923
    csy5 = (csx3) * (-3**0.5) + 20.7846
    csy6 = (csx3) * (3**0.5) - 20.7846

    dsx1 = np.linspace(-9.5, -4.75, 100)
    dsx2 = np.linspace(-4.75, 4.75, 100)
    dsx3 = np.linspace(4.75, 9.5, 100)
    dsy1 = (dsx1) * (3**0.5) + 16.45448
    dsy2 = (-3**0.5) * (dsx1) - 16.45448
    dsy3 = (dsx2) * 0 + 8.22724
    dsy4 = (dsx2) * 0 - 8.22724
    dsy5 = (dsx3) * (-3**0.5) + 16.45448
    dsy6 = (dsx3) * (3**0.5) - 16.45448

    esx1 = np.linspace(-7, -3.5, 100)
    esx2 = np.linspace(-3.5, 3.5, 100)
    esx3 = np.linspace(3.5, 7.0, 100)
    esy1 = (esx1) * (3**0.5) + 12.12435565
    esy2 = (-3**0.5) * (esx1) - 12.12435565
    esy3 = (esx2) * 0 + 6.0622
    esy4 = (esx2) * 0 - 6.0622
    esy5 = (esx3) * (-3**0.5) + 12.12435565
    esy6 = (esx3) * (3**0.5) - 12.12435565

    frame127 = [asx1, asx2, asx3, asy1, asy2, asy3, asy4, asy5, asy6]
    frame91 = [bsx1, bsx2, bsx3, bsy1, bsy2, bsy3, bsy4, bsy5, bsy6]
    frame61 = [csx1, csx2, csx3, csy1, csy2, csy3, csy4, csy5, csy6]
    frame37 = [dsx1, dsx2, dsx3, dsy1, dsy2, dsy3, dsy4, dsy5, dsy6]
    frame19 = [esx1, esx2, esx3, esy1, esy2, esy3, esy4, esy5, esy6]
    
    if '127' in str(ifu):
        frame = frame127
    elif '91' in str(ifu):
        frame = frame91
    elif '61' in str(ifu):
        frame = frame61
    elif '37' in str(ifu):
        frame = frame37
    elif '19' in str(ifu):
        frame = frame19
    
    for i in range(6):
        ax.plot(frame[int(i / 2)], frame[i + 3], color=color, linestyle=ls, linewidth=lw)
        

# TODO ppxf_miles_fitting_lines(spectra, wave, z,)
