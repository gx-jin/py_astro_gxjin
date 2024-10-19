# -*- coding: utf-8 -*-
"""
    Copyright (C) 2019-2024, Gaoxiang Jin

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
from requests.auth import HTTPBasicAuth
import requests

def download(url: str, local_filename: str,
             username='None', password='None',
             quiet=True):
    
    headers = {'User-Agent':
        'Mozilla/5.0 (Windows NT 10.0; Win64; x64) \
            AppleWebKit/537.36 (KHTML, like Gecko) \
                Chrome/130.0.0.0 Safari/537.36 Edg/130.0.0.0'}
    
    # Check parent directory
    parent_directory = os.path.dirname(local_filename)
    if not os.path.isdir(os.path.dirname(local_filename)):
        print(f"The parent directory '{parent_directory}' does not exist")
        return
    
    # Check if file exists
    if os.path.exists(local_filename):
        print(f"File exists: '{local_filename}'")
        return
    
    # Try download
    try:
        response = requests.head(url, 
                                 headers=headers, stream=True, timeout=540,
                                 auth = HTTPBasicAuth(username, password))
        if response.status_code == requests.codes.ok:
            with requests.get(url,
                              headers=headers, stream=True, timeout=540,
                              auth = HTTPBasicAuth(username, password)) as r:
                r.raise_for_status()
                with open(local_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            if not quiet:
                print(f"File '{local_filename}' downloaded")
            return        
        else:
            print(f'HTTP error occurred: {response.status_code}')
    except requests.RequestException:
        print('Invalid url')