"""
Fit distributions from the cosmos field
"""

import os
version='COSMOS_23.5_training_sample_sub'
def get_dir():
    d=os.environ['COSMOS_DIR']
    return os.path.join(d, version)

def get_cat_path():
    """
    Get the path to the main catalog (not the one with fit parameters)
    """
    d=get_dir()
    path=os.path.join(d, 'real_galaxy_catalog_23.5_sub.fits')
    return path

def get_fits_cat_path():
    """
    Get the path to the fit parameters
    """
    d=get_dir()
    path=os.path.join(d,'real_galaxy_catalog_23.5_sub_fits.fits')
    return path
