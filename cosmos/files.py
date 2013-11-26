import os

VERSION='23.5'
def get_dir():
    d=os.environ['COSMOS_DIR']
    vdir='COSMOS_%s_training_sample_sub' % VERSION
    return os.path.join(d, vdir)

def get_cat_path():
    """
    Get the path to the main catalog (not the one with fit parameters)
    """
    d=get_dir()
    path=os.path.join(d, 'real_galaxy_catalog_%s_sub.fits' % VERSION)
    return path

def get_fits_cat_path():
    """
    Get the path to the fit parameters
    """
    d=get_dir()
    path=os.path.join(d,'real_galaxy_catalog_%s_sub_fits.fits' % VERSION)
    return path

def read_cat():
    import fitsio
    path=get_cat_path()
    print 'reading:',path
    return fitsio.read(path,lower=True)

def read_fits_cat():
    import fitsio
    path=get_fits_cat_path()
    print 'reading:',path
    return fitsio.read(path,lower=True)

def get_galsim_catalog():
    import galsim
    path=get_cat_path()
    return galsim.RealGalaxyCatalog(path)
