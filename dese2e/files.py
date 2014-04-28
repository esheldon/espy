from __future__ import print_function
import os
from os import path

def get_input_dir():
    """
    Mike's input directory
    """
    return "/direct/astro+astronfs03/workarea/mjarvis/DES0436-5748"

def get_meds_file():
    """
    My meds file copied into Mike's dir
    """
    d=get_input_dir()
    return path.join(d, "test-end-to-end-cutouts.fits.fz")

def open_meds():
    """
    Return a MEDS object
    """
    import meds
    mfile=get_meds_file()
    return meds.MEDS(mfile)

def get_base_dir():
    """
    The base directory for my stuff
    """
    return path.join(os.environ['LENSDIR'], 'test-end-to-end')

def get_output_dir(**keys):
    """
    parameters
    ----------
    run: keyword, string
        The run identifier
    """
    d=get_base_dir()

    run=keys['run']
    d=path.join(d, "processing", run, "output")
    return d

def get_output_file(**keys):
    """
    parameters
    ----------
    run: keyword, string
        The run identifier
    first: keyword, int, optional
        First object processed in this file
    last: keyword, int, optional
        Last object processed in this file
    """

    d=get_output_dir(**keys)
    if 'last' in keys:
        fname="test-end-to-end-nfit-%(run)s-%(first)06d-%(last)06d.fits"
    else:
        fname="test-end-to-end-nfit-%(run)s.fits"

    fname = fname % keys
    return path.join(d, fname)

def read_output(**keys):
    """
    parameters
    ----------
    run: keyword, string
        The run identifier
    first: keyword, int, optional
        First object processed in this file
    last: keyword, int, optional
        Last object processed in this file
    """
    import fitsio 
    fname=get_output_file(**keys)

    return fitsio.read(fname)

def get_plot_dir(**keys):
    """
    parameters
    ----------
    run: keyword, string
        The run identifier
    """
    d=get_base_dir()

    run=keys['run']
    d=path.join(d, "processing", run, "plots")
    return d


