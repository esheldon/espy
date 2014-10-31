"""

code to work with output files

"""

from __future__ import print_function
from .files_common import *

def get_output_basedir():
    """
    lensdir/output
    """
    d=get_des_lensdir()
    return os.path.join(d, 'output')

def get_output_dir(run):
    """
    lensdir/run/run_name
    """
    d=get_output_basedir()
    return os.path.join(d, run)

def get_output_file(run, lens_chunk, source_tilename):
    """
    the xshear output file for a lens chunk and source tilename
    """
    d=get_output_dir(run)
    fname="%(run)s-lens-%(lens_chunk)06d-src-%(source_tilename)s.dat"
    fname=fname % {'run':run,
                   'lens_chunk':lens_chunk,
                   'source_tilename':source_tilename}

    return os.path.join(d, fname)

def get_reduced_file(run, lens_chunk):
    """
    File reduced across sources
    """
    d=get_output_dir(run)
    fname="%(run)s-lens-%(lens_chunk)06d-reduced.dat"
    fname=fname % {'run':run,
                   'lens_chunk':lens_chunk}

    return os.path.join(d, fname)

def get_combined_file(run):
    """
    File combined for all lenses
    """
    d=get_output_dir(run)
    fname="%s-combined.dat" % run
    return os.path.join(d, fname)

