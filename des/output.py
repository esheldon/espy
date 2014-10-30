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
    the yaml wq file for a source tilename and lens chunk
    """
    d=get_output_dir(run)
    fname="%(run)s-lens-%(lens_chunk)06d-src-%(source_tilename)s.dat"
    fname=fname % {'run':run,
                   'lens_chunk':lens_chunk,
                   'source_tilename':source_tilename}

    return os.path.join(d, fname)


