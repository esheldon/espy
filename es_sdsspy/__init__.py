"""
    Package:
        sdsspy
    Modules:
        See these modules for more documentation.
        files:
            A set of functions for dealing with SDSS files, including reading
            files, finding file lists, getting run metadata, converting SDSS
            id numbers to strings, etc.

        transform:
            A set of transformation functions for use with SDSS coordinate 
            systems.

        window:
            Tools for working with the SDSS window function, as well as getting
            lists of runs that have certain properties.

        sdss_stomp:
            Tools for working with sdss stomp maps.

        yanny:
            tools for working with yanny parameter files.

        flags:
            Get SDSS flag values by name, e.g. target flags, etc.

        sweeps:
            tools for generic processing fo the data sweeps.

    Modification History:
        Erin Sheldon, NYU, BNL
"""

import sys


try:
    from . import transform
except:
    sys.stderr.write('sdsspy.transform not loaded\n')

try:
    from . import files
except:
    sys.stderr.write('sdsspy.files module not loaded\n')

try:
    from . import window
except:
    sys.stderr.write('sdsspy.window module not loaded\n')

try:
    from . import sweeps
except:
    sys.stderr.write('sdsspy.sweeps module not loaded\n')

try:
    import yanny
except:
    sys.stderr.write('sdsspy.yanny module not loaded\n')

try:
    import flags
    from flags import flagval
except:
    sys.stderr.write('sdsspy.flags module not loaded\n')


try:
    import util
    from util import *
except:
    sys.stderr.write('sdsspy.util module not loaded\n')

import stomp_maps
import cas

try:
    import pg
except:
    pass

import target
