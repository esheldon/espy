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

from . import sg

from . import select
from . import transform
from . import files
from . import window
from . import sweeps
from . import sweeps_collate
from . import starmask
from . import yanny

from . import flags
from .flags import flagval


from .import util
from .util import *

from . import stomp_maps
from . import cas

from . import pg

from . import target
