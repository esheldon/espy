"""
code to make randoms
"""
from __future__ import print_function
import numpy
from .files import *

class DESRandoms(dict):
    """
    match randoms by z or selection

    example
    -------
    rm = DESRandoms(rconfig)
    ra,dec,z=rm.sample()
    """
    def __init__(self, rconf):
        self['random_name']=rconf
        conf=read_config(rconf)
        self.update(conf)

    def sample(self, nrand):
        """
        sample random points
        """

        ra,dec=self.sample_radec()
        z=self.sample_z()
        return ra,dec,z

    def sample_radec(self):
        """
        genreate the random ra,dec points
        """
        pass

    def sample_z(self):
        """
        sample random z points with constant comoving density
        """
        pass

    def _load_mask(self):
        """
        load the healpix mask
        """
        pass


