"""

Since the bright star mask is simple disks centered on the stars,
it is faster to use a point matching code than a mangle mask.

"""
import os,sys
import time
import esutil as eu
from numpy import rad2deg, arccos, ones, where

class StarMask:
    """
    This is very simple; just use htm_match to match points to the star locations
    """
    def __init__(self, depth=8):
        self.data=None
        self.radii=None
        self.depth=depth

    def contains(self, ra, dec):
        """
        Returns 1 if the point did not match a bright star.
        """
        cont = ones(ra.size, dtype='u1')
        self.load_data()

        h=eu.htm.HTM(self.depth)

        mstars, minput, dis = h.match(self.data['ra'],self.data['dec'],
                                      ra,dec,
                                      self.radii,
                                      maxmatch=0)

        # it is a veto situation, so matches are *not* contained
        if minput.size > 0:
            cont[minput] = 0

        return cont

    def load_data(self):
        if self.data is None:
            dir=os.getenv('MASK_DIR')
            if dir is None:
                raise ValueError("MASK_DIR is not defined")
            dir = os.path.join(dir, 'mangle')
            f=os.path.join(dir,'bright_star_mask.fits')

            #print 'loading star defs:',f
            self.data = eu.io.read(f, lower=True)

            #print 'setting radii'
            self.radii = rad2deg(arccos(1.-self.data['cmcaps']))


