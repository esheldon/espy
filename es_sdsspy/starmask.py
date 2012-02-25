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
        cont = ones(ra.size, dtype='u1')
        self.load_data()

        h=eu.htm.HTM(self.depth)

        mstars, minput, dis = h.match(self.data['ra'],self.data['dec'],
                                      ra,dec,
                                      self.radii,
                                      maxmatch=0)

        if minput.size > 0:
            cont[minput] = 0

        return cont

    def load_data(self):
        if self.data is None:
            dir=os.getenv('MASK_DIR')
            if dir is None:
                raise ValueError("MASK_DIR is not defined")
            dir = os.path.join(dir, 'stomp-sdss')
            f=os.path.join(dir,'bright_star_mask.fits')

            print 'loading star defs:',f
            self.data = eu.io.read(f, lower=True)

            print 'setting radii'
            self.radii = rad2deg(arccos(1.-self.data['cmcaps']))

def compare_speed(n=100000):
    from . import sweeps_collate
    c=sweeps_collate.open_columns('primgal')

    ingood=c['ingood_8192'][:]

    w,=where(ingood==1)

    ra=c['ra'][w[0:n]]
    dec=c['dec'][w[0:n]]

    m=StarMask()
    m.load_data()
    for depth in [5,6,7,8,9,10]:
        m.depth=depth
        print 'depth:',depth,' ',
        cont=m.contains(ra,dec)

def compare_stomp(maxres=8192,n=100000):
    from . import sweeps_collate

    c=sweeps_collate.open_columns('primgal')

    ingood=c['ingood_%s' % maxres][:]

    w,=where(ingood==1)

    ra=c['ra'][w[0:n]]
    dec=c['dec'][w[0:n]]
    intycho_stomp=c['intycho_%s' % maxres][w[0:n]]

    m=StarMask()
    cont=m.contains(ra,dec)

    stomp_good,=where(intycho_stomp==1)
    this_good,=where(cont==1)

    print 'stomp:',stomp_good.size
    print 'this: ',this_good.size
