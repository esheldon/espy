"""
code to make randoms
"""
from __future__ import print_function
import numpy
from .files import *

def make_randoms(rconf_name, nrand, cat_name):
    """
    write a random catalog

    parameters
    ----------
    rconf_name: string
        name of randoms config, e.g. "randoms-sva1-spte-neff"
    nrand: int
        number of randoms to write
    cat_name:
        name of catalog to write
    """
    import fitsio

    r=DESRandoms(rconf_name)

    print("generating",nrand,"randoms")
    ra,dec,z=r.sample(nrand)

    data=numpy.zeros(nrand, dtype=[('rand_index','i8'),
                                   ('ra','f8'),
                                   ('dec','f8'),
                                   ('z','f8')])
    data['rand_index']=numpy.arange(nrand)
    data['ra']=ra
    data['dec']=dec
    data['z']=z

    d=get_cat_dir(cat_name)
    fname=get_lcat_original_file(cat_name)
    if not os.path.exists(d):
        print("making dir:",d)
        os.makedirs(d)

    print("writing:",fname)
    fitsio.write(fname,data,clobber=True)

class DESRandoms(dict):
    """
    match randoms by z or selection

    example
    -------
    rm = DESRandoms(rconfig)
    ra,dec,z=rm.sample()
    """
    def __init__(self, rconf_name):

        self._load_config(rconf_name)
        self._load_cosmo()
        self._load_mask()
        self._load_z_generator()


    def sample(self, nrand):
        """
        sample random points
        """

        ra,dec=self.sample_radec(nrand)
        z=self.sample_z(nrand)
        return ra,dec,z

    def sample_radec(self, nrand):
        """
        genreate the random ra,dec points
        """
        ra_range=self['limits']['ra_range']
        dec_range=self['limits']['dec_range']
        ra,dec=self.hmap.genrand(nrand,
                                 ra_range=ra_range,
                                 dec_range=dec_range)
        return ra,dec

    def sample_z(self, nrand):
        """
        sample random z points with constant comoving density
        """

        z=self.zgen.genrand(nrand)
        return z

    def _load_config(self, rconf_name):
        """
        load the various configuration files
        """

        self['random_name']=rconf_name
        conf=read_config(rconf_name)
        self.update(conf)

        c=read_config('constants')
        self['limits']=c['regions'][self['region']]

    def _load_mask(self):
        """
        load the healpix mask
        """
        import healpix_util as hu
        print("loading mask:",self['mask_file'])
        self.hmap=hu.readDensityMap(self['mask_file'])


    def _load_cosmo(self):
        """
        load cosmology calculator
        """
        import cosmology
        cconf=read_config(self['cosmo_vers'])
        self.cosmo = cosmology.Cosmo(omega_m=cconf['omega_m'],
                                     H0=cconf['H0'])

    def _load_z_generator(self):
        """
        load the random z generator
        """
        import esutil as eu
        self.zgen = eu.random.Generator(self.cosmo.dV, 
                                        xrange=[self['zmin'],self['zmax']], 
                                        nx=1000, 
                                        method='cut')

