from __future__ import print_function
import esutil as eu
from numpy import where, zeros
import es_sdsspy
import os

from . import files

class Random:
    def __init__(self, name, version, nsamples, nper):
        """
        name: e.g. rm for redmapper
        version: e.g. dr8-v3
        """

        self.name=name
        self.version=version
        self.nsamples = nsamples
        self.nper=nper

        self.good_mask= es_sdsspy.mangle_masks.load('boss','good')

        # veto masks
        self.star_mask = es_sdsspy.starmask.StarMask()
        self.badfield_mask = es_sdsspy.mangle_masks.load('boss','badfield')

    def make_randoms(self, nrand):
        """
        Generate randoms in the bood mask then check the veto masks
        """

        data = self.struct(nrand)

        print("Generating points in boss good mask")
        ra,dec = self.good_mask.genrand(nrand)


        print("checking against star mask")
        cstar = self.star_mask.contains(ra,dec)

        print("checking against badfield mask")
        cbadfield = self.badfield_mask.contains(ra,dec)

        data['ra'] = ra
        data['dec'] = dec
        data['instar'] = cstar.astype('i1')
        data['inbadfield'] = cbadfield.astype('i1')

        return data

    def plot_randoms(self, data):

        plt=eu.plotting.bscatter(data['ra'],data['dec'],
                                 type='dot',show=False)

        ws, = where(data['cstar']==0)
        wb, = where(data['cbadfield']==1)
        plt2 = eu.plotting.bscatter(data['ra'][ws], 
                                    data['dec'][ws], 
                                    show=False, type='filled circle', 
                                    color='red', plt=plt)
        plt3 = eu.plotting.bscatter(data['ra'][wb], 
                                    data['dec'][wb], 
                                    show=False, type='filled diamond', 
                                    color='blue', plt=plt2)

        plt3.show()

    def write_randoms(self):
        """
        Write nsample files with nper randoms each
        """

        d= files.random_dir(self.name, self.version)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        for i in xrange(self.nsamples):
            f=files.random_file(self.name, self.version, i)

            data = self.make_randoms(self.nper)
            
            print(f)
            eu.io.write(f, data, clobber=True)

    def struct(self, n):
        dt = [('ra','f8'),('dec','f8'),
              ('instar','i1'),('inbadfield','i1')]
        data = zeros(n, dtype=dt)
        return data
