"""

Generate random points

"""

from __future__ import print_function
import numpy
from  numpy import pi as PI
import os,sys
from sys import stdout

import lensing

import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1

import es_sdsspy
from es_sdsspy import stomp_maps

import cosmology

def instantiate_sample(sample):
    conf = lensing.files.read_config('rcat',sample)
    if conf['catalog'] == 'sdss':
        return SDSSRandom(sample)
    else:
        raise ValueError("Don't know about catalog: '%s'" % conf['catalog'])


def create_input(sample):
    """
    e.g.  create_input('01')
    """

    c = instantiate_sample(sample)
    c.create_objshear_input()



class SDSSRandom(lensing.lcat.LcatBase):
    def __init__(self, sample, **keys):

        self['lcat_type'] = 'rcat'

        conf = lensing.files.read_config('rcat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['catalog'] not in ['sdss']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

        cconf = lensing.files.read_config('cosmo',self['cosmo_sample'])
        self.cosmo = cosmology.Cosmo(omega_m=cconf['omega_m'], H0=cconf['H0'])

        self.zgen = eu.random.Generator(self.cosmo.dV, 
                                        xrange=[self['zmin'],self['zmax']], 
                                        nx=1000, 
                                        method='cut')
        self['mapname'] = 'boss'
        self['maptype'] = 'basic'
        self.map = es_sdsspy.stomp_maps.load(self['mapname'],self['maptype'])

    def create_objshear_input(self):
        fname = self.file()
        nrand = self['nrand']
        #nrand = 10000

        print("Generating",self['nrand'],"random points "
              "with z in [%0.2f,%0.2f]" % (self['zmin'],self['zmax']))
        n=0

        dt = lensing.files.lcat_dtype()
        output = numpy.zeros(nrand, dtype=dt)
        while n < nrand:
            if n == 0:
                print("  generating",nrand-n," ",end='')
            else:
                print("  re-generating",nrand-n," ",end='')
            print("z",end='')
            z = self.zgen.genrand(nrand-n)
            print(" -> ra,dec ",end='')
            ra,dec = self.map.GenerateRandomEq(nrand-n)

            print(" -> maskflags ", end='')
            maskflags = self.get_maskflags(ra,dec,z)

            wgood = es_sdsspy.stomp_maps.quad_check(maskflags)
            print(" -> good ones:",wgood.size)
            if wgood.size > 0:
                output['zindex'][n:n+wgood.size] = numpy.arange(n,n+wgood.size,dtype='i8')
                output['ra'][n:n+wgood.size] = ra[wgood]
                output['dec'][n:n+wgood.size] = dec[wgood]
                output['z'][n:n+wgood.size] = z[wgood]
                output['maskflags'][n:n+wgood.size] = maskflags[wgood]
                n += wgood.size

        lensing.files.lcat_write(self['sample'], output, type='rcat')

    def get_maskflags(self, ra, dec, z):
        """

        Run the stomp edge checking code.

        """


        # Da is in Mpc
        Da = self.cosmo.Da(0.0, z)

        # radius in *degrees*
        radius = self['rmax']/Da*180./PI
        
        maskflags = self.map.Contains(ra, dec, "eq", radius)

        return numpy.array(maskflags, dtype='i8')

    def plot_coverage(self, region='both', show=True, dops=True):
        """

        Plot a random subset of the randoms along with the
        boss survey geometry area as bounding boxes

        """
        import biggles
        from biggles import FramedPlot, Points, PlotKey

        symsize=0.5

        l = self.read()
        llam,leta = eu.coords.eq2sdss(l['ra'],l['dec'])

        lammin,lammax = (-70.,70.)

        if region=='ngc':
            # a bit high to make room for the legend
            emin,emax=(-40,60)

            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)
            width=1.5

        elif region=='sgc':
            emin,emax=(125,165)
            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)
            width=2
        else:
            emin,emax=(-40,165)
            biggles.configure('screen','width', 1140)
            biggles.configure('screen','height', 1140)
            width=2

        wl=where1((leta > emin) & (leta < emax))
        llam=llam[wl]
        leta=leta[wl]


        plt=FramedPlot()
        plt.xlabel=r'$\lambda$'
        plt.ylabel=r'$\eta$'
        xrng = (lammin, lammax)
        yrng = (emin, emax)
        plt.xrange = xrng
        plt.yrange = yrng


        print("adding random subset of randoms")

        ii = eu.numpy_util.random_subset(llam.size, 500000)
        allp = Points(llam[ii],leta[ii],type='dot',size=symsize)
        plt.add(allp)

        plt.aspect_ratio = (yrng[1]-yrng[0])/float(xrng[1]-xrng[0])

        #es_sdsspy.stomp_maps.plot_boss_geometry(color='blue',plt=plt,show=False)
        es_sdsspy.stomp_maps.plot_boss_geometry(plt=plt,show=False,width=width)

        if show:
            plt.show()

        if dops:
            d = lensing.files.sample_dir(self['lcat_type'],self['sample'])
            d = os.path.join(d,'plots')
            if not os.path.exists(d):
                os.makedirs(d)
            epsfile = os.path.join(d, '%s-%s-coverage.eps' % (self['lcat_type'],self['sample']))
            print("Writing to eps file:",epsfile)
            plt.write_eps(epsfile)
        return plt


