"""

Classes for each catalog type, e.g. desmocks

Also create_input() function to make the input files for objshear

"""
import numpy
from  numpy import pi as PI
import os,sys
from sys import stdout

import lensing

import esutil
from esutil.ostools import path_join
from esutil.numpy_util import where1

from es_sdsspy import stomp_maps

import cosmology

def instantiate_sample(sample):
    conf = lensing.files.read_config('lcat',sample)
    if conf['catalog'] == 'maxbcg-full':
        return MaxBCG(sample)
    else:
        return DESMockLensCatalog(sample)

def create_input(sample):
    """
    e.g.  create_input('01')
    """

    c = instantiate_sample(sample)
    c.create_objshear_input()

def original_file(sample):
    c = instantiate_sample(sample)
    return c.original_file()

def read_original(sample):
    c = instantiate_sample(sample)
    return c.read_original()


def output_array(num):
    dt = lensing.files.lcat_dtype()
    output = numpy.zeros(num, dtype=dt)
    return output




class MaxBCG(dict):
    def __init__(self, sample, **keys):
        conf = lensing.files.read_config('lcat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['catalog'] not in ['maxbcg-full']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

        self['min_ngals_r200'] = 3

    def create_objshear_input(self):
        fname = self.file()

        data = self.read_original()

        # trim not well understood low ngals stuff
        good_ngals = self.ngals_cuts(data['ngals_r200'])
        data = data[good_ngals]

        # trim z for speed
        good_z = self.zcuts(data['photoz_cts'])
        data = data[good_z]

        # for now removing all edge problems
        good_edge = self.edge_cuts(data['ra'], data['dec'], data['photoz_cts'])
        data = data[good_edge]

        print 'creating output array'
        output = output_array(data.size)

        print 'copying data'
        output['zindex'] = numpy.arange(data.size,dtype='i8')
        output['ra'] = data['ra']
        output['dec'] = data['dec']
        output['z'] = data['photoz_cts']

        lensing.files.lcat_write(self['sample'], output)

    def ngals_cuts(self, ngals):
        print("Cutting ngals >= %d" % self['min_ngals_r200'])
        w=where1(ngals >= self['min_ngals_r200'])

        print("Keeping %d/%d" % (w.size,ngals.size))
        if w.size == 0:
            raise ValueError("No objects passed z cut")

        return w


    def zcuts(self, z):
        print("Cutting z to [%f, %f]" % (self['zmin'],self['zmax']))
        w=where1( (z > self['zmin']) & (z < self['zmax']) )

        print("Keeping %d/%d" % (w.size,z.size))
        if w.size == 0:
            raise ValueError("No objects passed z cut")

        return w

    def edge_cuts(self, ra, dec, z):
        """

        This will remove any objects that intersect an edge, so
        this is not as nice as doing the quad check.

        """

        print("Doing edge check")
        # get radius for edge check
        cconf = lensing.files.read_config('cosmo',self['cosmo_sample'])
        print(cconf)
        
        c = cosmology.Cosmo(H0=cconf['H0'], omega_m=cconf['omega_m'])

        print("Getting radius")
        # Da is in Mpc
        Da = c.Da(0.0, z)

        # radius in *degrees*
        radius = self['rmax']/Da*180./PI

        print("radii are in range [%f,%f]" % (radius.min(), radius.max()))

        map = stomp_maps.load('boss','basic')
        
        print("Getting maskflags...")
        maskflags = map.Contains(ra, dec, "eq", radius)

        good = where1(maskflags == 0)
        print("Keeping %d/%d" % (good.size,ra.size))

        if good.size == 0:
            raise ValueError("No objects passed edge cut")

        #print("Getting lenses with adjacent quadrants")
        #good = stomp_maps.quad_check(maskflags)

        return good

    def file(self):
        fname = lensing.files.sample_file('lcat',self['sample'])
        return fname

    def read(self, split=None):
        return lensing.files.lcat_read(sample=self['sample'], split=split)


    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return str(d)

    def original_file(self):
        d = self.original_dir()
        f='catalog_full_bcgs_orig.fit'
        infile = path_join(d, f)
        return infile

    def read_original(self):
        infile = self.original_file()
        stdout.write("Reading lens catalog: %s\n" % infile)
        data = esutil.io.read(infile, lower=True, ensure_native=True)
        return data



class DESMockLensCatalog(dict):
    """
    Provides the interface needed by CatalogConverter
    """

    def __init__(self, sample, **keys):
        conf = lensing.files.read_config('lcat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['catalog'] not in ['desmocks-2.13']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

    def file(self):
        fname = lensing.files.sample_file('lcat',self['sample'])
        return fname

    def read(self, split=None):
        return lensing.files.lcat_read(sample=self['sample'], split=split)

    def create_objshear_input(self):
        fname = self.file()

        data = self.read_original()

        print 'creating output array'
        output = output_array(data.size)

        print 'copying data'
        output['zindex'] = numpy.arange(data.size,dtype='i8')
        output['ra'] = data['ra']
        output['dec'] = data['dec']
        output['z'] = data['z']
        #output['dc'] = -9999.0

        lensing.files.lcat_write(self['sample'], output)


    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return d

    def original_file(self):
        d = self.original_dir()
        f='%s-halos.fit' % self['catalog']
        infile = path_join(d, f)
        return infile

    def read_original(self):
        infile = self.original_file()
        if not os.path.exists(infile):
            raise ValueError("File not found: %s\n" % infile)

        stdout.write("Reading lens catalog: %s\n" % infile)
        data = esutil.io.read(infile, lower=True, verbose=True, 
                              ensure_native=True)
        return data


