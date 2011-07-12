"""

Classes for each catalog type, e.g. desmocks

Also create_input() function to make the input files for objshear

"""
import numpy
import os,sys
from sys import stdout

import lensing

import esutil
from esutil.ostools import path_join

def instantiate_sample(sample):
    if sample in ['04']:
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

    def create_objshear_input(self):
        fname = self.file()

        data = self.read_original()

        print 'creating output array'
        output = output_array(data.size)

        print 'copying data'
        output['zindex'] = numpy.arange(data.size,dtype='i8')
        output['ra'] = data['ra']
        output['dec'] = data['dec']
        output['z'] = data['photoz_cts']

        lensing.files.lcat_write(self['sample'], output)




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


