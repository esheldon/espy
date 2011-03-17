import numpy
import os,sys
from sys import stdout

import esutil
from esutil.ostools import path_join

import lensing

def create_input(sample):
    """
        create_input('01')
    """

    c = DESMockSrcCatalog(sample)
    c.create_objshear_input()


class DESMockSrcCatalog(dict):
    """
    This reads the mock catalog and creates an input
    catalog for objshear


    Need to have config files that tell us what kind
    of sigma crit inv we are using.
    """

    def __init__(self, sample):
        conf = lensing.files.json_read('scat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['catalog'] not in ['desmocks-2.13o-f8']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

        self['cosmo'] = lensing.files.json_read('cosmo',self['cosmo_sample'])

    def file(self):
        fname = lensing.files.sample_file('scat',self['sample'])
        return fname
    def read(self):
        return lensing.files.scat_read(sample=self['sample'])

    def create_objshear_input(self):
        outfile = self.file()

        data = self.read_original()

        print 'creating output array'
        output = self.output_array(data.size)

        print 'copying data'
        output['ra'] = data['ora']
        output['dec'] = data['odec']

        print 'not changing sign of gamma1'
        output['g1'] = data['gamma1']

        output['g2'] = data['gamma2']
        output['err'] = 0.0
        output['z'] = data['z']
        output['dc'] = -9999
        output['hpixid'] = -9999

        lensing.files.scat_write(output, file=outfile)


    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return d

    def original_file(self):
        d = self.original_dir()
        f='%s-sources.rec' % self['catalog']
        infile = path_join(d, f)
        return infile

    def read_original(self):
        infile = self.original_file()
        if not os.path.exists(infile):
            raise ValueError("File not found: %s\n" % infile)

        data = esutil.io.read(infile, lower=True, verbose=True)
        return data

    def output_array(self, num):
        dt = lensing.files.scat_dtype()
        output = numpy.zeros(num, dtype=dt)
        return output


