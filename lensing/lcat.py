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

def read_catalog(catalog, version):
    if catalog == 'desmocks':
        c = DESMockLensCatalog(version)
        return c.read_original()
    else:
        raise ValueError("don't know about catalog: '%s'" % catalog)



def create_input(catalog, version, sample, nsplit=0):
    """
    e.g.
        create_input('desmocks','2.13','desmocks-2.13')
        create_input('desmocks','2.13','desmocks-2.13-redo')

        create_input('desmocks','2.13','desmocks-2.13-redo', nsplit=100)
    """
    scc = lensing.convert.CatalogConverter('lcat')

    if catalog == 'desmocks':
        scc.convert(DESMockLensCatalog, version, sample)
        if nsplit != 0:
            fname = scc.sample_file(sample)
            split_cat(fname, nsplit)
    else:
        raise ValueError("don't know about catalog: '%s'" % catalog)

def split_cat(fname, nsplit):
    """
    Split the lens file into nsplit parts
    """

    data = lensing.files.lcat_read(file=fname)

    ntot = data.size
    nper = ntot/nsplit
    nleft = ntot % nsplit


    for i in xrange(nsplit):
        sstr = '%03d' % i
        beg = i*nper
        end = (i+1)*nper
        if i == (nsplit-1):
            end += nleft
        sdata = data[beg:end]
        sfile = fname.replace('.bin','-'+sstr+'.bin') 

        lensing.files.lcat_write(sdata, file=sfile)



def read_mocks(version):
    dm=DESMockLensCatalog(version)
    return dm.read_original()

class DESMockLensCatalog:
    """
    Provides the interface needed by CatalogConverter
    """

    def __init__(self, version, **keys):
        self.version = version
        self.lam = keys.get('lam',False)

    def create_objshear_input(self, filename):
        data = self.read_original()

        print 'creating output array'
        output = self.output_array(data.size)

        print 'copying data'
        output['zindex'] = numpy.arange(data.size,dtype='i4')
        output['ra'] = data['ra']
        output['dec'] = data['dec']
        output['z'] = data['z']
        output['dc'] = -9999.0
        output['padding'] = -9999

        lensing.files.lcat_write(output, file=filename)

    def type(self):
        return 'desmocks'

    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self.type()+'-'+self.version)
        return d

    def original_file(self):
        d = self.original_dir()
        extra=''
        if self.lam:
            extra='_lambda'
        f='DES_Mock_v{version}_halos{extra}.fit'.format(version=self.version,
                                                        extra=extra)
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

    def output_array(self, num):
        dt = lensing.files.lcat_dtype()
        output = numpy.zeros(num, dtype=dt)
        return output



