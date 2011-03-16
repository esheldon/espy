import numpy
import os,sys
from sys import stdout

import esutil
from esutil.ostools import path_join

import lensing

def create_input(run):
    """
        create_input('05')
    """

    conf = lensing.files.json_read(run)
    cat=conf['src_catalog']
    version=conf['src_version']
    sample=conf['src_sample']

    scc = lensing.convert.CatalogConverter('scat')

    if catalog == 'desmocks':
        scc.convert(DESMockSrcCatalog, version, sample)
    else:
        raise ValueError("don't know about catalog: '%s'" % catalog)




def read_mocks(version):
    dc = DESMockSrcCatalog(version)
    return dc.read_original()
class DESMockSrcCatalog:
    """
    This reads the mock catalog and creates an input
    catalog for objshear


    Need to have config files that tell us what kind
    of sigma crit inv we are using.
    """

    def __init__(self, version):

        self.version = version

    def create_objshear_input(self, outfile):
        data = self.read_original()

        print 'creating output array'
        output = self.output_array(data.size)

        print 'copying data'
        output['ra'] = data['ora']
        output['dec'] = data['odec']
        if self.version == '2.13':
            print 'changing sign of gamma1'
            output['g1'] = -data['gamma1']
        elif self.version in ['2.13o','2.13o-f8']:
            print 'not changing sign of gamma1'
            output['g1'] = data['gamma1']
        else:
            raise ValueError("Dont' know shear convention "
                             "for version '%s'" % self.version)

        output['g2'] = data['gamma2']
        output['err'] = 0.0
        output['z'] = data['z']
        output['dc'] = -9999
        output['hpixid'] = -9999

        lensing.files.scat_write(output, file=outfile)

    def type(self):
        return 'desmocks'

    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self.type()+'-'+self.version)
        return d

    def original_file(self):
        d = self.original_dir()
        extra=''
        if self.version == '2.13':
            extra='-fix'
            ext='.fits'
        elif self.version == '2.13o':
            extra='-origshear'
            ext='.rec'
        elif self.version == '2.13o-f8':
            extra='-origshear'
            ext='.rec'
        else:
            ext='.fits'
        f = 'galaxies-full-{version}{extra}{ext}'.format(version=self.version,
                                                         extra=extra,ext=ext)
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


