import os
import sys
from sys import stdout

import columns
import numpy
from numpy import where
import esutil
from esutil.ostools import path_join, expand_path, getenv_check

import sdsspy
import sdssgal

# select the maxbcg input data from a columns database
class MaxbcgColumnSelector:
    """

    mcs = MaxbcgColumnSelector('prim03')
    mcs.select()
    mcs.write('fits')
    mcs.write('rec')

    """
    def __init__(self, procrun, rmax=22.0, imax=21.0):
        """
        E.g.  procrun='prim03'
        """

        self.basedir = expand_path('~/lensinputs/maxbcg_input')
        self.prefix='maxbcg'
        self.proctype='sdssgal'
        self.procrun=procrun
        self.rmax=rmax
        self.imax=imax

        self.masktype = 'tycho'

        self.cols = sdssgal.open_columns(procrun)

        self.keep_indices=None

    def write(self, filetype):
        if self.keep_indices is None:
            raise ValueError("run select() first")

        dir = self.dir()
        if not os.path.exists(dir):
            os.makedirs(dir)

        outfile=self.filename(filetype)
        stdout.write("Will write to file: %s\n" % outfile)

        data = self.read_input_columns()

        if filetype == 'fits':
            stdout.write("Converting to big endian in place\n")
            esutil.numpy_util.to_big_endian(data, inplace=True)

        stdout.write("Writing to file: %s\n" % outfile)
        stdout.flush()

        # clobber is actually the default for sfile, so
        # the keyword will be ignored
        esutil.io.write(outfile, data, clobber=True)

        del data



    def select(self):
        # for shorthand notation
        c = self.cols

        stdout.write("Reading flags\n")
        flags = c['flags'][:]
        stdout.write('Getting binned logic\n')
        binned_logic = sdssgal.binned_logic_ri(flags)
        w,=where(binned_logic)
        stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
        del flags
        
        stdout.write("Reading objc_flags\n")
        objc_flags = c['objc_flags'][:]
        stdout.write('Getting objc_flags1 logic\n')
        objc_flags1_logic = sdssgal.object1_logic(objc_flags)
        w,=where(objc_flags1_logic)
        stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
        del objc_flags

        stdout.write("Reading cmodelmag_dered\n")
        mag = c['cmodelmag_dered'][:]
        stdout.write('Getting mag logic\n')
        mag_logic = (mag[:,2] < self.rmax) & (mag[:,3] < self.imax)
        w,=where(mag_logic)
        stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
        del mag


        mask_colname = 'in'+self.masktype
        stdout.write("Reading %s (and good)\n" % mask_colname)
        inmask = c[mask_colname][:]
        stdout.write("Getting %s mask logic\n" % self.masktype)
        mask_logic = (inmask == 1)
        w,=where(mask_logic)
        stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
        del inmask

        logic = binned_logic & objc_flags1_logic & mag_logic & mask_logic
        
        w,=where(logic)
        stdout.write("A total of %s/%s passed\n" % (w.size,c['flags'].size))

        self.keep_indices = w

    def dir(self):
        outdir=path_join(self.basedir, self.proctype+'-'+self.procrun)
        return outdir

    def filename(self, filetype):
        if filetype == 'fits':
            ext='fits'
        elif filetype == 'rec':
            ext='rec'
        else:
            raise ValueError("Choose filetype 'fits' or 'rec'")

        outdir = self.dir()

        printpars=(self.proctype,
                   self.procrun,
                   self.masktype,
                   self.rmax,
                   self.imax,
                   ext)
        outfile=self.prefix+'-%s-%s-%s-r%0.1f-i%0.1f.%s' % printpars
        outfile = path_join(outdir,outfile) 
        return outfile

    def read_input_columns(self):

        if self.keep_indices is None:
            raise ValueError("run select() first")

        colnames = ['photoid',
                    'ra',
                    'dec',
                    'modelmag_dered',
                    'modelmag_dered_err',
                    'cmodelmag_dered',
                    'cmodelmag_dered_err']

        stdout.write("Reading columns: %s\n" % colnames)
        data = self.cols.read_columns(colnames, 
                                      rows=self.keep_indices, 
                                      verbose=True)

        return data



def create_input_old(imax=21.0, rec=False):
    rmax=22.0
    proctype='sdssgal'
    procrun='03'
    masktype='tycho'

    if rec:
        ext='rec'
    else:
        ext='fits'

    outdir=expand_path('~/lensinputs/maxbcg_input')
    outdir=path_join(outdir, proctype+'-'+procrun)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outfile='maxbcg-%s-%s-%s-r%0.1f-i%0.1f.%s' % \
            (proctype,procrun,masktype,rmax,imax,ext)
    outfile = path_join(outdir,outfile) 

    stdout.write("Will write to file: %s\n" % outfile)

    c=sdssgal.open_columns(procrun)

    # get good galaxies


    stdout.write("Reading intycho, cmodelmag, extinction\n")
    intycho = c['intycho'][:]
    cmodelmag = c['cmodelmag'][:]
    ext = c['extinction'][:]
    imag = cmodelmag[:,3] - ext[:,3]
    del cmodelmag
    del ext


    # select things that pass the bound mask
    stdout.write("Getting unmasked objects\n")
    w, = where(intycho == 1)

    stdout.write("Found %s/%s passed mask\n" % (w.size,intycho.size))

    stdout.write("Getting dereddened i mag < %s\n" % imax)
    w, = where( (intycho == 1) & (imag < imax) )

    stdout.write("Found %s/%s passed mask and flux cut\n" % (w.size,intycho.size))

    del intycho
    del imag



    # now read in all objects passing the mask

    colnames = ['photoid',
                'ra',
                'dec',
                'modelmag',
                'modelmag_err',
                'cmodelmag',
                'cmodelmag_err']

    newnames= ['photoid',
               'ra',
               'dec',
               'modelmag_dered',
               'modelmag_dered_err',
               'cmodelmag_dered',
               'cmodelmag_dered_err']


    stdout.write("Reading columns: %s\n" % colnames)
    data = c.read_columns(colnames, rows=w, verbose=True)

    # get extinction separately
    stdout.write("Reading extinction\n")
    ext = c['extinction'].read(rows=w)
    stdout.write("Applying extinction correction\n")
    data['modelmag'] -= ext
    data['cmodelmag'] -= ext

    del ext
    del w

    # use the _dered name
    stdout.write("Renaming columns: %s\n" % newnames)
    data.dtype.names = newnames

    if not rec:
        stdout.write("Converting to big endian in place\n")
        esutil.numpy_util.to_big_endian(data, inplace=True)

    stdout.write("Writing to file: %s\n" % outfile)
    stdout.flush()
    if rec:
        esutil.sfile.write(data, outfile)
    else:
        esutil.io.write(outfile, data, clobber=True)

    del data

