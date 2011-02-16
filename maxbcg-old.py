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



def _create_input_select(c, rmax, imax):
    stdout.write("Reading flags\n")
    flags = c['flags'][:]
    stdout.write('Getting binned logic\n')
    binned_logic = sdssgal.binned_logic_ri(flags)
    w,=where(binned_logic)
    stdout.write("Found %s/%s passed\n" % (w.size,c['flags'].size))
    del flags
    
    stdout.write("Reading objc_flags\n")
    objc_flags = c['objc_flags'][:]
    stdout.write('Getting objc_flags1 logic\n')
    objc_flags1_logic = sdssgal.object1_logic(objc_flags)
    w,=where(objc_flags1_logic)
    stdout.write("Found %s/%s passed\n" % (w.size,c['flags'].size))
    del objc_flags

    stdout.write("Reading cmodelmag, extinction\n")
    cmodelmag = c['cmodelmag'][:]
    ext = c['extinction'][:]
    mag = cmodelmag - ext
    del ext
    stdout.write('Getting mag logic\n')
    mag_logic = (mag[:,2] < rmax) & (mag[:,3] < imax)
    w,=where(mag_logic)
    stdout.write("Found %s/%s passed\n" % (w.size,c['flags'].size))
    del mag


    stdout.write("Reading intycho (and good)\n")
    intycho = c['intycho'][:]
    stdout.write("Getting tycho mask logic\n")
    tycho_logic = (intycho == 1)
    w,=where(tycho_logic)
    stdout.write("Found %s/%s passed\n" % (w.size,c['flags'].size))
    del intycho

    logic = binned_logic & objc_flags1_logic & mag_logic & tycho_logic
    
    w,=where(logic)
    stdout.write("Found a total of %s/%s passed\n" % (w.size,c['flags'].size))

    return w


def create_input(imax=21.0, rec=False):
    rmax=22.0
    proctype='sdssgal'
    procrun='prim03'
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

    w=_create_input_select(c,rmax,imax)

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

