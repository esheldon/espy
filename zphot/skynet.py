from __future__ import print_function
import numpy
import esutil as eu
import fitsio

def read_ascii(fname):
    dt=[('coadd_objects_id','i8'),
        ('z_mode','f8'),
        ('z_mean','f8'),
        ('pofz','f8',30)]

    print("reading:",fname)

    with eu.recfile.Open(fname, mode='r', delim=' ', dtype=dt, skiplines=1) as fobj:
        data=fobj.read()
    return data

def get_zbins():
    nz=30
    zbinsize=0.05

    zmin=numpy.zeros(nz)*zbinsize
    zmax=numpy.zeros(nz)*(zbinsize + 1)

    ztbl = numpy.zeros(nz, dtype=[('zmin','f8'),
                                  ('zmax','f8'),
                                  ('zmean','f8')])
    ztbl['zmin'] = zmin
    ztbl['zmax'] = zmax
    ztbl['zmean'] = (zmin+zmax)/2.

    return ztbl

def convert2fits(fname_in, fname_out):

    with fitsio.FITS(fname_out,'rw',clobber=True) as fobj:
        data=read_ascii(fname_in)
        zbins = get_zbins()

        print("writing to:",fname_out)
        fobj.write(data, extname="pofz")
        fobj.write(zbins, extname="zbins")
