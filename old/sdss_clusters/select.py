from __future__ import print_function
import os
import sys

import columns
import numpy
from numpy import where, zeros
import esutil
from esutil.ostools import path_join, expand_path, getenv_check
from esutil.numpy_util import where1, to_big_endian

import es_sdsspy
import sdssgal

from . import files

from columns import Columns

def do_radcut(cat, rmpc=1.3):
    """
    apply the the BOSS basic mask with quadrant cuts
    """
    
    map = es_sdsspy.stomp_maps.load('boss','basic')

    c = esutil.cosmology.Cosmo()
    distmpc = c.Da(0.0, cat['z'])

    # radius must be in degrees
    radius=rmpc/distmpc*180.0/numpy.pi
    print("getting contained withing %0.2f Mpc" % rmpc)
    maskflags = map.Contains(cat['ra'], cat['dec'], "eq", radius)

    # we want to make sure all quadrants are contained
    w=es_sdsspy.stomp_maps.quad_check(maskflags, strict=True)
    print("Kept %s/%s" % (w.size,cat.size))
    return w

# select the input data from a columns database
# and into a columns database
class Selector:
    """
    mcs = Selector()
    mcs.select()
    mcs.write_columns()
    """
    def __init__(self, name, version, type='primgal', rmax=22.0):
        """
        E.g.  Selector('rm','dr8-v2')
        """

        self.version = version
        self.name=name
        self.rmax=rmax

        self.cols = es_sdsspy.sweeps_collate.open_columns(type)

        self.logic = None
        self.mag_and_mask_logic = None

    def write_columns(self):

        d = files.input_coldir(self.name, self.version)
        if os.path.exists(d):
            raise ValueError("coldir exists, start fresh: '%s'" % d)
        outcols = Columns(d)
        outcols.create()

        if self.logic is None:
            self.select()

        w=where1(self.mag_and_mask_logic)
        print("\nkeeping those that pass mag and mask logic: %s/%s" % (w.size, self.logic.size))


        colnames = ['photoid',
                    'ra',
                    'dec',
                    'objc_flags',
                    'objc_flags2',
                    'ingood',
                    'instar',
                    'inbadfield']

        for col in colnames:
            print('Creating:',col)
            data = self.cols[col][:]
            data = data[w]
            outcols.write_column(col, data)

            del data


        combcols = ['flags','flags2',
                    'modelflux','modelflux_ivar',
                    'cmodelflux','cmodelflux_ivar',
                    'extinction']

        for ccol in combcols:
            print('Creating:',ccol)
            colnames = [ccol+'_'+f for f in ['u','g','r','i','z'] ]
            data = self.cols.read_columns(colnames, rows=w, verbose=True)

            dts = data[ccol+'_'+f].dtype.descr[0][1]
            dt = [(ccol, dts, 5)]
            data = data.view(dt)

            # don't want it to show up as a .rec
            rawdata = data[ccol]
            outcols.write_column(ccol, rawdata)

            del data


        print('Adding more restrictive flags to "keep"')
        keep = zeros(w.size, dtype='u1')
        wrest = where1( self.logic[w] )
        keep[wrest] = 1

        outcols.write_column('keep', keep)


    def add_inbadfield(self):
        import es_sdsspy
        m=es_sdsspy.mangle_masks.load('boss','badfield')
        d = files.input_coldir(self.name, self.version)
        c = Columns(d)

        print("reading ra")
        ra=c['ra'][:]
        print("reading dec")
        dec=c['dec'][:]
 
        print("checking mask")
        cont = m.contains(ra,dec).astype('i1')

        c.write_column('inbadfield', cont)

    def add_colc_nmgypercount(self):
        """
        self.cols is the sweeps
        """
        import es_sdsspy
        d = files.input_coldir(self.name, self.version)

        ccols = Columns(d)
        print("cluster input cols:",ccols.dir)
        scols = self.cols
        print("sweep columns:",scols.dir)


        print("reading sweep photoid")
        spid = scols['photoid'][:]
        print("reading cluster input photoid")
        cpid = ccols['photoid'][:]

        print("matching")
        ms,mc = esutil.numpy_util.match(spid, cpid)

        n=cpid.size
        if mc.size != n:
            raise ValueError("not all matched: %d/%d" % (mt.size/n))

        # problem is the columns code always sorts the requested
        # indices so we will have lost order information.  Need
        # to pre-sort both sets of indices according to the sweep
        # row number
        print("ordering matches with sweep cols")
        s=ms.argsort()
        ms = ms[s]
        mc = mc[s]

        for colname in ['colc','nmgypercount']:
            data=numpy.zeros((n,5),dtype='f4')
            for i,band in enumerate(['u','g','r','i','z']):

                scolname='%s_%s' % (colname,band)
                print("    %s" % scolname)

                sdata=scols[scolname][ms]
                data[mc,i] = sdata

                del sdata

            print("writing: %s" % colname)
            ccols.write_column(colname, data)
            del data



    def select(self):
        # for shorthand notation
        c = self.cols

        mask_colname = 'inbasic'
        print("Reading mask info:",mask_colname)
        inmask = c[mask_colname][:]
        ntot = inmask.size

        print("Getting mask logic")
        mask_logic = (inmask == 1)
        w=where1(mask_logic)
        print("    %s/%s passed" % (w.size,ntot))


        # now read columns needed for binned_logic, object1_logic
        print("Reading flags, objc_flags")
        colnames=['objc_flags']
        colnames += ['flags_'+f  for f in ['u','g','r','i','z'] ]
        data = c.read_columns(colnames)

        selector = es_sdsspy.select.Selector(data)

        print("binned logic")
        binned_logic = selector.binned_logic()

        w=where1(binned_logic)
        print("    %s/%s passed" % (w.size,ntot))

        print("object1 logic")
        object1_logic = selector.object1_logic()
        
        w=where1(object1_logic)


        print("Reading cmodelmag_dered_r")
        cmag_r = c['cmodelmag_dered_r'][:]
        print('mag logic')
        cmag_r_logic = (cmag_r < self.rmax)
        w=where1(mask_logic & cmag_r_logic)
        print("    %s/%s passed" % (w.size,ntot))

        w=where1(cmag_r_logic)
        print("    %s/%s passed mag+mask" % (w.size,ntot))

        logic = binned_logic & object1_logic & cmag_r_logic & mask_logic
        
        w=where1(logic)
        print("A total of %s/%s passed" % (w.size,ntot))

        self.logic = logic
        self.mag_and_mask_logic = cmag_r_logic & mask_logic


