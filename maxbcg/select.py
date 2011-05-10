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

import maxbcg

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

# select the maxbcg input data from a columns database
# and into a columns database
class Selector:
    """

    mcs = Selector()
    mcs.select()
    mcs.write_columns()

    """
    def __init__(self, rmax=22.0):
        """
        E.g.  procrun='prim03'
        """

        self.basedir = os.environ['MAXBCG_INPUT']
        self.rmax=rmax

        self.cols = es_sdsspy.sweeps_collate.open_columns('primgal')

        self.logic = None
        self.mag_and_mask_logic = None

    def write_columns(self):

        d = maxbcg.files.input_coldir()
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
                    'intycho']

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


 

    def select(self):
        # for shorthand notation
        c = self.cols

        mask_colname = 'inbasic'
        print("Reading",mask_colname)
        inmask = c[mask_colname][:]
        print("Getting mask logic")
        mask_logic = (inmask == 1)
        w=where1(mask_logic)
        print("    %s/%s passed" % (w.size,c[mask_colname].size))


        # now read columns needed for binned_logic, object1_logic
        print("Reading flags, objc_flags")
        colnames=['objc_flags']
        colnames += ['flags_'+f  for f in ['u','g','r','i','z'] ]
        data = c.read_columns(colnames)

        selector = es_sdsspy.select.Selector(data)

        print("binned logic")
        binned_logic = selector.binned_logic()

        w=where1(binned_logic)
        print("    %s/%s passed" % (w.size,c[mask_colname].size))

        print("object1 logic")
        object1_logic = selector.object1_logic()
        
        w=where1(object1_logic)


        print("Reading cmodelmag_dered_r")
        cmag_r = c['cmodelmag_dered_r'][:]
        print('mag logic')
        cmag_r_logic = (cmag_r < self.rmax)
        w=where1(mask_logic & cmag_r_logic)
        print("    %s/%s passed" % (w.size,c[mask_colname].size))

        w=where1(cmag_r_logic)
        print("    %s/%s passed mag+mask" % (w.size,c[mask_colname].size))

        logic = binned_logic & object1_logic & cmag_r_logic & mask_logic
        
        w=where1(logic)
        print("A total of %s/%s passed" % (w.size,c[mask_colname].size))

        self.logic = logic
        self.mag_and_mask_logic = cmag_r_logic & mask_logic

