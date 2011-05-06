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

from columns import Columns

# select the maxbcg input data from a columns database
class MaxbcgColumnSelector:
    """

    mcs = MaxbcgColumnSelector('prim03')
    mcs.select()
    mcs.write('fits')
    mcs.write('rec')

    """
    def __init__(self, rmax=22.0):
        """
        E.g.  procrun='prim03'
        """

        self.basedir = expand_path('~/oh/maxbcg-input')
        self.rmax=rmax

        self.masktype = 'tycho'

        self.cols = es_sdsspy.sweeps_collate.open_columns('primgal')

        self.logic = None
        self.mag_and_mask_logic = None

    def write_columns(self):

        outcols = Columns( self.coldir() )
        if os.path.exists(outcols.dir):
            raise ValueError("coldir exists, start fresh: '%s'" % outcols.dir)
        outcols.create()

        if self.logic is None:
            self.select()

        w=where1(self.mag_and_mask_logic)
        print("\nkeeping those that pass mag and tycho logic: %s/%s" % (w.size, self.logic.size))


        colnames = ['photoid',
                    'ra',
                    'dec',
                    'objc_flags',
                    'objc_flags2']

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


 
    def write(self, filetype):
        if self.logic is None:
            self.select()

        dir = self.dir()
        if not os.path.exists(dir):
            os.makedirs(dir)

        outfile=self.filename(filetype)
        print("Will write to file:",outfile)

        data = self.read_input_columns()

        if filetype == 'fits':
            print("Converting to big endian in place")
            to_big_endian(data, inplace=True)

        print("Writing to file:",outfile)

        # clobber is actually the default for sfile, so
        # the keyword will be ignored
        esutil.io.write(outfile, data, clobber=True)

        del data



    def select(self):
        # for shorthand notation
        c = self.cols

        mask_colname = 'in'+self.masktype
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
        w=where1(cmag_r_logic)
        print("    %s/%s passed" % (w.size,c[mask_colname].size))


        logic = binned_logic & object1_logic & cmag_r_logic & mask_logic
        
        w=where1(logic)
        print("A total of %s/%s passed" % (w.size,c[mask_colname].size))

        self.logic = logic
        self.mag_and_mask_logic = cmag_r_logic & mask_logic

    def prefix(self):
        prefix = self.cols.dir.split('/')[-2].replace('_','-')
        if prefix.find('dr8') != -1:
            prefix='dr8'
        else:
            raise ValueError("only support dr8 sweeps for now")
        return prefix


    def coldir(self):
        prefix=self.prefix()
        coldir='maxbcg-input-%s.cols' % prefix
        coldir = path_join(self.basedir, coldir)
        return coldir


    def dir(self):
        outdir=path_join(self.basedir, self.prefix())
        return outdir

    def filename(self, filetype):
        if filetype == 'fits':
            ext='fits'
        elif filetype == 'rec':
            ext='rec'
        else:
            raise ValueError("Choose filetype 'fits' or 'rec'")

        outdir = self.dir()

        prefix = self.prefix()
        outfile='maxbcg-input-%s-%s-r%0.1f.%s' % (prefix,self.masktype,self.rmax,ext)
        outfile = path_join(outdir,outfile) 
        return outfile

    def read_input_columns(self):

        if self.logic is None:
            self.select()

        # column here is actually a dummy 'i2' which we'lll place the
        # "keep" flag
        colnames = ['thing_id',
                    'column',
                    'ra',
                    'dec',
                    'objc_flags']
        colnames += ['flags_'+f           for f in ['u','g','r','i','z'] ]
        colnames += ['modelflux_'+f       for f in ['u','g','r','i','z'] ]
        colnames += ['modelflux_ivar_'+f  for f in ['u','g','r','i','z'] ]
        colnames += ['cmodelflux_'+f      for f in ['u','g','r','i','z'] ]
        colnames += ['cmodelflux_ivar_'+f for f in ['u','g','r','i','z'] ]
        colnames += ['extinction_'+f      for f in ['u','g','r','i','z'] ]


        print("Reading columns:", colnames)
        w=where1(self.mag_and_mask_logic)

        d = self.cols.read_columns(colnames, 
                                   rows=w, 
                                   verbose=True)

        # re-interpret the modelflux, etc. as 5 element arrays
        dtype = [('thing_id',       'i4'),
                 ('keep',           'i2'),
                 ('ra',             'f8'),
                 ('dec',            'f8'),
                 ('objc_flags',     'i4'),
                 ('flags',          'i4', 5),
                 ('modelflux',      'f4', 5),
                 ('modelflux_ivar', 'f4', 5),
                 ('cmodelflux',     'f4', 5),
                 ('cmodelflux_ivar','f4', 5),
                 ('extinction',     'f4', 5)]

        data = d.view(dtype)

        data['keep'] = 0

        wrest = where1( self.logic[w] )
        data['keep'][wrest] = 1
                 
        return data


