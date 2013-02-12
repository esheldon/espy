from sys import stderr
import os
import numpy
import columns
import esutil as eu
from . import files


index_cols=['photoid','run','camcol','field',
            'am_s2n','s2n','Ts2n','sratio','model']

def open_columns(gmix_run=None):
    d=files.get_columns_dir(gmix_run=gmix_run)
    print 'opening:',d
    return columns.Columns(d)

class ColumnsMaker(dict):
    def __init__(self, **keys):
        self.update(keys)

        self._cdir=files.get_columns_dir(gmix_run=self['gmix_run'])
        self._start=keys.get('start',None)

    def make_columns(self):
        cols=self._open_for_writing()
        flist=files.read_field_cache(gmix_run=self['gmix_run'])

        runs=numpy.unique(flist['run'])
        nrun=len(runs)

        start=self._start
        if start is not None:
            skip=True
        else:
            skip=False

        for i,run in enumerate(runs):

            if start is not None and run < start['run']:
                continue

            print >>stderr,'%d/%d  run: %s' % (i+1,nrun,run),
            w,=numpy.where(flist['run']==run)
            camcols=numpy.unique(flist['camcol'][w])

            for camcol in camcols:
                if start is not None and camcol == start['camcol']:
                    skip=False

                if not skip:
                    print >>stderr,camcol,
                    data=files.read_sweep(gmix_run=self['gmix_run'],
                                          run=run,
                                          camcol=camcol)
                    if data.size > 0: # empty files
                        eu.numpy_util.to_native(data,inplace=True)
                        cols.write_columns(data) 

            print >>stderr,""

    def _open_for_writing(self):
        print >>stderr,'opening columns for writing:',self._cdir
        if os.path.exists(self._cdir) and self._start is None:
            raise ValueError("columns dir already exists: %s" % self._cdir)

        cols=columns.Columns(self._cdir)
        if not os.path.exists(self._cdir):
            cols.create()

        return cols

    def _open_for_reading(self):
        return columns.Columns(self._cdir)

