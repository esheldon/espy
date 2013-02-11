from sys import stderr
import os
import numpy
import columns
from . import files

def open_columns(gmix_run=None):
    d=files.get_columns_dir(gmix_run=gmix_run)
    return columns.Columns(d)

class ColumnsMaker(dict):
    def __init__(self, **keys):
        self.update(keys)

        self._cdir=files.get_columns_dir(gmix_run=self['gmix_run'])


    def make_columns(self):
        cols=self._open_for_writing()
        flist=files.read_field_cache(gmix_run=self['gmix_run'])

        runs=numpy.unique(flist['run'])
        nrun=len(runs)

        for i,run in enumerate(runs):
            print >>stderr,'%d/%d  run: %s' % (i+1,nrun,run),
            w,=numpy.where(flist['run']==run)
            camcols=numpy.unique(flist['camcol'][w])

            for camcol in camcols:
                print >>stderr,camcol,
                
                data=files.read_sweep(gmix_run=self['gmix_run'],
                                      run=run,
                                      camcol=camcol)
                if data.size > 0: # empty files
                    cols.write_columns(data) 
            print >>stderr,""

    def _open_for_writing(self):
        print >>stderr,'opening columns for writing:',self._cdir
        if os.path.exists(self._cdir):
            raise ValueError("columns dir already exists: %s" % self._cdir)

        cols=columns.Columns(self._cdir)
        cols.create()

        return cols

    def _open_for_reading(self):
        return columns.Columns(self._cdir)

