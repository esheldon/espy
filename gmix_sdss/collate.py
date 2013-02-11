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

        flist=self.flist
        runs=numpy.unique(flist['run'])

        for run in runs:
            w,=numpy.where(flist['run']==run)
            camcols=numpy.unique(flist['camcol'][w])

            for camcol in camcols:
                
                data=files.read_sweep(gmix_run=self['gmix_run'],
                                      run=run,
                                      camcol=camcol)
                cols.write_columns(data) 

    def _open_for_writing(self):
        if os.path.exists(self._cdir):
            raise ValueError("columns dir already exists: %s" % dir)

        return columns.Columns(self._cdir)

    def _open_for_reading(self):
        return columns.Columns(self._cdir)

