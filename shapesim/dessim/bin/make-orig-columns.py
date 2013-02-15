"""
    %prog [options] vers

Convert the original fits files for the aardvark sims to an indexed columns
database

vers is the originating version, e.g. aardvark_v0.5d  

Will read from {vers}_truth and write to {vers}_truth.cols
"""

import os
import sys
import glob
import numpy

import esutil as eu
from shapesim.dessim import files

from optparse import OptionParser
parser=OptionParser(__doc__)

def get_flist(vers):
    d=files.get_original_dir(vers)
    flist=glob.glob(os.path.join(d, '*.fit'))
    return flist

def open_columns(vers):
    import columns
    cdir=files.get_columns_dir(vers)
    if os.path.exists(cdir):
        raise ValueError("directory exists: %s" % cdir)
    print 'creating columns dir:',cdir
    cols=columns.Columns(cdir)
    cols.create()
    return cols

"""
def add_htm(cols, data, depth=12):
    h=eu.htm.HTM(depth)
    htmid=h.lookup_id(data['ra'], data['dec'])

    cols.write_column('htmid%s' % depth, htmid)
"""

def add_indexes(cols):
    icols=['ra','dec','fid']

    for icol in icols:
        print 'index:',icol
        cols[icol].add_index(create=True)

def get_fid(fname):
    fid = (fname.split('.'))[-2]
    return int( fid )

def add_fid(cols, fname, nrows):
    fid = get_fid(fname)
    fids=numpy.zeros(nrows, dtype='i1') + fid
    cols.write_column('fid', fids)

def add_data(cols, fname):
    data=eu.io.read(fname, ensure_native=True, lower=True)
    cols.write_columns(data)

    add_fid(cols, fname, data.size)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    vers=args[0]

    cols=open_columns(vers)
    flist=get_flist(vers)

    nfiles=len(flist)
    for i,fname in enumerate(flist):
        print '%d/%d  %s' % (i+1,nfiles,fname)
        add_data(cols, fname)

    add_indexes(cols)

main()
