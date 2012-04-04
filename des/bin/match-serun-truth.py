#!/usr/bin/env python
"""
    %prog [options] run mock_catalog

Set up the files for matching and write out a bash script that
will run smatch.  

Run this script again to put the match info into the columns database

mock_catalog should be e.g. 'desmocks-3.02'
"""
import os
import sys
from sys import stderr
import deswl
import des
import lensing
import des
import esutil as eu
import numpy

from des.compare_truth import get_match_files

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-r','--remake',action="store_true",
                  help=("force a remake of the ascii input catalogs "
                        "and rerun of smatch"))


def make_script(run, mock_catalog):
    fd=get_match_files(run, mock_catalog)
    rad = 1.0
    nside=4096
    script = """
runf=%(runf)s
mockf=%(mockf)s
matchf=%(matchf)s
rad=%(rad)0.2f
nside=%(nside)s
pv $runf | smatch -v -r $rad -n $nside $mockf > $matchf
    """.strip()

    fd['rad'] = rad
    fd['nside'] = nside
    script = script % fd

    with open(fd['scriptf'],'w') as fobj:
        fobj.write(script)
        fobj.write('\n')
    print >>stderr,'run this script',fd['scriptf']

def make_inputs(run, mock_catalog, remake=False):

    fd=get_match_files(run, mock_catalog)
    if not os.path.exists(fd['runf']) or remake:

        print 'reading run data'
        c=des.collate.open_columns(run)
        ra = c['ra'][:]
        dec = c['dec'][:]

        print 'writing run radec:',fd['runf']
        eu.misc.colprint(ra, dec, format='%.16g', file=fd['runf'])

    if not os.path.exists(fd['mockf']) or remake:
        dmc=lensing.scat.DESMockCatalog(mock_catalog)
        mock=dmc.read()
        print 'writing mock radec:',fd['mockf']
        mock_ra = eu.coords.shiftra(mock['ra'],shift=45.)
        mock_dec = -mock['dec']
        eu.misc.colprint(mock_ra, mock_dec, format='%.16g', file=fd['mockf'])

def make_match_column(run, mock_catalog):
    fd = get_match_files(run,mock_catalog)
    print >>stderr,"reading:",fd['matchf']
    rf=eu.recfile.Recfile(fd['matchf'],delim=' ',
                          dtype=[('seid','i8'),('mockid','i8')])
    matches=rf.read()

    c=des.collate.open_columns(run)
    coldata=-9999+numpy.zeros(c['ra'].size,dtype='i8')
    coldata[matches['seid']] = matches['mockid']

    colname = des.compare_truth.match_colname(mock_catalog)

    print >>stderr,'writing column:',colname
    c.write_column(colname, coldata, create=True)

    cname=c[colname].filename
    d=os.path.dirname(cname)
    d=os.path.dirname(d)
    fitsf=os.path.basename(cname)
    fitsf=fitsf.replace('.col','.fits')
    fitsf=os.path.join(d,run+'-fits',fitsf)
    print >>stderr,'writing fits:',fitsf
    eu.io.write(fitsf,coldata,clobber=True)


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) != 2:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    mock_catalog=args[1]

    make_inputs(run, mock_catalog, remake=options.remake)

    # if matches aren't made yet, or remaking, then make
    # the script and exit
    fd = get_match_files(run,mock_catalog)
    if not os.path.exists(fd['matchf']) or options.remake:
        make_script(run, mock_catalog)
    else:
        make_match_column(run, mock_catalog)

main()
