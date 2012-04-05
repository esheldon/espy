#!/usr/bin/env python
"""
    %prog [options] serun merun

Set up the files for matching and write out a bash script that
will run smatch.  
"""
import os
import sys
from sys import stderr
import deswl
import des
from des.compare_serun_merun import get_match_files, get_match_colname

import lensing
import esutil as eu
import numpy

import recfile

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-r','--remake',action="store_true",
                  help=("force a remake of the ascii input catalogs "
                        "and rerun of smatch"))


def make_script(serun, merun):

    fd = get_match_files(serun,merun)

    rad = 1.0
    nside=4096
    script = """
sef=%(sef)s
mef=%(mef)s
matchf=%(matchf)s
rad=%(rad)0.2f
nside=%(nside)s
pv $sef | smatch -v -r $rad -n $nside $mef > $matchf
    """.strip()

    fd['rad'] = rad
    fd['nside'] = nside
    script = script % fd

    if not os.path.exists(fd['dir']):
        os.makedirs(fd['dir'])
    with open(fd['scriptf'],'w') as fobj:
        fobj.write(script)
        fobj.write('\n')
    print >>stderr,'run this script',fd['scriptf']

def make_inputs(serun, merun, remake=False):

    fd=get_match_files(serun,merun)
    if not os.path.exists(fd['sef']) or remake:

        print 'reading run data for',serun
        c=des.collate.open_columns(serun)
        data = c.read_columns(['ra','dec'])

        print 'writing se radec:',fd['sef']
        with recfile.Recfile(fd['sef'],'w',delim=' ') as fobj:
            fobj.write(data)

    if not os.path.exists(fd['mef']) or remake:

        print 'reading run data for',merun
        c=des.collate.open_columns(merun)
        data = c.read_columns(['ra','dec'])

        print 'writing me radec:',fd['mef']
        with recfile.Recfile(fd['mef'],'w',delim=' ') as fobj:
            fobj.write(data)

def make_match_column(serun, merun):
    fd = get_match_files(serun, merun)
    print >>stderr,"reading:",fd['matchf']
    rf=eu.recfile.Recfile(fd['matchf'],delim=' ',
                          dtype=[('seid','i8'),('meid','i8')])
    matches=rf.read()

    c=des.collate.open_columns(serun)
    coldata=-9999+numpy.zeros(c['ra'].size,dtype='i8')
    coldata[matches['seid']] = matches['meid']

    colname = get_match_colname(merun)

    print >>stderr,'writing column:',colname
    c.write_column(colname, coldata, create=True)

    cname=c[colname].filename
    d=os.path.dirname(cname)
    d=os.path.dirname(d)
    fitsf=os.path.basename(cname)
    fitsf=fitsf.replace('.col','.fits')
    fitsf=os.path.join(d,serun+'-fits',fitsf)
    print >>stderr,'writing fits:',fitsf
    eu.io.write(fitsf,coldata,clobber=True)



def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) != 2:
        parser.print_help()
        sys.exit(45)

    serun=args[0]
    merun=args[1]

    make_inputs(serun, merun, remake=options.remake)

    # if matches aren't made yet, or remaking, then make
    # the script and exit
    fd = get_match_files(serun,merun)
    if not os.path.exists(fd['matchf']) or options.remake:
        make_script(serun,merun)
    else:
        make_match_column(serun,merun)

main()
