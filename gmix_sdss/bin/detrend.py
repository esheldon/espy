"""
    %prog [options] gmix_run

Description:

    Write detrended shears.  run regress.py first
"""

from __future__ import print_function
import sys
import numpy
import gmix_sdss
from gmix_sdss.cuts import SRATIO_MIN

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    gmix_run = args[0]

    cols=gmix_sdss.collate.open_columns(gmix_run)


    print("reading g")
    g=cols['g'][:] 

    print("reading s2n")
    s2n = cols['s2n'][:]
    print("reading camcol")
    camcol = cols['camcol'][:]

    selector=gmix_sdss.cuts.Selector(cols)
    selector.do_select()

    w=selector.indices

    gdt = g.copy()
    gdt[:,:] = -9999.
    for col in [1,2,3,4,5,6]:
        wcol,=numpy.where(camcol[w] == col)
        wcol = w[wcol]

        g1,g2 = gmix_sdss.regress.detrend_g_vs_s2n(g1=g[wcol,0],
                                                   g2=g[wcol,1],
                                                   gmix_run=gmix_run,
                                                   camcol=col,
                                                   s2n=s2n[wcol],
                                                   sratio=SRATIO_MIN)
        gdt[wcol,0] = g1
        gdt[wcol,1] = g2

    colname='g_dt'
    print("writing:",colname)
    cols.write_column(colname, gdt, create=True)
main()
