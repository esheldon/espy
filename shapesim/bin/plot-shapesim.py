"""
    %prog run

Description

    Compare the results with truth.

"""

import sys
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-y','--yrange',default='-0.05,0.05',
                  help='yrange, default %default')
parser.add_option('--s2max',default=None,
                  help="Max (spsf/sobj)**2 to plot, default %default")
parser.add_option('--noshow',action="store_true",
                  help="don't show")
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)


run=args[0]
s2max=options.s2max
if s2max is not None:
    s2max=float(s2max)

if options.noshow:
    show=False
else:
    show=True

yrng=options.yrange
if yrng is not None:
    yrng = yrng.split(',')
    if len(yrng) != 2:
        raise ValueError("expected yrange min,max")
    yrng = [float(yr) for yr in yrng]

p=shapesim.plotting.SimPlotter(run)
p.doplots(yrange=yrng,s2max=s2max,show=show)
