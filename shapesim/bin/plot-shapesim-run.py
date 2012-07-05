"""
    %prog run

Description

    Make plots comparing to truth for a single run, as a function
    of s2 and total ellipticity

"""

import sys
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-y','--yrange',default='-0.05,0.05',
                  help='yrange, default %default')
parser.add_option('-t','--type',default='diff',
                  help='yrange, default %default')
parser.add_option('--s2max',default=None,
                  help="Max (spsf/sobj)**2 to plot, default %default")
parser.add_option('--s2meas',action='store_true',
                  help="Use the measured s2")
parser.add_option('--noshow',action="store_true",
                  help="don't show")

parser.add_option('--skip1',default=None,
                  help="elements in index 1 to skip")
parser.add_option('--skip2',default=None,
                  help="elements in index 2 to skip")

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

skip1=options.skip1
if skip1 is None:
    skip1=[]
else:
    skip1 = [int(v) for v in skip1.split(',')]
skip2=options.skip2
if skip2 is None:
    skip2=[]
else:
    skip2 = [int(v) for v in skip2.split(',')]

yrng=options.yrange
if yrng is not None:
    yrng = yrng.split(',')
    if len(yrng) != 2:
        raise ValueError("expected yrange min,max")
    yrng = [float(yr) for yr in yrng]

c=shapesim.read_config(run)
runtype=c.get('runtype','byellip')

p=shapesim.plotting.SimPlotter(run)

if runtype == 'byellip':
    p.doplots_vs_e(yrange=yrng,s2max=s2max,s2meas=options.s2meas,
                   show=show,type=options.type,skip1=skip1,skip2=skip2)
else:
    p.doplots_vs_s2n(yrange=yrng, show=show, skip1=skip1,skip2=skip2)
