"""
    %prog run

Description

    Make plots comparing to truth.

    If a run is sent, can plot vs s2n or e depending
    on the run type.

    If something like "set-edg01" is sent, then a set
    of runs is plotted vs the shear
"""

import sys
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-x','--xrange',default=None,
                  help='xrange, default %default')
parser.add_option('-y','--yrange',default='-0.05,0.05',
                  help='yrange, default %default')
parser.add_option('--yrange2',default='-0.05,0.05',
                  help='yrange2, default %default')
parser.add_option('-t','--type',default='diff',
                  help='yrange, default %default')

parser.add_option('--title',default=None,
                  help='add a plot title, default %default')
parser.add_option('--maketitle',action="store_true",
                  help="make a title from the run")

parser.add_option('--s2min',default=None,
                  help="min value in s2")
parser.add_option('--noshow',action="store_true",
                  help="don't show")
parser.add_option('--etot',action='store_true',
                  help=('plot the difference between average '
                        'measured total ellip and true total ellip. '
                        'only for "byellip"'))
parser.add_option('--avg',action='store_true',
                  help="over plot the average over s2")
parser.add_option('--cum',action='store_true',
                  help="Accumulate as a function of s2")


parser.add_option('--skip1',default=None,
                  help="elements in index 1 to skip")
parser.add_option('--skip2',default=None,
                  help="elements in index 2 to skip")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)


run=args[0]

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

s2min=options.s2min
if s2min is not None:
    s2min=float(s2min)

xrng=options.xrange
if xrng is not None:
    xrng = xrng.split(',')
    if len(xrng) != 2:
        raise ValueError("expected xrange min,max")
    xrng = [float(yr) for yr in xrng]


yrng=options.yrange
if yrng is not None:
    yrng = yrng.split(',')
    if len(yrng) != 2:
        raise ValueError("expected yrange min,max")
    yrng = [float(yr) for yr in yrng]


if run[0:3] == 'set':
    if run == 'set-edg01':
        runs = ['gmix-fit-edg03r02',
                'gmix-fit-edg04r01',
                'gmix-fit-edg05r01',
                'gmix-fit-edg06r01',
                'gmix-fit-edg07r01',
                'gmix-fit-edg08r01']
    elif run == 'set-edg02':
        runs = ['gmix-fit-edg09r01',
                'gmix-fit-edg02r02',
                'gmix-fit-edg10r01',
                'gmix-fit-edg11r01',
                'gmix-fit-edg12r01'] #,
                #'gmix-fit-edg13r01']

    else:
        raise ValueError("don't know about set-edg01")
    p=shapesim.plotting.SimPlotterVsShear(runs, maketitle=options.maketitle,
                                          s2min=s2min, skip1=skip1, skip2=skip2,
                                          yrange=yrng,
                                          docum=options.cum)
    p.doplots()
else:

    c=shapesim.read_config(run)
    runtype=c.get('runtype','byellip')
    p=shapesim.plotting.SimPlotter(run, maketitle=options.maketitle,
                                   s2min=s2min, skip1=skip1, skip2=skip2,
                                   docum=options.cum)

    if runtype == 'byellip':
        if options.etot:
            yrng2 = options.yrange2
            if yrng2 is not None:
                yrng2 = yrng2.split(',')
                yrng2 = [float(yrng2[0]),float(yrng2[1])]
            p.plot_ediff_Rshear_vs_e(yrng=yrng, yrng2=yrng2, show=show,
                                     title=options.title)
        else:
            p.plot_shear_vs_e(yrng=yrng,
                              show=show,type=options.type,
                              title=options.title,
                              doavg=options.avg)
    else:
        p.plots_shear_vs_s2n(yrng=yrng, xrng=xrng, type=options.type, 
                             doavg=options.avg,
                             title=options.title,
                             show=show)
