"""
    %prog run type

Description

    type is 'regauss', 'am+', 'noweight' and can be comma separated list

"""

import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-y','--yrange',default='-0.1,0.1',
                  help='yrange, default %default')
parser.add_option('-R','--Rmin',default=0.33,
                  help="Min R value to plot, default %default")
parser.add_option('--reduce-plot-key',action='store_true')
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(45)


run=args[0]
types=args[1].split(',')

yrng=options.yrange
if yrng is not None:
    yrng = yrng.split(',')
    if len(yrng) != 2:
        raise ValueError("expected yrange min,max")
    yrng = [float(yr) for yr in yrng]

Rmin=0.0
if options.Rmin is not None:
    Rmin=float(options.Rmin)

p=lensing.regauss_sim.RegaussSimPlotter(run)
for type in types:
    p.plot(type,Rmin=Rmin,yrange=yrng,reduce_plot_key=options.reduce_plot_key)
