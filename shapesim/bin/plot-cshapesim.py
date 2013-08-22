"""
    %prog run shear_true
"""

import sys
import biggles
import esutil as eu
import shapesim

import numpy
from numpy import sqrt

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-y','--yrange',default='-0.05,0.05',
                  help='y range')

parser.add_option('--s2n-field',default='s2n_matched',
                  help="field for s2n")

parser.add_option('--eps',default=None,
                  help="eps file to write")
parser.add_option('--png',default=None,
                  help="png file to write")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    shear_true=float(args[1])

    c = shapesim.read_config(run)

    url=shapesim.get_averaged_url(run, 0)
    data=eu.io.read(url)

    s2n_vals = data[options.s2n_field]
    xrange=[0.5*s2n_vals[0], 1.5*s2n_vals[-1]]

    if options.s2n_field=='s2n_matched':
        xlabel=r'$(S/N)_{matched}$'
    elif options.s2n_field=='flux_s2n':
        xlabel=r'$(S/N)_{flux}$'
    else:
        xlabel=options.s2n_field

    yrange=options.yrange.split(',')
    yrange=[float(yr) for yr in yrange]

    err=sqrt(data['shear_cov'][:,0,0])

    plt0=biggles.FramedPlot()
    plt0.xlog=True
    plt0.add( biggles.FillBetween([1.e-6,5000], [0.004,0.004], 
                                  [1.e-6,5000], [-0.004,-0.004],
                                  color='grey80') )
    plt0.add( biggles.Curve([1.e-6,5000],[0,0]) )

    plt=eu.plotting.bscatter(s2n_vals,
                             data['shear'][:,0]/shear_true-1,
                             yerr=err/shear_true,
                             xlabel=xlabel,
                             ylabel=r'$\Delta \gamma/\gamma$',
                             xrange=xrange,
                             yrange=yrange,
                             color='darkblue',
                             show=False,
                             size=2,
                             plt=plt0)


    plt=eu.plotting.bscatter(s2n_vals,
                             data['shear'][:,0]/shear_true-1,
                             color='darkblue',
                             type='solid',
                             width=2,
                             plt=plt)

    if options.png:
        print options.png
        plt.write_img(800,800,options.png)
    if options.eps:
        print options.eps
        plt.write_eps(options.eps)
main()
