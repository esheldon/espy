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

    s2n_vals    = c['s2n_vals']

    xrange=[0.5*s2n_vals[0], 1.5*s2n_vals[-1]]
    yrange=options.yrange.split(',')
    yrange=[float(yr) for yr in yrange]

    err=sqrt(data['shear_cov'][:,0,0])

    plt0=biggles.FramedPlot()
    plt0.xlog=True
    plt0.add( biggles.FillBetween([1.e-6,5000], [0.004,0.004], 
                                  [1.e-6,5000], [-0.004,-0.004],
                                  color='grey80') )
    plt0.add( biggles.Curve([1.e-6,5000],[0,0]) )

    plt=eu.plotting.bscatter(data['s2n_matched'],
                             data['shear'][:,0]/shear_true-1,
                             yerr=err/shear_true,
                             xlabel=r'$(S/N)_{matched}$',
                             ylabel=r'$\Delta \gamma/\gamma$',
                             xrange=xrange,
                             yrange=yrange,
                             color='darkblue',
                             show=False,
                             size=2,
                             plt=plt0)


    plt=eu.plotting.bscatter(data['s2n_matched'],
                             data['shear'][:,0]/shear_true-1,
                             color='darkblue',
                             type='solid',
                             width=2,
                             plt=plt)

main()
