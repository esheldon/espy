"""
    %prog run shear_true

run can be a csv
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
parser.add_option('-x','--xrange',default=None,
                  help='x range')

parser.add_option('--s2n-field',default='s2n_matched',
                  help="field for s2n")

parser.add_option('--eps',default=None,
                  help="eps file to write")
parser.add_option('--png',default=None,
                  help="png file to write")


def plot_run(plt, run, shear_true, symbol, color, linestyle, options,
             with_points=True, with_curve=True):
    c = shapesim.read_config(run)
    url=shapesim.get_averaged_url(run, 0)
    data=eu.io.read(url)

    s2n_vals = data[options.s2n_field]
    err=sqrt(data['shear_cov'][:,0,0])

    if options.xrange is None:
        xrng=[0.75*s2n_vals[0], 1.25*s2n_vals[-1]]
    else:
        xrng=options.xrange.split(',')
        xrng=[float(xr) for xr in xrng]


    pts = biggles.Points(s2n_vals, 
                         data['shear'][:,0]/shear_true-1,
                         type=symbol,
                         size=2,
                         color=color)
    ep = biggles.SymmetricErrorBarsY(s2n_vals,
                                     data['shear'][:,0]/shear_true-1,
                                     err/shear_true,
                                     width=2,
                                     color=color)
    crv = biggles.Curve(s2n_vals, 
                        data['shear'][:,0]/shear_true-1,
                        type=linestyle,
                        width=2,
                        color=color)

    crv.label=run
    pts.label=run
    if with_points:
        plt.add(pts,ep)
    if with_curve:
        plt.add(crv)

    """
    plt=eu.plotting.bscatter(s2n_vals,
                             data['shear'][:,0]/shear_true-1,
                             yerr=err/shear_true,
                             xlabel=xlabel,
                             ylabel=r'$\Delta \gamma/\gamma$',
                             color='darkblue',
                             show=False,
                             size=2,
                             plt=plt0)


    plt=eu.plotting.bscatter(s2n_vals,
                             data['shear'][:,0]/shear_true-1,
                             color='darkblue',
                             type='solid',
                             width=2,
                             plt=plt,
                             show=False)
    """

    return xrng, crv, pts

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    runs=args[0].split(',')
    shear_true=float(args[1])

    colors=['darkblue','red','darkgreen','magenta']
    linestyles=['solid','dotdashed','shortdashed','dotted']
    symbols=['filled circle','filled triangle','filled square','filled diamond']

    ylabel=r'$\Delta \gamma/\gamma$'
    if options.s2n_field=='s2n_matched':
        xlabel=r'$(S/N)_{matched}$'
    elif options.s2n_field=='flux_s2n':
        xlabel=r'$(S/N)_{flux}$'
    elif options.s2n_field=='T_s2n':
        xlabel=r'$(S/N)_{size}$'
    else:
        xlabel=options.s2n_field

    yrng=options.yrange.split(',')
    yrng=[float(yr) for yr in yrng]


    plt=biggles.FramedPlot()

    plt.xlabel=xlabel
    plt.ylabel=ylabel
    plt.aspect_ratio=1
    plt.xlog=True
    plt.add( biggles.FillBetween([1.e-6,5000], [0.004,0.004], 
                                  [1.e-6,5000], [-0.004,-0.004],
                                  color='grey80') )
    plt.add( biggles.Curve([1.e-6,5000],[0,0]) )

    
    plt.yrange=yrng

    cobj=[]
    for irun,run in enumerate(runs):
        with_curve=False
        xrng_run, crv_run, pts_run = plot_run(plt, run, shear_true, 
                                              symbols[irun],
                                              colors[irun],
                                              linestyles[irun],
                                              options,
                                              with_curve=with_curve)
        if with_curve:
            cobj.append(crv_run)
        else:
            cobj.append(pts_run)

        if irun==0:
            xrng=xrng_run
        else:
            xrng[0] = min(xrng[0],xrng_run[0])
            xrng[1] = max(xrng[1],xrng_run[1])
    
    key=biggles.PlotKey(0.9,0.9,cobj,halign='right')
    plt.add(key)
    plt.xrange=xrng

    plt.show()

    if options.png:
        print options.png
        plt.write_img(800,800,options.png)
    if options.eps:
        print options.eps
        plt.write_eps(options.eps)
main()
