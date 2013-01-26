"""
    %prog [options]
"""

import sys
import os
import cluster_step
from cluster_step import files
from cluster_step.fitters import BiasFitter

from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')
parser.add_option('--s2n',None, help="s2n range as csv, %default")
parser.add_option('--sratio',default='1.0,1.e6', help='sratio range, %default')

parser.add_option('-p','--psfnums',default='1,2,3,4,5,6',
                  help='restrict to these PSFs, comma separated')

parser.add_option('-s','--set',default='use1',
                  help="selection set")

parser.add_option('-t','--type',default=None,
                  help="limit to objects best fit by this model")
parser.add_option('-f','--field',default='s2n_w',
                  help="field for S/N, default %default")

parser.add_option('--fitter',default='mcmc',
                  help="type of fitter, mcmc or lm, default %default")


parser.add_option('-y','--yrange',default=None,
                  help="yrange for plot")
parser.add_option('--show',action='store_true',
                  help="show the plot on the screen")


def doplot(bf, yrange=None, fitter_type='mcmc'):
    import biggles
    plt=biggles.FramedPlot()

    color1='red'
    color2='blue'
    g1pts=biggles.Points(bf.g1true, bf.g1diff, color=color1,type='filled circle')
    gerr1pts=biggles.SymmetricErrorBarsY(bf.g1true, bf.g1diff, bf.g1err, color=color1)
    g1pts.label=r'$\gamma_1$'

    g2pts=biggles.Points(bf.g2true, bf.g2diff, color=color2,type='filled circle')
    gerr2pts=biggles.SymmetricErrorBarsY(bf.g2true, bf.g2diff, bf.g2err, color=color2)
    g2pts.label=r'$\gamma_2$'

    ply1=bf.get_g1poly()
    ply2=bf.get_g2poly()
    g1c=biggles.Curve(bf.g1true, ply1(bf.g1true),color=color1)
    g2c=biggles.Curve(bf.g2true, ply2(bf.g2true),color=color2)

    res1=bf.g1fit.get_result()
    res2=bf.g2fit.get_result()
    m1_mess=r'$m1: %.4g \pm %.4g$' % (res1['pars'][0],res1['perr'][0])
    c1_mess=r'$c1: %.4g \pm %.4g$' % (res1['pars'][1],res1['perr'][1])
    m2_mess=r'$m2: %.4g \pm %.4g$' % (res2['pars'][0],res2['perr'][0])
    c2_mess=r'$c2: %.4g \pm %.4g$' % (res2['pars'][1],res2['perr'][1])

    y=0.3
    dec=-0.05
    m1lab=biggles.PlotLabel(0.1,y,m1_mess,halign='left')
    y += dec
    c1lab=biggles.PlotLabel(0.1,y,c1_mess,halign='left')
    y += dec
    m2lab=biggles.PlotLabel(0.1,y,m2_mess,halign='left')
    y += dec
    c2lab=biggles.PlotLabel(0.1,y,c2_mess,halign='left')

    key=biggles.PlotKey(0.9,0.15,[g1pts,g2pts],halign='right')

    plt.add(g1pts, gerr1pts, g1c, g2pts, gerr2pts, g2c, key)
    plt.add(m1lab,c1lab,m2lab,c2lab)

    plt.xlabel=r'$\gamma_{true}$'
    plt.ylabel=r'$\Delta \gamma$'
    plt.aspect_ratio=1
    if yrange is not None:
        plt.yrange=yrange
    plt.show()

    if fitter_type=='mcmc':
        import mcmc
        names=['m','c']
        plt1=mcmc.plot_results(bf.g1fit.trials,names=names,show=False)
        plt2=mcmc.plot_results(bf.g2fit.trials,names=names,show=False)

        plt1.title='g1'
        plt2.title='g2'

        plt1.show()
        plt2.show()

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if options.run is None or options.s2n is None:
        parser.print_help()
        sys.exit(1)

    run=options.run
    doshow  = options.show
    yrange=options.yrange
    if yrange is not None:
        yrange=[float(y) for y in yrange.split(',')]

    s2n_range=[float(s) for s in options.s2n.split(',')]
    sratio_range=[float(s) for s in options.sratio.split(',')]
    psfnums=[int(s) for s in options.psfnums.split(',')]

    s2n_field=options.field

    reader=files.Reader(run=run, 
                        psfnums=psfnums,
                        s2n_range=s2n_range,
                        sratio_range=sratio_range,
                        setname=options.set,
                        progress=True)

    data=reader.get_data()

    bf=BiasFitter(data, order=1, s2n_field=options.field, 
                  fitter=options.fitter)

    print bf.g1fit
    print bf.g2fit

    if doshow:
        doplot(bf,yrange=yrange, fitter_type=options.fitter)

main()
