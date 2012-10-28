"""
    %prog run
"""

import sys
import shapesim
import biggles
import lensing
from numpy import zeros

from optparse import OptionParser
parser=OptionParser(__doc__)

profmap={1:'gauss', 2:'exp', 3:'dev', 4:'Bulge Disk conc',
         5:'Bulge Disk ortho', 6:'Bulge Disk offcen'}

g1true=0.03
g2true=0.02

ecolors=['blue','red','magenta']
evals=[0.0, 0.3, 0.6]
hlrvals=[0.4,0.8,1.2,1.6,2.]

yrange=[-0.014,0.014]
def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]

    biggles.configure("default","fontsize_min",1.0)
    biggles.configure("screen","width",1800)
    biggles.configure("screen","height",1100)

    nrows=2
    ncols=3
    #arr=biggles.FramedArray(3,2)
    tab=biggles.Table(nrows,ncols)

    for ip,p in enumerate([1,2,3,4,5,6]):
        prow = ip / ncols
        pcol = ip % ncols
        arr=biggles.FramedArray(2,1)
        arr.xlabel='half light radius'
        arr.ylabel=r'$\Delta \gamma/\gamma$'
        arr.title=profmap[p]

        profile='p%02i' % p
        pname=profmap[p]

        curves=[]
        for ie,ellip in enumerate([0,3,6]):

            g1fracdiff=zeros(len(hlrvals))
            g2fracdiff=zeros(len(hlrvals))
            g1errs=zeros(len(hlrvals))
            g2errs=zeros(len(hlrvals))
            for ires,res in enumerate([4,8,12,16,20]):
                g=shapesim.gmix_fit_sim.GMixGalSim(run,profile,ellip,res)
                t=g.read_output()
                g1,g2,R,g1err,g2err=lensing.util.average_shear(t['e1'],t['e2'],
                                                               e1err=t['e1err'],
                                                               e2err=t['e2err'])

                g1fracdiff[ires] = (g1-g1true)/g1true
                g2fracdiff[ires] = (g2-g2true)/g2true
                g1errs[ires] = g1err/g1true
                g2errs[ires] = g2err/g2true

            c1=biggles.Curve(hlrvals,g1fracdiff,color=ecolors[ie])
            c1err=biggles.SymmetricErrorBarsY(hlrvals,g1fracdiff,g1errs,color=ecolors[ie])
            c2=biggles.Curve(hlrvals,g2fracdiff,color=ecolors[ie])
            c2err=biggles.SymmetricErrorBarsY(hlrvals,g2fracdiff,g2errs,color=ecolors[ie])

            c1.label='e: %.2f' % evals[ie]
            curves.append(c1)

            arr[0,0].add(c1,c1err)
            arr[1,0].add(c2,c2err)
            arr[0,0].add(biggles.Curve(hlrvals,zeros(len(hlrvals))))
            arr[1,0].add(biggles.Curve(hlrvals,zeros(len(hlrvals))))

            arr[0,0].add(biggles.PlotLabel(0.9,0.1,r'$\gamma_1$',halign='right'))
            arr[1,0].add(biggles.PlotLabel(0.9,0.1,r'$\gamma_2$',halign='right'))

        key=biggles.PlotKey(0.9,0.9,curves,halign='right')
        arr[0,0].add(key)
        arr.yrange=yrange
        #arr.show()

        tab[prow,pcol] = arr

    tab.show()

 
main()
