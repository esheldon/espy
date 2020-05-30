
"""
    %prog [options] run1 run2 bintype
"""
from __future__ import print_function
import sys
import numpy
import lensing
from lensing.plotting import plot2dsig_over
import esutil as eu
import biggles
from biggles import Table, FramedPlot, Points, Curve,\
        SymmetricErrorBarsY as SymErrY, PlotLabel
from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-t",'--type',default='corrected',
                  help="Should be binned, corrected, corrected-ssh, jackknife.  Default %default")
parser.add_option("--rrange",
                  help="range over which to compute mean ratio")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(1)

    run1 = args[0]
    run2 = args[1]
    bintype = args[2]

    rrange=options.rrange
    if rrange is not None:
        rrange=[float(r) for r in rrange.split(',')]

    biggles.configure('screen','width', 1400)
    biggles.configure('default','fontsize_min',1.0)

    b = lensing.binning.instantiate_binner(bintype)

    name=b.get_name()
    data1 = lensing.files.sample_read(type=options.type, sample=run1, name=name)
    data2 = lensing.files.sample_read(type=options.type, sample=run2, name=name)

    nbin=data1.size

    tab=Table(1,2)
    for i in xrange(nbin):
        tab[0,0] = plot2dsig_over(
            data1['r'][i],data1['dsig'][i],data1['dsigerr'][i],
            data2['r'][i],data2['dsig'][i],data2['dsigerr'][i],
            label1=run1,label2=run2,label=b.bin_label(i),
            show=False)
        tab[0,0].aspect_ratio=1

        pbot = FramedPlot()

        rat=data1['dsig'][i]/data2['dsig'][i]
        raterr = rat*numpy.sqrt( (data1['dsigerr'][i]/data1['dsig'][i])**2
                                +(data2['dsigerr'][i]/data2['dsig'][i])**2 )

        rpts=Points(data1['r'][i], rat, type='filled circle',color='blue')
        rerr=SymErrY(data1['r'][i], rat, raterr,color='blue')

        if rrange is None:
            rat_use,raterr_use=rat,raterr
        else:
            w,=numpy.where( (data1['r'][i] > rrange[0]) & (data1['r'][i] < rrange[1]))
            rat_use,raterr_use=rat[w],raterr[w]
        mn,err = eu.stat.wmom(rat_use, 1/raterr_use**2, calcerr=True)

        plot_rrange=[0.5*data1['r'][i].min(),data1['r'][i].max()*2]
        z=Curve(plot_rrange, [1,1])
        pbot.add(rpts,rerr,z)
        pbot.xlog=True
        pbot.xlabel = lensing.plotting.labels['rproj']
        pbot.ylabel = '%s/%s' % (run1,run2)
        pbot.yrange=[0,2]
        pbot.xrange = plot_rrange
        pbot.aspect_ratio=1

        lab='<ratio>'
        if rrange is not None:
            lab = r'$%s r \in [%.2f,%.2f]$' % (lab,rrange[0],rrange[1])

        lab='%s: %0.2f +/- %0.2f' % (lab,mn,err)
        pbot.add(PlotLabel(0.9,0.9,lab,halign='right'))

        tab[0,1] = pbot
        tab.show()

        epsfile=lensing.files.sample_file(type='binned-plots',
                                          sample=run1,
                                          name=name,
                                          extra='compare-%s-%02d' % (run2,i),
                                          ext='eps')
        print("writing:",epsfile)
        tab.write_eps(epsfile)
        if nbin > 1:
            key=raw_input('hit a key: ')
            if key.lower() == 'q':
                return

main()
