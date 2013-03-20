
"""
    %prog [options] run1 run2 bintype nbin
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

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 4:
        parser.print_help()
        sys.exit(1)

    run1 = args[0]
    run2 = args[1]
    bintype = args[2]
    nbin = int(args[3])

    biggles.configure('screen','width', 1400)

    b = lensing.binning.instantiate_binner(bintype, nbin)

    data1 = lensing.files.sample_read(type='corrected', sample=run1, name=b.name())
    data2 = lensing.files.sample_read(type='corrected', sample=run2, name=b.name())

    tab=Table(1,2)
    for i in xrange(nbin):
        tab[0,0] = plot2dsig_over(
            data1['r'][i],data1['dsig'][i],data1['dsigerr'][i],
            data2['r'][i],data2['dsig'][i],data2['dsigerr'][i],
            label1=run1,label2=run2,label=b.bin_label(i),
            show=False)
        tab[0,0].aspect_ratio=1

        pbot = FramedPlot()
        diff = data1['dsig'][i] - data2['dsig'][i]
        differr = numpy.sqrt(data1['dsigerr'][i]**2 + data2['dsigerr'][i]**2)
        dpts=Points(data1['r'][i], diff, type='filled circle',color='blue')
        derr=SymErrY(data1['r'][i], diff, differr,color='blue')

        mn,err = eu.stat.wmom(diff, 1/differr**2)

        rrange=[0.5*data1['r'][i].min(),data1['r'][i].max()*2]
        z=Curve(rrange, [0,0])
        pbot.add(dpts,derr,z)
        pbot.xlog=True
        pbot.xlabel = lensing.plotting.labels['rproj']
        pbot.ylabel = r'$\Delta\Sigma_1 \Delta\Sigma_2 ~ [M_{sun} pc^{-2}]$'
        pbot.yrange=[-10,10]
        pbot.xrange = rrange
        pbot.aspect_ratio=1
        lab='<diff> = %0.2f +/- %0.2f' % (mn,err)
        pbot.add(PlotLabel(0.9,0.9,lab,halign='right'))

        tab[0,1] = pbot
        tab.show()
        key=raw_input('hit a key: ')
        if key.lower() == 'q':
            return

main()
