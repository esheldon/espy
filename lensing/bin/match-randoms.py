
"""
    %prog [options] lensrun randrun

Description:

    Match redshift histograms between the lenses
    and random.  Writes out a file holding the lens averages."""

from __future__ import print_function

import os
import sys
import lensing
import weighting
import converter

import esutil as eu

import numpy
from numpy import where, zeros

from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")
parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")
parser.add_option("-b",dest="binsize",default=0.01,
                  help="Binsize to use in z for hist matching, default %default")
parser.add_option("-s",dest="show",default=False,
                  help="Show histogram comparisons on the screen. default %default")


def makedir_fromfile(path):
    d=os.path.dirname(path)
    if not os.path.exists(d):
        print("Making dir:",d)
        os.makedirs(d)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    lensrun=args[0]
    randrun=args[1]

    bintype=options.bintype
    nbin=int(options.nbin)
    binsize=float(options.binsize)
    show=options.show

    if bintype is None or nbin is None:
        raise ValueError("currently demand some kind of binning")

    b = lensing.binning.instantiate_binner(bintype, nbin)

    # read collated lens catalog and select lenses with the
    # bin criteria.  Then match the randoms redshift histogram
    data = lensing.files.sample_read('collated', lensrun)
    rand = lensing.files.sample_read('collated', randrun)

    output = lensing.binning.lensbin_struct(data['rsum'][0].size, n=nbin)

    outextra='randmatch-%s' % randrun
    weights_file=lensing.files.sample_file('weights',lensrun,name=b.name(),extra=outextra)
    print("opening weights file for writing:",weights_file)
    wobj=eu.sfile.Open(weights_file,'w')
    for binnum in xrange(nbin):

        print("-"*70)
        print(b.bin_label(binnum))

        eps_extra='%02d-randmatch-%s' % (binnum,randrun)
        epsfile=lensing.files.sample_file('binned-plots',
                                          lensrun,
                                          name=b.name(),
                                          extra=eps_extra, ext='eps')
        makedir_fromfile(epsfile)


        w = b.select_bin(data, binnum)

        weights = weighting.hist_match(rand['z'], data['z'][w], binsize)


        effnum = weights.sum()
        effperc = effnum/rand.size
        print("effective number: %d/%d = %0.2f" % (effnum,rand.size, effperc))


        tit=b.bin_label(binnum)
        tit+=' rand: '+randrun
        weighting.plot_results1d(rand['z'], data['z'][w], weights, binsize, 
                                 epsfile=epsfile, title=tit, show=show)
        converter.convert(epsfile, dpi=120, verbose=True)

        print("combining randoms with weights")
        comb = lensing.outputs.average_lensums(rand, weights=weights)

        wstruct=zeros(1, dtype=[('weights','f8',weights.size)])
        wstruct['weights'] = weights
        wobj.write(wstruct)

        # copy all common tags
        for n in comb.dtype.names:
            output[n][binnum] = comb[n][0]

    wobj.close()
    lensing.files.sample_write(output,'binned',lensrun,name=b.name(),extra=outextra)

main()
