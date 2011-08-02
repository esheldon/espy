
"""
    %prog [options] lensrun randrun

Description:

    Match redshift histograms between the lenses
    and random.  Writes out a file holding the lens averages."""

from __future__ import print_function

import sys
import lensing
import weighting
import converter

from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")
parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")
parser.add_option("--n-near",dest="n_near",default=100,
                  help="The number of nearest neighbors to use, default %default")

def get_n_near(zmin, zmax, npoints, resolution):
    nbin = float(zmax-zmin)/resolution
    x = npoints/nbin
    if x > 100:
        n_near=100
    elif x > npoints:
        n_near=npoints
    else:
        n_near=x

    return n_near


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    lensrun=args[0]
    randrun=args[1]

    bintype=options.bintype
    nbin=int(options.nbin)
    n_near=int(options.n_near)

    if bintype is None or nbin is None:
        raise ValueError("currently demand some kind of binning")

    b = lensing.binning.instantiate_binner(bintype, nbin)

    # read collated lens catalog and select lenses with the
    # bin criteria.  Then match the randoms redshift histogram
    data = lensing.files.sample_read('collated', lensrun)
    rand = lensing.files.sample_read('collated', randrun)

    # resolution gives is a rough number for our desired resolution
    resolution=0.01

    output = lensing.binning.lensbin_struct(data['rsum'][0].size, n=nbin)

    for binnum in xrange(nbin):

        eps_extra='%02d-randmatch-%s' % (binnum,randrun)
        epsfile=lensing.files.sample_file('binned-plots',
                                          lensrun,
                                          name=b.name(),
                                          extra=eps_extra, ext='eps')

        tit=b.bin_label(binnum)
        tit+=' rand: '+randrun

        w = b.select_bin(data, binnum)
        wc = weighting.WeightCalculator(rand['z'], data['z'][w])

        wc.calc_1pass(n_near)

        weighting.plot_results1d(wc, resolution, epsfile=epsfile, title=tit)
        converter.convert(epsfile, dpi=100, verbose=True)

        print("combining randoms with weights")
        comb = lensing.outputs.average_lensums(rand, weights=wc.weight_data['weight'])

        # copy all common tags
        for n in comb.dtype.names:
            output[n][binnum] = comb[n][0]

    outextra='randmatch-%s' % randrun
    lensing.files.sample_write(output,'binned',lensrun,name=b.name(),extra=outextra)

main()
