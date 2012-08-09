
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

import fitsio

import esutil as eu

import numpy
from numpy import where, zeros, ones

from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")
parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")
parser.add_option("-b",dest="binsize",default=0.01,
                  help="Binsize to use in z for hist matching, default %default")
parser.add_option("-s",dest="show",action="store_true", default=False,
                  help="Show histogram comparisons on the screen. default %default")
parser.add_option("-r",dest="remove",action="store_true", default=False,
                  help="Instead of using weighting, actually remove "
                       "randoms to match z hist")


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
    remove=options.remove

    if bintype is None or nbin is None:
        raise ValueError("currently demand some kind of binning")

    b = lensing.binning.instantiate_binner(bintype, nbin)

    # read collated lens catalog and select lenses with the
    # bin criteria.  Then match the randoms redshift histogram
    conf=lensing.files.cascade_config(lensrun)
    # this is where z is, may be a different name in the collated data
    lcat = lensing.files.lcat_read(sample=conf['lens_sample'])
    data = lensing.files.sample_read(type='collated', sample=lensrun, fs='hdfs')
    rand = lensing.files.sample_read(type='collated', sample=randrun, fs='hdfs')

    output = lensing.binning.lensbin_struct(data['rsum'][0].size, n=nbin)

    if not remove:
        outextra='randmatch-%s' % randrun
        weights_file=lensing.files.sample_file(type='weights',
                                               sample=lensrun,
                                               name=b.name(),
                                               extra=outextra)
        print("opening weights file for writing:",weights_file)
        wfits=fitsio.FITS(weights_file,'rw',clobber=True)
    else:
        outextra='randmatch-rm-%s' % randrun

    html_name=lensing.files.sample_file(type='binned-plots', 
                                        sample=lensrun, name=b.name(),
                                        extra=outextra, ext='html')
    makedir_fromfile(html_name)
    html_file=open(html_name,'w')
    html_file.write('<html>\n<body bgcolor=white>\n')

    for binnum in xrange(nbin):

        print("-"*70)
        print("%s/%s: " % (binnum+1,nbin), b.bin_label(binnum))

        if remove:
            eps_extra='%02d-randmatch-rm-%s' % (binnum,randrun)
        else:
            eps_extra='%02d-randmatch-%s' % (binnum,randrun)

        epsfile=lensing.files.sample_file(type='binned-plots',
                                          sample=lensrun,
                                          name=b.name(),
                                          extra=eps_extra, ext='eps')
        pngfile=epsfile.replace('.eps','.png')


        tit=b.bin_label(binnum)
        tit+=' rand: '+randrun

        w = b.select_bin(data, binnum)

        # simple histogram matching
        if remove:
            print("matching hist with removal")
            keep = weighting.hist_match_remove(rand['z'], lcat['z'][w], binsize)
            rkeep = rand[keep]
            perc = keep.size/(1.*rand.size)

            print("used number: %d/%d = %0.2f" % (keep.size,rand.size, perc))

            print("combining kept randoms")
            comb = lensing.outputs.average_lensums(rkeep)

            weights=ones(keep.size)
            weighting.plot_results1d(rkeep['z'], lcat['z'][w], weights, binsize, 
                                     epsfile=epsfile, title=tit, show=show)
        else:
            print("matching hist with weights")
            weights = weighting.hist_match(rand['z'], lcat['z'][w], binsize)


            effnum = weights.sum()
            effperc = effnum/rand.size
            print("effective number: %d/%d = %0.2f" % (effnum,rand.size, effperc))

            print("combining randoms with weights")
            comb = lensing.outputs.average_lensums(rand, weights=weights)

            weighting.plot_results1d(rand['z'], lcat['z'][w], weights, binsize, 
                                     epsfile=epsfile, title=tit, show=show)

            #wstruct=zeros(1, dtype=[('weights','f8',weights.size)])
            wstruct=zeros(weights.size, dtype=[('weights','f8')])
            wstruct['weights'] = weights

            # a new extension for each bin
            wfits.write(wstruct)

        converter.convert(epsfile, dpi=120, verbose=True)
        html_file.write('    <img src="%s"><p>\n' % os.path.basename(pngfile))


        # copy all common tags
        for n in comb.dtype.names:
            output[n][binnum] = comb[n][0]

    html_file.write('</body>\n</html>\n')
    html_file.close()

    if not remove:
        wfits.close()

    lensing.files.sample_write(data=output,type='binned',sample=lensrun,name=b.name(),extra=outextra)

main()
