
"""
    %prog [options] lensrun randrun

Description:

    Match redshift histograms between the lenses
    and random.  Writes out a file holding the lens averages.
    
    bintype and nbin are required
"""

from __future__ import print_function

import os
import sys
import lensing
import weighting
import converter

import fitsio

from numpy import zeros

from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option("-w","--extra-weights",default=None,
                  help="field in randoms to use as extra weight")

parser.add_option("--erf-mean-field",default=None,
                  help="field in randoms to use as mean for erf bin weighting")
parser.add_option("--erf-sigma-field",default=None,
                  help="field in randoms to use as sigma for erf bin weighting")


parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")

parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")

parser.add_option("-b",dest="binsize",default=0.01,
    help="Binsize to use in z for hist matching, default %default")

parser.add_option("-s",dest="show",action="store_true", default=False,
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

    conf=lensing.files.cascade_config(lensrun)
    z_field=conf['lens_config']['z_field']
    print('z_field:',z_field)

    b = lensing.binning.instantiate_binner(bintype, nbin)

    # read collated lens catalog and select lenses with the
    # bin criteria.  Then match the randoms redshift histogram

    # this is where z is, may be a different name in the collated data
    data = lensing.files.collated_read(sample=lensrun)
    rand = lensing.files.collated_read(sample=randrun)

    z=data[z_field]

    output = lensing.binning.lensbin_struct(data['rsum'][0].size, n=nbin)

    outextra = 'randmatch-%s' % randrun

    weights_file=lensing.files.sample_file(type='weights',
                                           sample=lensrun,
                                           name=b.name(),
                                           extra=outextra)

    print("opening weights file for writing:",weights_file)
    wfits=fitsio.FITS(weights_file,'rw',clobber=True)

    html_name=lensing.files.sample_file(type='binned-plots', 
                                        sample=lensrun, name=b.name(),
                                        extra=outextra, ext='html')

    makedir_fromfile(html_name)
    html_file=open(html_name,'w')
    html_file.write('<html>\n<body bgcolor=white>\n')

    for binnum in xrange(nbin):

        print("-"*70)
        print("%s/%s: " % (binnum+1,nbin), b.bin_label(binnum))

        png_extra='%02d-%s-%s' % (binnum,outextra,randrun)

        pngfile=lensing.files.sample_file(type='binned-plots',
                                          sample=lensrun,
                                          name=b.name(),
                                          extra=png_extra, ext='png')


        tit=b.bin_label(binnum)
        tit+=' rand: '+randrun

        w = b.select_bin(data, binnum)

        print("matching hist with weights")

        extra_weights=None
        if options.erf_mean_field is not None:
            print("erf weights from:",options.erf_mean_field) 
            if options.erf_sigma_field is None:
                raise ValueError("send both --erf-mean-field and --erf-sigma-field")

            x = rand[options.erf_mean_field]
            sigma=rand[options.erf_sigma_field]
            weights=b.get_bin_erf_weights_1d(binnum, x, sigma)
        else:

            if options.extra_weights is not None:
                extra_weights=rand[extra_weights]

            weights = weighting.hist_match(rand['z'], z[w], binsize,
                                           extra_weights1=extra_weights)

        weights *= ( 1.0/weights.max() )
        effnum = weights.sum()
        effperc = effnum/rand.size
        print("effective number: %d/%d = %0.2f" % (effnum,rand.size, effperc))

        print("combining randoms with weights")
        comb = lensing.outputs.average_lensums(rand, weights=weights)

        weighting.plot_results1d(rand['z'], z[w], weights, binsize, 
                                 pngfile=pngfile, title=tit, show=show)

        wstruct=zeros(weights.size, dtype=[('weights','f8')])
        wstruct['weights'] = weights

        # a new extension for each bin
        wfits.write(wstruct)

        html_file.write('    <img src="%s"><p>\n' % os.path.basename(pngfile))


        # copy all common tags
        for n in comb.dtype.names:
            output[n][binnum] = comb[n][0]

    html_file.write('</body>\n</html>\n')
    html_file.close()

    wfits.close()

    lensing.files.sample_write(data=output,type='binned',
                               sample=lensrun,name=b.name(),extra=outextra)

main()
