
"""
    %prog [options] lensrun

Description:

    Correct the delta sigma for ssh only, no randoms
"""

from __future__ import print_function
import sys, os
import lensing
import numpy
from numpy import where
from optparse import OptionParser
import biggles
import esutil as eu
import converter

parser=OptionParser(__doc__)
parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")
parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")
parser.add_option("--noylog",action='store_true',
                  help="don't use log on y axis")

parser.add_option("-x","--xrange",default=None,
                  help="range for x axis")
parser.add_option("-y","--yrange",default=None,
                  help="range for y axis")

parser.add_option("-p",dest="prompt",action='store_true',default=False,
                  help="Prompt between each plot.  default %default")
parser.add_option("--show",action='store_true',default=False,
                  help="Show plots on screen.  default %default")

def doplot(binned_data, corr_data, label, ylog=True, xrng=None, yrng=None,
           show=False):
    plt=lensing.plotting.plot2dsig_over(binned_data['r'],binned_data['dsig'],binned_data['dsigerr'],
                                        corr_data['r'],corr_data['dsig'],corr_data['dsigerr'],
                                        label=label,
                                        label1='orig', 
                                        label2='ssh corrected',
                                        ylog=ylog,
                                        xrange=xrng,
                                        yrange=yrng,
                                        show=show)
    return plt

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    lensrun=args[0]

    bintype=options.bintype
    nbin=int(options.nbin)

    if bintype is None or nbin is None:
        raise ValueError("currently demand some kind of binning")

    ylog=True
    if options.noylog:
        ylog=False
    xrng=options.xrange
    if xrng is not None:
        xrng = [float(xr) for xr in xrng.split(',')]
    yrng=options.yrange
    if yrng is not None:
        yrng = [float(yr) for yr in yrng.split(',')]


    b = lensing.binning.instantiate_binner(bintype, nbin)

    binned_data = lensing.files.sample_read(type='binned', sample=lensrun, name=b.name())

    alldata = lensing.correct.correct(binned_data, None)

    lensing.files.sample_write(data=alldata,
                               type='corrected',
                               sample=lensrun,
                               name=b.name()) 

    # now some plots
    biggles.configure('screen','width', 1100)
    biggles.configure('screen','height', 1100)
    
    range4var = [0.1,100]
    for binnum in xrange(nbin):
        eps_dsigcorr_extra='dsigcorr-%02d' % binnum
        eps_dsigcorr=lensing.files.sample_file(type='corrected-plots',
                                               sample=lensrun,
                                               name=b.name(),
                                               extra=eps_dsigcorr_extra, ext='eps')


        lensing.files.make_dir_from_path(eps_dsigcorr)

        data = alldata[binnum]

        label=b.bin_label(binnum)

        plt=doplot(binned_data[binnum], data, label, show=options.show, 
                   xrng=xrng, yrng=yrng, ylog=ylog)

        # the corr(r)-1 plot
        plt.write_eps(eps_dsigcorr)
        converter.convert(eps_dsigcorr, dpi=90, verbose=True)

        if options.prompt:
            key=raw_input("hit a key (q to quit): ")
            if key == 'q':
                return

    d=os.path.dirname(eps_dsigcorr)
    print("chdir:",d)
    os.chdir(d)

    outfile = 'dsigcorr.html'
    pattern = eps_dsigcorr.replace('%02d.eps' % (nbin-1,), '*.png')
    pattern=os.path.basename(pattern)
    print("making disg corr html file:",outfile)
    os.system('im2html -p '+pattern+' > '+outfile)


main()
