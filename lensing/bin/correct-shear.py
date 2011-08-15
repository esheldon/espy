
"""
    %prog [options] lensrun randrun

Description:

    Correct the delta sigma for 

        - subtract the signal from randoms
        - 1/ssh  the shear polarizability correction
        - clustering of sources with lenses
        
    Note the subtraction of randoms should happen before the other two.  You
    need to pick a radius outside of which the subtraction occurs if you don't
    like the default.

    The inputs to this program are the binned outputs from /bin/bin-lenses.py
    and the matched random outputs from /bin/match-randoms.py"""

import sys
import lensing
import numpy
from numpy import where
from optparse import OptionParser
import biggles
import esutil as eu

parser=OptionParser(__doc__)
parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")
parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")

parser.add_option("-s",dest="subtract_rand",action='store_true',default=False,
                  help="Subtract the randoms.  default %default")
parser.add_option("-r",dest="minrad",default=1.,
                  help="The minimum radius for the subtraction of random signal. "
                       "Set to a large number to turn off subtraction. Default %default Mpc")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    lensrun=args[0]
    randrun=args[1]

    bintype=options.bintype
    nbin=int(options.nbin)
    minrad=float(options.minrad)
    subtract_rand = bool(options.subtract_rand)

    if bintype is None or nbin is None:
        raise ValueError("currently demand some kind of binning")

    b = lensing.binning.instantiate_binner(bintype, nbin)

    alldata_precorr = lensing.files.sample_read('binned', lensrun, name=b.name())
    extra='randmatch-%s' % randrun
    allrand = lensing.files.sample_read('binned', lensrun, name=b.name(), extra=extra)

    alldata = lensing.correct.correct(alldata_precorr, allrand, 
                                      subtract_rand=subtract_rand, minrad=minrad)
    
    biggles.configure('screen','width', 1100)
    biggles.configure('screen','height', 1100)

    
    range4var = [0.1,100]
    for binnum in xrange(nbin):

        data = alldata[binnum]
        rand = allrand[binnum]
        

        tab = biggles.Table(1,2)

        label=b.bin_label(binnum)

        arr=lensing.plotting.plot2dsig(data['r'], 
                                       data['dsig'], data['dsigerr'],
                                       rand['dsig'], rand['dsigerr'],
                                       plot_label=label, label1='data', label2='random',
                                       range4var=[0.1,100],show=False)

        
        tab[0,0] = arr


        cplt = eu.plotting.bscatter(data['r'], data['clust_corr']-1, yerr=data['clust_corr_err'],
                                    xlabel=lensing.plotting.labels['rproj'], 
                                    ylabel='Corr-1', show=False, xlog=True, ylog=True)

        cplt.aspect_ratio=1
        tab[0,1] = cplt

        tab.show()

        key=raw_input("hit a key (q to quit): ")
        if key == 'q':
            return


main()
