
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

parser=OptionParser(__doc__)
parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")
parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")
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

    if bintype is None or nbin is None:
        raise ValueError("currently demand some kind of binning")

    b = lensing.binning.instantiate_binner(bintype, nbin)

    data = lensing.files.sample_read('binned', lensrun, name=b.name())
    extra='randmatch-%s' % randrun
    rand = lensing.files.sample_read('binned', lensrun, name=b.name(), extra=extra)

    biggles.configure('screen','width', 1000)
    biggles.configure('screen','height', 1000)
    for binnum in xrange(nbin):
        w,=where(data['r'][binnum] > minrad)
        wr,=where(rand['r'][binnum] > minrad)
        if w.size != wr.size:
            raise ValueError("Found > minrad %d from data but %d "
                             "from rand" % (w.size,wr.size))
        

        label=b.bin_label(binnum)

        lensing.plotting.plot2dsig_new(data['r'][binnum], 
                                   data['dsig'][binnum], data['dsigerr'][binnum],
                                   rand['dsig'][binnum], rand['dsigerr'][binnum],
                                   plot_label=label, label1='data', label2='random',
                                   range4var=[0.1,100])


        key=raw_input("hit a key (q to quit): ")
        if key == 'q':
            return


main()
