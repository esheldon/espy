"""
    %prog [options] procrun sweeptype field plot_type

Description:

    type is either 'meane' or 'residual'.  residual is only supported
    for star sweeptype, as it is e_star-e_psf.

    Plot the mean ellipticity vs the indicated field.  procrun should
    be a sweeps regauss run or princeton

    field can be comma-separated.  
    
    For R/amflags/e1/e2/uncer/e1_psf/e2_psf/I*/I*psf you don't need to specify
    the band.  Use rmag for either cmodelmag_dered_r or psfmag_dered_r

    Some pre-calculated fields are 'sigma' and 'sigma_psf','fwhm_psf'

    Other things we want to implement besides meane:
        for sweeptype 'star' we want to plot the residual between
        e* and e*_psf, and sigma/sigma_psf as a function of various things


"""

import sys
import lensing
import numpy

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-f","--filter", default=None, help="the band for rg runs")
parser.add_option("-r","--run", default='any', 
                  help=("Only plot the input run. "
                        "'any' means combine all runs.  'all' means doall "
                        "individual runs.  default %default"))
parser.add_option("-c","--camcol", default='any', 
                  help=("Only plot the input camcol. 'any' means combine all camcols.  'all' means "
                        "do all individual camcols.  default %default"))
parser.add_option("-d","--coldir", default=None, help="Use this columns dir location")
parser.add_option("--rmag-max", default=None, help="Max mag in r")
parser.add_option("--rmag-min", default=None, help="Min mag in r")
parser.add_option("-y", "--yrng", default=None, help="y range on plots.")
parser.add_option("-x", "--xrng", default=None, help="x range on plots.")
parser.add_option("--nperbin", default=None, help="number per bin")


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 4:
        parser.print_help()
        sys.exit(45)

    procrun = args[0]
    sweeptype = args[1]
    fieldstring = args[2]
    plot_type = args[3]

    filter = options.filter
    runs = options.run
    camcols=options.camcol
    nperbin=options.nperbin
    rmag_max=options.rmag_max
    rmag_min=options.rmag_min

    if rmag_min is not None:
        rmag_min=float(rmag_min)
    if rmag_max is not None:
        rmag_max=float(rmag_max)

    xrng = options.xrng
    yrng = options.yrng

    if camcols != 'any':
        if camcols != 'all':
            camcols=camcols.split(',')
            camcols=int(camcols)
        else:
            camcols=[1,2,3,4,5,6]
    else:
        camcols = ['any']

    if runs != 'any':
        if runs != 'all':
            runs=runs.split(',')
            runs = [int(r) for r in runs]
        else:
            tmp = lensing.regauss_test.Tester(procrun, sweeptype, filter, run='any', coldir=options.coldir)
            c=tmp.open_columns()
            runs = c['run'][:]
            runs = numpy.unique(runs)

        if nperbin is None:
            if plot_type == 'residual':
                nperbin=10000
            else:
                # gals need more per bin
                nperbin=50000
    else:
        runs=['any']
        if nperbin is None:
            if camcols[0] != 'any':
                nperbin = 100000
            else:
                nperbin = 500000

    if xrng is not None:
        xrng=xrng.split(',')
        xrng = [float(xr) for xr in xrng]

    if yrng is None:
        if plot_type == 'residual':
            yrng = [-0.004,0.004]
        else:
            if run == 'any':
                yrng=[-0.016,0.016]
            else:
                yrng = [-0.02,0.02]
    else:
        yrng = yrng.split(',')
        yrng = [float(yr) for yr in yrng]

    fields = fieldstring.split(',')

    print 'fields:',fields
    print 'xrng:',xrng
    print 'yrng:',yrng

    if procrun == 'princeton':
        t = lensing.princeton.Tester(run)
        for field in fields:
            t.plot_vs_field(field, plot_type, nperbin=nperbin, 
                            yrng=yrng, xrng=xrng, show=False, 
                            rmag_min=rmag_min, rmag_max=rmag_max)
    else:
        if filter is None:
            raise ValueError("You must enter -f filter for regauss runs")

    
        print runs
        for run in runs:
            # note camcol could be 'any'
            for camcol in camcols:
                t = lensing.regauss_test.Tester(procrun, sweeptype, filter, 
                                                run=run, camcol=camcol, coldir=options.coldir)

                for field in fields:
                    t.plot_vs_field(field, plot_type, nperbin=nperbin, 
                                    yrng=yrng, xrng=xrng, show=False, 
                                    rmag_min=rmag_min, rmag_max=rmag_max)

main()
