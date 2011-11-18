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

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-f","--filter", default=None, help="the band for rg runs")
parser.add_option("-r","--run", default='any', help="Only plot the input run")
parser.add_option("-d","--coldir", default=None, help="Use this columns dir location")
parser.add_option("--rmag-max", default=None, help="Max mag in r")
parser.add_option("--rmag-min", default=None, help="Min mag in r")
parser.add_option("-y", "--yrange", default=None, 
                  help="y range on plots. default +/-0.01 for meane, +/-0.005 for residual")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 4:
    parser.print_help()
    sys.exit(45)

procrun = args[0]
sweeptype = args[1]
fieldstring = args[2]
plot_type = args[3]

filter = options.filter
run = options.run
if run != 'any':
    run = int(run)
    nperbin=50000
else:
    nperbin = 500000

if options.yrange is None:
    if plot_type == 'residual':
        #yrange = [-0.005,0.005]
        yrange = None
    else:
        yrange = [-0.01,0.01]
else:
    yrange = options.yrange.split(',')
    yrange = [float(yr) for yr in yrange]

fields = fieldstring.split(',')

print 'fields:',fields
print 'yrange:',yrange

if procrun == 'princeton':
    t = lensing.princeton.Tester(run)
else:
    if filter is None:
        raise ValueError("You must enter -f filter for regauss runs")
    t = lensing.regauss_test.Tester(procrun, sweeptype, filter, run, coldir=options.coldir)

for field in fields:
    t.plot_vs_field(field, plot_type, nperbin=nperbin, yrange=yrange, show=False, 
                    rmag_min=options.rmag_min, rmag_max=options.rmag_max)
