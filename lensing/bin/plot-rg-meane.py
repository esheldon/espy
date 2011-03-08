"""
    %prog [options] procrun field

Description:
    Plot the mean ellipticity vs the indicated field.  procrun should
    be a sweeps regauss run or princeton

    field can be comma-separated
        %prog 04
"""

import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-f","--filter", default=None, help="the band for rg runs")
parser.add_option("-r","--run", default='any', help="Only plot the input run")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(45)

procrun = args[0]
fieldstring = args[1]

filter = options.filter
run = options.run
if run != 'any':
    run = int(run)
    nperbin=50000
else:
    nperbin = 500000

fields = fieldstring.split(',')

if procrun == 'princeton':
    t = lensing.princeton.Tester(run)
else:
    if filter is None:
        raise ValueError("You must enter -f filter for regauss runs")
    t = lensing.regauss_test.Tester(procrun, filter, run)

for field in fields:
    t.plot_ellip_vs_field(field, nperbin=nperbin, yrange=[-0.04,0.04], show=False)
