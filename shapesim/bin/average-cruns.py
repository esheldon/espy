"""
    %prog new_run_name run1 run2 ..

Description

    Average the set of runs
"""

import sys
import shapesim

import esutil as eu

from optparse import OptionParser
parser=OptionParser(__doc__)


def do_sums(data, data0):
    data['nsum'] += data0['nsum']
    data['Q_sum'] += data0['Q_sum']
    data['Cinv_sum'] += data0['Cinv_sum']
    data['shear_cov_inv_sum'] += data0['shear_cov_inv_sum']

def do_avg(data):

    for is2n in xrange(data.size):
        C = numpy.linalg.inv(data['Cinv_sum'][is2n])
        shear = numpy.dot(C,data['Q_sum'][is2n])

        shear_cov = numpy.linalg.inv(data['shear_cov_inv_sum'][is2n])

        data['shear'][is2n] = shear
        data['shear_cov'][is2n] = shear_cov

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(45)


    new_run_name=args[0]
    runs2average = args[1:]

    print 'new run name:',new_run_name
    print 'runs2average:',runs2average

    for i,run in enumerate(runs2average):
        fname=shapesim.get_averaged_url(run, 0)
        print fname
        data0=eu.io.read(fname)

        if i==0:
            data=data0
        else:
            do_sums(data, data0)

    do_avg(data)

    fout=shapesim.get_averaged_url(new_run_name, 0)
    print 'writing run average:',fout
    eu.io.write(fout, data, clobber=True)

main()
