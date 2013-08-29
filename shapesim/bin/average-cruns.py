"""
    %prog new_run_name run1 run2 ..

Description

    Average the set of runs
"""

import sys
import numpy
from numpy import sqrt

import esutil as eu
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)


def do_sums(data, data0):
    data['nsum'] += data0['nsum']
    data['Q_sum'] += data0['Q_sum']
    data['Cinv_sum'] += data0['Cinv_sum']
    data['shear_cov_inv_sum'] += data0['shear_cov_inv_sum']

    if 'flux_sum' in data.dtype.names:
        data['flux_sum'] += data0['flux_sum']
        data['flux_err2invsum'] += data0['flux_err2invsum']
        data['flux_s2n_sum'] += data0['flux_s2n_sum']

        data['T_sum'] += data0['T_sum']
        data['T_err2invsum'] += data0['T_err2invsum']
        data['T_s2n_sum'] += data0['T_s2n_sum']




def do_avg(data):

    if 'flux_sum' in data.dtype.names:
        data['flux'] = data['flux_sum']/data['nsum']
        data['flux_err'] = sqrt(1.0/data['flux_err2invsum'])

        data['T'] = data['T_sum']/data['nsum']
        data['T_err'] = sqrt(1.0/data['T_err2invsum'])

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
    eu.ostools.makedirs_fromfile(fout)
    print 'writing run average:',fout
    eu.io.write(fout, data, clobber=True)

main()
