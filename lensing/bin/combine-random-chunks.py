"""
    %prog sample nchunk

Description:

    Make condor scripts for generation of random lcat in chunks.

"""
from __future__ import print_function
import sys, os
import lensing
import esutil as eu
import numpy
from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

sample = args[0]
nchunk = int(args[1])


alldata=[]
for i in xrange(nchunk):
    extra_name='chunk%02d' % i

    t=lensing.files.lcat_read(sample,extra=extra_name)

    alldata.append(t)

data=eu.numpy_util.combine_arrlist(alldata)

data['zindex'] = numpy.arange(data.size, dtype='i8')
lensing.files.lcat_write(sample, data)
