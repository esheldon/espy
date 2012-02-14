"""
    %prog [options] run

Description:
    Check all expected outputs exist for the input run
"""
import sys
from lensing import files
import esutil as eu
from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])


if len(args) < 1:
    parser.print_help()
    sys.exit(1)


run=args[0]

conf = files.cascade_config(run)
nsplit=conf['src_config']['nsplit']

pattern = files.sample_file('lensout',run, split=0,fs='hdfs')
pattern = pattern.replace('-000.dat','-*.dat')

dulist = eu.hdfs.du(pattern, dict=True)

for split in xrange(nsplit):
    f = files.sample_file('lensout',run, split=split,fs='hdfs')
    if f not in dulist:
        print 'File not found:',f



