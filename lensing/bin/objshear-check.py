"""
    %prog [options] run

Description:
    Check all expected outputs exist for the input run

    Also verified they have the same number of lines (nlens)
"""
import sys
from sys import stderr
import os
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

pattern = files.sample_file(type='lensout',sample=run, split=0,fs='hdfs')
pattern = pattern.replace('-000.dat','-*.dat')

dulist = eu.hdfs.du(pattern, dict=True)

nlines_expected =-1
for split in xrange(nsplit):
    f = files.sample_file(type='lensout',sample=run, split=split,fs='hdfs')
    stderr.write('.')
    if f not in dulist:
        print >>stderr,'\nFile not found:',f
        continue

    # now count lines
    with eu.hdfs.HDFSFile(f) as fobj:
        fobj.stage()

        nlines = os.popen('cat %s | wc -l' % fobj.localfile).read().strip()
        nlines = int(nlines)

        if nlines_expected < 0:
            nlines_expected = nlines

        if nlines != nlines_expected:
            print >>stderr,'\nfound file with different line count:',f,"(%s/%s)" % (nlines,nlines_expected)

print '\ndone'
