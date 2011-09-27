"""
    %prog [options] run

Description:

    Reduce the lenses into a single lensout file.

"""

from __future__ import print_function

import os
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("--fs",default="hdfs",
                  help="File system to read from, 'nfs','hdfs','local'.  Default '%default'")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    run = args[0]
    fs=options.fs


    conf = lensing.files.cascade_config(run)
    nsplit = conf['src_config']['nsplit']

    outfile = lensing.files.sample_file('reduced', run)
    print("Will combine into file:",outfile)
    print("Combining",nsplit,"splits from run",run)

    d=lensing.files.sample_dir('lensout',run)
    if not os.path.exists(d):
        print("Making output dir:",d)
        os.makedirs(d)

    for i in xrange(nsplit):
        tdata = lensing.files.sample_read('lensout',run, split=i, fs=fs)
        if i == 0:
            data = tdata
        else:
            print("summing")
            lensing.outputs.add_lensums(data, tdata)
        del tdata


    lensing.files.sample_write(data, 'reduced', run, clobber=True)
    

main()

