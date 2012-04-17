#!/usr/bin/env python
"""
    %prog [options] imfile catfile amfile
"""
import os
import sys
from sys import stderr
import deswl
import lensing
import des
import esutil as eu
import numpy

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) != 3:
        parser.print_help()
        sys.exit(45)

    imfile  = args[0]
    catfile = args[1]
    amfile  = args[2]
    out, meta = des.admom_des.process_image_cat(imfile,
                                                catfile,
                                                get_meta=True)

    eu.io.write(amfile, out, header=meta)

main()
