#!/usr/bin/env python
"""
    %prog [options] run mock_catalog

Plot comparisons of shears

mock_catalog should be e.g. 'desmocks-3.02'
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

    if len(args) != 2:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    mock_catalog=args[1]

    c=des.compare_truth.Comparator(run,mock_catalog) 

    c.plot_sheardiff()

main()
