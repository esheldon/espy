"""
    %prog simname
"""

import sys
import os
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    simname=args[0]
    
    shapesim.dessim.pointings.make_pointings(simname)
    

main()
