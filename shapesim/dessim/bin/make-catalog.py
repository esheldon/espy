"""
    %prog simname pointing_id
"""

import sys
import os
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-a','--ascii',
                  action='store_true',
                  help="write an ascii file")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    simname=args[0]
    pointing_id=int(args[1])
    
    maker=shapesim.dessim.catmaker.SimpleCatalogMaker(simname, pointing_id)
    maker.go()
    if options.ascii:
        maker.write_ascii()
    else:
        maker.write()
     

main()
