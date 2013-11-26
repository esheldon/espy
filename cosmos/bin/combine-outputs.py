"""
    %prog version
"""
import sys, os
from sys import stderr
import cosmos

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    version=args[0]

    cosmos.files.combine_outputs(version)

main()
