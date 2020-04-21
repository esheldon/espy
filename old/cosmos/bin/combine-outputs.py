"""
    %prog config_file
"""
import sys, os
from sys import stderr
import cosmos
import yaml

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    config_file=args[0]
    conf=yaml.load(open(config_file))
    version=conf['version']

    cosmos.files.combine_outputs(version)

main()
