"""
    %prog version model1 model2 ...
"""
import sys, os
from sys import stderr
import cosmos

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    version=args[0]
    models=args[1:]

    output_dir=cosmos.files.get_gmix_output_dir(version)
    if not os.path.exists(output_dir):
        print >>stderr,'making dir:',output_dir
        os.makedirs(output_dir)

    cosmos.files.make_master(version,models)
    cosmos.files.make_condor_script(version)

main()
