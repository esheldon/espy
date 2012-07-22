"""
    %prog new_run_name run1 run2 ..

Description

    Average the set of runs
"""

import sys
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(45)


    new_run_name=args[0]
    runs2average = args[1:]

    print 'new run name:',new_run_name
    print 'runs2average:',runs2average
    shapesim.shapesim.average_runs(runs2average, new_run_name)

main()
