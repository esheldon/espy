"""
    %prog [options] type
"""
import sys
import sdsspy

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-r","--runs",
                  dest="runs",
                  default=None,
                  help="Comma separated list of runs to stuff. Default is all")

parser.add_option("-i","--indices",
                 action="store_true",
                 dest="create_indices",
                 default=False,
                 help="make the indices. Default %default")

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 1:
        parser.print_help()
        return
    type = args[0]
    runs = options.runs
    create_indices = options.create_indices


    if runs is not None:
        runs = runs.split(',')
        runs = [int(run) for run in runs]


    # put it on the same array where the server is
    tmpdir = '/export/tutti1/esheldon/tmp'
    ss = sdsspy.pg.SweepStuffer(type, tmpdir=tmpdir)
    ss.stuff_runs(runs)


    if create_indices:
        ss.create_indices()

if __name__=="__main__":
    main()
