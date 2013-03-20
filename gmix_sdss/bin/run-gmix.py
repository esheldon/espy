"""
    %prog [options] gmix_run sdss_run camcol field
"""
import sys, os
import gmix_sdss

from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option('--start-field',default=None,
                  help="start full column processing at this field")

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 3:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]
    run=int(args[1])
    camcol = int(args[2])
    if len(args) > 3:
        field=int(args[3])
        pipe=gmix_sdss.pipe.GMixField(gmix_run, run, camcol, field)
        pipe.go()
    else:
        gmix_sdss.pipe.process_camcol(gmix_run=gmix_run,
                                      run=run,
                                      camcol=camcol,
                                      start_field=options.start_field)

main()
