"""
    %prog run psfnum shnum ccd

Run the adaptive moment code and write the results as well as a size-magnitude
diagram.

"""
import sys
import cluster_step

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-s','--show',action='store_true',
                  help="show the plot on the screen")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 4:
        parser.print_help()
        sys.exit(1)

    run=args[0]
    psfnum=int(args[1])
    shnum=int(args[2])
    ccd=int(args[3])

    pipe=cluster_step.pipe.Pipe(run=run,
                                psfnum=psfnum,
                                shnum=shnum,
                                ccd=ccd)
    pipe.run_admom()
    pipe.plot_admom_sizemag(show=options.show)

main()
