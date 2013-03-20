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
parser.add_option('--run-admom',action='store_true',
                  help="Force a re-run of adaptive moments")
parser.add_option('--run-psf',action='store_true',
                  help="Force a re-run of psf measurement")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 4:
        parser.print_help()
        sys.exit(1)

    run=args[0]
    psfnum=int(args[1])
    shnum=int(args[2])
    ccd=int(args[3])

    if 'umom' in run:
        pipe=cluster_step.pipe.MomentPipe(run=run,
                                          psfnum=psfnum,
                                          shnum=shnum,
                                          ccd=ccd)

        pipe.run()

    elif 'fftstack' in run:
        pipe=cluster_step.pipe.FFTStackPipe(run=run,
                                            psfnum=psfnum,
                                            shnum=shnum,
                                            ccd=ccd)

        pipe.run()

    elif 'stack' in run:
        pipe=cluster_step.pipe.StackPipe(run=run,
                                         psfnum=psfnum,
                                         shnum=shnum,
                                         ccd=ccd)

        pipe.run()

    else:
        pipe=cluster_step.pipe.Pipe(run=run,
                                    psfnum=psfnum,
                                    shnum=shnum,
                                    ccd=ccd)

        pipe.run_shear(run_admom=options.run_admom,
                       run_psf=options.run_psf)

main()