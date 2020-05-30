"""
    %prog run

Makes wq scripts for all profiles, ellip and resolutions
"""

import sys
import os
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-g','--groups',default=None,
                  help='groups for wq, csv. Default is unconstrained')
parser.add_option('-p','--priority',default='med',
                  help='priority for queue')

#module unload wl && module load wl/work
_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/shapesim/bin/run-galsim-fit.py %(run)s %(profile)s %(ellip)d %(res)d

%(groups)s
job_name: %(job_name)s
priority: %(pri)s\n"""


def get_wq_dir(run):
    return os.path.expanduser('~/galsim-outputs/wq/%s' % run)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]

    wqd = get_wq_dir(run)
    if not os.path.exists(wqd):
        os.makedirs(wqd)

    extra=''

    groups=''
    if options.groups is not None:
        groups = 'group: [%s]' % options.groups

    for p in [1,2,3,4,5,6]:
        profile='p%02i' % p
        for ellip in [0,3,6]:
            for res in [4,8,12,16,20]:

                job_name='%s-%s-%03i-%03i' % (run,profile,ellip,res)
                wqurl=os.path.join(wqd, job_name+'.yaml')

                print "writing wq script:",wqurl
                with open(wqurl,'w') as fobj:
                    wqscript=_wqtemplate % {'job_name':job_name,
                                            'run':run, 
                                            'profile':profile,
                                            'ellip':ellip,
                                            'res':res,
                                            'groups':groups,
                                            'pri':options.priority}
                    fobj.write(wqscript)


main()
