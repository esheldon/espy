"""
    %prog run
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-g','--groups',default=None,
                  help='groups for wq, csv. Default is unconstrained')
parser.add_option('-p','--priority',default='med',
                  help='priority for queue')
parser.add_option('--bynode',action='store_true',
                  help='select full nodes')
parser.add_option('--ncores',default=None,
                  help='select this many cores on single nodes')

parser.add_option('--bytrial',action='store_true',
                  help=('Write a wq script for each trial '
                        'separately (with possible repeats)'))
_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload wl && module load wl/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/shapesim/bin/run-shapesim.py %(run)s %(is2)d %(ie)d

%(extra)s
%(groups)s
job_name: %(job_name)s
priority: %(pri)s\n"""

_wqtemplate_bytrial="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload wl && module load wl/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/shapesim/bin/run-shapesim.py %(run)s %(is2)d %(ie)d %(itrial)d

%(extra)s
%(groups)s
job_name: %(job_name)s
priority: %(pri)s\n"""


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]

    c = shapesim.read_config(run)
    cs = shapesim.read_config(c['sim'])

    if options.bytrial:
        orient=cs.get('orient','rand')
        if orient == 'ring':
            ntrial = cs['nring']
        else:
            ntrial = c['ntrial']

    groups=options.groups
    if groups is None:
        groups=''
    else:
        groups = 'group: [%s]' % groups

    wqd = shapesim.get_wq_dir(run, bytrial=options.bytrial)
    if not os.path.exists(wqd):
        os.makedirs(wqd)

    extra=''
    if options.bynode:
        extra='mode: bynode\nN: 1'
    elif options.ncores is not None:
        ncores=int(options.ncores)
        extra='mode: bycore1\nN: %d' % ncores

    if run[0:8] == 'gmix-fit':
        rstr=run.replace('gmix-fit','gmix')
    else:
        rstr=run

    for is2 in xrange(cs['nums2']):
        for ie in xrange(cs['nume']):
            job_name='%s-%03i-%03i' % (rstr,is2,ie)

            if options.bytrial:
                for itrial in xrange(ntrial):
                    wqurl = shapesim.get_wq_url(run,is2,ie,itrial=itrial)
                    wlog("writing wq script:",wqurl)
                    with open(wqurl,'w') as fobj:
                        d={'job_name':job_name,
                           'run':run, 
                           'is2':is2,
                           'ie':ie,
                           'itrial':itrial,
                           'groups':groups,
                           'extra':extra,
                           'pri':options.priority}
                        wqscript=_wqtemplate_bytrial % d
                        fobj.write(wqscript)


            else:
                wqurl = shapesim.get_wq_url(run,is2,ie)

                wlog("writing wq script:",wqurl)
                with open(wqurl,'w') as fobj:
                    wqscript=_wqtemplate % {'job_name':job_name,
                                            'run':run, 
                                            'is2':is2,
                                            'ie':ie,
                                            'groups':groups,
                                            'extra':extra,
                                            'pri':options.priority}
                    fobj.write(wqscript)


main()
