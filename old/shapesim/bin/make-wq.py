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

parser.add_option('--i2new',default=None,
                  help='use new,new2 groups for i2 in [0,value]')

parser.add_option('--i2new1',default=None,
                  help='use new group for i2 in [0,value]')
parser.add_option('--i2new2',default=None,
                  help='use new2 group for i2 in [0,value]')

parser.add_option('--bytrial',action='store_true',
                  help=('Write a wq script for each trial '
                        'separately (with possible repeats)'))
#module unload wl && module load wl/work
_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/shapesim/bin/run-shapesim.py %(run)s %(i1)d %(i2)d

%(extra)s
%(groups)s
job_name: %(job_name)s
priority: %(pri)s\n"""

_wqtemplate_bytrial="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/shapesim/bin/run-shapesim.py %(run)s %(i1)d %(i2)d %(itrial)d

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
        if 'stack-' in run:
            ntrial=c['nsplit']
        else:
            orient=cs.get('orient','rand')
            if orient == 'ring':
                ntrial = cs['nsplit']
            else:
                ntrial = cs['ntrial']


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

    n1 = shapesim.get_numT(cs)

    runtype = c['runtype']
    if runtype == 'byellip':
        n2 = cs['nume']
    else:
        n2 = shapesim.get_nums2n(c)



    for i1 in xrange(n1):
        for i2 in xrange(n2):
            groups=''
            if options.i2new is not None and i2 <= int(options.i2new):
                groups='group: [new,new2]'
            elif options.i2new1 is not None and i2 <= int(options.i2new1):
                groups='group: [new]'
            elif options.i2new2 is not None and i2 <= int(options.i2new2):
                groups='group: [new2]'
            elif options.groups is not None:
                groups = 'group: [%s]' % options.groups

            if options.bytrial:
                for itrial in xrange(ntrial):
                    job_name='%s-%03i-%03i-%02i' % (rstr,i1,i2,itrial)
                    wqurl = shapesim.get_wq_url(run,i1,i2,itrial=itrial)
                    wlog("writing wq script:",wqurl)
                    with open(wqurl,'w') as fobj:
                        d={'job_name':job_name,
                           'run':run, 
                           'i1':i1,
                           'i2':i2,
                           'itrial':itrial,
                           'groups':groups,
                           'extra':extra,
                           'pri':options.priority}
                        wqscript=_wqtemplate_bytrial % d
                        fobj.write(wqscript)


            else:
                job_name='%s-%03i-%03i' % (rstr,i1,i2)
                wqurl = shapesim.get_wq_url(run,i1,i2)

                wlog("writing wq script:",wqurl)
                with open(wqurl,'w') as fobj:
                    wqscript=_wqtemplate % {'job_name':job_name,
                                            'run':run, 
                                            'i1':i1,
                                            'i2':i2,
                                            'groups':groups,
                                            'extra':extra,
                                            'pri':options.priority}
                    fobj.write(wqscript)


main()
