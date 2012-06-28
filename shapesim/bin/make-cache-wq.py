"""
    %prog simname numper
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim
import esutil as eu

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

_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload wl && module load wl/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/shapesim/bin/cache-sim.py %(simname)s %(is2)d %(ie)d %(numper)d

%(extra)s
%(groups)s
job_name: %(job_name)s
priority: %(pri)s\n"""

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    simname=args[0]
    numper=int(args[1])

    ss=shapesim.ShapeSim(simname)
    if ss.fs != 'hdfs':
        raise ValueError("This only works for HDFS right now "
                         "would need to worry about making dirs")

    c = shapesim.read_config(simname)

    groups=options.groups
    if groups is None:
        groups=''
    else:
        groups = 'group: [%s]' % groups

    wqd = shapesim.get_cache_wq_dir(simname)
    if not os.path.exists(wqd):
        os.makedirs(wqd)

    extra=''
    if options.bynode:
        extra='mode: bynode\nN: 1'
    elif options.ncores is not None:
        ncores=int(options.ncores)
        extra='mode: bycore1\nN: %d' % ncores

    for is2 in xrange(c['nums2']):
        for ie in xrange(c['nume']):
            job_name='%s-%i-%i' % (simname,is2,ie)

            wqurl = shapesim.get_cache_wq_url(simname,is2,ie)

            wlog("writing wq script:",wqurl)
            with open(wqurl,'w') as fobj:
                wqscript=_wqtemplate % {'job_name':job_name,
                                        'simname':simname, 
                                        'is2':is2,
                                        'ie':ie,
                                        'numper':numper,
                                        'groups':groups,
                                        'extra':extra,
                                        'pri':options.priority}
                fobj.write(wqscript)


main()
