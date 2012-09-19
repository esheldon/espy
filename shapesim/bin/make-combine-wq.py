"""
    %prog run

Write wq scripts to run combine-trials.py to combine the "bytrial"
files.
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-p','--priority',default='med',
                  help='priority for queue')
parser.add_option('--bynode',action='store_true',
                  help="grab a full node")
parser.add_option('-g','--group',default=None,
                  help='groups/groups')
parser.add_option('-e','--extra',default=None,
                  help='; separated')

#module unload wl && module load wl/work
_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/shapesim/bin/combine-trials.py %(run)s %(i1)d %(i2)d

%(group)s
%(mode)s
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

    group=options.group
    if group is not None:
        group = 'group: ['+group+']'
    else:
        group=''

    mode=''
    if options.bynode:
        mode='mode: bynode'

    extra=''
    if options.extra is not None:
        extra='\n'.join( options.extra.split(';') )

    wqd = shapesim.get_wq_dir(run, combine=True)
    if not os.path.exists(wqd):
        os.makedirs(wqd)

    if run[0:8] == 'gmix-fit':
        rstr=run.replace('gmix-fit','gmix')
    else:
        rstr=run
 
    n1 = cs['nums2']

    runtype = c['runtype']
    if runtype == 'byellip':
        n2 = cs['nume']
    else:
        n2 = shapesim.get_nums2n(c)

    for i1 in xrange(n1):
        for i2 in xrange(n2):
            job_name='%s-combine-%03i-%03i' % (rstr,i1,i2)

            wqurl = shapesim.get_wq_url(run,i1,i2,combine=True)

            wlog("writing wq script:",wqurl)
            with open(wqurl,'w') as fobj:
                wqscript=_wqtemplate % {'job_name':job_name,
                                        'run':run, 
                                        'i1':i1,
                                        'i2':i2,
                                        'group':group,
                                        'mode':mode,
                                        'extra':extra,
                                        'pri':options.priority}
                fobj.write(wqscript)


main()
