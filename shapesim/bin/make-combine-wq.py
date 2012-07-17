"""
    %prog run
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-p','--priority',default='med',
                  help='priority for queue')

_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload wl && module load wl/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/shapesim/bin/combine-trials.py %(run)s %(i1)d %(i2)d

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

    wqd = shapesim.get_wq_dir(run, combine=True)
    if not os.path.exists(wqd):
        os.makedirs(wqd)

    n1 = cs['nums2']

    runtype = c['runtype']
    if runtype == 'byellip':
        n2 = cs['nume']
    else:
        n2 = c['nums2n']

    for i1 in xrange(n1):
        for i2 in xrange(n2):
            job_name='%s-combine-%03i-%03i' % (run,i1,i2)

            wqurl = shapesim.get_wq_url(run,i1,i2,combine=True)

            wlog("writing wq script:",wqurl)
            with open(wqurl,'w') as fobj:
                wqscript=_wqtemplate % {'job_name':job_name,
                                        'run':run, 
                                        'i1':i1,
                                        'i2':i2,
                                        'pri':options.priority}
                fobj.write(wqscript)


main()
