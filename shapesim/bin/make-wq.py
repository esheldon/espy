"""
    %prog run
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-g','--groups',default='gen345',help='groups for wq, csv')

_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload fimage && module load fimage/work
    module unload wl && module load wl/work
    python $ESPY_DIR/shapesim/bin/run-shapesim.py %(run)s %(is2)d %(ie)d

group: [%(groups)s]
job_name: %(job_name)s
priority: med\n"""

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]

    c = shapesim.read_config(run)
    groups=options.groups

    wqd = shapesim.get_wq_dir(run)
    if not os.path.exists(wqd):
        os.makedirs(wqd)
    # make this to avoid race conditions later
    od = shapesim.get_output_dir(run)
    if not os.path.exists(od):
        os.makedirs(od)

    for is2 in xrange(c['nums2']):
        for ie in xrange(c['nume']):
            job_name='%s-%i-%i' % (run,is2,ie)

            wqurl = shapesim.get_wq_url(run,is2,ie)

            wlog("writing wq script:",wqurl)
            with open(wqurl,'w') as fobj:
                wqscript=_wqtemplate % {'job_name':job_name,
                                        'run':run, 
                                        'is2':is2,
                                        'ie':ie,
                                        'groups':groups}
                fobj.write(wqscript)


main()
