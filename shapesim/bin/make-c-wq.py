"""
    %prog run

wq for the c code
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


_wqtemplate="""
command: |
    source ~/.bashrc
    module unload gsim_ring && module load gsim_ring/work
    gsim-ring-mcmc %(sim_config)s %(mcmc_config)s %(s2n)g %(npair)d

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

    run         = c['run']
    sim_config  = c['sim_config']
    mcmc_config = c['mcmc_config']
    nsplit      = c['nsplit']
    s2n_vals    = c['s2n_vals']
    s2n_fac     = c['s2n_fac']

    ns2n = len(s2n_vals)

    wqd = shapesim.get_wq_dir(run, bytrial=True)
    if not os.path.exists(wqd):
        os.makedirs(wqd)

    groups=''
    if options.groups is not None:
        groups = 'group: [%s]' % options.groups

    for is2n in xrange(ns2n):

        s2n = s2n_vals[is2n]
        npair = shapesim.get_s2n_nrepeat(s2n, fac=s2n_fac)

        for isplit in xrange(nsplit):
            job_name='%s-%02i-%03i' % (run,is2n,isplit)
            wqurl = shapesim.get_wq_url(run,0,is2n,itrial=isplit)
            wlog("writing wq script:",wqurl)
            with open(wqurl,'w') as fobj:
                d={'job_name':job_name,
                   'groups':groups,
                   'pri':options.priority,
                   'sim_config':sim_config,
                   'mcmc_config':mcmc_config,
                   's2n':s2n,
                   'npair':npair}
                wqscript=_wqtemplate % d
                fobj.write(wqscript)

main()
