"""
    %prog run

wq for the c code
"""

import sys
import os
import numpy

import esutil as eu
from esutil.misc import wlog
from shapesim import shapesim

from shapesim.shapesim import get_npair_nsplit

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-g','--groups',default=None,
                  help='groups for wq, csv. Default is unconstrained')
parser.add_option('-p','--priority',default='med',
                  help='priority for queue')
parser.add_option('-v','--version',default='work',
                  help='priority for queue')


_wqtemplate="""
command: |
    source ~/.bashrc
    module unload gsim_ring && module load gsim_ring/%(version)s
    output=%(output)s

    rm -f ${output}

    tmpdir=/data/esheldon/tmp
    mkdir -p $tmpdir
    tmp_output=$tmpdir/gsim-ring-mcmc-$RANDOM-$RANDOM.rec
    while [[ -e $tmp_output ]]; do
        tmp_output=$tmpdir/gsim-ring-mcmc-$RANDOM-$RANDOM.rec
    done

    s2n=%(s2n)g
    npair=%(npair)d
    seed=%(seed)d
    gsim-ring-mcmc %(sim_config)s %(mcmc_config)s ${s2n} ${npair} ${seed} > ${tmp_output}

    err="$?"
    if [[ $err != "0" ]]; then
        echo "error running gsim-ring-mcmc.  code: $err"
    else

        # try a couple of times
        echo "copying to output: ${output}"
        cp ${tmp_output} ${output}
        if [[ $? != "0" ]]; then
            echo "failed to copy, will try again in a few seconds...."
            sleep 10
            cp ${tmp_output} ${output}
        fi
    fi

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
    s2n_vals    = c['s2n_vals']

    ns2n = len(s2n_vals)

    d = shapesim.get_wq_dir(run, bytrial=True, fs='local')
    if not os.path.exists(d):
        os.makedirs(d)
    d = shapesim.get_output_dir(run, sub='bytrial')
    if not os.path.exists(d):
        os.makedirs(d)


    groups=''
    if options.groups is not None:
        groups = 'group: [%s]' % options.groups

    smax=numpy.iinfo('i8').max
    for is2n in xrange(ns2n):

        s2n = s2n_vals[is2n]

        npair, nsplit = get_npair_nsplit(c, is2n)

        for isplit in xrange(nsplit):
            job_name='%s-%03i-%03i' % (run,is2n,isplit)

            # note the wq logs are local
            wqurl = shapesim.get_wq_url(run,0,is2n,itrial=isplit,fs='local')
            output = shapesim.get_output_url(run, 0, is2n, itrial=isplit)

            seed = numpy.random.randint(smax)

            wlog("writing wq script:",wqurl)
            with open(wqurl,'w') as fobj:
                d={'job_name':job_name,
                   'version':options.version,
                   'groups':groups,
                   'pri':options.priority,
                   'sim_config':sim_config,
                   'mcmc_config':mcmc_config,
                   's2n':s2n,
                   'npair':npair,
                   'seed':seed,
                   'output':output}
                wqscript=_wqtemplate % d
                fobj.write(wqscript)

main()
