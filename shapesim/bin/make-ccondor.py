"""
    %prog run
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-v','--version',default='work',
                  help='priority for queue')



_condor_template_head="""
Universe        = vanilla

Notification    = Never 

# Run this exe with these args
Executable      = {master_script}


# Estimate of init. image size
Image_Size      = 14M


# You could require (CPU_Experiment == "phenix") or something like that...
#Requirements = ( CPU_Type == "mapr" )

kill_sig        = SIGINT

+Experiment     = "astro"
"""

_queue_template="""
Output          = ./logs/{job_name}.$(cluster).out
Error           = ./logs/{job_name}.$(cluster).err
Log             = ./logs/{job_name}.$(cluster).log

Arguments       = {s2n} {npair} {output}
Queue
"""

_master_template="""#!/bin/bash

source ~/.bashrc
module unload gsim_ring && module load gsim_ring/%(version)s

sim_config=%(sim_config)s
mcmc_config=%(mcmc_config)s

s2n=$1
npair=$2
output=$3

rm -vf ${output}

gsim-ring-mcmc ${sim_config} ${mcmc_config} ${s2n} ${npair}  > ${output}
status=$?

if [[ $status != "0" ]]; then
    echo "error running gsim-ring-mcmc: $status"
fi

echo "time: $SECONDS"
exit $status
"""


def get_seconds_per_pair():
    """
    need to figure this out since we don't know how to read the c config files
    """
    return 2.0

def write_master(c):
    master_url=shapesim.get_condor_master_url(c['run'])
    eu.ostools.makedirs_fromfile(master_url)

    print master_url
    with open(master_url,'w') as fobj:
        fobj.write(_master_template % c)
    print

    return master_url

def make_some_dirs(run):
    pbsd = shapesim.get_pbs_dir(run)
    logdir=os.path.join(pbsd,'logs')
    outd = shapesim.get_output_dir(run, sub='bytrial')

    if not os.path.exists(pbsd):
        os.makedirs(pbsd)
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    if not os.path.exists(outd):
        os.makedirs(outd)

def write_condor_file(c, master_script):
    run=c['run']
    job_name = '-'.join( (run.split('-'))[1:] )

    condor_job_url=shapesim.get_condor_job_url(run)

    print condor_job_url
    with open(condor_job_url,'w') as fobj:

        text = _condor_template_head.format(master_script=master_script)
        fobj.write(text)

        for is2n in xrange(ns2n):

            s2n = c['s2n_vals'][is2n]
            npair = shapesim.get_s2n_nrepeat(s2n, fac=c['s2n_fac'])
            if npair < c['min_npair']:
                npair = c['min_npair']

            time_hours = npair*seconds_per/(np-1)/60.0
            if time_hours > 2:
                raise ValueError("time is greater than two hours")

            print 'time (hours):',time_hours

            for isplit in xrange(c['nsplit']):
                output = shapesim.get_output_url(run, 0, is2n, itrial=isplit)

                this_job_name='%s-%03d'
                qdata=_queue_template.format(job_name=this_job_name,
                                             s2n=s2n,
                                             npair=npair,
                                             output=output)
                fobj.write(qdata)


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]

    seconds_per=get_seconds_per_pair()

    c = shapesim.read_config(run)
    c['version'] = options.version

    ns2n = len(c['s2n_vals'])

    make_some_dirs(run)
    master_script=write_master(c)

    write_condor_file(c, master_script)

main()
