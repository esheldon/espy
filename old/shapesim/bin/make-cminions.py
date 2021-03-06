"""
    %prog run nodes ppn
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-v','--version',default='work',
                  help='priority for queue')



_pbs_template="""#!/bin/bash -l
#PBS -N {job_name}
#PBS -j oe
#PBS -l nodes={nodes}:ppn={ppn},walltime={hours}:00:00
#PBS -q regular
#PBS -o {pbslog}
#PBS -A des

if [[ "Y${{PBS_O_WORKDIR}}" != "Y" ]]; then
    echo "moving do directory: $PBS_O_WORKDIR"
    cd $PBS_O_WORKDIR
fi


module load openmpi-gnu
mpirun -np {np} {minions_cmd} < {commands_file}

echo "done minions"
"""

_script_template="""#!/bin/bash

nsetup_ess

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



_cmd_template="""
bash ./{script_file} {s2n} {npair} {output} &> logs/{logfile}
""".strip()


def get_seconds_per_pair():
    """
    need to figure this out since we don't know how to read the c config files
    """
    return 2.0

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    nodes=int(args[1])
    ppn=int(args[2])

    np=nodes*ppn
    # seconds per ellip ring pair. this is for exp...

    seconds_per=get_seconds_per_pair()

    c = shapesim.read_config(run)

    sim_config  = c['sim_config']
    mcmc_config = c['mcmc_config']
    nsplit      = c['nsplit']
    s2n_vals    = c['s2n_vals']
    s2n_fac     = c['s2n_fac']
    min_npair   = c['min_npair']

    ns2n = len(s2n_vals)

    pbsd = shapesim.get_pbs_dir(run)
    logdir=os.path.join(pbsd,'logs')
    outd = shapesim.get_output_dir(run, sub='bytrial')

    if not os.path.exists(pbsd):
        os.makedirs(pbsd)
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    if not os.path.exists(outd):
        os.makedirs(outd)

    rstr=run

    job_name = '-'.join( (rstr.split('-'))[1:] )

    script_url=shapesim.get_minions_script_url(run)
    script_base=os.path.basename(script_url)

    print script_url
    with open(script_url,'w') as fobj:
        c['version'] = options.version
        fobj.write(_script_template % c)
    print

    pbsf=shapesim.get_minions_url(run,0)
    cmdf=shapesim.get_commands_url(run,0)

    print cmdf
    ntot=0
    with open(cmdf,'w') as fobj:

        for is2n in xrange(ns2n):

            s2n = s2n_vals[is2n]
            npair = shapesim.get_s2n_nrepeat(s2n, fac=s2n_fac)
            if npair < min_npair:
                npair = min_npair

            ntot += npair*nsplit

            for isplit in xrange(nsplit):
                output = shapesim.get_output_url(run, 0, is2n, itrial=isplit)
                outbname=os.path.basename(output)

                logf=outbname.replace('.rec','.log')
                cmd=_cmd_template.format(script_file=script_base,
                                         s2n=s2n,
                                         npair=npair,
                                         output=output,
                                         logfile=logf)
                fobj.write(cmd)
                fobj.write('\n')


    time_seconds = ntot*seconds_per/(np-1)

    print pbsf
    hours_raw = time_seconds/3600.
    hours = int(round(hours_raw)) + 1
    print '    raw hours:',hours_raw,'hours used:',hours

    if options.version=='hopper':
        minions_cmd='minions-hopper'
    else:
        minions_cmd='minions'

    with open(pbsf,'w') as fobj:
        cmd_bname=os.path.basename(cmdf)
        pbslog=os.path.basename(pbsf)+'.pbslog'
        pbs_text=_pbs_template.format(minions_cmd=minions_cmd,
                                      job_name=job_name,
                                      pbslog=pbslog,
                                      nodes=nodes,
                                      ppn=ppn,
                                      np=np,
                                      hours=hours,
                                      commands_file=cmd_bname)
        fobj.write(pbs_text)
main()
