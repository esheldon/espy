"""
    %prog run nodes
"""

import sys
import os
from esutil.misc import wlog
from shapesim import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-p','--ppn',default=8, help='processors per node')


_pbs_template="""
#!/bin/bash -l
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
mpirun -np {np} minions < {commands_file}

echo "done minions"
"""

_script_template="""
#!/bin/bash

run=$1
i1=$2
i2=$3
itrial=$4

nsetup_ess

module unload espy && module load espy/work
module unload fimage && module load fimage/work
module unload gmix_image && module load gmix_image/work

python $ESPY_DIR/shapesim/bin/run-shapesim.py $run $i1 $i2 $itrial
"""


#_cmd_template="""
#source ~/.bashrc && nsetup_ess && module unload espy && module load espy/work && module unload fimage && module load fimage/work && module unload gmix_image && module load gmix_image/work && python $ESPY_DIR/shapesim/bin/run-shapesim.py {run} {i1} {i2} {itrial} &> logs/{logfile}
#""".strip()

_cmd_template="""
bash ./{script_file} {run} {i1} {i2} {itrial} &> logs/{logfile}
""".strip()


def get_nellip(conf, is2n):
    s2n = shapesim.get_s2n(conf, is2n)
    s2n_fac = conf['s2n_fac']
    nellip = shapesim.get_s2n_nrepeat(s2n, fac=s2n_fac)

    if nellip < conf['min_gcount']:
        nellip=conf['min_gcount']
    return nellip


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    nodes=int(args[1])
    ppn=int(options.ppn)

    np=nodes*ppn
    # seconds per ellip ring pair. this is for exp...

    c = shapesim.read_config(run)
    cs = shapesim.read_config(c['sim'])

    seconds_per=0.0
    for model in c['fitmodel']:
        if model=='gexp':
            seconds_per += 3.0
        else:
            seconds_per += 3.2
    nsplit = cs['nsplit']

    pbsd = shapesim.get_pbs_dir(run)
    logdir=os.path.join(pbsd,'logs')
    if not os.path.exists(pbsd):
        os.makedirs(pbsd)
    if not os.path.exists(logdir):
        os.makedirs(logdir)

    if run[0:8] == 'gmix-fit':
        rstr=run.replace('gmix-fit','gmix')
    else:
        rstr=run

    job_name = '-'.join( (rstr.split('-'))[1:] )

    n1 = cs['nums2']
    n2 = shapesim.get_nums2n(c)


    script_url=shapesim.get_minions_script_url(run)
    script_base=os.path.basename(script_url)

    print script_url
    with open(script_url,'w') as fobj:
        fobj.write(_script_template)
    print

    for i1 in xrange(n1):
        pbsf=shapesim.get_minions_url(run,i1)
        cmdf=shapesim.get_commands_url(run,i1)

        print cmdf
        ntot=0
        with open(cmdf,'w') as fobj:
            for i2 in xrange(n2):
                nellip=get_nellip(c,i2)

                for itrial in xrange(nsplit):
                    outf=shapesim.get_output_url(run, i1, i2, itrial=itrial)
                    outbname=os.path.basename(outf)
                    logf=outbname.replace('.rec','.log')
                    cmd=_cmd_template.format(script_file=script_base,
                                             run=run,
                                             i1=i1,
                                             i2=i2,
                                             itrial=itrial,
                                             logfile=logf)

                    fobj.write(cmd)
                    fobj.write('\n')

                    ntot += nellip

        time_seconds = ntot*seconds_per/np

        print pbsf
        hours_raw = time_seconds/3600.
        hours = int(round(hours_raw)) + 1
        print '    raw hours:',hours_raw,'hours used:',hours

        with open(pbsf,'w') as fobj:
            cmd_bname=os.path.basename(cmdf)
            pbslog=os.path.basename(pbsf)+'.pbslog'
            pbs_text=_pbs_template.format(job_name=job_name,
                                          pbslog=pbslog,
                                          nodes=nodes,
                                          ppn=ppn,
                                          np=np,
                                          hours=hours,
                                          commands_file=cmd_bname)
            fobj.write(pbs_text)
main()
