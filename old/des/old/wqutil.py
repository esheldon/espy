import os
from sys import stderr
import deswl

def parallelize_by_exposure(serun, type, cmd, groups=None, priority='med'):
    """
    create a job for each exposure

    parameters
    ----------
    serun:
        SE run
    type:
        the subdirectory name under ~/des-wq/{serun}-{type}/
        also job names are {expname}-{type}
    cmd:
        The command to run for each exposure.  You can expect the variables
        {expname}, {serun} to be substituded

            python $ESPY_DIR/des/bin/plot_size_mag.py --fs 'hdfs' -e {expname} {serun}

    groups: optional
        csv list of groups, default use all machines
    priority:
        default med
    """
    outdir=os.path.expanduser('~/des-wq/{serun}-{type}'.format(serun=serun,
                                                               type=type))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    sf=deswl.files.ShearFiles(serun)
    expnames = sf.get_expnames()
    for expname in expnames:

        command=cmd.format(serun=serun,
                           expname=expname)
        job_name=expname+'-'+type

        if groups is None:
            grps=''
        else:
            grps='group: [' + groups +']'

        text="""
command: |
    #source /opt/astro/SL53/bin/setup.hadoop.sh
    #source ~astrodat/setup/setup.sh
    #source ~/.dotfiles/bash/astro.bnl.gov/modules.sh
    source ~/.bashrc
    module unload esutil && module load esutil/work
    module unload tmv    && module load tmv/work
    module unload wl     && module load wl/work

    {command} &> {expname}-{type}.out

{groups}
priority: {priority}
job_name: {job_name}\n""".format(command=command,
                                 expname=expname,
                                 type=type,
                                 groups=grps, 
                                 priority=priority,
                                 job_name=job_name)


        outfile='{expname}-{type}.yaml'.format(expname=expname,
                                               type=type)
        outfile=os.path.join(outdir,outfile)
        print >>stderr,'Writing submit file:',outfile
        with open(outfile,'w') as fobj:
            fobj.write(text)
