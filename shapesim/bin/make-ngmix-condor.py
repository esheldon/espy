"""
    %prog run
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

parser.add_option('-v','--version',default='work',
                  help='priority for queue')

MAXTIME_HOURS=1.5

_condor_template_head="""
Universe        = vanilla

Notification    = Never 

# Run this exe with these args
Executable      = {master_script}


# Estimate of init. image size.  This is actually high
# condor always raises it to 60-90 Meg
Image_Size      = 90000

GetEnv = True

kill_sig        = SIGINT

+Experiment     = "astro"

Output          = ./{overall_name}.$(cluster).out
Error           = ./{overall_name}.$(cluster).err
Log             = /data/esheldon/tmp/{overall_name}.$(cluster).log

"""

_queue_template="""
Arguments       = {s2n} {npair} {output} {logfile}
Queue
"""


_master_template="""#!/bin/bash
function runsim {
    host=$(hostname)


    echo "host: $host"
    echo "writing to temp file: $tmpfile"

    python $ESPY_DIR/shapesim/bin/run-ngmix-sim.py ${run} ${s2n} ${npair} ${output}
    status=$?

    echo "time: $SECONDS"

    if [[ $status != "0" ]]; then
        echo "error running sim: $status"
    else
        cp -v "$tmpfile" "$output"
        status=$?
        if [[ $status != "0" ]]; then
            echo "error copying to output: $output"
        fi
    fi

    return $status
}

source ~/.bashrc
module unload gsim_ring && module load gsim_ring/%(version)s
module unload espy && module load espy/%(version)s


s2n=$1
npair=$2
output=$3
logfile=$4

run=%(run)s

bname=$(basename $output)
log_bname=$(basename $logfile)
tmpfile="$_CONDOR_SCRATCH_DIR/$bname"
tmplog="$_CONDOR_SCRATCH_DIR/$log_bname"

rm -vf ${logfile}
rm -vf ${output}

runsim &> $tmplog
status=$?

echo "copying log file $tmplog -> $logfile" >> $tmplog

# any errors will go to the jobs stderr
cp "$tmplog" "$logfile" 1>&2
status2=$?

if [[ $status != "0" ]]; then
    # this error message will go to main error file
    echo "error copying to output: $output" 1>&2

    status=$status2
fi

exit $status
"""


def get_seconds_per_pair():
    """
    Boosting this for the slowest machines
    """
    return 4.0

def write_master(c):
    master_url=shapesim.get_condor_master_url(c['run'])
    eu.ostools.makedirs_fromfile(master_url)

    print master_url
    with open(master_url,'w') as fobj:
        fobj.write(_master_template % c)
    print

    os.system('chmod 755 %s' % master_url)
    return master_url

def make_some_dirs(run):
    d = shapesim.get_condor_dir(run)
    outd = shapesim.get_output_dir(run, sub='bytrial')

    if not os.path.exists(d):
        os.makedirs(d)

    if not os.path.exists(outd):
        os.makedirs(outd)

    tmpdir='/data/esheldon/tmp'
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

def write_condor_file(c, master_script, equal_time=False):
    run=c['run']
    job_name = '-'.join( (run.split('-'))[1:] )

    condor_job_url=shapesim.get_condor_job_url(run)
    ns2n=len(c['s2n_vals'])
    seconds_per=get_seconds_per_pair()

    if 'desired_err' in c:
        do_by_noise=True
    else:
        do_by_noise=False

    print condor_job_url
    njobs=0

    smax=numpy.iinfo('i8').max
    with open(condor_job_url,'w') as fobj:

        text = _condor_template_head.format(master_script=master_script,
                                            overall_name=job_name)
        fobj.write(text)

        for is2n in xrange(ns2n):
            s2n=c['s2n_vals'][is2n]

            npair, nsplit = get_npair_nsplit(c, is2n)

            time_hours = npair*seconds_per/3600.0
            if time_hours > MAXTIME_HOURS:
                raise ValueError("time is greater than %.2f "
                                 "hours: %d*%.2f/3600.0 = %s" % (MAXTIME_HOURS,npair,seconds_per,time_hours))

            print 'nsplit:',nsplit,'npair:',npair,'time (hours):',time_hours
            njobs += nsplit

            for isplit in xrange(nsplit):
                output = shapesim.get_output_url(run, 0, is2n, itrial=isplit)
                logfile = output.replace('.rec','.log')

                this_job_name='%s-%03d-%03d' % (job_name,is2n,isplit)
                qdata=_queue_template.format(job_name=this_job_name,
                                             s2n=s2n,
                                             npair=npair,
                                             output=output,
                                             logfile=logfile)
                fobj.write(qdata)


    print 'total jobs:',njobs

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]


    c = shapesim.read_config(run)
    c['version'] = options.version

    make_some_dirs(run)
    master_script=write_master(c)

    write_condor_file(c, master_script)

main()
