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
parser.add_option('--missing',action='store_true',
                  help='write a condor file for the missing files')
parser.add_option('--max-jobs',default=10000,
                  help='max jobs per condor file, as split on run')


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
+job_name       = "{job_name}"
Arguments       = {s2n} {npair} {output} {logfile}
Queue
"""


_master_template="""#!/bin/bash
function runsim {
    echo "host: $(hostname)"
    echo "will write to file: $output"

    python $ESPY_DIR/shapesim/bin/run-ngmix-sim.py ${run} ${s2n} ${npair} ${output}
    status=$?

    echo "time: $SECONDS"

    if [[ $status != "0" ]]; then
        echo "error running sim: $status see log ${logfile}"
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

# start with a clean log file
rm -f ${logfile}

# temporary log file, to be transferred later
log_bname=$(basename $logfile)
tmplog="$_CONDOR_SCRATCH_DIR/$log_bname"

runsim &> ${tmplog}
status=$?

echo "copying log file ${tmplog} -> ${logfile}" >> ${tmplog}

# errors go to the jobs stderr
cp "${tmplog}" "${logfile}" 1>&2
status2=$?

if [[ $status2 != "0" ]]; then
    # this error message will go to main error file
    echo "error ${status2} copying to log: ${logfile}" 1>&2

    status=$status2
fi

exit $status
"""


def get_seconds_per_pair(c):
    """
    Boosting this for the slowest machines
    """
    nstep=c['nwalkers']*(c['nstep'] + c['burnin'])
    nstep_cal=float(20*(200+400))
    if c['fit_model']=='exp':
        secper=4.0
    else:
        secper=6.0

    secper *= nstep/nstep_cal
    return secper

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

def get_flist(run):
    import glob
    fs=shapesim.get_default_fs()
    f=shapesim.get_output_url(run, 0, 0, itrial=0, fs=fs)
    d=os.path.dirname(f)

    if 'ngmix' in run:
        flist=glob.glob(d+'/*.fits')
    else:
        flist=glob.glob(d+'/*.rec')
    return flist



def write_condor_file(c, master_script, equal_time=False, missing=False, max_jobs=10000):
    run=c['run']
    overall_name = '-'.join( (run.split('-'))[1:] )

    if missing:
        flist=get_flist(run)
        overall_name += 'missing'
    else:
        flist=None


    ns2n=len(c['s2n_vals'])
    seconds_per=get_seconds_per_pair(c)
    print 'seconds per:',seconds_per

    if 'desired_err' in c:
        do_by_noise=True
    else:
        do_by_noise=False

    njobs_thisfile=0
    njobs=0

    old_filenum=-1
    filenum=0

    smax=numpy.iinfo('i8').max

    fobj=start_new_file(run, filenum, master_script, 
                        overall_name, missing=missing)
    for is2n in xrange(ns2n):

        s2n=c['s2n_vals'][is2n]

        npair, nsplit = get_npair_nsplit(c, is2n)

        time_hours = npair*seconds_per/3600.0
        if time_hours > MAXTIME_HOURS:
            raise ValueError("time is greater than %.2f "
                             "hours: %d*%.2f/3600.0 = %s" % (MAXTIME_HOURS,npair,seconds_per,time_hours))

        print '    nsplit:',nsplit,'npair:',npair,'time (hours):',time_hours

        for isplit in xrange(nsplit):
            output = shapesim.get_output_url(run, 0, is2n, itrial=isplit)
            logfile = output.replace('.fits','.log')

            this_job_name='%s-%03d-%05d' % (overall_name,is2n,isplit)
            qdata=_queue_template.format(job_name=this_job_name,
                                         s2n=s2n,
                                         npair=npair,
                                         output=output,
                                         logfile=logfile)
            do_write=True
            if missing and output in flist:
                do_write=False
            if do_write:
                njobs += 1
                njobs_thisfile += 1
                fobj.write(qdata)

            if njobs_thisfile >= max_jobs:
                filenum += 1
                njobs_thisfile=0
                fobj.close()
                fobj=start_new_file(run, filenum, master_script, 
                                    overall_name, missing=missing)


    print 'total jobs:',njobs

def start_new_file(run, filenum, master_script, overall_name, missing=False):
    condor_job_url=shapesim.get_condor_job_url(run,
                                               filenum=filenum,
                                               missing=missing)
    oname='%s-%03d' % (overall_name,filenum)
    print 'staring new job file:'
    print condor_job_url
    fobj=open(condor_job_url,'w')
    text = _condor_template_head.format(master_script=master_script,
                                        overall_name=oname)
    fobj.write(text)

    return fobj

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]

    max_jobs=int(options.max_jobs)


    c = shapesim.read_config(run)
    c['version'] = options.version

    make_some_dirs(run)
    master_script=write_master(c)

    write_condor_file(c, master_script, missing=options.missing,
                      max_jobs=max_jobs)

main()
