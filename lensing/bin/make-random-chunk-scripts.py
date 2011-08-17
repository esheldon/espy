"""
    %prog sample nchunk

Description:

    Make condor scripts for generation of random lcat in chunks.

"""
from __future__ import print_function
import sys, os
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

sample = args[0]
nchunk = int(args[1])

c = lensing.lcat.instantiate_sample(sample)


nperchunk = c['nrand']/nchunk
nleft = c['nrand'] % nchunk

dir = lensing.files.sample_dir('lcat',sample)
dir=os.path.join(dir, 'condor')
if not os.path.exists(dir):
    os.makedirs(dir)

for i in xrange(nchunk):
    extra_name='chunk%02d' % i
    script_base='gen-%s' % extra_name
    script_fname=os.path.join(dir, script_base+'.sh')
    condor_fname=os.path.join(dir, script_base+'.condor')
    print(condor_fname)

    nrand = nperchunk
    if (i == (nchunk-1)):
        nrand += nleft

    script="""#!/bin/bash
script_dir=~esheldon/python/lensing/bin
python -u $script_dir/make-objshear-input.py -n {nrand} -e {extra_name} lcat {sample} {nchunk}
echo Done
    \n""".format(sample=sample,
                 nchunk=nchunk,
                 nrand=nrand,
                 extra_name=extra_name)

    condor_script="""
Universe        = vanilla
Notification    = Error
GetEnv          = True
Notify_user     = esheldon@bnl.gov
Requirements    = (CPU_Experiment == "astro") && (TotalSlots == 12 || TotalSlots == 8) && (Machine  != "astro0033.rcf.bnl.gov")
+Experiment     = "astro"
Initialdir      = {dir}

Executable      = {script_base}.sh
Output          = {script_base}.out
Error           = {script_base}.err
Log             = {script_base}.log

Queue
    
    \n""".format(dir=dir, script_base=script_base)



    with open(script_fname,'w') as fobj:
        fobj.write(script)
    with open(condor_fname,'w') as fobj:
        fobj.write(condor_script)

    os.system('chmod 755 '+script_fname)
