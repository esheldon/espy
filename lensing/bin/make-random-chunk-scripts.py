"""
    %prog sample

Description:

    Make wq scripts for generation of random lcat in chunks.
"""
from __future__ import print_function
import sys, os
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option("--gen-radec", action='store_true', 
                  help=("For this sample we generate the ra,dec points"))

parser.add_option("-g", "--group", default=None, help=("Groups for wq"))

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(1)

sample = args[0]

c = lensing.lcat.instantiate_sample(sample)

nsplit = c['nsplit']
nperchunk = c['nrand']/nsplit
nleft = c['nrand'] % nsplit

dir='/data/esheldon/lensing/lcat/%s/wq' % sample
if not os.path.exists(dir):
    os.makedirs(dir)

group=options.group
if group:
    group='group: ['+group+']'
else:
    group=''

for i in xrange(nsplit):
    chunk_name='chunk%03d' % i
    script_base='gen-%s' % chunk_name
    wq_fname=os.path.join(dir, script_base+'.yaml')
    print(wq_fname)

    if options.gen_radec:
        nrand = nperchunk
        if (i == (nsplit-1)):
            nrand += nleft
        opts='-n {nrand} -e {extra_name}'.format(nrand=nrand,
                                                    extra_name=chunk_name)
    else:
        opts='--split {split}'.format(split=i)

    text="""
command: |
    source ~/.bashrc
    script_dir=${{ESPY_DIR}}/lensing/bin
    python -u ${{script_dir}}/make-objshear-input.py {options} lcat {sample}
    echo Done
    
job_name: {job_name}
{group}
\n""".format(sample=sample, options=opts, job_name=script_base, group=group)



    with open(wq_fname,'w') as fobj:
        fobj.write(text)
