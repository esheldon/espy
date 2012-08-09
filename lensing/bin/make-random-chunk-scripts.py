"""
    %prog sample nsplit

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

parser.add_option("--groups", default=None, help=("Groups for wq"))

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

sample = args[0]
nsplit = int(args[1])

c = lensing.lcat.instantiate_sample(sample)

nperchunk = c['nrand']/nsplit
nleft = c['nrand'] % nsplit

dir = lensing.files.sample_dir(type='lcat',sample=sample)
dir=os.path.join(dir, 'wq')
if not os.path.exists(dir):
    os.makedirs(dir)

groups=options.groups
if groups:
    groups='groups: ['+groups+']'
else:
    groups=''

for i in xrange(nsplit):
    chunk_name='chunk%03d' % i
    script_base='gen-%s' % chunk_name
    wq_fname=os.path.join(dir, script_base+'.yaml')
    print(wq_fname)

    if options.gen_radec:
        opts='-n {nrand} -e {extra_name}'.format(nrand=nrand,
                                                    extra_name=chunk_name)
    else:
        nrand = nperchunk
        if (i == (nsplit-1)):
            nrand += nleft
        opts='--nsplit {nsplit} --split {split}'.format(nsplit=nsplit,
                                                           split=i)

    text="""
command: |
    source ~/.bashrc
    script_dir=${{ESPY_DIR}}/lensing/bin
    python -u ${{script_dir}}/make-objshear-input.py {options} lcat {sample}
    echo Done
    
job_name: {job_name}
{groups}
\n""".format(sample=sample, options=opts, job_name=script_base, groups=groups)



    with open(wq_fname,'w') as fobj:
        fobj.write(text)
