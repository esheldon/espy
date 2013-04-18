"""
    %prog sample

Description:

    Make wq scripts for generation of random lcat in chunks.

    Note if gen-radec is not sent, the randoms must have been
    pre-generated.
"""
from __future__ import print_function
import sys, os
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option("--gen-radec", action='store_true', 
                  help=("For this sample we generate the ra,dec points"))

parser.add_option("-g", "--group", default=None, help=("Groups for wq"))

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    sample = args[0]

    c = lensing.lcat.instantiate_sample(sample=sample)

    nsplit = c['nsplit']
    if 'nrand' not in c:
        nrand=c.get_nrand()
    else:
        nrand=c['nrand']

    nperchunk = nrand/nsplit
    nleft = nrand % nsplit

    f0=lensing.files.sample_file(type='wq-genrand',
                                 sample=sample,
                                 lens_split=0)
    dir=os.path.dirname(f0)
    if not os.path.exists(dir):
        os.makedirs(dir)

    group=options.group
    if group:
        group='group: ['+group+']'
    else:
        group=''

    for i in xrange(nsplit):

        wq_fname=lensing.files.sample_file(type='wq-genrand',
                                           sample=sample,
                                           lens_split=i)


        chunk_name='chunk%03d' % i
        script_base='gen-%s' % chunk_name
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
    module unload espy && module load espy/work
    dir=${{ESPY_DIR}}/lensing/bin
    python -u ${{dir}}/make-objshear-input.py {options} -t lcat -s {sample}
    echo Done
        
job_name: {job_name}
{group}\n"""
        text=text.format(sample=sample,
                         options=opts,
                         job_name=script_base,
                         group=group)



        with open(wq_fname,'w') as fobj:
            fobj.write(text)

main()
