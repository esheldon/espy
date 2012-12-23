"""
    %prog run ftype

ftype admom, psf, shear
"""

import sys
import os
import cluster_step
from cluster_step import files

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-g','--groups',default=None,
                  help='groups for wq, csv. Default is unconstrained')
parser.add_option('-p','--priority',default='med',
                  help='priority for queue')


parser.add_option('--run-admom',action='store_true',
                  help="Force a re-run of adaptive moments")
parser.add_option('--run-psf',action='store_true',
                  help="Force a re-run of psf")

_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/work
    module unload gmix_image && module load gmix_image/work
    python $ESPY_DIR/cluster_step/bin/%(cmd)s %(run)s %(psfnum)s %(shnum)s %(ccd)s

%(groups)s
job_name: %(job_name)s
priority: %(pri)s\n"""

def get_cmd(ftype, options):
    if ftype=='admom':
        cmd='run-admom.py'
    elif ftype=='psf':
        cmd='run-psf.py'
    elif ftype=='shear':
        cmd='run-shear.py'
    else:
        raise ValueError("bad ftype: '%s'" % ftype)

    if options.run_admom:
        cmd = '%s --run-admom' % cmd
    elif options.run_psf:
        cmd = '%s --run-psf' % cmd

    return cmd

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    ftype=args[1]

    conf = files.read_config(run)

    cmd=get_cmd(ftype,options)

    groups=''
    if options.groups is not None:
        groups = 'group: [%s]' % options.groups


    wqd = files.get_wq_dir(run=run, ftype=ftype)
    if not os.path.exists(wqd):
        os.makedirs(wqd)

    for psfnum in files.psfnums:
        for shnum in files.shnums:
            for ccd in files.ccds:
                wqfile=files.get_wq_path(run=run,psfnum=psfnum,shnum=shnum,
                                         ccd=ccd,ftype=ftype)
                job_name='%s-%s-%s-%s-%s' % (run,psfnum,shnum,ccd,ftype)

                print 'writing:',wqfile
                with open(wqfile,'w') as fobj:
                    d={'job_name':job_name,
                       'cmd':cmd,
                       'run':run, 
                       'psfnum':psfnum,
                       'shnum':shnum,
                       'ccd':ccd,
                       'groups':groups,
                       'pri':options.priority}

                    text=_wqtemplate % d
                    fobj.write(text)


main()
