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

parser.add_option('-n','--nprocess',default=128,
                  help="number or parallel processes")
parser.add_option('--run-admom',action='store_true',
                  help="Force a re-run of adaptive moments")
parser.add_option('--run-psf',action='store_true',
                  help="Force a re-run of psf")

parser.add_option('--pbs1',action='store_true',
                  help="instead of a script do a pbs script for each ccd")
parser.add_option('--walltime1',default='8:00:00',
                  help="walltime for pbs1")

# don't need to source bashrc because pbs will load environment
_template="""#!/bin/bash
module unload espy && module load espy/work
module unload gmix_image && module load gmix_image/work

logfile="%(logfile)s"
echo "host: `hostname`" > "$logfile"
python -u $ESPY_DIR/cluster_step/bin/%(cmd)s %(run)s %(psfnum)s %(shnum)s %(ccd)s 2>&1 > "$logfile"
"""


_pbs1_template="""#!/bin/bash -l
#PBS -N %(run)sp%(psfnum)ss%(shnum)sc%(ccd)02d
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=%(walltime)s
#PBS -q serial
#PBS -o %(base)s.pbsout
#PBS -A des

if [[ "Y${PBS_O_WORKDIR}" != "Y" ]]; then
    cd $PBS_O_WORKDIR
fi

module unload espy && module load espy/work
module unload gmix_image && module load gmix_image/work

logfile="%(base)s.out"
echo "host: `hostname`" > "$logfile"
python -u $ESPY_DIR/cluster_step/bin/%(cmd)s %(run)s %(psfnum)s %(shnum)s %(ccd)s 2>&1 > "$logfile"
"""


_pbs_template_all="""#!/bin/bash -l
#PBS -N %(run)s
#PBS -j oe
#PBS -l nodes=%(nodes)s:ppn=%(ppn)s,walltime=36:00:00
#PBS -q regular
#PBS -o %(base)s.pbsout
#PBS -A des

if [[ "Y${PBS_O_WORKDIR}" != "Y" ]]; then
    cd $PBS_O_WORKDIR
fi

find ./%(ftype)s -name "*%(ftype)s*.sh" | sort | mpirun -np %(np)s minions
"""

# each is 62 ccds plus master: round 63 to 64 so we get full nodes.
# 8 hours because some jobs take just over 4 and the short queue stops at 4
_pbs_template="""#!/bin/bash -l
#PBS -N %(run)sp%(psfnum)ss%(shnum)s
#PBS -j oe
#PBS -l nodes=8:ppn=8,walltime=8:00:00
#PBS -q regular
#PBS -o %(base)s.pbsout
#PBS -A des

if [[ "Y${PBS_O_WORKDIR}" != "Y" ]]; then
    cd $PBS_O_WORKDIR
fi

find ./%(ftype)s -name "*-p%(psfnum)s-s%(shnum)s*%(ftype)s*.sh" | sort | mpirun -np 64 minions
"""


def get_specification(nprocess):
    ppn=8
    if ( (nprocess % ppn) != 0):
        raise ValueError("nprocess=%d is not a multiple of ppn=%d" % (nprocess,ppn))

    nodes=nprocess/ppn
    return nodes,ppn

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
    if options.run_psf:
        cmd = '%s --run-psf' % cmd

    return cmd

def write_pbs_all(run,ftype,nodes,ppn):
    f=files.get_pbs_all_path(run=run,ftype=ftype)
    print 'writing:',f
    base=os.path.basename(f).replace('.pbs','')
    with open(f,'w') as fobj:
        text=_pbs_template_all % {'run':run,
                                  'ftype':ftype,
                                  'nodes':nodes,
                                  'ppn':ppn,
                                  'np':nodes*ppn,
                                  'base':base}
        fobj.write(text)

def write_pbs(run,ftype,psfnum,shnum):
    f=files.get_pbs_path(run=run,ftype=ftype,psfnum=psfnum,shnum=shnum)
    print 'writing:',f
    base=os.path.basename(f).replace('.pbs','')
    with open(f,'w') as fobj:
        text=_pbs_template % {'run':run,
                              'ftype':ftype,
                              'psfnum':psfnum,
                              'shnum':shnum,
                              'base':base}
        fobj.write(text)


def write_script(run,ftype,psfnum,shnum,ccd,cmd):
    sfile=files.get_script_path(run=run,psfnum=psfnum,shnum=shnum,
                                ccd=ccd,ftype=ftype)

    logfile=sfile.replace('.sh','.out')

    print 'writing:',sfile
    with open(sfile,'w') as fobj:
        d={'cmd':cmd,
           'run':run, 
           'psfnum':psfnum,
           'shnum':shnum,
           'ccd':ccd,
           'logfile':logfile}

        text=_template % d
        fobj.write(text)

    os.system('chmod u+x "%s"' % sfile)

def write_pbs1(run,ftype,psfnum,shnum,ccd,cmd,walltime):
    sfile=files.get_script_path(run=run,psfnum=psfnum,shnum=shnum,
                                ccd=ccd,ftype=ftype)
    pbsfile=sfile.replace(".sh",".pbs")

    base=os.path.basename(pbsfile).replace('.pbs','')

    print 'writing:',pbsfile
    with open(pbsfile,'w') as fobj:
        d={'cmd':cmd,
           'run':run, 
           'psfnum':psfnum,
           'shnum':shnum,
           'ccd':ccd,
           'base':base,
           'walltime':walltime}

        text=_pbs1_template % d
        fobj.write(text)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    ftype=args[1]
    nprocess=int(options.nprocess)


    nodes,ppn=get_specification(nprocess)

    conf = files.read_config(run)

    cmd=get_cmd(ftype,options)

    sdir = files.get_script_dir(run=run, ftype=ftype)
    if not os.path.exists(sdir):
        os.makedirs(sdir)

    write_pbs_all(run,ftype,nodes,ppn)

    for psfnum in files.PSFNUMS:
        for shnum in files.SHNUMS:
            write_pbs(run,ftype,psfnum,shnum)
            for ccd in files.CCDS:
                if options.pbs1:
                    write_pbs1(run,ftype,psfnum,shnum,ccd,cmd,
                               options.walltime1)
                else:
                    write_script(run,ftype,psfnum,shnum,ccd,cmd)

main()
