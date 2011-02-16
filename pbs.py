from __future__ import print_function
import esutil
from esutil.misc import isstring
from esutil.ostools import expand_path

import os

class PBS:
    def __init__(self, 
                 filename,
                 command, 
                 job_name=None, 
                 queue='fast', 
                 walltime=None, # e.g. '24:00:00'
                 nodes=1, 
                 ppn=1, 
                 setups=None,
                 no_redirect=False):

        filename = expand_path(filename)
        self.filename = filename

        pbslog=filename+'.pbslog'
        logf = filename+'.log'

        # don't put job name there at all if not sent
        if job_name is not None:
            jobstring = "#PBS -N %s" % job_name
        else:
            jobstring=""
        if walltime is not None:
            timestring = '#PBS -l walltime=' % walltime
        else:
            timestring = ''

        if setups is not None:
            if not isstring(setups):
                setups = '\n'.join(setups) 
        else:
            setups=""


        script="""# vim: set ft=sh :
#PBS -S /bin/bash
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -q {queue}
#PBS -j oe
#PBS -o {pbslog}
#PBS -m a
#PBS -V
#PBS -r n
#PBS -W umask=0022
{timestring}
{jobstring}

{setups}

logf="{logf}"

{command}"""
        if not no_redirect:
            command += ''' &> "$logf"'''

        script += "\n"

        script=script.format(nodes=nodes,
                             ppn=ppn,
                             queue=queue,
                             pbslog=pbslog,
                             timestring=timestring,
                             jobstring=jobstring,
                             setups=setups,
                             logf=logf,
                             command=command)

        self.script=script


    def write(self):
        fobj = open(self.filename, 'w')

        fobj.write(self.script)

        fobj.close()

def create_many(outdir, name, commands, parlist, python=False, **keys):
    """
    Generate pbs scripts for the input command and the parameter list.

    The parlist is a list of dictionaries.

    The commands in the string will have formats {something} and
    these must be keys of the dictionaries.

    files will be

        outdir/name-{jobnum}.pbs

    And the job name will be name-{jobnum}

    e.g.
        parlist=[{'x':3}, {'x':5}]
        pbs_python_many('/some/path', 'job', 'x={x}', parlist)
    """

    outdir=os.path.expandvars(outdir)
    outdir=os.path.expanduser(outdir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    subfile=open( os.path.join(outdir, 'submit-%s.sh' % name),'w')

    njob = len(parlist)
    if njob < 10:
        jobformat = name+'-%i'
    elif njob < 100:
        jobformat = name+'-%02i'
    elif njob < 1000:
        jobformat = name+'-%03i'
    else:
        jobformat = name+'-%04i'

    i=0
    for pars in parlist:
        jobname = jobformat % i
        bfname = jobname+'.pbs'
        fname = os.path.join(outdir, bfname)

        comm=commands.format(**pars)

        print(fname) 

        if python:
            p = PBSPython(fname, comm, 
                          job_name=jobname,
                          **keys)
        else:
            p = PBS(fname, comm, 
                    job_name=jobname,
                    **keys)

        p.write()

        subfile.write('echo -n "%s "; qsub %s\n' % (bfname,bfname))
        i+= 1
    subfile.close()
    print("submit file:",subfile.name)

def test_many():

    setups="""
setup esutil -r ~/exports/esutil-work
setup admom -r ~/exports/admom-work
setup sdsspy -r ~/exports/sdsspy-work
setup espy -r ~/exports/espy-work
    """
    commands="""
x={x}
print 'x=',x
    """

    parlist=[{'x':x} for x in xrange(15)]
    outdir='~/tmp/test-pbsmany'
    name='tpbs'
    queue='batch'

    create_many(outdir, name, commands, parlist, python=True,
                queue=queue, setups=setups) 

class PBSPython():
    def __init__(self, 
                 filename,
                 commands, 
                 job_name=None, 
                 queue='fast', 
                 walltime=None, # e.g. '24:00:00', 
                 nodes=1, 
                 ppn=1, 
                 setups=None,
                 buffer=False):

        filename = expand_path(filename)
        self.filename = filename

        pbslog=filename+'.pbslog'
        logf = filename+'.log'

        # don't put job name there at all if not sent
        if job_name is not None:
            jobstring = "#PBS -N %s" % job_name
        else:
            jobstring=""

        if walltime is not None:
            timestring = '#PBS -l walltime=' % walltime
        else:
            timestring = ''

        if setups is not None:
            if not isstring(setups):
                setups = '\n'.join(setups) 
        else:
            setups=""

        # now the python commands
        if isinstance(commands, basestring):
            commands = [commands]

        commands = [c.strip() for c in commands]

        commands = '\n'.join(commands)

        if buffer:
            ex='python'
        else:
            ex='python -u'



        self.script="""
#PBS -S /bin/bash
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -q {queue}
#PBS -j oe
#PBS -o {pbslog}
#PBS -m a
#PBS -V
#PBS -r n
#PBS -W umask=0022
{timestring}
{jobstring}

{setups}

logf="{logf}"

{ex} &> "$logf" <<EOF
{commands}
EOF

        \n""".format(nodes=nodes,
                     ppn=ppn,
                     queue=queue,
                     pbslog=pbslog,
                     timestring=timestring,
                     jobstring=jobstring,
                     setups=setups,
                     logf=logf,
                     ex=ex,
                     commands=commands)




    def write(self):
        fobj = open(self.filename, 'w')

        fobj.write(self.script)

        fobj.close()
