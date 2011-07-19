from __future__ import print_function

import os
from . import files

def write_submit_scripts(run):
    sc=LensCondor(run)
    sc.write_scripts()
    sc.write_allone()

_script_head="""
Universe        = vanilla
Notification    = Error
GetEnv          = True
Notify_user     = esheldon@bnl.gov
Requirements    = (CPU_Experiment == "astro")
+Experiment     = "astro"
Initialdir      = {proc_dir}
"""
_script_text="""
Executable      = {script_base}.sh
Output          = {script_base}.out
Error           = {script_base}.err
Log             = {script_base}.log
"""



class LensCondor(dict):
    def __init__(self, run):
        conf = files.cascade_config(run)
        for key in conf:
            self[key] = conf[key]

        if self['run'] != run:
            raise ValueError("mismatch between run ids: %s %s" \
                             % (run,self['run']))


    def write_allone(self):
        d=files.sample_dir('condor',self['run'])
        if not os.path.exists(d):
            os.makedirs(d)

        d = files.sample_dir('condor',self['run'])
        script_file=files.sample_file('condor',self['run'])

        print("writing all in one script:",script_file)
        with open(script_file,'w') as fobj:

            text=_script_head.format(proc_dir=d)
            fobj.write(text)

            nsplit=self['src_config']['nsplit']
            for split in xrange(nsplit):
                s=files.sample_file('script',self['run'],split=split)

                sbase = os.path.basename(s).replace('.sh','')

                text=_script_text.format(script_base=sbase) + '\nQueue\n'
                fobj.write(text)

    def write_scripts(self):
        d=files.sample_dir('condor',self['run'])
        if not os.path.exists(d):
            os.makedirs(d)

        nsplit=self['src_config']['nsplit']
        for split in xrange(nsplit):
            cfile=files.sample_file('condor',self['run'],split=split)
            text=self.script_text(split)

            print("writing condor submit script:",cfile)
            with open(cfile,'w') as f:
                f.write(text)
                f.write('\n')

    def script_text(self, split):
        d = files.sample_dir('condor',self['run'])
        script_file=files.sample_file('script',self['run'],split=split)
        script_base = os.path.basename(script_file).replace('.sh','')

        text=[_script_head.format(proc_dir=d)] 
        text += [_script_text.format(script_base=script_base)]
        text += ['Queue']

        text = '\n'.join(text)
        return text
