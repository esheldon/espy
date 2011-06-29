from __future__ import print_function

import os
from . import files

def write_submit_scripts(run):
    sc=LensCondor(run)
    sc.write_scripts()

_script_text="""
Universe        = vanilla
Notification    = Always
Initialdir      = {proc_dir}
Executable      = {script_base}.sh
Output          = {script_base}.out
Error           = {script_base}.err
Log             = {script_base}.log
GetEnv          = False
Notify_user     = esheldon@bnl.gov
+Experiment     = "astro"

Queue\n"""



class LensCondor(dict):
    def __init__(self, run):
        conf = files.cascade_config(run)
        for key in conf:
            self[key] = conf[key]

        if self['run'] != run:
            raise ValueError("mismatch between run ids: %s %s" \
                             % (run,self['run']))


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

    def script_text(self, split):
        d = files.sample_dir('condor',self['run'])
        script_file=files.sample_file('script',self['run'],split=split)
        script_base = os.path.basename(script_file).replace('.sh','')

        text=_script_text.format(proc_dir=d, script_base=script_base)
        return text
