from __future__ import print_function

import os
from . import files

def write_scripts(run):
    sc=Scripts(run)
    sc.write_scripts()

_script_text="""#!/bin/bash

source /astro/u/astrodat/products/eups/bin/setups.sh
setup objshear -r /astro/u/esheldon/exports/objshear-work

/usr/bin/time objshear {config_file}\n"""



class Scripts(dict):
    def __init__(self, run):
        conf = files.cascade_config(run)
        for key in conf:
            self[key] = conf[key]

        if self['run'] != run:
            raise ValueError("mismatch between run ids: %s %s" \
                             % (run,self['run']))


    def write_scripts(self):
        d=files.sample_dir('script',self['run'])
        if not os.path.exists(d):
            os.makedirs(d)

        nsplit=self['src_config']['nsplit']
        for split in xrange(nsplit):
            sfile=files.sample_file('script',self['run'],split=split)
            text=self.script_text(split)

            print("writing script:",sfile)
            with open(sfile,'w') as f:
                f.write(text)

    def script_text(self, split):
        config_file=files.sample_file('config',self['run'],split=split)

        text=_script_text.format(config_file=config_file)
        return text
