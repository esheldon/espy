# this is obsolete
from __future__ import print_function

import os
from . import files

class Scripts(dict):
    def __init__(self, run, manager='modules'):
        conf = files.cascade_config(run)
        self.manager=manager
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
            os.popen('chmod 755 '+sfile)

    def script_text(self, split):
        config_file=files.sample_file('config',self['run'],split=split)

        local_file=files.sample_file('lensout',self['run'],split=split,fs='local')
        local_dir=files.sample_dir('lensout',self['run'],fs='local')

        hadoop_file=files.sample_file('lensout',self['run'],split=split,fs='hdfs')
        hadoop_dir=files.sample_dir('lensout',self['run'],fs='hdfs')

        masktype=self['lens_config'].get('masktype',None)
        if masktype == 'sdss':
            objshear_version='sdssmask-work'
        else:
            objshear_version='work'

        if self.manager == 'eups':
            objshear_dir = 'objshear-'+objshear_version
            text=_script_text_eups.format(config_file=config_file, 
                                          objshear_dir=objshear_dir,
                                          local_file=local_file,
                                          local_dir=local_dir,
                                          hadoop_file=hadoop_file,
                                          hadoop_dir=hadoop_dir)
        elif self.manager == 'modules':
            text=_script_text_modules.format(config_file=config_file, 
                                             objshear_version=objshear_version,
                                             local_file=local_file,
                                             local_dir=local_dir,
                                             hadoop_file=hadoop_file,
                                             hadoop_dir=hadoop_dir)

        return text
