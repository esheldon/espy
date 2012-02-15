from __future__ import print_function

import os
from . import files


def write_reduce_script(run):
    pattern = files.sample_file('lensout', run, split=0,fs='hdfs')
    pattern = pattern.replace('-000.dat','-*.dat')
    outf = files.sample_file('lensout', run, fs='hdfs')

    _reduce_script="""
#!/bin/bash

module load sobjshear/sdssmask-work
pattern="{pattern}"
outf={oufile}

tmp_dir=/data/esheldon/redshear-$RANDOM-$RANDOM
mkdir -vp $tmp_dir
tmp_outf=$tmp_dir/$(basename $outf)

echo "reducing files with pattern $pattern"
hadoop fs -cat "$pattern" | redshear > $tmp_outf

err=$?
if [[ $err != "0" ]]; then
    echo "Error running redshear: $err"
fi
if [[ -e $tmp_outf ]]; then
    echo -e "pushing temp file to\\n  $outf"
    hadoop fs -rm $outf 2> /dev/null 1 > /dev/null
    hadoop fs -put $tmp_outf $outf
fi
rm -rvf $tmp_dir

echo `date`
    """.format(pattern=pattern,
               outf=outf)



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
