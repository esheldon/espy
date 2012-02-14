from __future__ import print_function
import os
from . import files

class WQLens(dict):
    def __init__(self, run, groups, priority):
        conf = files.cascade_config(run)
        for key in conf:
            self[key] = conf[key]

        self['groups'] = groups
        self['priority'] = priority

        if self['run'] != run:
            raise ValueError("mismatch between run ids: %s %s" \
                             % (run,self['run']))

    def write_scripts(self):
        d=files.sample_dir('wq',self['run'])
        if not os.path.exists(d):
            os.makedirs(d)

        nsplit=self['src_config']['nsplit']
        for split in xrange(nsplit):
            fname=files.sample_file('wq',self['run'],split=split)

            text=self.wq_text(split)

            print("writing wq submit script:",fname)
            with open(fname,'w') as f:
                f.write(text)

    def wq_text(self, split):
        fs='hdfs'
        config_file=files.sample_file('config',self['run'],fs=fs)

        scat=files.sample_file('scat',self['src_config']['sample'],split=split, fs=fs)
        out_file = files.sample_file('lensout',self['run'], split=split,fs=fs)

        log_file=files.sample_file('log',self['run'],split=split)

        groups = self['groups']
        if groups != '':
            groups='group: [' + groups +']'

        job_name = '%s-%03i' % (self['run'],split)

        s=_wqscript % {'config_file':config_file,
                       'scat':scat,
                       'out_file':out_file,
                       'log_file':log_file,
                       'groups':groups,
                       'priority':self['priority'],
                       'job_name':job_name}
        return s

_wqscript="""
command: |
    source ~esheldon/.bashrc

    module load sobjshear/sdssmask-work

    config_file=%(config_file)s
    scat=%(scat)s

    out_file=%(out_file)s
    log_file=%(log_file)s

    tmp_dir=/data/esheldon/sobjshear
    mkdir -p $tmp_dir
    tmp_file=$tmp_dir/$(basename $out_file)
    tmp_file=${tmp_file/\.dat/-$RANDOM\.dat}

    hadoop fs -cat $scat | sobjshear $config_file 2> $log_file 1> $tmp_file
    err=$?
    if [[ $err != "0" ]]; then
        echo "Error running sobjshear: $err" >> $log_file
    fi
    if [[ -e $tmp_file ]]; then
        echo "pushing temp file to hdfs $out_file" >> $log_file
        hadoop fs -test -e $out_file
        if [[ $? == "0" ]]; then
            echo "    removing existing file $out_file" >> $log_file
            hadoop fs -rm $out_file
        fi

        hadoop fs -put $tmp_file $out_file
        echo "cleaning up temp file: $tmp_file" >> $log_file
        rm -f $tmp_file
    else
        echo "temporary output file was not found: $tmp_file" >> $log_file
    fi

%(groups)s
priority: %(priority)s
job_name: %(job_name)s\n""" 


