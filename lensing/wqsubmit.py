from __future__ import print_function
import os
from . import files

class WQLens(dict):
    def __init__(self, run, groups, priority):
        conf = files.cascade_config(run)
        for key in conf:
            self[key] = conf[key]

        if groups != '':
            groups='group: [' + groups +']'
        self['groups'] = groups
        self['priority'] = priority

        if self['run'] != run:
            raise ValueError("mismatch between run ids: %s %s" \
                             % (run,self['run']))

        d=files.sample_dir('wq',self['run'])
        if not os.path.exists(d):
            os.makedirs(d)

    def write_shear_scripts(self):
        nsplit=self['src_config']['nsplit']
        for split in xrange(nsplit):
            fname=files.sample_file('wq',self['run'],split=split)

            text=self.shear_text(split)

            print("writing wq submit script:",fname)
            with open(fname,'w') as f:
                f.write(text)


    def write_reduce_script(self):
        fname=files.sample_file('wq-reduce',self['run'])

        text=self.reduce_text()

        print("writing wq reduce submit script:",fname)
        with open(fname,'w') as f:
            f.write(text)

    def shear_text(self, split):
        fs='hdfs'
        config_file=files.sample_file('config',self['run'],fs=fs)

        scat=files.sample_file('scat',self['src_config']['sample'],split=split, fs=fs)
        lcat=files.sample_file('lcat',self['lens_config']['sample'], fs=fs)
        out_file = files.sample_file('lensout',self['run'], split=split,fs=fs)

        log_file=files.sample_file('log',self['run'],split=split)

        groups = self['groups']

        job_name = '%s-%03i' % (self['run'],split)

        s=_shear_script % {'config_file':config_file,
                           'scat':scat,
                           'lcat':lcat,
                           'out_file':out_file,
                           'log_file':log_file,
                           'groups':groups,
                           'priority':self['priority'],
                           'job_name':job_name}
        return s

    def reduce_text(self):

        pattern = files.sample_file('lensout', self['run'], split=0,fs='hdfs')
        pattern = pattern.replace('-000.dat','-*.dat')
        out_file = files.sample_file('reduced', self['run'], fs='hdfs')
        log_file=files.sample_file('log-reduce',self['run'])
        config_file=files.sample_file('config',self['run'],fs='hdfs')

        groups = self['groups']

        job_name = 'reduce-%s' % self['run']

        s=_reduce_script % {'pattern':pattern,
                            'config_file':config_file,
                            'out_file':out_file,
                            'log_file':log_file,
                            'groups':groups,
                            'priority':self['priority'],
                            'job_name':job_name}
        return s


_shear_script="""
command: |
    source ~esheldon/.bashrc

    module load sobjshear/sdssmask-work

    config=%(config_file)s
    scat=%(scat)s
    lcat=%(lcat)s

    outf=%(out_file)s
    logf=%(log_file)s

    echo `hostname` > $logf

    tmp_dir=/data/esheldon/sobjshear-$RANDOM-$RANDOM
    mkdir -vp $tmp_dir 2>> $logf

    tmp_outf=$tmp_dir/$(basename $outf)
    tmp_lcat=$tmp_dir/$(basename $lcat)
    tmp_config=$tmp_dir/$(basename $config)

    echo -e "staging lcat file to local disk:\\n  $lcat\\n  $tmp_lcat" >> $logf
    hadoop fs -cat $lcat > $tmp_lcat
    echo -e "staging config file to local disk:\\n  $config\\n  $tmp_config" >> $logf
    hadoop fs -cat $config > $tmp_config

    hadoop fs -cat $scat | sobjshear $tmp_config $tmp_lcat 2>> $logf 1> $tmp_outf

    err=$?
    if [[ $err != "0" ]]; then
        echo "Error running sobjshear: $err" >> $logf
    fi
    if [[ -e $tmp_outf ]]; then
        echo -e "pushing temp file to\\n  $outf" >> $logf
        hadoop fs -rm $outf 2> /dev/null 1 > /dev/null
        hadoop fs -put $tmp_outf $outf
    fi
    rm -rvf $tmp_dir 2>&1 >> $logf

    echo `date` >> $logf

%(groups)s
priority: %(priority)s
job_name: %(job_name)s
""" 

_reduce_script="""
command: |
    source ~esheldon/.bashrc

    module load sobjshear/sdssmask-work

    config=%(config_file)s
    pattern="%(pattern)s"
    outf=%(out_file)s
    logf=%(log_file)s

    echo `hostname` > $logf

    tmp_dir=/data/esheldon/redshear-$RANDOM-$RANDOM
    mkdir -vp $tmp_dir 2>> $logf

    tmp_outf=$tmp_dir/$(basename $outf)
    tmp_config=$tmp_dir/$(basename $config)

    echo -e "staging config file to local disk:\\n  $config\\n  $tmp_config" >> $logf
    hadoop fs -cat $config > $tmp_config

    echo "reducing files with pattern $pattern" >> $logf
    hadoop fs -cat "$pattern" | redshear $tmp_config 2>> $logf 1> $tmp_outf

    err=$?
    if [[ $err != "0" ]]; then
        echo "Error running redshear: $err" >> $logf
    fi
    if [[ -e $tmp_outf ]]; then
        echo -e "pushing temp file to\\n  $outf" >> $logf
        hadoop fs -rm $outf 2> /dev/null 1 > /dev/null
        hadoop fs -put $tmp_outf $outf
    fi
    rm -rvf $tmp_dir 2>&1 >> $logf

    echo `date` >> $logf

%(groups)s
priority: %(priority)s
job_name: %(job_name)s
""" 


