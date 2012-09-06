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

        d=files.sample_dir(type='wq',sample=self['run'])
        if not os.path.exists(d):
            os.makedirs(d)

    def write_shear_scripts(self):
        src_nsplit=self['src_config']['nsplit']
        lens_nsplit=self['lens_config']['nsplit']

        for src_split in xrange(src_nsplit):
            for lens_split in xrange(lens_nsplit):
                fname=files.sample_file(type='wq',
                                        sample=self['run'],
                                        lens_split=lens_split,
                                        src_split=src_split)
                text=self.shear_text(lens_split=lens_split, src_split=src_split)
                print("writing wq submit script:",fname)
                with open(fname,'w') as f:
                    f.write(text)

    def shear_text(self, lens_split=None, src_split=None):
        """
        Each of the scripts to run the lens and source splits
        """
        fs='hdfs'
        config_file=files.sample_file(type='config',sample=self['run'],fs=fs)

        scat=files.sample_file(type='scat',
                               sample=self['src_config']['sample'],
                               src_split=src_split,
                               fs=fs)
        lcat=files.sample_file(type='lcat',
                               sample=self['lens_config']['sample'], 
                               lens_split=lens_split,
                               fs=fs)
        out_file = files.sample_file(type='lensout',
                                     sample=self['run'], 
                                     src_split=src_split,
                                     lens_split=lens_split,
                                     fs=fs)


        groups = self['groups']

        job_name = self['run']
        if lens_split is not None:
            job_name = '%s-%03i' % (job_name, lens_split)
        if src_split is not None:
            job_name = '%s-%03i' % (job_name, src_split)

        s=_shear_script % {'config_file':config_file,
                           'scat':scat,
                           'lcat':lcat,
                           'out_file':out_file,
                           'groups':groups,
                           'priority':self['priority'],
                           'job_name':job_name}
        return s


    def write_all_reduce_script(self):
        """
        reduce all
        """
        fname=files.sample_file(type='wq-reduce',sample=self['run'])
        text=self.reduce_all_text()

        print("writing wq reduce all script:",fname)
        with open(fname,'w') as f:
            f.write(text)

    def write_src_reduce_scripts(self):
        """
        reduce across sources at fixed lens split
        """
        lens_nsplit=self['lens_config']['nsplit']

        for i in xrange(lens_nsplit):
            text=self.reduce_src_text(i)
            fname=files.sample_file(type='wq-src-reduce-split',sample=self['run'], 
                                    lens_split=i)
            print("writing wq src reduce script:",fname)
            with open(fname,'w') as f:
                f.write(text)

    def write_lens_collate_scripts(self):
        """
        Collate each of the lens splits
        """
        lens_nsplit=self['lens_config']['nsplit']

        for i in xrange(lens_nsplit):
            text=self.collate_text(i)
            fname=files.sample_file(type='wq-collate-split',sample=self['run'], 
                                    lens_split=i)
            print("writing wq collate script:",fname)
            with open(fname,'w') as f:
                f.write(text)


    def write_lens_concat_script(self):
        """
        Reduce the outputs from the lens reduction across sources.
        """
        fname=files.sample_file(type='wq-lens-concat',sample=self['run'])
        text=self.concat_lenses_text()

        print("writing wq lens concat script:",fname)
        with open(fname,'w') as f:
            f.write(text)



    def reduce_all_text(self):
        """
        A script for reducing every file, even if there are both lens
        and source splits
        """
        src_nsplit=self['src_config']['nsplit']
        lens_nsplit=self['lens_config']['nsplit']
        nfiles=src_nsplit*lens_nsplit

        pattern = files.sample_file(type='lensout', sample=self['run'], 
                                    fs='hdfs')
        pattern = pattern.replace('.dat','-*.dat')
        out_file = files.sample_file(type='reduced', sample=self['run'], 
                                     fs='hdfs')
        config_file=files.sample_file(type='config',sample=self['run'],
                                      fs='hdfs')

        groups = self['groups']

        job_name = 'reduce-%s' % self['run']

        extra='mode: bynode\nN: 1\nmin_mem: 25'
        s=_reduce_script % {'pattern':pattern,
                            'config_file':config_file,
                            'out_file':out_file,
                            'nfiles':nfiles,
                            'groups':groups,
                            'priority':self['priority'],
                            'job_name':job_name,
                            'extra':''}
        return s

    def reduce_src_text(self, lens_split):
        """
        A script for reducing across sources at fixed lens split.
        """
        nfiles=self['src_config']['nsplit']

        pattern = files.sample_file(type='lensout', sample=self['run'], fs='hdfs')
        # just the 3 digits for all of the lens splits followed by the src
        # split
        pattern = pattern.replace('.dat','-%03i-[0-9][0-9][0-9].dat' % lens_split)

        # this part the same as all reduce above
        out_file = files.sample_file(type='src-reduced-split', 
                                     sample=self['run'], 
                                     lens_split=lens_split, 
                                     fs='hdfs')
        config_file=files.sample_file(type='config',sample=self['run'],fs='hdfs')

        groups = self['groups']

        job_name = 'srcred-%s-%03i' % (self['run'],lens_split)

        extra=''
        s=_reduce_script % {'pattern':pattern,
                            'config_file':config_file,
                            'out_file':out_file,
                            'nfiles':nfiles,
                            'groups':groups,
                            'priority':self['priority'],
                            'extra':extra,
                            'job_name':job_name}
        return s



    def concat_lenses_text(self):
        """
        A script for reducing across lenses when we have already produced the
        reductions across sources.  If you just want to reduce all at once, e.g.
        if you don't have a large number of lens splits (or none) use the reduce
        all script.
        """
        nfiles=self['src_config']['nsplit']

        pattern = files.sample_file(type='src-reduced-split', 
                                    lens_split=0, 
                                    sample=self['run'], 
                                    fs='hdfs')
        # just the 3 digit ending rather than both we would have before lens
        # reduction
        pattern = pattern.replace('-000.dat','-[0-9][0-9][0-9].dat')

        # this part the same as all reduce above
        out_file = files.sample_file(type='reduced', sample=self['run'], fs='hdfs')
        config_file=files.sample_file(type='config',sample=self['run'],fs='hdfs')

        groups = self['groups']

        job_name = 'lens-concat-%s' % self['run']

        extra=''
        s=_concat_script % {'pattern':pattern,
                            'out_file':out_file,
                            'nfiles':nfiles,
                            'groups':groups,
                            'priority':self['priority'],
                            'job_name':job_name,
                            'extra':extra}
        return s


    def collate_text(self, lens_split):
        extra=''
        job_name = 'collate-%s-%03i' % (self['run'],lens_split)
        return _collate_script % {'run':self['run'],
                                  'lens_split':lens_split,
                                  'job_name':job_name,
                                  'groups':self['groups'],
                                  'priority':self['priority'],
                                  'extra':extra}

_shear_script="""
command: |
    source ~esheldon/.bashrc

    module load sobjshear/work

    config=%(config_file)s
    scat=%(scat)s
    lcat=%(lcat)s

    outf=%(out_file)s

    echo `hostname`

    tmp_dir=/data/esheldon/sobjshear-$RANDOM-$RANDOM
    mkdir -vp $tmp_dir

    tmp_outf=$tmp_dir/$(basename $outf)
    tmp_lcat=$tmp_dir/$(basename $lcat)
    tmp_config=$tmp_dir/$(basename $config)

    echo -e "staging lcat file to local disk:\\n  $lcat\\n  $tmp_lcat"
    hadoop fs -cat $lcat > $tmp_lcat
    echo -e "staging config file to local disk:\\n  $config\\n  $tmp_config"
    hadoop fs -cat $config > $tmp_config

    hadoop fs -cat $scat | sobjshear $tmp_config $tmp_lcat 1> $tmp_outf

    err=$?
    if [[ $err != "0" ]]; then
        echo "Error running sobjshear: $err"
    fi
    if [[ -e $tmp_outf ]]; then
        echo -e "pushing temp file to\\n  $outf"
        hadoop fs -rm $outf 2> /dev/null 1 > /dev/null
        hadoop fs -put $tmp_outf $outf
    fi
    rm -rvf $tmp_dir 2>&1

    echo `date`

%(groups)s
priority: %(priority)s
job_name: %(job_name)s
""" 

_reduce_script="""
command: |
    source ~esheldon/.bashrc

    module load sobjshear/work

    config="%(config_file)s"
    pattern="%(pattern)s"
    outf="%(out_file)s"

    echo `hostname`

    nfiles=$(hadoop fs -ls "$pattern" | wc -l)
    if [[ $nfiles != "%(nfiles)s" ]]; then
        echo "Error: expected %(nfiles)s files but found $nfiles"
        exit 1
    fi

    tmp_dir=/data/esheldon/redshear-$RANDOM-$RANDOM
    mkdir -vp $tmp_dir

    tmp_outf=$tmp_dir/$(basename $outf)
    tmp_config=$tmp_dir/$(basename $config)

    echo -e "staging config file to local disk:\\n  $config\\n  $tmp_config"
    hadoop fs -cat $config > $tmp_config

    echo "reducing files with pattern $pattern"
    hadoop fs -cat "$pattern" | redshear $tmp_config 1> $tmp_outf

    err=$?
    if [[ $err != "0" ]]; then
        echo "Error running redshear: $err"
    fi
    if [[ -e $tmp_outf ]]; then
        echo -e "pushing temp file to\\n  $outf"
        hadoop fs -rm $outf 2> /dev/null 1 > /dev/null
        hadoop fs -put $tmp_outf $outf
    fi
    rm -rvf $tmp_dir 2>&1

    echo `date`

%(groups)s
priority: %(priority)s
job_name: %(job_name)s
%(extra)s
""" 

_concat_script="""
command: |
    echo `hostname`

    source ~esheldon/.bashrc

    pattern="%(pattern)s"
    outf="%(out_file)s"

    nfiles=$(hadoop fs -ls "$pattern" | wc -l)
    if [[ $nfiles != "%(nfiles)s" ]]; then
        echo "Error: expected %(nfiles)s files but found $nfiles"
        exit 1
    fi

    tmp_dir=/data/esheldon/lens-concat-$RANDOM-$RANDOM
    mkdir -vp $tmp_dir

    tmp_outf=$tmp_dir/$(basename $outf)

    echo "concating files with pattern $pattern"
    hadoop fs -cat "$pattern" > $tmp_outf

    err=$?
    if [[ $err != "0" ]]; then
        echo "Error concatenating files: $err"
    else
        if [[ -e $tmp_outf ]]; then
            echo -e "pushing temp file to\\n  $outf"
            hadoop fs -rm $outf 2> /dev/null 1 > /dev/null
            hadoop fs -put $tmp_outf $outf
        else
            echo "Error: tmp file missing: $tmp_outf"
        fi
    fi
    rm -rvf $tmp_dir 2>&1

    echo `date`

%(groups)s
priority: %(priority)s
job_name: %(job_name)s
%(extra)s
""" 


_collate_script="""
command: |
    source ~esheldon/.bashrc
    python ${ESPY_DIR}/lensing/bin/collate-reduced.py -s %(lens_split)s %(run)s
%(groups)s
priority: %(priority)s
job_name: %(job_name)s
%(extra)s
"""
