"""
classes for creating wq files
"""

from __future__ import print_function
from .files_common import *

class MakeLcatWQJob(dict):
    """
    class to write wq scripts to make lcats
    """
    def __init__(self, lcat_vers, chunk):
        self['lcat_vers']=lcat_vers
        self['chunk']=chunk

    def write(self):
        """
        write the wq script
        """
        fname=get_make_lcat_wq_file(self['lcat_vers'], self['chunk'])
        d=get_make_lcat_wq_dir(self['lcat_vers'])

        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        text=self.get_text()

        print("writing:",fname)
        with open(fname,'w') as fobj:
            fobj.write(text)

    def get_text(self):
        """
        text for the wq file
        """
        job_name="%s-%06d" % (self['lcat_vers'], self['chunk'])

        text=_make_lcat_wq_template.format(lcat_vers=self['lcat_vers'],
                                           chunk=self['chunk'],
                                           job_name=job_name)

        return text

class XShearWQJob(dict):
    """
    For writing xshear wq files
    """
    def __init__(self, run, lens_chunk, source_tilename):
        conf=cascade_config(run)
        self.update(conf)

        self['lens_chunk']=lens_chunk
        self['source_tilename']=source_tilename
    
    def write(self):
        """
        write the wq file
        """

        d=get_xshear_wq_dir(self['run'])
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        fname=get_xshear_wq_file(self['run'],
                                 self['lens_chunk'],
                                 self['source_tilename'])
        text=self.get_text()
        print("writing:",fname)
        with open(fname,'w') as fobj:
            fobj.write(text)

    def get_text(self):
        """
        get the wq job text
        """

        c={}
        c['config_file']=get_xshear_config_file(self['run'])

        c['scat_file']=get_scat_file(self['scat_vers'],
                                     self['source_tilename'])
        c['lcat_file']=get_lcat_file(self['lcat_vers'],
                                     self['lens_chunk'])
        c['output_file']=get_output_file(self['run'], self['lens_chunk'],
                                         self['source_tilename'])
        job_name='%(run)s-%(lens_chunk)06d-%(source_tilename)s'
        c['job_name']=job_name % self

        text=_xshear_wq_template % c

        return text

class RedshearWQJob(dict):
    """
    make wq for reduced files
    """
    def __init__(self, run, lens_chunk):
        conf=cascade_config(run)
        self.update(conf)
        self['lens_chunk']=lens_chunk

    def write(self):
        """
        make the wq file to reduce across sources
        """

        wqfile=get_redshear_wq_file(self['run'],self['lens_chunk'])
        d=get_redshear_wq_dir(self['run'])

        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        text=self.get_text()

        print("writing:",wqfile)
        with open(wqfile,'w') as fobj:
            fobj.write(text)

    def get_text(self):

        self['config_file']=get_xshear_config_file(self['run'])
        self['pattern']=get_output_file(self['run'],self['lens_chunk'],'*')
        self['reduced_dir']=get_reduced_dir(self['run'])
        self['reduced_file']=get_reduced_file(self['run'],self['lens_chunk'])
        self['job_name']='redshear-%s-%06d' % (self['run'],self['lens_chunk'])

        text=_redshear_template % self
        return text

class CombineWQJob(dict):
    """
    make wq for reduced files
    """
    def __init__(self, run):
        conf=cascade_config(run)
        self.update(conf)

    def write(self):
        """
        make the wq file to reduce across sources
        """

        wqfile=get_combine_wq_file(self['run'])
        d=get_combine_wq_dir(self['run'])

        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        text=self.get_text()

        print("writing:",wqfile)
        with open(wqfile,'w') as fobj:
            fobj.write(text)

    def get_text(self):

        self['pattern']=get_reduced_file(self['run'],'*')
        self['combined_dir']=get_combined_dir(self['run'])
        self['combined_file']=get_combined_file(self['run'])
        self['job_name']='combine-%s' % self['run']

        text=_combine_template % self
        return text

class CollateWQJob(dict):
    """
    make wq for collating
    """
    def __init__(self, run):
        conf=cascade_config(run)
        self.update(conf)

    def write(self):
        """
        make the wq file to reduce across sources
        """

        wqfile=get_collate_wq_file(self['run'])
        d=get_collate_wq_dir(self['run'])

        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        text=self.get_text()

        print("writing:",wqfile)
        with open(wqfile,'w') as fobj:
            fobj.write(text)

    def get_text(self):
        self['job_name']='collate-%s' % self['run']
        text=_collate_template % self
        return text



def get_make_lcat_wq_dir(lcat_vers):
    """
    dir holding wq files to make lens catalog
    """
    dir=get_lcat_dir(lcat_vers)
    return os.path.join(dir, 'wq')

def get_make_lcat_wq_file(lcat_vers, chunk):
    dir=get_make_lcat_wq_dir(lcat_vers)
    name='%s-%06d.yaml' % (lcat_vers, chunk)
    return os.path.join(dir, name)


def get_xshear_wq_dir(run):
    """
    dir holding xshear wq files
    """
    d=get_run_dir(run)
    return os.path.join(d, 'wq-xshear')

def get_xshear_wq_file(run, lens_chunk, source_tilename):
    """
    the yaml wq file for a source tilename and lens chunk
    """
    d=get_xshear_wq_dir(run)
    fname="%(run)s-lens-%(lens_chunk)06d-src-%(source_tilename)s.yaml"
    fname=fname % {'run':run,
                   'lens_chunk':lens_chunk,
                   'source_tilename':source_tilename}

    return os.path.join(d, fname)

def get_redshear_wq_dir(run):
    """
    dir holding reduce wq scripts
    """
    d=get_run_dir(run)
    return os.path.join(d, 'wq-redshear')

def get_redshear_wq_file(run, lens_chunk):
    """
    the yaml wq file for a lens chunk
    """
    d=get_redshear_wq_dir(run)
    fname="%(run)s-lens-%(lens_chunk)06d.yaml"
    fname=fname % {'run':run,
                   'lens_chunk':lens_chunk}

    return os.path.join(d, fname)

def get_combine_wq_dir(run):
    """
    dir holding reduce wq scripts
    """
    d=get_run_dir(run)
    return os.path.join(d, 'wq-combine')

def get_combine_wq_file(run):
    """
    the yaml wq to combine all lens chunks
    """
    d=get_combine_wq_dir(run)
    fname="%(run)s-combine.yaml"
    fname=fname % {'run':run}

    return os.path.join(d, fname)

def get_collate_wq_dir(run):
    """
    dir holding reduce wq scripts
    """
    d=get_run_dir(run)
    return os.path.join(d, 'wq-collate')

def get_collate_wq_file(run):
    """
    the yaml wq to collate all lens chunks
    """
    d=get_collate_wq_dir(run)
    fname="%(run)s-collate.yaml"
    fname=fname % {'run':run}

    return os.path.join(d, fname)

_make_lcat_wq_template="""
command: |
    lcat_vers={lcat_vers}
    chunk={chunk}
    $ESPY_DIR/des/bin/make-xshear-lcat --chunk $chunk $lcat_vers

job_name: {job_name}
"""


_xshear_wq_template="""
command: |
    source ~esheldon/.bashrc

    hostname
    module load xshear/work

    config=%(config_file)s
    scat=%(scat_file)s
    lcat=%(lcat_file)s

    outf=%(output_file)s

    tmp_dir=$TMPDIR/sobjshear-$RANDOM-$RANDOM
    mkdir -vp $tmp_dir

    dname=$(dirname $outf)
    mkdir -vp $dname

    tmp_outf=$tmp_dir/$(basename $outf)
    tmp_lcat=$tmp_dir/$(basename $lcat)
    tmp_scat=$tmp_dir/$(basename $scat)
    tmp_config=$tmp_dir/$(basename $config)

    echo -e "staging lcat file to local disk:\\n  $lcat\\n  $tmp_lcat"
    cp $lcat $tmp_lcat
    echo -e "staging scat file to local disk:\\n  $scat\\n  $tmp_scat"
    cp $scat $tmp_scat

    echo -e "staging config file to local disk:\\n  $config\\n  $tmp_config"
    cp $config $tmp_config

    xshear $tmp_config $tmp_lcat < $tmp_scat 1> $tmp_outf

    err=$?
    if [[ $err != "0" ]]; then
        echo "Error running xshear: $err"
    fi
    if [[ -e $tmp_outf ]]; then
        echo -e "pushing temp file to\\n  $outf"
        mv -fv $tmp_outf $outf
    fi
    rm -rvf $tmp_dir 2>&1

    date

job_name: "%(job_name)s"
"""

_redshear_template="""
command: |
    source ~/.bashrc
    module load xshear/work

    dir=%(reduced_dir)s
    conf=%(config_file)s
    outf=%(reduced_file)s

    mkdir -p $dir
    cat %(pattern)s | redshear $conf > $outf

job_name: "%(job_name)s"
"""

_combine_template="""
command: |
    dir=%(combined_dir)s
    outf=%(combined_file)s
    mkdir -p $dir

    echo "concatenating to $outf"
    cat %(pattern)s > $outf
    echo
    echo done
    echo

job_name: "%(job_name)s"
"""

_collate_template="""
command: |
    source ~/.bashrc
    $ESPY_DIR/des/bin/collate %(run)s

    echo
    echo done
    echo

job_name: "%(job_name)s"
"""


