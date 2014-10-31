"""
classes for creating wq files
"""

from __future__ import print_function
from .files_common import *

class MakeLcatWQJob(dict):
    def __init__(self, lcat_vers):
        pass

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
        self['len_chunk']=lens_chunk

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

        self['output_dir']=get_output_dir(self['run'])
        self['config_file']=get_xshear_config_file(self['run'])
        self['pattern']=get_output_file(self['run'],self['lens_chunk'],'*')
        self['reduced_file']=get_reduced_file(self['run'],self['lens_chunk'])
        self['job_name']='redshear-%s-%06d' % (self['run'],self['lens_chunk'])

        text=_redshear_template % self
        return text

def get_xshear_wq_dir(run):
    """
    dir holding xshear wq files
    """
    d=get_run_dir(run)
    return os.path.join(d, 'xshear-wq')

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

def get_resdshear_wq_dir(run):
    """
    dir holding reduce wq scripts
    """
    d=get_run_dir(run)
    return os.path.join(d, 'redshear-wq')

def get_redshear_wq_file(run, lens_chunk):
    """
    the yaml wq file for a lens chunk
    """
    d=get_redshear_wq_dir(run)
    fname="%(run)s-lens-%(lens_chunk)06d.yaml"
    fname=fname % {'run':run,
                   'lens_chunk':lens_chunk}

    return os.path.join(d, fname)



