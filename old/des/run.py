"""

runs are a combination of source and lens sample with configuration details

"""

from __future__ import print_function
from .files import *
from .import scat
from . import lcat

from .wqscripts import XShearWQJob, RedshearWQJob, CombineWQJob, CollateWQJob
from .xshear_config import XShearConfig

def write_run(run):
    """
    write the xshear config, and wq submit scripts for xshear, redshear,
    combination, and collation
    """
    r=Run(run)
    r.write_all()
class Run(dict):
    """
    for writing run configs and wq submit files
    """
    def __init__(self, run):
        conf=cascade_config(run)
        self.update(conf)

    def write_all(self):
        """
        write config and yaml files
        """
        self.write_config()
        self.write_wq()

    def write_config(self):
        """
        write the xshear config
        """
        xshear_conf=XShearConfig(self['run'])
        xshear_conf.write()

    def write_wq(self):
        """
        write the cfg file
        """
        self.write_xshear_wq()
        self.write_redshear_wq()
        self.write_combine_wq()
        self.write_collate_wq()

    def write_xshear_wq(self):
        """
        write all the chunks and tiles
        """
        lens_nchunk=self['lens_conf']['nchunk']
        tilenames=scat.get_tilenames(self['source_conf']['scat_table'])

        ntile=len(tilenames)
        for lens_chunk in xrange(lens_nchunk):
            for i,tilename in enumerate(tilenames):
                print("    %d/%d: %s" % (i,ntile,tilename))
                # first check if this source catalog exists
                if self._scat_exists(tilename):
                    job=XShearWQJob(self['run'],
                                    lens_chunk,
                                    tilename)
                    job.write()
                else:
                    print("    skipping due to missing scat")

    def check_xshear_output(self):
        """
        write all the chunks and tiles
        """
        lens_nchunk=self['lens_conf']['nchunk']
        tilenames=scat.get_tilenames(self['source_conf']['scat_table'])

        ntile=len(tilenames)
        for lens_chunk in xrange(lens_nchunk):
            print("    checking chunk: %d/%d" % (lens_chunk+1, lens_nchunk))
            for i,tilename in enumerate(tilenames):
                # first check if this source catalog exists
                if self._scat_exists(tilename):
                    job=XShearWQJob(self['run'],
                                    lens_chunk,
                                    tilename)
                    info=job.get_info()
                    if not os.path.exists(info['output_file']):
                        print("missing output:",info['output_file'])


    def write_redshear_wq(self):
        """
        write all the chunks
        """
        lens_nchunk=self['lens_conf']['nchunk']

        for lens_chunk in xrange(lens_nchunk):
            job=RedshearWQJob(self['run'], lens_chunk)
            job.write()

    def write_combine_wq(self):
        """
        write all the chunks
        """
        job=CombineWQJob(self['run'])
        job.write()

    def write_collate_wq(self):
        """
        write all the chunks
        """
        job=CollateWQJob(self['run'])
        job.write()


    def _scat_exists(self, tilename):
        """
        check if the source catalog exists
        """
        fname=scat.get_scat_file(self['scat_vers'],
                                 tilename)
        return os.path.exists(fname)
