"""

runs are a combination of source and lens sample with configuration details

"""

from __future__ import print_function
from .files_common import *
from .import scat
from . import lcat
from . import output

from .wqscripts import XShearWQJob, RedshearWQJob
from .xshear_config import XShearConfig

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
        return
        self.write_collate_wq()
    
    def write_xshear_wq(self):
        """
        write all the chunks and tiles
        """
        lens_nchunk=self['lens_conf']['nchunk']
        tilenames=scat.get_tilenames(self['source_conf']['scat_name'])

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

    def write_redshear_wq(self):
        """
        write all the chunks
        """
        lens_nchunk=self['lens_conf']['nchunk']

        for lens_chunk in xrange(lens_nchunk):
            job=RedShearWQJob(self['run'], lens_chunk)
            job.write()

    def _scat_exists(self, tilename):
        """
        check if the source catalog exists
        """
        fname=scat.get_scat_file(self['scat_vers'],
                                 tilename)
        return os.path.exists(fname)

                                    

