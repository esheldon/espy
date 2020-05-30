import lensing
import os
import esutil as eu
from esutil.ostools import path_join
import numpy

def write_configs(run):
    rc = ObjshearRunConfig(run)
    rc.write_config()

MASK_STYLES = {None:1, 'sdss':2, 'healpix':3}

class ObjshearRunConfig(dict):
    """

    This is reads the lcat, scat, run etc. config files and writes out the
    simple config files for input to objshear, as well as the pbs files to run
    the code.

    Note the output is sent to file system fs='local', the script must make
    the copy to nfs or hadoop

    The lensing.files config reading routines read the original hand written
    files, the "config" routines write the objshear inputs.

    """

    def __init__(self, run, fs='nfs'):

        self.fs=fs
        conf = lensing.files.cascade_config(run)
        for key in conf:
            self[key] = conf[key]

        if self['run'] != run:
            raise ValueError("mismatch between run ids: %s %s" \
                             % (run,self['run']))
        
    def make_zlvals(self):
        import sigmacrit
        src_conf=self['src_config']
        if 'zlmin' in src_conf:
            zlmin=src_conf['zlmin']
            zlmax=src_conf['zlmax']
            nzl=src_conf['nzl']
        else:
            sconf = self['scinv_config']
            zlmin=sconf['zlmin']
            zlmax=sconf['zlmax']
            nzl=sconf['nzl']
        return numpy.linspace(zlmin,zlmax,nzl)

    # inputs to objshear
    def write_config(self):
        config_file = lensing.files.sample_file(type='config',
                                                sample=self['run'],
                                                fs=self.fs)

        print 'Writing config file:',config_file
        # should automate this type of thing; maybe an
        # "auto" file class for hdfs that returns the open
        # file handle for the local file?

        mask_type=self['lens_config']['mask_type']
        mask_style=MASK_STYLES[mask_type]

        if self.fs=='hdfs':
            hdfs_file= eu.hdfs.HDFSFile(config_file)
            local_fname=hdfs_file.localfile
        else:
            local_fname=config_file
            eu.ostools.makedirs_fromfile(local_fname)

        with open(local_fname,'w') as local_file:

            fmt='%-17s = %s\n'

            for key in ['H0','omega_m','npts']:
                local_file.write(fmt % (key, self['cosmo_config'][key]))
            for key in ['nside']:
                local_file.write(fmt % (key, self[key]))

            local_file.write(fmt % ("mask_style",mask_style))

            for key in ['sigmacrit_style']:
                local_file.write(fmt % (key,self['src_config'][key]))
            for key in ['nbin','rmin','rmax']:
                local_file.write(fmt % (key,self['lens_config'][key]))

            if self['src_config']['sigmacrit_style'] == 2:
                zlvals=self.make_zlvals()

                local_file.write('zlvals = [')
                zlvals.tofile(local_file, sep=' ')
                local_file.write(']\n')

            if 'zmin' in self['lens_config']:
                st=fmt % ('min_zlens_interp',self['lens_config']['zmin'])
                local_file.write(st)

            for key in ['mag_range','R_range']:
                if key in self:
                    ls=['%s' % v for v in self[key]]
                    vstr = '[' + ' '.join(ls) + ']'
                    local_file.write(fmt % (key,vstr))

        if self.fs=='hdfs':
            hdfs_file.put(clobber=True)
            hdfs_file.cleanup()

      
