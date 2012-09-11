import lensing
import os
import esutil as eu
from esutil.ostools import path_join

def write_configs(run):
    rc = ObjshearRunConfig(run)
    rc.write_config()


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

    def __init__(self, run):

        conf = lensing.files.cascade_config(run)
        for key in conf:
            self[key] = conf[key]

        if self['run'] != run:
            raise ValueError("mismatch between run ids: %s %s" \
                             % (run,self['run']))

        # make sure the cosmology is consistent
        l = self['lens_config']
        s = self['src_config']
        if (l['cosmo_sample'] != s['cosmo_sample']):
            err= """
            cosmo sample mismatch:
                lcat config: %s
                scat config: %s
            """ % (l['cosmo_sample'],s['cosmo_sample'])
            raise ValueError(err)
        
    def make_zlvals(self):
        import sigmacrit
        sconf = self['src_config']
        return sigmacrit.make_zlvals(sconf['dzl'], sconf['zlmin'], sconf['zlmax'])

    # inputs to objshear
    def write_config(self):
        fs='hdfs'
        #lens_file = lensing.files.sample_file(type='lcat',sample=self['lens_sample'],fs=fs)
        config_file = lensing.files.sample_file(type='config',sample=self['run'], fs=fs)

        print 'Writing config file:',config_file
        # should automate this type of thing; maybe an
        # "auto" file class for hdfs that returns the open
        # file handle for the local file?
        with eu.hdfs.HDFSFile(config_file) as hdfs_file:
            with open(hdfs_file.localfile,'w') as local_file:

                fmt='%-17s = %s\n'

                for key in ['H0','omega_m','npts']:
                    local_file.write(fmt % (key, self['cosmo_config'][key]))
                for key in ['nside']:
                    local_file.write(fmt % (key, self[key]))

                for key in ['mask_style']:
                    local_file.write(fmt % (key,self[key]))

                for key in ['sigmacrit_style']:
                    local_file.write(fmt % (key,self['src_config'][key]))
                for key in ['nbin','rmin','rmax']:
                    local_file.write(fmt % (key,self['lens_config'][key]))

                if self['src_config']['sigmacrit_style'] == 2:
                    zlvals=self.make_zlvals()
                    #local_file.write(fmt % ('nzl',zlvals.size))

                    local_file.write('zlvals = [')
                    zlvals.tofile(local_file, sep=' ')
                    local_file.write(']\n')

                if 'zmin' in self['lens_config']:
                    local_file.write(fmt % ('min_zlens_interp',self['lens_config']['zmin']))

                for key in ['mag_range','R_range']:
                    if key in self:
                        vstr = '[' + ' '.join( ['%s' % v for v in self[key]] ) + ']'
                        local_file.write(fmt % (key,vstr))

            hdfs_file.put(clobber=True)

      
