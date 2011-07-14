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
        
    # inputs to objshear
    def write_config(self):

        lf = lensing.files.sample_file('lcat',self['lens_sample'])
        if self['src_config']['nsplit'] == 0:
            cf = lensing.files.sample_file('config',self['run'])
            sf = lensing.files.sample_file('scat',self['src_sample'])
            of = lensing.files.sample_file('lensout',self['run'])

            self._write_config(cf,lf,sf,of)
        else:
            for i in xrange(self['src_config']['nsplit']):
                cf = lensing.files.sample_file('config',self['run'], split=i)
                sf = lensing.files.sample_file('scat',self['src_sample'], split=i)
                of = lensing.files.sample_file('lensout',self['run'], split=i)
                self._write_config(cf,lf,sf,of)


    def _write_config(self, config_file, lfile, sfile, ofile):
        d = os.path.dirname(config_file)
        if not os.path.exists(d):
            print "Making config dir:",d
            os.makedirs(d)
        d = os.path.dirname(ofile)
        if not os.path.exists(d):
            print "Making lensout output dir:",d
            os.makedirs(d)

        print 'Writing config file:',config_file
        fobj = open(config_file,'w')

        fmt='%-17s %s\n'
        fobj.write(fmt % ("lens_file",lfile))
        fobj.write(fmt % ("source_file",sfile))
        fobj.write(fmt % ("output_file",ofile))

        for key in ['H0','omega_m','npts']:
            fobj.write(fmt % (key, self['cosmo_config'][key]))
        for key in ['nside']:
            fobj.write(fmt % (key, self[key]))

        for key in ['sigmacrit_style']:
            fobj.write(fmt % (key,self['src_config'][key]))
        for key in ['nbin','rmin','rmax']:
            fobj.write(fmt % (key,self['lens_config'][key]))

        fobj.close()

       
