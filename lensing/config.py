import lensing
import os
import esutil as eu
from esutil.ostools import path_join

def write_config(run):
    rc = RunConfig(run)
    rc.write_config()

class RunConfig(dict):
    """

    This is reads the json config files and writes out the simple config files
    for input to objshear, as well as the pbs files to run the code.

    The "json" routines read the original hand written json files, 
    the "config" routines write the objshear inputs.
    """

    def __init__(self, run):
        conf = lensing.files.json_read('run',run)
        for key in conf:
            self[key] = conf[key]

        if self['run'] != run:
            raise ValueError("mismatch between run ids: %s %s" \
                             % (run,self['run']))

        # make sure the cosmology is consistent
        l = lensing.files.json_read('lcat',self['lens_sample'])
        s = lensing.files.json_read('scat',self['src_sample'])

        csample = self['cosmo_sample']
        if (csample != l['cosmo_sample'] or csample != s['cosmo_sample']):
            err= """
            cosmo sample ismatch:
                run config:  %s
                lcat config: %s
                scat config: %s
            """ % (self['cosmo_sample'],l['cosmo_sample'],s['cosmo_sample'])
            raise ValueError(err)

        self['lens_config'] = l
        self['src_config'] = s
        cosmo = lensing.files.json_read('cosmo',csample)
        for key in cosmo:
            self[key] = cosmo[key]
        self['nsplit'] = l['nsplit']
        self['sigmacrit_style'] = s['sigmacrit_style']
        
    # inputs to objshear
    def write_config(self):

        sf = lensing.files.sample_file('scat',self['src_sample'])
        if self['nsplit'] == 0:
            cf = lensing.files.sample_file('config',self['run'])

            lf = lensing.files.sample_file('lcat',self['lens_sample'])
            of = lensing.files.sample_file('lensout',self['run'])

            self._write_config(cf,lf,sf,of)
        else:
            for i in xrange(self['nsplit']):
                cf = lensing.files.sample_file('config',self['run'],
                                               split=i)
                lf = lensing.files.sample_file('lcat',self['lens_sample'],
                                               split=i)
                of = lensing.files.sample_file('lensout',self['run'],
                                               split=i)
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

        fobj.write('lens_file\n')
        fobj.write(lfile+'\n')
        fobj.write('source_file\n')
        fobj.write(sfile+'\n')
        fobj.write('output_file\n')
        fobj.write(ofile+'\n')

        for key in ['H0','omega_m','npts','nside',
                    'sigmacrit_style','nbin','rmin','rmax']:
            fobj.write(key+'\n')
            fobj.write('%s\n' % self[key])

        fobj.close()

       
