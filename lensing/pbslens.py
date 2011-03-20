import os
from sys import stdout

import lensing
import pbs

def write_pbs(run):
    rp = RunPbs(run)
    rp.write_pbs()


class RunPbs:
    def __init__(self, run):
        #self.conf = lensing.files.json_read(run)
        self.conf = lensing.config.RunConfig(run)

    def write_pbs(self):
        if self.conf['nsplit'] == 0:
            cf = lensing.files.sample_file('config',self.conf['run'])
            pbsf = lensing.files.sample_file('pbs',self.conf['run'])
            self._write_pbs(pbsf,cf)
        else:
            for i in xrange(self.conf['nsplit']):
                cf = lensing.files.sample_file('config',self.conf['run'],
                                               split=i)
                pbsf = lensing.files.sample_file('pbslens',self.conf['run'],
                                                 split=i)
                self._write_pbs(pbsf,cf,i)
    def _write_pbs(self, pbsf, conf_file,split=None):

        d = os.path.dirname(pbsf)
        if not os.path.exists(d):
            print "Making pbs dir:",d
            os.makedirs(d)

        jobname='shear-%s' % self.conf['run']
        if split is not None:
            jobname += '-%03i' % split

        setups=['source /global/data/products/eups/bin/setups.sh',
                'setup objshear -r ~/exports/objshear-work']
        command='time -p OMP_NUM_THREADS=2 objshear %s' % conf_file

        stdout.write('Writing pbs file: %s\n' % pbsf)
        p = pbs.PBS(pbsf,command,setups=setups,job_name=jobname,ppn=2)
        p.write()


