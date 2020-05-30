import os
from sys import stdout,stderr
import subprocess
import numpy
from numpy import where
import esutil
import esutil as eu
from esutil import recfile

from esutil.ostools import path_join, expand_path
from esutil import numpy_util

import sdssgal


# this is imported into __init__
class CasZphot:
    """

    Work with the ascii files from the CAS.  These are the p(z) files created
    by carlos.

    """
    def __init__(self, origin='uchicago', type='dr7pofz'):
        self.good_origins = ['uchicago']
        if origin not in self.good_origins:
            raise ValueError("origin should be in: %s" % self.good_origins)

        self.good_types = ['dr7pofz']
        if type not in self.good_types:
            raise ValueError("type should be in: %s" % self.good_types)

        self.origin = origin
        self.type=type

        self.basedir = path_join('~esheldon','photoz',origin,type)
        self.basedir = os.path.expanduser(self.basedir)


    def columns_dir(self):
        return self.basedir + '.cols'
    def open_columns(self):
        import columns
        coldir=self.columns_dir()
        stdout.write("Opening zphot columns dir: '%s'\n" % coldir)
        return columns.Columns(coldir)

    def output_dir(self):
        return path_join(self.basedir, 'outputs')

    def raw_dtype(self):
        if self.type == 'dr7pofz':
            # note flags is not important, only
            # means rmag > 20, so discard it
            dt=[('run','i2'),
                ('rerun','i2'),
                ('camcol','i1'),
                ('field','i2'),
                ('id','i2'),
                ('photoz_cc2','f4'),
                ('photoz_cc2_err','f4'),
                ('photoz_d1','f4'),
                ('photoz_d1_err','f4'),
                ('flags','i1'),
                ('casid','u8'),
                ('ra','f8'),
                ('dec','f8'),
                ('pz','100f4')]

        return dt

    def read_raw(self, fname):
        path = os.path.expanduser(fname)  
        stdout.write("Reading file: %s\n" % path)

        if not os.path.exists(path):
            raise ValueError("file not found: %s" % path)
        
        stdout.write("Getting nlines.... ")
        stdout.flush()
        nlines = file_len(path)
        stdout.write("%d\n" % nlines)
        stdout.flush()

        if nlines == 0:
            raise ValueError("File is empty: %s" % path)

        # hmm, see if delim ' ' works
        dt = self.raw_dtype()
        rec=recfile.Open(path, 'r', delim=' ', dtype=dt, nrows=nlines)

        data = rec.read()
        rec.close()

        return data


    def dr7pofz_zvals(self):
        arr = numpy.arange(100, dtype='f4')
        return numpy_util.arrscl( arr, 0.03, 1.47 )



    def scinv_dtype(self, n_zlens):
        return [('ra','f8'),('dec','f8'),
                ('mean_scinv','f4',n_zlens)]

    def scinv_array(self, n_zlens, nrows):
        dt = self.scinv_dtype(n_zlens)
        arr = numpy.zeros(nrows, dtype=dt)
        return arr

    def scinv_name(self, rawname):
        newname = rawname.replace('.dat','-scinv.rec')
        return newname

    def add_scinv_to_raw(self, fname, clobber=True):
        """

        Read in the raw file and write out a new file with just ra,dec and the
        mean scinv as a function of lens redshift

        """

        import lensing

        outfile = self.scinv_name(fname)
        if os.path.exists(outfile) and not clobber:
            stdout.write("file exists, skipping\n")
            return
        
        if self.type != 'dr7pofz':
            raise ValueError("only support dr4cc2 for now")

        data = self.read_raw(fname)

        zsvals = self.dr7pofz_zvals()
        zsmax = zsvals.max()
        zsmin = zsvals.min()
        scalc = lensing.tools.ScinvCalculator(zsmin, zsmax)

        output = self.scinv_array(scalc.n_zlens, data.size)

        stdout.write("Copying in common columns\n")
        stdout.flush()
        numpy_util.copy_fields(data, output)

        for i in range(data.size):
            if ((i+1) % 10000) == 0:
                stdout.write("%d/%d\n" % ((i+1), data.size))
                stdout.flush()
            pz = data['pz'][i]

            output['mean_scinv'][i] = scalc.calc_mean_sigmacritinv(zsvals, pz)

        hdr = {'zlvals': list(scalc.zlvals)}

        stdout.write("Writing scinv file: %s\n" % outfile)
        esutil.sfile.write(output, outfile, header=hdr)
        stdout.write("done\n")

    def combined_file(self):
        dir=self.output_dir()
        fname = path_join(dir, self.type+'-combined-scinv.rec')
        return fname

    def scinv_file_pattern(self):
        if self.type == 'dr7pofz':
            return 'pofz.ra*-scinv.rec'
        else:
            raise ValueError("need to define type: %s" % self.type)

    def combine_scinv(self):
        from glob import glob

        outfile=self.combined_file()
        stdout.write("Will write to file: %s\n" % outfile)

        dir=self.output_dir()
        pattern = self.scinv_file_pattern()
        pattern=path_join(dir, pattern)
        flist = glob(pattern)

        datalist = []
        idmin = 0
        for f in flist:
            print f

            tdata = esutil.io.read(f)
            data = numpy_util.add_fields(tdata, [('zid','i4')])
            data['zid'] = idmin + numpy.arange(data.size)
            idmin += data.size

            print data['zid'].min(), data['zid'].max()
            datalist.append(data)

        hdr = esutil.sfile.read_header(flist[0])
        print 'combining data'
        data = numpy_util.combine_arrlist(datalist)

        print 'writing file: %s' % outfile
        esutil.sfile.write(data, outfile, header=hdr)

    def make_columns(self):
        """
        Make a columns database from all available data
        """

        from glob import glob
        import columns
        import sdsspy
        
        coldir = self.columns_dir()

        stdout.write("Will create columns: %s\n" % coldir)

        if os.path.exists(coldir):
            raise ValueError("coldir exists, please start from scratch")

        c=columns.Columns(coldir, verbose=False)


        dir=self.output_dir()
        pattern = self.scinv_file_pattern()
        pattern=path_join(dir, pattern)
        flist = glob(pattern)

        # z values in p(z)
        zvals = list( self.dr7pofz_zvals() )
        pofz_meta = {'zvals': zvals}
        i=1
        ntot=len(flist)
        for scinv_file in sorted(flist):
            stdout.write('-'*70+'\n')
            stdout.write('%s/%s\n' % (i,ntot))
            raw_file = scinv_file.replace('-scinv.rec','.dat')

            # read both scinv and raw, should line up
            stdout.write("Reading scinv data from: %s\n" % scinv_file)
            scinv, scinv_meta = esutil.sfile.read(scinv_file, header=True)
            scinv_meta = {'zlvals': scinv_meta['zlvals']}

            raw = self.read_raw(raw_file)

            if scinv.size != raw.size:
                raise ValueError("scinv and raw file are different sizes")

            photoid = sdsspy.photoid(raw)

            # write data from the raw file
            c.write_column('photoid',photoid, type='rec')
            c.write_column('casid',raw['casid'], type='rec')
            c.write_column('run',raw['run'], type='rec')
            c.write_column('rerun',raw['rerun'], type='rec')
            c.write_column('camcol',raw['camcol'], type='rec')
            c.write_column('field',raw['field'], type='rec')
            c.write_column('id',raw['id'], type='rec')
            c.write_column('ra',raw['ra'], type='rec')
            c.write_column('dec',raw['dec'], type='rec')
            c.write_column('pofz',raw['pz'], type='rec',
                           meta=pofz_meta)

            # now the scinv data
            c.write_column('mean_scinv',scinv['mean_scinv'], type='rec',
                           meta=scinv_meta)

            i+=1





def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return long(result.strip().split()[0])




