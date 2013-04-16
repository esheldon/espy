from __future__ import print_function
import numpy
import os,sys
from sys import stdout

import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1
from esutil import cosmology
#import cosmology

import lensing
from . import files
from . import sigmacrit
from . import regauss

import copy

try:
    import zphot
except:
    pass

def instantiate_sample(sample):
    conf = lensing.files.read_config('scat', sample)
    if 'dr8regauss' in conf['catalog']:
        c = DR8RegaussCatalog(sample)
    else:
        c = DESMockSrcCatalog(sample)
    return c

def create_input(sample):
    """
        create_input('01')
    """
    c = instantiate_sample(sample)
    c.create_objshear_input()

def original_file(sample):
    c = instantiate_sample(sample)
    return c.original_file()

def read_original(sample):
    c = instantiate_sample(sample)
    return c.read_original()


class GenericSrcCatalog(dict):
    def __init__(self, sample, fs='hdfs'):
        conf = lensing.files.read_config('scat', sample)
        for k in conf:
            self[k] = conf[k]

        # not all configs have this option
        self['detrend']=self.get('detrend',False)

        self.open_all_columns()
        if self['sigmacrit_style'] != 2:
            raise ValueError("Expected sigmacrit_style 2")

    def scinv_colname(self):
        return 'scinv%s' % self['sample']


class DR8RegaussCatalog(GenericSrcCatalog):
    """

    Before using, make sure you have matched the regauss cols with your chosen
    photoz sample using lensing.regauss.zphot_match()

    Run add_scinv() before running create_objshear_input()
    (or /bin/add-scinv.py before /bin/make-objshear-input.py)

    """
    def __init__(self, sample, fs='hdfs'):
        super(DR8RegaussCatalog,self).__init__(sample, fs=fs)


        if 'dr8regauss' not in self['catalog']:
            raise ValueError("Expected dr8regauss as catalog")

        self.open_all_columns()


    def get_colnames(self):
        if self['detrend']:
            rmstr='%0.1f' % self['rmag_max']
            rmstr = rmstr.replace('.','')
            e1name = 'e1_rg_dt'+rmstr+'_eq'
            e2name = 'e2_rg_dt'+rmstr+'_eq'
            flagname = 'dt'+rmstr+'_flag'

        else:
            e1name = 'e1_rg_eq'
            e2name = 'e2_rg_eq'
            flagname = 'corrflags_rg'

        magname='cmodelmag_dered_r'
        Rname = 'R_rg_'+self['filter']

        e1name += '_'+self['filter']
        e2name += '_'+self['filter']
        errname = 'uncer_rg_'+self['filter']
        flagname += '_'+self['filter']
        return e1name,e2name,errname,flagname,magname,Rname
    

    def create_objshear_input(self):
        filter=self['filter']
        zlvals=lensing.sigmacrit.make_zlvals(self['dzl'], 
                                             self['zlmin'], 
                                             self['zlmax'])
        nzl = zlvals.size

        keep,zphot_matches = self.select()
        self._keep=keep
        self._zphot_matches=zphot_matches

        print("creating output struct")
        dt = lensing.files.scat_dtype(self['sigmacrit_style'], nzl=nzl)
        output = numpy.zeros(keep.size, dtype=dt)

        print("extracting ra,dec,g1,g2,err and copying into output")

        print("  reading ra,dec,g1,g2,err")
        output['ra'][:] = self.scols['ra'][keep]
        output['dec'][:] = self.scols['dec'][keep]

        e1name,e2name,errname,flagname,magname,Rname=self.get_colnames()
        # the code requires no minus sign for matt's simulations, but
        # does seem to require one for my shapes converted to equatorial
        print("Adding minus sign to e1")
        print(e1name)
        output['g1'][:] = -self.scols[e1name][keep]/2
        print(e2name)
        output['g2'][:] =  self.scols[e2name][keep]/2
        output['err'][:] = self.scols[errname][keep]/2

        output['mag'][:] = self.scols[magname][keep]
        output['R'][:] = self.scols[Rname][keep]

        scinvcol = self.scinv_colname()
        print("  reading",scinvcol)
        output['scinv'][:] = self.scols[scinvcol][keep]
        
        if self['nsplit'] > 0:
            self.split(data=output)
        else:
            lensing.files.scat_write_ascii(sample=self['sample'], data=output)

    def split(self, data=None):
        """
        Split the source file into nsplit parts
        """
        
        nsplit = self['nsplit']
        if nsplit == 0:
            return

        print('splitting into:',self['nsplit'])

        if data is None:
            data = self.read()


        ntot = data.size
        nper = ntot/nsplit
        nleft = ntot % nsplit

        for i in xrange(nsplit):

            beg = i*nper
            end = (i+1)*nper
            if i == (nsplit-1):
                end += nleft
            sdata = data[beg:end]

            lensing.files.scat_write_ascii(sample=self['sample'], data=sdata, src_split=i)


    def read(self, split=None):
        return lensing.files.scat_read_ascii(sample=self['sample'], split=split)

    def select(self):
        """

        Apply selections to the sources, such as regauss processing flags,
        resolution cuts and selecting those with good photoz matches.

        """

        scinvcol = self.scinv_colname()
        if scinvcol not in self.scols:
            raise ValueError("you need to run add-scinv for this sample")

        # first make sure scinv column exists

        filter = self['filter']

        match_column = 'match_zphot%s' % self['pzrun']
        if match_column not in self.scols:
            raise ValueError("First use regauss.zphot_match() to match zphot and regauss")

        e1name,e2name,errname,flagname,magname,Rname=self.get_colnames()

        print("reading:",match_column)
        m = self.scols[match_column][:]
        print("Reading R")
        R = self.scols[Rname][:]


        print("reading e1:",e1name)
        e1 = self.scols[e1name][:]
        print("reading e2:",e2name)
        e2 = self.scols[e2name][:]

        print("Reading corrflags:",flagname)
        flags = self.scols[flagname][:]


        # note this is *not* in a where statement!
        match_logic = (m >= 0)
        wmatch=where1(match_logic)
        print("Found %s/%s zphot matches" % (wmatch.size, m.size))
        

        # sample 5 didn't have range
        print("Getting R logic")
        if 'Rcut' in self:
            R_logic = (R > self['Rcut'])
            wR = where1(match_logic & R_logic)
            print("R > %s: %i/%i" % (self['Rcut'], wR.size, wmatch.size))
        elif 'Rrange' in self:
            Rrange = self['Rrange']
            R_logic = ((R > Rrange[0]) & (R < Rrange[1]))
            wR = where1(match_logic & R_logic)
            print("R in %s: %i/%i" % (Rrange, wR.size, wmatch.size))
        else:
            raise ValueError("need Rcut or Rrange")


        # sample 5 didn't have this cut
        if 'erange' in self:
            print("Getting erange logic")
            erange = self.get('erange',[-4,4])
            erange_logic = ((e1 > erange[0]) 
                            & (e1 < erange[1])
                            & (e2 > erange[0])
                            & (e2 < erange[1]) )
            we=where1(match_logic & erange_logic)
            print("erange %s: %i/%i" % (erange,wR.size, wmatch.size))


        print("Getting flag logic")
        flag_logic = (flags == 0)
        wflag = where1(match_logic & flag_logic)
        print("flags == 0: %i/%i" % (wflag.size, wmatch.size))


        # combined logic
        if 'erange' in self:
            wgood=where1(match_logic & R_logic & erange_logic & flag_logic)
        else:
            wgood=where1(match_logic & R_logic & flag_logic)

        print("Found %s/%s good objects" % (wgood.size, m.size))

        wgood.sort()
        # wgood is index into overall set
        matches = m[wgood]
        return wgood, matches

    def create_pofz_files(self, fs='nfs'):
        """
        For joseph clampett.  Note nfs
        """
        filter=self['filter']
        zlvals=lensing.sigmacrit.make_zlvals(self['dzl'], 
                                             self['zlmin'], 
                                             self['zlmax'])
        nzl = zlvals.size

        pzcols = self.pzcols
        print(pzcols)

        keep,zphot_matches = self.select()

        zphot_matchname = 'match_zphot%s' % self['pzrun']

        nsplit = self['nsplit']
        print('splitting into:',self['nsplit'])

        ntot = zphot_matches.size
        nper = ntot/nsplit
        nleft = ntot % nsplit

        dir=lensing.files.sample_dir(sample=self['sample'], type='pofz', fs='nfs')
        if not os.path.exists(dir):
            print('making dir:',dir)
            os.makedirs(dir)

        tmp=pzcols['pofz'][0:10]
        dt=[('photoid','i8'),
            ('ra','f8'),
            ('dec','f8'),
            ('pofz','f8',tmp.shape[1])]
        for i in xrange(nsplit):

            beg = i*nper
            end = (i+1)*nper
            if i == (nsplit-1):
                end += nleft

            pzind=zphot_matches[beg:end]
            rgind=keep[beg:end]

            pofz = pzcols['pofz'][pzind]
            photoid = pzcols['photoid'][pzind]

            ra=self.scols['ra'][rgind]
            dec=self.scols['dec'][rgind]

            print(photoid.size)
            fname=lensing.files.sample_file(sample=self['sample'], type='pofz', fs='nfs', src_split=i)
            print(fname)

            if os.path.exists(fname):
                os.remove(fname)

            out=numpy.zeros(photoid.size, dtype=dt)
            out['photoid'] = photoid
            out['ra']      = ra
            out['dec']     = dec
            out['pofz']    = pofz

            with eu.recfile.Open(fname, mode='w', delim=' ') as fobj:
                fobj.write(out)

    def add_scinv(self, clobber=False, chunksize=100000):
        """

        Add a column to the REGAUSS db that is the mean scinv as a function of
        lens redshift.

        scinv are created using the *corrected* p(z).  This is separate so I
        can run them in parallel.

        This will be particular to the source sample, and will have name
        scinv{sample}

        The zlvals will be in the meta data, although they are
        easily generated using sigmcrit.make_zlvals(dzl, zlmin, zlmax)

        Procedure
        ---------
        - get the correction factor N(z)/sumpofz(z)
            calls zphot.weighting.pofz_correction(pzrun)

        - for each p(z) multiply by that factor
        - generate mean inverse critical density.


        """

        print("\nopening columns '%s'" % self['procrun'])
        cols = regauss.open_columns(self['procrun'],self['sweeptype'])

        # the column which will hold the inverse critical density.
        # depending on keywords, we might want to raise an error
        colname = 'scinv%s' % self['sample']
        print("Writing to column:\n",colname)
        if colname in cols:
            if not clobber:
                raise ValueError("Column already exists")
            else:
                print("  removing column")
                cols[colname].delete()

        # get the matches
        zphot_matchname = 'match_zphot%s' % self['pzrun']
        if zphot_matchname not in cols:
            raise ValueError("zphot match column not found: '%s'" % zphot_matchname)

        # get the source z vals
        pzconf = zphot.cascade_config(self['pzrun'])
        cosmo = files.read_config('cosmo', self['cosmo_sample'])

        """
        wo = zphot.weighting.WeightedOutputs()
        pofz_file = wo.zhist_file(self['pzrun'])
        print("reading summed p(z)")
        sumpofz_struct = eu.io.read(pofz_file)
        zs = (sumpofz_struct['zmin']+sumpofz_struct['zmax'])/2.
        """

        #
        # correction factor to apply to all p(z)
        #
        print("getting p(z) correction function\n")
        corrstruct = eu.io.read(zphot.weighting.pofz_correction_file(self['pzrun']))
        zs = (corrstruct['zmax']+corrstruct['zmin'])/2.
        corr = corrstruct['corr']


        print("")
        scalc = sigmacrit.ScinvCalculator(zs, self['dzl'], self['zlmin'], self['zlmax'],
                                          H0=cosmo['H0'], omega_m=cosmo['omega_m'])

        zlvals = scalc.zlvals



        print("opening corresponding p(z) columns: '%s'\n" % self['pzrun'])
        pzcols = zphot.weighting.open_pofz_columns(self['pzrun'])

        # work on chunks
        print("using chunksize:",chunksize)
        ntot = cols['run'].size
        nchunk = ntot/chunksize
        if (ntot % chunksize) > 0:
            nchunk += 1

        print("nchunk:",nchunk)
        for i in xrange(nchunk):
            imin = i*chunksize
            imax = (i+1)*chunksize
            print("  ",imin,imax)
            match = cols[zphot_matchname][imin:imax]
            
            nrows = match.size
            scinv = numpy.zeros((nrows, zlvals.size), dtype='f8') -9999.

            w=where1(match >= 0)
            print("    ",w.size,"with good matches")
            if w.size > 0:
                # read the p(z) for the matches
                pofz = pzcols['pofz'][match[w]]

                pofz = pofz[:]*corr[:]

                for j in xrange(w.size):
                    scinv[w[j], :] = scalc.calc_mean_scinv(pofz[j])
            # metadata only gets written once
            cols.write_column(colname, scinv, meta={'zlvals':zlvals})

    def open_all_columns(self):
        print("opening regauss columns for procrun:",self['procrun'],"sweeptype:",self['sweeptype'])
        self.scols = lensing.regauss.open_columns(str(self['procrun']), str(self['sweeptype']))
        print("  #rows:",self.scols['photoid'].size)
        print("opening zphot columns for pzrun:",self['pzrun'])
        self.pzcols = zphot.weighting.open_pofz_columns(str(self['pzrun']))
        print("  #rows:",self.pzcols['photoid'].size)


class DESMockCatalog(dict):
    """
    This is just for reading the original catalog
    """
    def __init__(self, catalog):
        self['catalog'] = catalog

    def dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return d

    def fname(self, scinv=False):
        d = self.dir()
        f='%s-sources.fits' % self['catalog']
        if scinv:
            f=f.replace('sources','sources-scinv')
        infile = path_join(d, f)
        return infile

    def read(self, scinv=False):
        infile = self.fname(scinv=scinv)
        if not os.path.exists(infile):
            raise ValueError("File not found: %s\n" % infile)

        data = eu.io.read(infile, lower=True, verbose=True)
        return data

class DESMockSrcCatalog(dict):
    """
    This reads the mock catalog and creates an input
    catalog for objshear

    Need to have config files that tell us what kind
    of sigma crit inv we are using.
    """

    def __init__(self, sample):
        conf = lensing.files.read_config('scat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['catalog'] not in ['desmocks-2.13o-f8','desmocks-3.02']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

        self['cosmo'] = lensing.files.read_config('cosmo',self['cosmo_sample'])

    def read(self, split=None):
        return lensing.files.scat_read_ascii(sample=self['sample'], split=split)

    def add_scinv(self):
        from . import sigmacrit
        import zphot
        if 'zlmin' not in self or 'zlmax' not in self or 'dzl' not in self:
            raise ValueError("You must have zlmin,zlmax,dzl in config")


        """
        note we don't really need this because we have perfect zs, but this
        is what you would do with a real sample

        wo = zphot.WeightedOutputs()
        z_file = wo.z_file(self['pzrun'], chunk=0)
        zscols = zphot.weighting.read_z(z_file)
        zsvals = (zscols['zmax']+zscols['zmin'])/2.0
        
        sc = sigmacrit.ScinvCalculator(zsvals, 
                                       self['dzl'],
                                       self['zlmin'],
                                       self['zlmax'])
        """

        # ok, now read the original file and add to it
        # the man scinv as a function of zlens
        # for now just using exact z

        zlvals = sigmacrit.make_zlvals(self['dzl'],
                                       self['zlmin'],
                                       self['zlmax'])
        n_zlens=zlvals.size
        
        data = self.read_original()
        dt=data.dtype.descr
        new_dtype = dt + [('scinv','f8',n_zlens)]
        out = numpy.zeros(data.size, dtype=new_dtype)
        eu.numpy_util.copy_fields(data, out)


        # actually not using sc since redshifts are exact
        cosmo=self['cosmo']
        c = cosmology.Cosmo(omega_m=cosmo['omega_m'], H0=cosmo['H0'])

        print('adding scinv to each')

        for i in xrange(data.size):
            if ((i+1) % 1000) == 0:
                print("%s/%s  %s%%" % (i+1,data.size,float(i+1)/data.size*100.))

            zs = data['z'][i]
            out['scinv'][i,:] = c.sigmacritinv(zlvals, zs)

        h=copy.deepcopy(cosmo)
        h['zlvals'] = zlvals

        outf=self.original_file(scinv=True)
        print("writing scinv file:",outf)

        eu.io.write(outf, out, header=h)

    def create_objshear_input(self):

        print("sigmacrit_style:",self['sigmacrit_style'])
        if self['sigmacrit_style'] == 2:
            data = self.read_original(scinv=True)
            print('creating output array')
            output = self.output_array(data.size, nzl=data['scinv'][0,:].size)
        else:
            data = self.read_original()
            print('creating output array')
            output = self.output_array(data.size)


        print('copying data')
        output['ra'] = data['ora']
        output['dec'] = data['odec']

        print('not changing sign of gamma1')
        output['g1'] = data['gamma1']

        output['g2'] = data['gamma2']
        output['err'] = 0.0
        #output['hpixid'] = -9999

        if self['sigmacrit_style'] == 1:
            output['z'] = data['z']
            #output['dc'] = -9999
        else:
            output['scinv'] = data['scinv']

        lensing.files.scat_write_ascii(sample=self['sample'], data=output)

        if self['nsplit'] > 0:
            self.split(data=output)

    def split(self, data=None):
        """
        Split the source file into nsplit parts
        """
        
        nsplit = self['nsplit']
        if nsplit == 0:
            return
        print('splitting into:',self['nsplit'])

        data = self.read()

        fname=self.file()

        ntot = data.size
        nper = ntot/nsplit
        nleft = ntot % nsplit

        for i in xrange(nsplit):

            beg = i*nper
            end = (i+1)*nper
            if i == (nsplit-1):
                end += nleft
            sdata = data[beg:end]

            lensing.files.scat_write_ascii(sample=self['sample'], data=sdata, src_split=i)


    def original_file(self, scinv=False):
        dmc=DESMockCatalog(self['catalog'])
        return dmc.fname(scinv=scinv)

    def read_original(self, scinv=False):
        dmc=DESMockCatalog(self['catalog'])
        return dmc.read(scinv=scinv)

    def output_array(self, num, nzl=None):
        dt = lensing.files.scat_dtype(self['sigmacrit_style'], nzl=nzl)
        output = numpy.zeros(num, dtype=dt)
        return output


