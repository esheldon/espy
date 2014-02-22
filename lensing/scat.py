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
    elif 'im3' in conf['catalog']:
        c=IM3ShapePointz(sample)
    elif 'gmix' in conf['catalog']:
        c = DR8GMixCatalog(sample)
    elif 'desmocks' in conf['catalog']:
        c = DESMockSrcCatalog(sample)
    else:
        raise ValueError('bad catalog: %s' % conf['catalog'])
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

def scinv_colname(scinv_sample):
    return 'scinv%s' % scinv_sample

def zphot_matchname(pzrun):
    zphot_matchname = 'match_zphot%s' % pzrun
    return zphot_matchname 

class GenericSrcCatalog(dict):
    def __init__(self, sample, fs='hdfs'):
        conf = lensing.files.read_config('scat', sample)
        for k in conf:
            self[k] = conf[k]

        self['scinv_config'] = \
            lensing.files.read_config('scinv', conf['scinv_sample'])
        # not all configs have this option
        self['detrend']=self.get('detrend',False)

        self.open_all_columns()
        if self['sigmacrit_style'] != 2:
            raise ValueError("Expected sigmacrit_style 2")

        self.fs=fs

    def dir_original(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return d

    def open_all_columns(self):
        raise RuntimeError("over-ride open_all_columns")
    def create_objshear_input(self):
        raise RuntimeError("over-ride create_objshear_input")

    def read(self, split=None):
        return lensing.files.scat_read_ascii(sample=self['sample'], split=split)

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

            lensing.files.scat_write_ascii(sample=self['sample'],
                                           data=sdata,
                                           src_split=i,
                                           fs=self.fs)


class DR8GMixCatalog(GenericSrcCatalog):
    def __init__(self, sample, fs='hdfs'):
        super(DR8GMixCatalog,self).__init__(sample, fs=fs)

        if 'gmix' not in self['catalog']:
            raise ValueError("Expected gmix as catalog")

    def create_objshear_input(self):
        """
        - Trim to selected
        - extract columns
        - write to file or split files
        """
        keep = self.select()

        sconf=self['scinv_config']
        zlvals=lensing.sigmacrit.make_zlvals(sconf['dzl'], 
                                             sconf['zlmin'], 
                                             sconf['zlmax'])
        nzl = zlvals.size

        print("creating output struct")
        dt = lensing.files.scat_gmix_dtype(self['sigmacrit_style'], nzl=nzl)
        output = numpy.zeros(keep.size, dtype=dt)

        scinvcol = scinv_colname(self['scinv_sample'])
        print("  reading",scinvcol)
        scinv = self.scols[scinvcol][keep]
        output['scinv'][:] = scinv
        del scinv

        print("  ra,dec")
        ra = self.scols['ra'][keep]
        output['ra'][:] = ra
        del ra

        dec = self.scols['dec'][keep]
        output['dec'][:] = dec
        del dec

        print("reading g")
        g=self.scols['g_dt_eq'][keep]
        output['g1'] = g[:,0]
        output['g2'] = g[:,1]
        del g

        print("reading gcov")
        gcov=self.scols['gcov'][keep]
        output['gcov11'] = gcov[:,0,0]
        output['gcov12'] = gcov[:,0,1]
        output['gcov22'] = gcov[:,1,1]

        del gcov

        if 'gsens' in self.scols:
            print("reading gsens")
            gsens=self.scols['gsens'][keep]
            output['gsens1'][:] = gsens[:,0]
            output['gsens2'][:] = gsens[:,1]
            del gsens
        else:
            print('gsens column not found, using 1.0')
            output['gsens1'][:] = 1.0
            output['gsens2'][:] = 1.0


        # the code requires no minus sign for matt's simulations, but
        # does seem to require one for my shapes converted to equatorial
        print("Adding minus sign to g1")
        output['g1'] *= (-1)

        
        if self['nsplit'] > 0:
            self.split(data=output)
        else:
            lensing.files.scat_write_ascii(sample=self['sample'], data=output)

    def select(self):
        import gmix_sdss
        selector=gmix_sdss.cuts.Selector(self.scols,self['psf_model'])
        selector.do_select(pzrun=self['pzrun'],
                           s2n_min=self['s2n_range'][0],
                           s2n_max=self['s2n_range'][1],
                           sratio_min=self['sratio_min'])

        print("kept",selector.indices.size)
        return selector.indices

    def open_all_columns(self):
        import gmix_sdss
        conf=gmix_sdss.files.read_config(self['catalog'])
        self['psf_model'] = conf['psf_model']

        self.scols = gmix_sdss.collate.open_columns(self['catalog'])
        print("  #rows:",self.scols['photoid'].size)

        print("opening zphot columns for pzrun:",self['pzrun'])
        self.pzcols = zphot.weighting.open_pofz_columns(str(self['pzrun']))
        print("  #rows:",self.pzcols['photoid'].size)


class DR8RegaussCatalog(GenericSrcCatalog):
    """

    Before using, make sure you have matched the regauss cols with your chosen
    photoz sample using lensing.regauss.zphot_match()
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

        scinvcol = scinv_colname(self['scinv_sample'])
        print("  reading",scinvcol)
        output['scinv'][:] = self.scols[scinvcol][keep]
        
        if self['nsplit'] > 0:
            self.split(data=output)
        else:
            lensing.files.scat_write_ascii(sample=self['sample'], data=output)



    def select(self):
        """

        Apply selections to the sources, such as regauss processing flags,
        resolution cuts and selecting those with good photoz matches.

        """

        scinvcol = scinv_colname(self['scinv_sample'])
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

class IM3ShapePointz(GenericSrcCatalog):
    def __init__(self, sample, fs='nfs'):
        conf = lensing.files.read_config('scat', sample)
        for k in conf:
            self[k] = conf[k]

        if 'scinv_sample' in self:
            self['scinv_config'] = \
                    lensing.files.read_config('scinv', conf['scinv_sample'])

        self.pixscale=0.263
        self.fs=fs

    def fname_original(self):
        d = self.dir_original()
        f='011-im3shape-3_%s_photoz_gold_Feb12.fits' % self['band']
        fname = path_join(d, f)
        return fname

    def read_original(self):
        fname=self.fname_original()
        print("reading:",fname)
        return eu.io.read(fname, ext=1, lower=True)

    def create_objshear_input(self):
        """
        - Trim to selected
        - write to file or split files
        """

        orig = self.read_original()
        output = self.create_output(orig)

        # the code requires no minus sign for matt's simulations, but
        # does seem to require one for my shapes converted to equatorial
        #print("Adding minus sign to g1")
        #output['g1'] *= (-1)
        print("Adding minus sign to g2")
        output['g2'] *= (-1)
        
        if self['nsplit'] > 0:
            self.split(data=output)
        else:
            lensing.files.scat_write_ascii(sample=self['sample'],
                                           data=output,
                                           fs=self.fs)

    def get_good(self, data, do_srat=True):
        from numpy import where, isnan

        ntot=data.size

        # remove size and snr problems entirely from cat so we don't divide by
        # zero

        ws,=where(  (data['psf_fwhm']*self.pixscale > 0) & (data['radius'] > 0) )
        print("sizes cut: %s/%s" % (ws.size,data.size))

        data=data[ws]

        w,=where( data['snr'] > self['min_s2n'] )
        print("s/n cut %s: %s/%s" % (self['min_s2n'],w.size,ntot))

        data = data[w]

        ntot=data.size

        nan_test=(isnan(data['snr']) == False)
        w,=where( nan_test )
        print("nan cut: %s/%s" % (w.size,ntot))


        esq = data['s1']**2 + data['s2']**2
        e_test = ( esq < ( self['max_ellip']**2) )
        w,=where( e_test )
        print("e cut %s: %s/%s" % (self['max_ellip'],w.size,ntot))

        flags=self.get_flags()
        flag_test=( (data['flag'] & flags) == 0 )
        w,=where( flag_test )
        print("flag cut %s: %s/%s" % (flags,w.size,ntot))

        if do_srat: 
            srat = 2*data['radius']/(data['psf_fwhm']*self.pixscale)
            srat_test=(srat > self['min_srat'])
            w,=where( srat_test )
            print("srat cut %s: %s/%s" % (self['min_srat'],w.size,ntot))

        sx_flags_test=(data['flags_r']==0) & (data['flags_i']==0)
        w,=where( sx_flags_test )
        print("sx flags cut: %s/%s" % (w.size,ntot))

        sg_test = self.get_sg_flags(data)
        w,=where( sg_test )
        print("sg cut: %s/%s" % (w.size,ntot))

        tests = (nan_test & flag_test & sx_flags_test & e_test & sg_test )
        if do_srat:
            tests = tests & srat_test
        w,=where(tests)


        print("finally keeping %s/%s" % (w.size, ntot))

        data = data[w]
        return data

    def get_flags(self):
        each_flag = [2**bit for bit in self['flags']]
        flags=sum(each_flag)
        return flags

    def get_sg_flags(self, data):
        class_star_i = data['class_star_i']
        spread_model_i = data['spread_model_i']
        spreaderr_model_i = data['spreaderr_model_i']
        mag_auto_i = data['mag_auto_i']
        #mag_psf_i = data['mag_psf_i']

        #sg_testp = (  ((class_star_i > 0.3) & (mag_auto_i < 18.0)) 
        #            | ((spread_model_i + 3*spreaderr_model_i) < 0.003) |  ( (mag_psf_i > 30.0) & (mag_auto_i < 21.0) )   )

        sg_testp = (  ((class_star_i > 0.3) & (mag_auto_i < 18.0)) 
                    | ((spread_model_i + 3*spreaderr_model_i) < 0.003)  )

        sg_test = (sg_testp == False) & (numpy.abs(spread_model_i) < 0.1)

        return sg_test



    def create_output(self, orig):
        data = self.get_good(orig)

        assert self['sigmacrit_style']==1
        dt=files.scat_im3shape_dtype(self['sigmacrit_style'])
        output = numpy.zeros(data.size, dtype=dt)

        weights = self.get_weights(data)

        output['ra']     = data['ra']
        output['dec']    = data['dec']
        output['g1']     = data['s1']
        output['g2']     = data['s2']
        output['weight'] = weights
        output['z']      = data['zp_2']

        return output

    def get_weights(self, data):
        snr = data['snr']
        eprox2 = (0.1/(snr/40.))**2
        weight = 1.0/(0.2**2 + eprox2)

        return weight

    def plot_size_mag(self, data, nplot=100000, **keys):
        import biggles

        mag_range=[15,25]
        rad_range=[0,6]

        r=numpy.random.random(data.size)
        s=r.argsort()
        ind=s[0:nplot]

        pts=biggles.Points(data['mag_auto_i'][ind], data['radius'][ind],
                           type='dot')

        plt=biggles.FramedPlot()
        plt.xtitle=r'mag_auto$_i$'
        plt.ytitle=r'$r_{1/2} [arcsec]$'
        plt.aspect_ratio=1
        plt.xrange=mag_range
        plt.yrange=rad_range

        title=keys.get('title',None)
        if title:
            plt.title=title

        plt.add( pts )

        write_plot(plt, **keys)

    def plot_sizes(self, data, binsize=0.01, **keys):
        import biggles

        rad_range=[0,5]

        wts=self.get_weights(data)
        wts *= (1./wts.max())
        wts1=numpy.ones(data.size)

        bs=eu.stat.Binner(data['radius'],
                          weights=wts1)
        bs.dohist(binsize=binsize,
                  min=rad_range[0],
                  max=rad_range[1])
        bs.calc_stats()

        wbs=eu.stat.Binner(data['radius'], weights=wts)
        wbs.dohist(binsize=binsize,
                   min=rad_range[0],
                   max=rad_range[1])
        wbs.calc_stats()

        hplt=biggles.Histogram(bs['whist'], 
                               x0=rad_range[0], binsize=binsize,
                               color='red')
        whplt=biggles.Histogram(wbs['whist'], 
                                x0=rad_range[0], binsize=binsize,
                                color='blue')

        hplt.label='unweighted'
        whplt.label='weighted'

        key=biggles.PlotKey(0.9, 0.9, [hplt,whplt],
                            halign='right')


        plt=biggles.FramedPlot()
        plt.add( hplt, whplt, key )


        plt.xtitle=r'$r_{1/2} [arcsec]$'
        plt.aspect_ratio=1

        title=keys.get('title',None)
        if title:
            plt.title=title

        write_plot(plt, **keys)

    def plot_srat(self, data, binsize=0.01, **keys):
        import biggles

        srat = 2*data['radius']/(data['psf_fwhm']*self.pixscale)

        srat_range=[0,10]

        wts=self.get_weights(data)
        wts *= (1./wts.max())
        wts1=numpy.ones(data.size)

        bs=eu.stat.Binner(srat, 
                          weights=wts1)
        bs.dohist(binsize=binsize,
                  min=srat_range[0],
                  max=srat_range[1])
        bs.calc_stats()

        wbs=eu.stat.Binner(srat, weights=wts)
        wbs.dohist(binsize=binsize,
                   min=srat_range[0],
                   max=srat_range[1])
        wbs.calc_stats()

        hplt=biggles.Histogram(bs['whist'], 
                               x0=srat_range[0], binsize=binsize,
                               color='red')
        whplt=biggles.Histogram(wbs['whist'], 
                                x0=srat_range[0], binsize=binsize,
                                color='blue')

        hplt.label='unweighted'
        whplt.label='weighted'

        key=biggles.PlotKey(0.9, 0.9, [hplt,whplt],
                            halign='right')


        plt=biggles.FramedPlot()
        plt.add( hplt, whplt, key )


        plt.xtitle=r'$FWHM_{gal}/FWHM_{PSF}$'
        plt.aspect_ratio=1

        title=keys.get('title',None)
        if title:
            plt.title=title

        write_plot(plt, **keys)




    def plot_ellip(self, data, binsize=0.005, **keys):
        import biggles

        e_range=[0,1]

        wts=self.get_weights(data)
        wts *= (1./wts.max())
        wts1=numpy.ones(data.size)

        e=numpy.sqrt( data['s1']**2 + data['s2']**2)
        bs=eu.stat.Binner(e, weights=wts1)
        bs.dohist(binsize=binsize,
                  min=e_range[0],
                  max=e_range[1])
        bs.calc_stats()

        wbs=eu.stat.Binner(e, weights=wts)
        wbs.dohist(binsize=binsize,
                   min=e_range[0],
                   max=e_range[1])
        wbs.calc_stats()

        hplt=biggles.Histogram(bs['whist'], 
                               x0=e_range[0], binsize=binsize,
                               color='red')
        whplt=biggles.Histogram(wbs['whist'], 
                                x0=e_range[0], binsize=binsize,
                                color='blue')

        hplt.label='unweighted'
        whplt.label='weighted'

        key=biggles.PlotKey(0.2, 0.1, [hplt,whplt])


        plt=biggles.FramedPlot()
        plt.add( hplt, whplt, key )

        plt.xtitle='|e|'
        plt.aspect_ratio=1

        title=keys.get('title',None)
        if title:
            plt.title=title

        write_plot(plt, **keys)


def write_plot(plt, show=False, eps=None,
               png=None, width=800, height=800,
               convert=False, dpi=90, **keys):
    """
    If convert=True and eps keyword is sent,
    then the converter is used.
    
    Otherwise send png
    """
    if show:
        plt.show()

    if eps is not None:
        print(eps)
        plt.write_eps(eps)

    if eps is not None and convert:
        import converter
        converter.convert(eps, dpi=dpi, verbose=True)
    elif png is not None:
        print(png)
        plt.write_img(width, height, png)

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
            raise ValueError("The config sample '%s' doesn't "
                             "match input '%s'" % (self['sample'],sample))

        self['cosmo'] = lensing.files.read_config('cosmo',self['cosmo_sample'])

    def read(self, split=None):
        return lensing.files.scat_read_ascii(sample=self['sample'], split=split)

    def add_scinv(self):
        raise RuntimeError("note using scinv configs not, need to update "
                           "this.  see gmix_sdss/bin/add-scinv.py")
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

        ntot = data.size
        nper = ntot/nsplit
        nleft = ntot % nsplit

        for i in xrange(nsplit):

            beg = i*nper
            end = (i+1)*nper
            if i == (nsplit-1):
                end += nleft
            sdata = data[beg:end]

            lensing.files.scat_write_ascii(sample=self['sample'],
                                           data=sdata,
                                           src_split=i)


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


