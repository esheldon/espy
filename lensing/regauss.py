"""

classes:

    RegaussSweep
    RegaussAtlas
    SweepCache
    Collator

    RunTester
    Tester

"""
from __future__ import print_function
import os
import glob
import sdsspy
import es_sdsspy
import columns
import numpy
import esutil as eu
from esutil.numpy_util import where1
from esutil.ostools import path_join

import biggles
from biggles import FramedPlot, PlotKey, Table, PlotLabel, Points, \
            SymmetricErrorBarsY as SymErrY, SymmetricErrorBarsX as SymErrX

import images
import fimage
import admom

def open_columns(procrun):
    coll = Collator(procrun)
    return coll.open_columns()

def output_dir(procrun):
    coll = Collator(procrun)
    return coll.output_dir()

class RegaussSweep:
    def __init__(self, objs, **keys):
        self.objs = objs
        self.ra = RegaussAtlas()

        self.rmax = keys.get('rmax',22.)

        self.verbosity = keys.get('verbosity',1)

    def select(self):
        print("Selecting objects")
        s = es_sdsspy.select.Selector(self.objs)
        print("  getting resolve logic")
        resolve_logic = s.resolve_logic()
        print("  getting tycho logic")
        tycho_logic = s.mask_logic('tycho')

        print("  getting flag logic")
        flag_logic = s.flag_logic()

        print("  getting rmag logic")
        rmag_logic = s.cmodelmag_logic("r", self.rmax)

        logic = \
            resolve_logic & flag_logic & rmag_logic

        keep = where1(logic)
        print("  keeping %i/%i" % (keep.size, self.objs.size))

        if keep.size > 0:
            # original objs will live on with the caller
            self.objs = self.objs[keep]
        else:
            self.objs = None

    def process(self):

        self.select()
        if self.objs == None:
            return None

        ra = self.ra

        run=self.objs['run'][0]
        camcol=self.objs['camcol'][0]

        objs = self.objs
        self.make_output()

        b=eu.stat.Binner(objs['field'])
        b.dohist(binsize=1, rev=True)
        
        rev = b['rev']
        print("Processing objects by field")
        for i in xrange(b['hist'].size):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                field = objs['field'][w[0]]
                print("  Processing %4i in %06i-%i-%04i" % (w.size,run,camcol,field))
                
                for index in w:

                    obj = objs[index]
                    id = obj['id']
                    if ra.has_atlas(run,camcol,field,id):

                        if self.verbosity > 1:
                            print("    Processing %06i-%i-%04i-%04i" % (run,camcol,field,id))
                        self.output['has_atlas'][index] = 1

                        for fnum in [0,1,2,3,4]:
                            if self.verbosity > 1:
                                print("        ",sdsspy.FILTERCHAR[fnum])

                            rg = ra.regauss_obj(obj, fnum)
                            if rg is not None:
                                self.copy2output(index, rg, fnum)
        return self.output

    def copy2output(self, index, rg, fnum):
        
        # get a reference to this row, modifications to this will
        # also affect self.output
        data = self.output[index]

        if 'imstats' in rg:
            s = rg['imstats']
            if s is not None:
                data['wrow'][fnum]    = s['wrow']
                data['wcol'][fnum]    = s['wcol']
                data['Irr'][fnum]     = s['Irr']
                data['Irc'][fnum]     = s['Irc']
                data['Icc'][fnum]     = s['Icc']
                data['a4'][fnum]      = s['a4']
                data['uncer'][fnum]   = s['uncer']
                data['amflags'][fnum] = s['whyflag']
                data['amflags_str'][fnum] = s['whystr']
                data['numiter'][fnum] = s['numiter']
                
        if 'psfstats' in rg:
            s = rg['psfstats']
            if s is not None:
                data['Irr_psf'][fnum]     = s['Irr']
                data['Irc_psf'][fnum]     = s['Irc']
                data['Icc_psf'][fnum]     = s['Icc']
                data['a4_psf'][fnum]      = s['a4']
                data['amflags_psf'][fnum] = s['whyflag']
                data['amflags_psf_str'][fnum] = s['whystr']
 
        if 'corrstats' in rg:
            s = rg['corrstats']
            if s is not None:
                data['e1_lin'][fnum]        = s['e1']
                data['e2_lin'][fnum]        = s['e2']
                data['R_lin'][fnum]         = s['R']
                data['corrflags_lin'][fnum] = s['flags']


        if 'rgstats' in rg:
            s = rg['rgstats']
            if s is not None:
                data['f0flags'][fnum]    = s['f0flags']
                if s['f0flags'] == 0:
                    # we only get here if we were able to make f0
                    data['Irr_rg'][fnum]     = s['Irr']
                    data['Irc_rg'][fnum]     = s['Irc']
                    data['Icc_rg'][fnum]     = s['Icc']
                    data['a4_rg'][fnum]      = s['a4']
                    data['uncer_rg'][fnum]   = s['uncer']
                    data['amflags_rg'][fnum] = s['whyflag']
                    data['amflags_rg_str'][fnum] = s['whystr']
                    data['numiter_rg'][fnum] = s['numiter']
     

        if 'rgcorrstats' in rg:
            s = rg['rgcorrstats']
            if s is not None:
                data['e1_rg'][fnum]        = s['e1']
                data['e2_rg'][fnum]        = s['e2']
                data['R_rg'][fnum]         = s['R']
                data['corrflags_rg'][fnum] = s['flags']



    def make_output(self):
        print("creating output")
        dtype = self.output_dtype()
        self.output = numpy.zeros(self.objs.size, dtype=dtype)

        print("    copying from objs")
        eu.numpy_util.copy_fields(self.objs, self.output)

        # set defaults
        output = self.output

        print("    setting defaults")
        # these get set to 9999
        flist = ['uncer','uncer_rg']
        for f in flist:
            for fnum in [0,1,2,3,4]:
                output[f][:,fnum] = 9999

        # these get set to 2**15 = -32768 which means not processed
        flist = ['amflags','amflags_psf','amflags_rg',
                 'corrflags_lin','corrflags_rg']
        for f in flist:
            for fnum in [0,1,2,3,4]:
                output[f][:,fnum] = 2**15

        # these get set to -9999
        flist = ['wrow','wcol',
                 'Irr','Irc','Icc','a4',
                 'Irr_psf','Irc_psf','Icc_psf','a4_psf',
                 'Irr_rg','Irc_rg','Icc_rg','a4_rg',
                 'e1_lin','e2_lin','R_lin',
                 'e1_rg','e2_rg','R_rg']
        for f in flist:
            for fnum in [0,1,2,3,4]:
                output[f][:,fnum] = -9999

    def output_dtype(self):
        dt=[('run','i2'),
            ('rerun','i2'),
            ('camcol','i1'),
            ('field','i2'),
            ('id','i2'),
            ('thing_id','i4'),
            ('ra','f8'),  # do we need these?
            ('dec','f8'),

            ('has_atlas','i1'), # zero for no atlas

            ('wrow','5f8'),
            ('wcol','5f8'),
            ('Irr','5f8'),
            ('Irc','5f8'),
            ('Icc','5f8'),
            ('a4','5f8'),
            ('uncer','5f8'),
            ('numiter','5i1'),
            ('amflags','5i2'),
            ('amflags_str','5S5'),

            ('Irr_psf','5f8'),
            ('Irc_psf','5f8'),
            ('Icc_psf','5f8'),
            ('a4_psf','5f8'),
            ('amflags_psf','5i2'),
            ('amflags_psf_str','5S5'),

            ('e1_lin','5f8'),
            ('e2_lin','5f8'),
            ('R_lin','5f8'),
            ('corrflags_lin','5i2'),  # flags during compea4 correction

            ('f0flags','5i1'),
            ('Irr_rg','5f8'),
            ('Irc_rg','5f8'),
            ('Icc_rg','5f8'),
            ('a4_rg','5f8'),
            ('uncer_rg','5f8'),
            ('numiter_rg','5i1'),
            ('amflags_rg','5i2'), # flags running admom on f0-epsilon in the rg process
            ('amflags_rg_str','5S5'),

            ('e1_rg','5f8'),    # corrected ellipticity
            ('e2_rg','5f8'),
            ('R_rg','5f8'),     # this is more stable than R_lin
            ('corrflags_rg','5i2')]  # flags during compea4 correction

        return dt


class RegaussAtlas:
    def __init__(self, verbose=False):
        self.verbose=verbose

        self.sweep_cache = SweepCache(verbose=self.verbose)
        self.atls_key = None
        self.field_key = None


    def regauss_id(self, type, run, camcol, field, id, filter):
        """
        Run regauss on a single object in a single band
        """

        obj = self.readobj(type, run, camcol, field, id)
        return self.regauss_obj(obj, filter)

    def regauss_obj(self, obj, filter):
        run = obj['run']
        camcol = obj['camcol']
        field = obj['field']
        id = obj['id']

        c = sdsspy.FILTERNUM[filter]

        self.cache_psf(run, camcol, field)

        try:
            self.cache_atlas(run, camcol, field, id)
        except NoAtlasImageError:
            # we can safely ignore such errors, but we return
            # nothing
            return None

        im = numpy.array( self.atls['images'][c], dtype='f8')
        im -= self.atls['SOFT_BIAS']


        rowc=obj['rowc'][c]
        colc=obj['colc'][c]
        row=rowc - self.atls['row0'][c] - 0.5
        col=colc - self.atls['col0'][c] - 0.5

        sigsky = self.psfield['skysig'][0,c]
        psf = self.psfKLs[c].rec(rowc, colc, trim=True)

        # sometimes we get very negative pixels, not sure why
        #im=im.clip(-3*sigsky, im.max())

        guess_psf=admom.fwhm2mom( obj['psf_fwhm'][c], pixscale=0.4)/2.
        guess = guess_psf
        if 'm_rr_cc' in obj.dtype.names:
            if obj['m_rr_cc'][c] > guess_psf:
                guess=obj['m_rr_cc'][c]/2

        rg = admom.ReGauss(im, row, col, psf, 
                           sigsky=sigsky, guess=guess, guess_psf=guess_psf,
                           verbose=self.verbose)
        rg.do_all() 

        # keep these things in case we want to call show()
        self.rg=rg
        self.im=im
        self.psf=psf

        return rg

    def show(self, min=0, **keys):
        plt=images.view(self.im, show=False, min=min, **keys)
        s=self.rg['imstats']

        if s['whyflag'] == 0:
            levels=7

            wrow = s['wrow']
            wcol = s['wcol']

            model = fimage.model_image('gauss',
                                       self.im.shape,
                                       [wrow,wcol],
                                       [s['Irr'],s['Irc'],s['Icc']],
                                       counts=self.im.sum())

            cmod = biggles.Contours(model.transpose(), color='grey')
            cmod.levels = levels
            plt.add(cmod)
        plt.show()


    def readobj(self, type, run, camcol, field, id):
        obj = self.sweep_cache.readobj(type, run, camcol, field, id)
        return obj

    def has_atlas(self, run, camcol, field, id):
        try:
            self.cache_atlas(run, camcol, field, id)
            return True
        except NoAtlasImageError:
            # we can safely ignore such errors
            return False


    def cache_atlas(self, run, camcol, field, id):
        """
        This will cache reads so it can be used for multiple
        filter processings
        """
        key = '%06i-%i-%04i-%05i' % (run,camcol,field,id)
        if self.atls_key != key:
            self.atls = sdsspy.read('fpAtlas',run,camcol,field,id, 
                                    trim=True,
                                    verbose=self.verbose)
            self.atls_key = key

    def cache_psf(self, run, camcol, field):
        """

        This will cache reads so it can be used for multiple filter
        processings.  So process all objects from a field in order to be most
        efficient.

        """
        key = '%06i-%i-%04i' % (run,camcol,field)
        if self.field_key != key:
            self.psfield = sdsspy.read('psfield',run,camcol,field, 
                                       verbose=self.verbose)
            self.psfKLs = []
            for filter in sdsspy.FILTERCHARS:
                kl = sdsspy.read('psField',run,camcol,field,filter=filter,
                                 verbose=self.verbose)
                self.psfKLs.append(kl)

            self.field_key = key

def regauss_atlas(type, run, camcol, field, id, filter, 
                  show=False, showmin=0, verbose=False, **keys):
    """

    Run regauss on a given atlas image.  This is not optimal, since only a
    single flter is processed at a time, but is more for interactive work

    """

    rga = RegaussAtlas(verbose=verbose)
    rga.regauss1(type, run, camcol, field, id, filter)

    if show:
        rga.show(min=showmin)

    return rga.rg


def admom_atlas(type, run, camcol, field, id, filter, 
                objs=None,
                show=False, showmin=0, **keys):
    """
    Objs must be for this camcol!
    """
    c = sdsspy.FILTERNUM[filter]

    if objs is None:
        sc=SweepCache()
        obj = sc.readobj(type, run, camcol, field, id)
    else:
        w=where1( (objs['field'] == field) & (objs['id'] == id) )
        if w.size == 0:
            raise ValueError("Object not found objs: "
                             "%06i-%i-%04i" % (run,camcol,id))
        obj = objs[w[0]]

    atls = sdsspy.read('fpAtlas',run,camcol,field,id, trim=True)
    im = numpy.array( atls['images'][c], dtype='f8') - atls['SOFT_BIAS']

    psfield = sdsspy.read('psField', run, camcol, field, lower=True)

    rowc=obj['rowc'][c]
    colc=obj['colc'][c]

    row=rowc - atls['row0'][c] - 0.5
    col=colc - atls['col0'][c] - 0.5

    sigsky = psfield['skysig'][0,c]
    Tguess_psf=admom.fwhm2mom( obj['psf_fwhm'][c], pixscale=0.4)
    if obj['m_rr_cc'][c] > Tguess_psf:
        Tguess=obj['m_rr_cc'][c]
    else:
        Tguess=Tguess_psf

    out = admom.admom(im,
                      row, 
                      col, 
                      sky=0.0,
                      sigsky=sigsky,
                      Tguess=Tguess)

    if show:
        plt=images.view(im, min=showmin, show=False, **keys)

        if out['whyflag'] == 0:
            levels=7

            wrow = out['wrow']
            wcol = out['wcol']

            model = fimage.model_image('gauss',
                                       im.shape,
                                       [wrow,wcol],
                                       [out['Irr'],out['Irc'],out['Icc']],
                                       counts=im.sum())


            cmod = biggles.Contours(model.transpose(), color='grey')
            cmod.levels = levels
            plt.add(cmod)
        plt.show()

    return out



class SweepCache:
    def __init__(self, type=None, data=None, verbose=False):
        self.verbose=verbose
        self.init(type=type, data=data)

    def init(self, type=None, data=None):
        if type is not None and data is not None:
            run=data['run'][0]
            camcol=data['camcol'][0]
            key = '%s-%06i-%d' % (type,run,camcol)
            self.key = key
            self.data = data
        else:
            self.key=None
            self.data=None


    def read(self,type,run,camcol, **keys):
        self.cache_column(type, run, camcol, **keys)
        return self.data

    def readfield(self, type, run, camcol, field, **keys):    
        self.cache_column(type, run, camcol, **keys)
        data = self.data
        w=where1(data['field'] == field)
        if w.size == 0:
            raise ValueError("field not found: %06i-%i-%04i" % (run,camcol,field))
        return data[w]

    def readobj(self, type, run, camcol, field, id, **keys):    
        self.cache_column(type, run, camcol, **keys)
        data = self.data
        w=where1( (data['field'] == field) & (data['id'] == id) )
        if w.size == 0:
            raise ValueError("object not found: %06i-%i-%04i-%05i" % (run,camcol,field,id))
        return data[w[0]]


    def cache_column(self, type,run,camcol, **keys):
        key = '%s-%06i-%d' % (type,run,camcol)
        if self.key != key:
            data = sdsspy.read('calibobj.%s' % type,
                               run,camcol, lower=True, verbose=self.verbose, 
                               **keys)
            self.key = key
            self.data = data


class RunTester(dict):
    """
    testing a single run
    """
    def __init__(self,procrun,band, run='any'):
        self.procrun = procrun
        self.run = run
        self.band = band


    def load_data(self):
        if 'r' not in self:
            print("Opening columns")
            c = open_columns(self.procrun)

            if self.run != 'any':
                print('Getting indices of run:',self.run)
                w=(c['run'].match(self.run))

                if w.size == 0:
                    raise ValueError("Run",self.run,"not found")
            else:
                w=None

            # basic tags to load
            cn_byband=['r','whyflag_rg','e1_rg','e2_rg','momerr',
                       'cmodelmag_dered','psf_fwhm','e1_psf','e2_psf']
            cn = ['field']

            for name in cn_byband:
                cname = name+'_'+self.band
                print("    loading col: ",cname)
                if w is None:
                    self[name] = c[cname][:]
                else:
                    self[name] = c[cname][w]
            for name in cn:
                print("    loading col: ",name)
                if w is None:
                    self[name] = c[name][:]
                else:
                    self[name] = c[name][w]



    def plot_ellip_vs_field(self, field, fmin=None, fmax=None, nbin=20, nperbin=50000):
        self.load_data()

        w=where1(  (self['r'] > 1.0/3.0) 
                 & (self['r'] < 1.0)
                 & (self['whyflag_rg'] == 0)
                 & (self['e1_rg'] < 4)
                 & (self['e1_rg'] > -4)
                 & (self['cmodelmag_dered'] > 18.0)
                 & (self['cmodelmag_dered'] < 21.5) )

        if w.size == 0:
            print("no good objects")
            return

        weights = 1.0/(0.32**2 + self['momerr'][w]**2)

        if field == 'psfsize':
            field_data = self['ixx_psf'] + self['iyy_psf']
            fstr = 'PSF Ixx+Iyy'
        else:
            field_data = self[field][w]
            fstr=field

        fstr = fstr.replace('_','\_')

        be1 = eu.stat.Binner(field_data, self['e1_rg'][w], weights=weights)
        be2 = eu.stat.Binner(field_data, self['e2_rg'][w], weights=weights)

        print("hist e1")
        be1.dohist(nperbin=nperbin, min=fmin, max=fmax)
        #be1.dohist(nbin=nbin, min=fmin, max=fmax)
        print("stats e1")
        be1.calc_stats()
        print("hist e2")
        be2.dohist(nperbin=nperbin, min=fmin, max=fmax)
        #be2.dohist(nbin=nbin, min=fmin, max=fmax)
        print("stats e2")
        be2.calc_stats()

        plt = FramedPlot()
        p1 = Points( be1['wxmean'], be1['wymean'], type='filled circle', color='blue')
        p1err = SymErrY( be1['wxmean'], be1['wymean'], be1['wyerr2'], color='blue')
        p1.label = r'$e_1$'

        p2 = Points( be2['wxmean'], be2['wymean'], type='filled circle', color='red')
        p2.label = r'$e_2$'
        p2err = SymErrY( be2['wxmean'], be2['wymean'], be2['wyerr2'], color='red')

        key = PlotKey(0.1,0.9, [p1,p2])
        plt.add(p1, p1err, p2, p2err, key)


        plt.xlabel = r'$'+fstr+'$'
        plt.ylabel = r'$<e>$'

        plt.show()

       

class Tester:
    def __init__(self, procrun):
        self.procrun=procrun

    def plot_mean_ellip(self, band, run='all', type='rg', data=None):
        """
        Plot the mean ellip as a function of r
        """
        c = open_columns(self.procrun)
        rmin = 1.0/3.0

        type_sup =''
        if type != '':
            type_sup='^{'+type+'}'
            type = '_'+type

        # always use regular r, not r returned after regaussianization
        r_name = 'r_'+band
        e1_name = 'e1'+type+'_'+band
        momerr_name = 'momerr'+type+'_'+band


        #print("Getting good values")
        #w = (c['whyflag_rg_'+band] == 0) & (c['r_rg_'+band].between(0,1)) & (c['e1_rg_'+band].between(-4.0,4.0))
        if data is None:
            print("Reading",r_name)
            r = c[r_name][:]
            print("Reading",e1_name)
            e1 = c[e1_name][:]
            print("Reading",momerr_name)
            momerr = c[momerr_name][:]

            print("Getting good values")
            w=where1( (r >= rmin) & (r <= 1) & (e1 > -4.0) & (e1 < 4.0))

            weights = 1.0/(0.32**2 + momerr[w]**2)

            print("Binning")
            b = eu.stat.Binner(r[w], e1[w], weights=weights)

            nperbin = 100000
            print("    histogram")
            b.dohist(nperbin=nperbin)
            print("    calc stats")
            b.calc_stats()

            data={'r':r, 'e1':e1, 'momerr':momerr, 'binner':b}

        else:
            r = data['r']
            e1 = data['e1']
            momerr = data['momerr']
            b = data['binner']


        plt = FramedPlot()
        p = Points( b['wxmean'], b['wymean'], type='filled circle')
        perr = SymErrY( b['wxmean'], b['wymean'], b['wyerr2'])

        plt.add(p, perr)

        r_label = 'r'+type+'_'+band
        e1_label = 'e1'+type+'_'+band

        plt.xlabel = r'$<r'+type_sup+'_'+band+'>$'
        plt.ylabel = r'$<e1'+type_sup+'_'+band+'>$'

        plt.show()

        return data


class Collator:
    """

    Collate the regauss outputs into columns.

    """
    
    def __init__(self, procrun):
        self._proctype='regauss'
        self._procrun=procrun

    def collate_as_columns_byband(self):
        """
        Collate the collated files from the run into a columns directory.
        Multi-band columns are split to, e.g., flux_r, flux_i, ...
        """

        c = self.open_columns()
        print("Will write in coldir:", c.dir)
        c.create()

        files = self.file_list()
        files.sort()

        # create the meta column from input data
        tmp=eu.io.read(files[0],ext=1,lower=True)
        meta = self.make_meta(tmp)
        c.write_column('meta', meta)

        for f in files:
            print("Reading:",f)
            
            tmp=eu.io.read(f,ext=2,lower=True,ensure_native=True)

            print("    Creating output")
            out_dict = self.create_output(tmp)

            print("    Writing columns")
            for name in out_dict:
                c.write_column(name, out_dict[name])
        self.match_primgal()

        if 'ra' not in c:
            add_radec_to_columns(self.procrun)
        if 'e1_psf_r' not in c:
            self.add_psf_ellip()

    def add_psf_ellip(self):
        c = self.open_columns()
        for band in ['u','g','r','i','z']:
            print("band:",band)
            print("    reading data")
            ixx=c['ixx_psf_'+band][:]
            ixy=c['ixy_psf_'+band][:]
            iyy=c['iyy_psf_'+band][:]
            print("    calculating e1")
            T=ixx+iyy
            e1 = (ixx-iyy)/T
            print("    writing e1")
            c.write_column('e1_psf_'+band, e1)
            print("    calculating e2")
            e2 = 2.0*ixy/T
            print("    writing e2")
            c.write_column('e2_psf_'+band, e2)

    def create_indexes(self):
        c = self.open_columns()

        cnames = ['photoid','run','thing_id',
                  'whyflag_rg_r','whyflag_rg_i',
                  'r_r','r_i',
                  'e1_rg_r','e1_rg_i',
                  'e2_rg_r','e2_rg_i',
                  'cmodelmag_dered_r']

        for n in cnames:
            if not c[n].has_index():
                print("Creating index for:",n)
                c[n].create_index()




    def collate_as_columns(self):
        """
        Collate the collated files from the run into a columns directory.
        """

        c = self.open_columns()
        print("Will write in coldir:", c.dir)
        c.create()

        files = self.file_list()
        files.sort()

        tmp=eu.io.read(files[0],ext=1,lower=True)
        meta = self.make_meta(tmp)
        c.write_column('meta', meta)

        for f in files:
            print("Reading:",f)
            
            tmp=eu.io.read(f,ext=2,lower=True,ensure_native=True)

            output = self.collate_struct(tmp.size) 
            eu.numpy_util.copy_fields(tmp, output)

            cmodelmag_dered = sdsspy.make_cmodelmag(tmp, doerr=False, dered=True)
            modelmag_dered = sdsspy.nmgy2mag(tmp['modelflux']) - tmp['extinction']

            output['cmodelmag_dered'] = cmodelmag_dered
            output['modelmag_dered'] = modelmag_dered
            output['photoid'] = sdsspy.photoid(tmp)

            if f == files[0]:
                eu.numpy_util.ahelp(output)

            c.write_columns(output)
        self.match_primgal()

        if 'ra' not in c:
            add_radec_to_columns(self.procrun)

    def match_primgal(self):
        """

        Find the matching indices, by thing_id, into the sweeps primgal columns
        and save as the column 'primgal_match'

        """

        rg_cols = self.open_columns()
        sweep_cols = es_sdsspy.sweeps.open_columns('gal')

        print("Getting thing_id from regauss")
        thing_id = rg_cols['thing_id'][:]
        print("Matching to the primgal thing_id column")
        mid = sweep_cols['thing_id'].match(thing_id)

        print("Writing match ids")
        rg_cols.write_column('primgal_id', mid)



    def open_columns(self):
        d = self.columns_dir()
        return columns.Columns(d)

    def columns_dir(self):
        p=es_sdsspy.sweeps.Proc(self._proctype,self._procrun) 
        return p.columns_dir()

    def output_dir(self):
        p=es_sdsspy.sweeps.Proc(self._proctype,self._procrun) 
        return p.output_dir()

    def file_list(self):
        indir = self.output_dir()
        pattern='sweepgal-%s-%s-collate-[0-9][0-9][0-9].fits' % (self._proctype,self._procrun)
        print("Searching for pattern:",pattern)
        pattern = os.path.join(indir,pattern)
        files = glob.glob(pattern)

        return files




    def make_meta(self, data):
        import datetime
        date = str( datetime.datetime.now() )
        out = {}
        tags=['proctype','procrun',
              'photoop_v','idlutils_v',
              'photo_sweep','photo_resolve','photo_calib']
        for n in tags:
            out[n] = str(data[n][0])
        out['date'] = date
  
        return out


    
    def collate_struct(self, n=1):
        descr = [('photoid','i8'), # generated
                 ('thing_id','i4'),
                 ('run', 'i2'), 
                 ('rerun', 'S3'), 
                 ('camcol', 'i1'), 
                 ('field', 'i2'), 
                 ('id', 'i2'), 

                 ('cmodelmag_dered','5f4'), # generated
                 ('modelmag_dered','5f4'),  # generated
                 ('extinction','5f4'),

                 ('colc','5f4'),
                 ('rowc','5f4'),
                 ('m_rr_cc','5f4'),

                 ('ixx', '5f8'), 
                 ('iyy', '5f8'), 
                 ('ixy', '5f8'), 
                 ('a4', '5f8'), 
                 ('s2', '5f8'), 
                 ('momerr', '5f8'), 
                 ('whyflag', '5i2'), 

                 ('ixx_psf', '5f8'), 
                 ('iyy_psf', '5f8'), 
                 ('ixy_psf', '5f8'), 
                 ('a4_psf', '5f8'), 
                 ('s2_psf', '5f8'), 
                 ('whyflag_psf', '5i2'), 

                 ('ixx_rg', '5f8'), 
                 ('iyy_rg', '5f8'), 
                 ('ixy_rg', '5f8'), 
                 ('a4_rg', '5f8'), 
                 ('s2_rg', '5f8'), 
                 ('momerr_rg', '5f8'), 
                 ('whyflag_rg', '5i2'), 

                 ('ixx_f0', '5f8'), 
                 ('ixy_f0', '5f8'), 
                 ('iyy_f0', '5f8'), 
                 ('detf0', '5f8'), 

                 ('r', '5f8'), 
                 ('e1', '5f8'), 
                 ('e2', '5f8'), 

                 ('r_rg', '5f8'), 
                 ('e1_rg', '5f8'), 
                 ('e2_rg', '5f8')]

        st = numpy.zeros(n, dtype=descr)
        return st

    def create_output(self, st):
        bands = ['u','g','r','i','z']

        out_dict = {}
        for d in st.dtype.descr:
            name = str( d[0] )

            if len(d) == 3:
                if d[2] == 5:
                    for bandi in xrange(5):
                        fname = name+'_'+bands[bandi]
                        out_dict[fname] = st[name][:, bandi]
            else:
                out_dict[name] = st[name]


        cmodelmag_dered = sdsspy.make_cmodelmag(st, doerr=False, dered=True)
        modelmag_dered = sdsspy.nmgy2mag(st['modelflux']) - st['extinction']

        out_dict['cmodelmag_dered_u'] = cmodelmag_dered[:,0]
        out_dict['cmodelmag_dered_g'] = cmodelmag_dered[:,1]
        out_dict['cmodelmag_dered_r'] = cmodelmag_dered[:,2]
        out_dict['cmodelmag_dered_i'] = cmodelmag_dered[:,3]
        out_dict['cmodelmag_dered_z'] = cmodelmag_dered[:,4]

        out_dict['modelmag_dered_u'] = modelmag_dered[:,0]
        out_dict['modelmag_dered_g'] = modelmag_dered[:,1]
        out_dict['modelmag_dered_r'] = modelmag_dered[:,2]
        out_dict['modelmag_dered_i'] = modelmag_dered[:,3]
        out_dict['modelmag_dered_z'] = modelmag_dered[:,4]

        out_dict['photoid'] = sdsspy.photoid(st)

        return out_dict

    def collate_dtype_gen(self, descr):
        new_descr = [('photoid', 'i8'),
                     ('modelmag_dered_u', 'f4'),
                     ('modelmag_dered_g', 'f4'),
                     ('modelmag_dered_r', 'f4'),
                     ('modelmag_dered_i', 'f4'),
                     ('modelmag_dered_z', 'f4'),
                     ('cmodelmag_dered_u', 'f4'),
                     ('cmodelmag_dered_g', 'f4'),
                     ('cmodelmag_dered_r', 'f4'),
                     ('cmodelmag_dered_i', 'f4'),
                     ('cmodelmag_dered_z', 'f4')]
        for d in descr:
            if d[0] == 'camcol':
                new_descr.append( ('camcol','i1') )
            else:
                if len(d) == 3:
                    if d[2] == 5:
                        for b in ['u','g','r','i','z']:
                            n = d[0]+'_'+b
                            new_descr.append( (n, d[1]) )
        return new_descr



def add_radec_to_columns(procrun):
    gals = es_sdsspy.sweeps.open_columns('gal')
    c = open_columns(procrun)

    if 'ra' in c or 'dec' in c:
        raise ValueError("ra or dec already in columns: "+c.dir)

    print("Reading primgal_id for correspondence")
    galid = c['primgal_id'][:]

    print("Reading ra")
    ra = gals['ra'][galid]
    print("Writing matching ra")
    c.write_column('ra', ra)

    print("Reading dec")
    dec = gals['dec'][galid]
    print("Writing matching dec")
    c.write_column('dec', dec)

def match_yale(procrun):
    dir = '~/lensing/yale'
    yfile=path_join(dir,'sdss_images1_iv.rec')
    outfile = path_join(dir, 'match-yale-%s.rec' % procrun)
    print("Will write to output file:",outfile)

    yale = eu.io.read(yfile, verbose=True)

    urun = numpy.unique(yale['run'])

    max_rmag=22

    c=open_columns(procrun)
    append=False
    for run in urun:
        
        wyale = where1(yale['run'] == run)

        print('-'*70)
        print("Matching run:",run)
        wrun = c['run'] == run

        print("Found",wrun.size,"in run")

        print("Getting whyflag== 0 and r <",max_rmag)
        whyflag = c['whyflag_rg_r'][wrun]
        mag = c['cmodelmag_dered_r'][wrun]

        w=where1( (whyflag == 0) & (mag < max_rmag) )

        w = wrun[w]

        print("Found:",w.size)

        print("Reading ra")
        ra=c['ra'][w]
        print("Reading dec")
        dec=c['dec'][w]

        h = eu.htm.HTM()

        matchrad = 1.0
        print("Matching ra/dec within",matchrad,"arcsec")
        myale, mrg, dist = h.match(yale['ra'][wyale], yale['dec'][wyale], ra, dec, matchrad/1000.)

        if myale.size == 0:
            raise ValueError("no matches!")
        print("Found %d/%d ra/dec matches" % (mrg.size, wyale.size))

        # index back into the original columns
        w = w[mrg]

        columns = ['run','rerun','camcol','field','id',
                   'ra','dec','e1_rg_r','e2_rg_r','momerr_r','r_r',
                   'cmodelmag_dered_r']
        data = c.read_columns(columns, rows=w)

        gc = es_sdsspy.sweeps.open_columns('gal')
        meta = gc['meta'].read()

        hdr={'regauss_procrun':procrun,
             'photo_sweep': meta['photo_sweep'],
             'matchrad':matchrad,
             'max_rmag':max_rmag}

        eu.io.write(outfile, data, header=hdr, append=append, delim=' ')
        append=True
