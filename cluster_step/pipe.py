import os
from sys import stdout, stderr
import pprint
from math import ceil
import numpy
from numpy import where, sqrt, random, zeros, arange, median, random, array

from . import files
from . import prior

import gmix_image
from gmix_image.gmix import GMix, GMixCoellip
from gmix_image.gmix_mcmc import MixMCStandAlone, MixMCPSF
from gmix_image.gmix_em import GMixEMPSF
from gmix_image.gmix_fit import GMixFitSimple

CLUSTERSTEP_GAL=2**0
CLUSTERSTEP_STAR=2**1
CLUSTERSTEP_PSF_STAR=2**2

UW_BAD_ADMOM=2**0
UW_TOO_FEW_PIXELS=2**1
UW_CROWDED=2**2

class Pipe(dict):
    def __init__(self, **keys):
        """
        parameters
        ----------
        run:
            run id
        psfnum:
            psf number
        shnum:
            shear number
        ccd:
            ccd number

        version: optional
            version id

        bugs
        ----

        the psf and shear row,col outputs are in the sub-image coordinates
        Currently need to add row_range[0],col_range[0] to get back to
        image coords

        """

        self._check_keys(**keys)

        for k in keys:
            self[k] = keys[k]

        conf=files.read_config(self['run'])

        for k in conf:
            self[k]=conf[k]
        
        self._set_gpriors()
        self._set_Tpriors()

        random.seed(self['seed'])
        self._load_data()

    def _check_keys(self, **keys):
        if ('run' not in keys
                or 'psfnum' not in keys
                or 'shnum' not in keys
                or 'ccd' not in keys):
            raise ValueError("send run=, psfnum=, "
                             "shnum=, ccd=")

    def _set_Tpriors(self):
        self['Tprior_type'] = self.get('Tprior_type',None)
        Tprior_type=self['Tprior_type']

        if Tprior_type is None:
            self.Tpriors=None
            return

        if Tprior_type == 'fits-vs-mag':
            self._set_Tpriors_vs_mag()
        else:
            raise ValueError("bad Tprior type: '%s'" % Tprior_type)

    def _set_gpriors(self):
        if self['gprior_type'] is None:
            self.gpriors=None
            return

        priors={}
        if self['gprior_type'] == 'old':
            priors['gexp']=gmix_image.priors.GPrior(A=12.25,
                                                    B=0.2,
                                                    C=1.05,
                                                    D=13.)
            priors['gdev'] = priors['gexp']
            self.gpriors=priors
        elif self['gprior_type'] == 'fits-vs-mag-nosplit':
            self.set_priors_vs_mag_exponly()
        elif self['gprior_type'] == 'fits-vs-mag':
            self._set_gpriors_vs_mag()
        else:
            priors['gexp']=prior.GPriorExp(self['gprior_pars_exp'])
            priors['gdev']=prior.GPriorDev(self['gprior_pars_dev'])
            self.gpriors=priors

    
    def _set_gpriors_vs_mag(self):
        """
        Note we use the GPriorExp for both dev and exp, but
        different galaxies were used to train it
        """
        exp_prior_pars=files.read_gprior(type='gexp')
        dev_prior_pars=files.read_gprior(type='gdev')

        gpriors={}

        exp_plist=[]
        for i in xrange(exp_prior_pars.size):
            exp_pdict={}

            pexp=prior.GPriorExp(exp_prior_pars['pars'][i])

            exp_pdict['gprior'] = pexp
            exp_pdict['minmag'] = exp_prior_pars['minmag'][i]
            exp_pdict['maxmag'] = exp_prior_pars['maxmag'][i]

            exp_plist.append(exp_pdict)


        dev_plist=[]
        for i in xrange(dev_prior_pars.size):
            dev_pdict={}

            pdev=prior.GPriorExp(dev_prior_pars['pars'][i])

            dev_pdict['gprior'] = pdev
            dev_pdict['minmag'] = dev_prior_pars['minmag'][i]
            dev_pdict['maxmag'] = dev_prior_pars['maxmag'][i]

            dev_plist.append(dev_pdict)


        gpriors['gexp']=exp_plist
        gpriors['gdev']=dev_plist

        self.gpriors=gpriors


    def _set_Tpriors_vs_mag(self):
        exp_pars=files.read_sprior(type='gexp')
        dev_pars=files.read_sprior(type='gdev')

        Tpriors={}

        exp_plist=[]
        for i in xrange(exp_pars.size):
            pdict={}

            sigma_mean=exp_pars['pars'][i,1]
            sigma_width=exp_pars['pars'][i,2]

            p=prior.TPrior(sigma_mean, sigma_width)

            pdict['Tprior'] = p
            pdict['minmag'] = exp_pars['minmag'][i]
            pdict['maxmag'] = exp_pars['maxmag'][i]

            exp_plist.append(pdict)


        dev_plist=[]
        for i in xrange(dev_pars.size):
            pdict={}

            sigma_mean=dev_pars['pars'][i,1]
            sigma_width=dev_pars['pars'][i,2]

            p=prior.TPrior(sigma_mean, sigma_width)

            pdict['Tprior'] = p
            pdict['minmag'] = dev_pars['minmag'][i]
            pdict['maxmag'] = dev_pars['maxmag'][i]

            dev_plist.append(pdict)


        Tpriors['gexp']=exp_plist
        Tpriors['gdev']=dev_plist

        self.Tpriors=Tpriors



    def set_priors_vs_mag_exponly(self):
        """
        deprecated
        """
        prior_pars=files.read_gprior(nosplit=True)

        gpriors={}
        plist=[]
        for i in xrange(prior_pars.size):
            pdict={}
            p=prior.GPriorExp(prior_pars['pars'][i])
            pdict['gprior'] = p
            pdict['minmag'] = prior_pars['minmag'][i]
            pdict['maxmag'] = prior_pars['maxmag'][i]

            plist.append(pdict)

        gpriors['gexp']=plist
        gpriors['gdev']=plist

        self.gpriors=gpriors

    def get_gprior(self, index, fitmodel):
        if self['gprior_type']==None:
            gprior=None
        elif self['gprior_type'] in ['fits-vs-mag','fits-vs-mag-nosplit']:
            gprior=self.get_gprior_vs_mag(index, fitmodel)
        else:
            gprior=self.gpriors[fitmodel]
        return gprior

    def get_gprior_vs_mag(self, index, fitmodel):
        mag=self.cat['mag_auto_r'][index]
        gprior=None
        for pdict in self.gpriors[fitmodel]:
            if mag >= pdict['minmag'] and mag <= pdict['maxmag']:
                gprior=pdict['gprior']
        if gprior is None:
            if mag < self.gpriors[fitmodel][0]['minmag']:
                gprior=self.gpriors[fitmodel][0]['gprior']
            elif mag > self.gpriors[fitmodel][-1]['maxmag']:
                gprior=self.gpriors[fitmodel][-1]['gprior']
            else:
                raise ValueError("not possible error finding mag: nan?")

        return gprior


    def get_Tprior(self, index, fitmodel):
        if self.Tpriors==None:
            Tprior=None

        mag=self.cat['mag_auto_r'][index]

        Tprior=None
        for pdict in self.Tpriors[fitmodel]:
            if mag >= pdict['minmag'] and mag <= pdict['maxmag']:
                Tprior=pdict['Tprior']

        if Tprior is None:
            if mag < self.Tpriors[fitmodel][0]['minmag']:
                Tprior=self.Tpriors[fitmodel][0]['Tprior']
            elif mag > self.Tpriors[fitmodel][-1]['maxmag']:
                Tprior=self.Tpriors[fitmodel][-1]['Tprior']
            else:
                raise ValueError("not possible error finding mag: nan?")

        return Tprior




    def _load_data(self):
        self.cat=files.read_cat(**self)
        image_orig, self.hdr=files.read_image(**self)

        admom_path=files.get_output_path(ftype='admom', **self)
        if os.path.exists(admom_path):
            self.ares=files.read_fits_output(ftype='admom',**self)

        psf_path=files.get_output_path(ftype='psf', **self)
        if os.path.exists(psf_path):
            psfres=files.read_fits_output(ftype='psf',**self)
            self._set_psfres(psfres)

        shear_path=files.get_output_path(ftype='shear', **self)
        if os.path.exists(shear_path):
            self.shear_res=files.read_fits_output(ftype='shear',**self)

        self.seg, self.seg_hdr=files.read_image(ftype='seg', **self)

        skypersec=self.hdr['sky']
        exptime=self.hdr['exptime']

        # assume gain same for both amplifiers
        gain=self.hdr['gaina']
        read_noise=self.hdr['rdnoisea']

        # note unit mismatch? Both check out
        self['sky'] = skypersec*exptime/gain
        self['skyvar'] = self.hdr['sky'] + read_noise
        self['skysig'] = sqrt(self['skyvar'])
        self['ivar'] = 1.0/self['skyvar']

        self.image = image_orig.astype('f8') - self['sky']
        del image_orig

    def _set_psfres(self, psfres):
        self.psfres=psfres
        self._set_best_psf_model()

    def get_cutout(self, index, size=None):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]
        id=self.cat['id'][index]

        padding=self.get('seg_padding',0)
        include_all_seg=self.get('include_all_seg',True)

        cutout=CutoutWithSeg(self.image, self.seg, cen, id, size,
                             include_all_seg=include_all_seg,
                             padding=padding)

        return cutout

    def get_full_cutout(self, index, size=None):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]

        cutout=Cutout(self.image, cen, size)

        return cutout


    def get_zerod_cutout(self, index, **keys):
        try:
            c=self.get_cutout(index, **keys)
        except NoSegMatches as excpt:
            print >>stderr,str(excpt)
            # sometimes there are no associated
            # segment pixels!
            c=self.get_full_cutout(index, **keys)
            return c

        im=c.subimage
        seg=c.seg_subimage

        id=self.cat['id'][index]
        self.zero_seg(im, seg, id)

        return c

    def zero_seg(self, im, seg, id):
        w=where( (seg != id) & (seg != 0) )
        if w[0].size != 0:
            im[w] = self['skysig']*random.randn(w[0].size)


    def get_stars(self, good=True):
        """
        Get things labeled as stars.  if good, also trim
        to those with good adaptive moments
        """
        logic=self.cat['class']==CLUSTERSTEP_STAR

        if good:
            logic = logic & (self.ares['whyflag']==0)

        w,=where(logic)
        return w

    def get_gals(self, good=True, s2n_min=None):
        """
        Get things labeled as galaxies.  if good, also trim
        to those with good adaptive moments
        """

        logic=self.cat['class']==CLUSTERSTEP_GAL

        if good:
            logic = logic & (self.ares['whyflag']==0)

        if s2n_min is not None:
            logic=logic & (self.ares['s2n'] > s2n_min)

        #logic = logic & (self.ares['s2n'] > 210) & (self.ares['s2n'] < 400)
        #logic = logic & (self.ares['s2n'] > 200) & (self.ares['s2n'] < 210)
        #logic = logic & (self.ares['s2n'] > 100)
        #print "FIX THIS!!!!!!"

        w,=where(logic)
        return w

    def get_psf_stars(self, good=True):
        """
        Get things labeled as stars in the psf star mag range.  if good, also
        trim to those with good adaptive moments
        """
        star_logic = self.cat['class']==CLUSTERSTEP_STAR
        mag_logic  = ( (self.cat['mag_auto_r'] > self['star_minmag'])
                       & (self.cat['mag_auto_r'] < self['star_maxmag']))

        logic = star_logic & mag_logic

        if good:
            logic = logic & (self.ares['whyflag']==0)

        w,=where(logic)
        return w

    def run_shear(self, run_admom=False, run_psf=False):
        """
        Run all things labelled as galaxies with good
        admom measurements through the shear code
        """
        if not hasattr(self,'ares') or run_admom:
            self.run_admom()
            try:
                self.plot_admom_sizemag()
            except:
                # some systems at nersc don't even have X installed
                # and biggles always trys to load it
                print 'failed to make sizemag plots'
        if not hasattr(self,'psfres') or run_psf:
            self.run_psf()

        print 'using psf model:',self._best_psf_model
        ares=self.ares

        wgal=self.get_gals(s2n_min=self['shear_s2n_min'])

        out=self.get_shear_struct(wgal.size)

        out['id'][:] = self.cat['id'][wgal]
        out['simid'][:] = self.cat['simid'][wgal]
        out['psfnum'][:] = self['psfnum']
        out['shnum'][:] = self['shnum']
        out['ccd'][:] = self['ccd']
        out['mag_auto_r'][:] = self.cat['mag_auto_r'][wgal]
        out['row_range'][:] = self.ares['row_range'][wgal]
        out['col_range'][:] = self.ares['col_range'][wgal]

        for igal in xrange(wgal.size):
            from esutil.numpy_util import aprint

            index=wgal[igal]

            print "%s/%s s2n: %s" % (igal+1,wgal.size,ares['s2n'][index])

            c=self.get_zerod_cutout(index)

            im=c.subimage

            ares0=ares[index].copy()
            if self['psf_interp']=='median-render':
                gmix_psf=self.fit_gmix_psf_rendered(im.shape, 
                                                    ares0['wrow']-ares0['row_range'][0],
                                                    ares0['wcol']-ares0['col_range'][0])

            else:
                gmix_psf=self.get_gmix_psf()

            if gmix_psf is None:
                out['flags'][index] = 2**16
                print 'gmix_psf is None'
                continue

            res=self.fit_shear_models(index, im, ares0, gmix_psf)
            self.copy_shear_results(out, res, gmix_psf, igal, ares0)

        self.shear_res=out
        files.write_fits_output(data=out, ftype='shear', **self)


    def fit_shear_models(self, index, im, ares, gmix_psf):
        """
        Fit all listed models, return the best fitting
        """
        aic=9999.e9
        fitmodels=self.get_fitmodels()
        for fitmodel in fitmodels:
            fitter=self.run_shear_model(index, im, ares, gmix_psf, fitmodel)

            res0 = fitter.get_result()
            print '  model: %s prob: %.6f aic: %.6f bic: %.6f Ts/n: %.6f ' % \
                    (fitmodel,res0['fit_prob'],res0['aic'],res0['bic'],res0['Ts2n'])
            if res0['aic'] < aic:
                res=res0
                aic = res0['aic']

        if False:
            import images
            images.multiview(im*0.01,nonlinear=0.4)
            key=raw_input('hit a key: ')
            if key=='q':
                stop


        if len(fitmodels) > 1:
            print '    best model:',res['model']
        return res

    def run_shear_model(self, index, im, ares0, gmix_psf, fitmodel):
        """
        Run the shear model though the mixmc code
        """


        ares={'wrow':ares0['wrow']-ares0['row_range'][0],
              'wcol':ares0['wcol']-ares0['col_range'][0],
              'Irr':ares0['Irr'],
              'Irc':ares0['Irc'],
              'Icc':ares0['Icc'],
              'whyflag':ares0['whyflag']}
 
        cen_width=self.get('cen_width',1.0)
        fitstyle=self.get('fitstyle','mcmc')
        if fitstyle=='lm':
            fitter=self.run_lm(index, im, self['ivar'], gmix_psf, fitmodel,
                               ares,cen_width)

        else:
            nsub=self.get('object_nsub',None)
            gprior=self.get_gprior(index, fitmodel)
            Tprior=self.get_Tprior(index, fitmodel)

            nwalkers=self['nwalkers']
            burnin=self['burnin']
            nstep=self['nstep']
            fitter=MixMCStandAlone(im, self['ivar'],
                                   gmix_psf, gprior, fitmodel,
                                   nwalkers=nwalkers,
                                   nstep=nstep,
                                   burnin=burnin,
                                   mca_a=self['mca_a'],
                                   iter=self.get('iter',False),
                                   draw_gprior=self['draw_gprior'],
                                   ares=ares,
                                   cen_width=cen_width,
                                   nsub=nsub,
                                   Tprior=Tprior,
                                   make_plots=False)
        return fitter

    def run_lm(self, index, im, ivar, gmix_psf, fitmodel, ares, cen_width):
        aic=9.999e9

        # can be None
        gprior=self.get_gprior(index, fitmodel)
        gprior_like=self['gprior_like']
        for i in xrange(4):
            fitter0=GMixFitSimple(im, ivar,
                                  gmix_psf, fitmodel,
                                  ares,cen_width=cen_width,
                                  gprior=gprior,
                                  gprior_like=gprior_like)
            res=fitter0.get_result()
            if res['aic'] < aic:
                aic=res['aic']
                fitter=fitter0
        return fitter

    def get_fitmodels(self):
        fitmodels=self['fitmodel']
        if not isinstance(fitmodels,list):
            fitmodels=[fitmodels]
        return fitmodels

    def get_gmix_psf(self):
        """
        Get a psf model

        ares is not necessarily used
        """
        if self['psf_interp']=='random':
            return self.get_random_gmix_psf()
        elif self['psf_interp'] in ['mean','median']:
            return self.get_mean_gmix_psf()
        else:
            raise ValueError("bad psf interp: '%s'" % self['psf_interp'])

    def fit_gmix_psf_rendered(self, shape, row, col):
        im=self.render_psf(shape, row, col)
        im *= (10000./im.sum())

        skysig=1.0
        im += skysig*(random.randn(im.size).reshape(im.shape))

        ngauss=self.get_model_ngauss(self['psf_render_fit_model'])

        print '  fitting rendered psf'
        ntry=1
        while ntry <= 20:
            res=gmix_image.gmix_fit.quick_fit_psf_coellip(im, skysig, ngauss,
                                                          cen=[row,col])
            if res is None:
                # admom failed
                continue
            if res.flags==0:
                break
            ntry += 1

        if res is None or res.flags != 0:
            print "could not fit psf"
            return None

        return res.get_gmix()


    def render_psf(self, shape, row, col):
        """
        Render the mean model at the indicated position
        """
        import images
        mpsf_gmix=self.get_mean_gmix_psf()
        mpsf_gmix.set_cen(row,col)


        nsub=self['psf_admom_nsub']
        print '  rendering psf at:',row,col,'nsub:',nsub
        im=gmix_image.render.gmix2image(mpsf_gmix, 
                                        shape,
                                        nsub=nsub)
        #print mpsf_gmix
        #print 'im sum:',im.sum()
        #print im.max(),im.min()
        #print shape
        #images.multiview(im)
        #stop
        return im
    def _set_best_psf_model(self):
        if self['psf_models'][0]=='admom':
            self._best_psf_model='admom'
            return

        best_aic=1.e9
        best_model='none'
        for model in self['psf_models']:
            w, = where(self.psfres[model+'_flags'] == 0)
            aic=median(self.psfres[model+'_aic'][w])
            if aic < best_aic:
                best_aic=aic
                best_model=model

        print 'best psf model:',best_model
        self._best_psf_model = best_model

    def get_random_gmix_psf(self):
        """
        Get a random psf model
        """
        model=self._best_psf_model
        wpsf,=where(self.psfres[model+'_flags']==0)

        irand=random.randint(wpsf.size)  
        irand=wpsf[irand]

        pars=self.psfres[model+'_pars'][irand]
        gmix=self._get_model_psf_gmix(pars,model)
        return gmix

    def get_mean_gmix_psf(self):
        """
        Get a the mean of the good psf models
        """

        if not hasattr(self,'_mean_gmix_psf'):
            print 'getting the',self['psf_interp'],'psf model'
            model=self._best_psf_model
            wpsf,=where(self.psfres[model+'_flags']==0)

            pars=self.psfres[model+'_pars'][wpsf[0]]
            allpars=self.psfres[model+'_pars'][wpsf]

            for i in xrange(len(pars)):
                if self['psf_interp'] in ['median','median-render']:
                    pars[i]=median(allpars[:,i])
                else:
                    pars[i]=allpars[:,i].mean()
            gmix=self._get_model_psf_gmix(pars,model)
            self._mean_gmix_psf=gmix
        return self._mean_gmix_psf

    def _get_model_psf_gmix(self, pars, model):
        if model=='gturb':
            # note the MixMCPSF uses e1,e2 not g1,g2
            gmix=GMix(pars, model='gturb')
        elif model in ['gmix1','gmix2','gmix3']:
            gmix=GMixCoellip(pars)
        else:
            # assum pars are full [pi,ri,ci,irri,irci,icci,...]
            gmix=GMix(pars)
        return gmix


    def run_psf(self, run_admom=False):
        """
        Run the PSF measurement on the psf stars

        Only run admom if not already run, or run_admom=True
        """

        #print 'running em on PSF stars with ngauss:',self['ngauss_psf']

        if not hasattr(self,'ares') or run_admom:
            self.run_admom()
            self.plot_admom_sizemag()

        if self['psf_models'][0] == 'admom':
            # just copy the admom results
            self.set_psfres_from_admom()
            self.write_psfres()
            return

        ares=self.ares

        wpsf=self.get_psf_stars()

        out=self.get_psf_struct(wpsf.size)
        out['id'][:] = self.cat['id'][wpsf]
        out['simid'][:] = self.cat['simid'][wpsf]
        out['row_range'][:] = self.ares['row_range'][wpsf]
        out['col_range'][:] = self.ares['col_range'][wpsf]

        for model in self['psf_models']:
            print 'model:',model
            for ipsf in xrange(wpsf.size):
                stdout.write(".")
                #stdout.flush()
                index=wpsf[ipsf]

                # put back in the sub-coordinate system
                # hoops because of endian bug in numpy
                wrow=ares['wrow'][index]-ares['row_range'][index,0]
                wcol=ares['wcol'][index]-ares['col_range'][index,0]
                Irr=ares['Irr'][index]
                Irc=ares['Irc'][index]
                Icc=ares['Icc'][index]
                whyflag=ares['whyflag'][index]
                aresi = {'wrow':wrow,'wcol':wcol,'Irr':Irr,'Irc':Irc,'Icc':Icc,
                         'e1':ares['e1'][index],'e2':ares['e2'][index],
                         'whyflag':whyflag}

                c=self.get_zerod_cutout(index)

                im=c.subimage

                if model in ['em1','em2','em2cocen']:
                    self.run_em_psf(im, aresi, model, out, ipsf)
                elif model in ['gmix1','gmix2','gmix3']:
                    self.run_gmix_fit_psf(im, aresi, model, out, ipsf)
                elif model=='gturb':
                    self.run_turb_psf(im, aresi, out, ipsf)
                else:
                    raise ValueError("bad model: %s" % model)
            print

        best_aic=1.e9
        best_model='none'
        for model in self['psf_models']:
            w, = where(out[model+'_flags'] == 0)
            wbad, = where(out[model+'_flags'] != 0)
            aic=median(out[model+'_aic'][w])
            if aic < best_aic:
                best_aic=aic
                best_model=model
            print 'psf model:',model
            print '  median prob:',median(out[model+'_prob'][w])
            print '  median aic:',aic
            print '  median bic:',median(out[model+'_bic'][w])
            print '  %d/%d bad   %s' % (wbad.size, wpsf.size, wbad.size/float(wpsf.size))

        self._set_psfres(out)
        self.write_psfres()

    def write_psfres(self):
        files.write_fits_output(data=self.psfres,
                                ftype='psf',
                                **self)


    def set_psfres_from_admom(self):
        print 'setting psfres from admom'

        # get psf stars with good admom measurements
        wpsf=self.get_psf_stars()

        psfres=self.get_psf_struct(wpsf.size)
        psfres['simid'][:]        = self.ares['simid'][wpsf]
        psfres['id'][:]           = self.ares['id'][wpsf]
        psfres['row_range'][:]    = self.ares['row_range'][wpsf]
        psfres['col_range'][:]    = self.ares['col_range'][wpsf]

        psfres['admom_s2n'][:]    = self.ares['s2n'][wpsf]
        psfres['admom_flags'][:]  = self.ares['whyflag'][wpsf]
        psfres['admom_pars'][:,0] = 1.0
        psfres['admom_pars'][:,1] = self.ares['wrow'][wpsf]
        psfres['admom_pars'][:,2] = self.ares['wcol'][wpsf]
        psfres['admom_pars'][:,3] = self.ares['Irr'][wpsf]
        psfres['admom_pars'][:,4] = self.ares['Irc'][wpsf]
        psfres['admom_pars'][:,5] = self.ares['Icc'][wpsf]

        self._set_psfres(psfres)

    def get_model_ngauss(self, model):
        if model in ['em1','gmix1']:
            ngauss=1
        elif model in ['em2','gmix2']:
            ngauss=2
        elif model=='em2cocen':
            ngauss=2
        elif model in ['em3','gmix3']:
            ngauss=3
        else:
            raise ValueError("bad psf model: '%s'" % model)

        return ngauss

    def run_em_psf(self, im, aresi, model, out, ipsf):
        ngauss=self.get_model_ngauss(model)
        if 'cocen' in model:
            cocenter=True
        else:
            cocenter=False
        gpsf=GMixEMPSF(im, self['ivar'], ngauss,
                       ares=aresi, 
                       maxiter=self['em_maxiter'], 
                       tol=self['em_tol'],
                       cocenter=cocenter)
        res=gpsf.get_result()
        
        out[model+'_flags'][ipsf] = res['flags']
        out[model+'_numiter'][ipsf] = res['numiter']
        out[model+'_fdiff'][ipsf] = res['fdiff']
        out[model+'_ntry'][ipsf] = res['ntry']

        out[model+'_pars'][ipsf,:] = res['gmix'].get_pars()

        if res['flags']==0:
            stats=gpsf.get_stats()
            out[model+'_prob'][ipsf] = stats['fit_prob']
            out[model+'_aic'][ipsf] = stats['aic']
            out[model+'_bic'][ipsf] = stats['bic']

        if False and res['flags'] != 0:
        #if True:
            print '\nipsf:',ipsf
            pprint.pprint(res)
            pprint.pprint(gpsf.get_stats())
            if True:
                gpsf.compare_model()
                key=raw_input("hit a key (q to quit): ")
                if key=='q':
                    stop

    def run_gmix_fit_psf(self, im, aresi, model, out, ipsf):
        cocenter=False
        if model=='gmix1':
            ngauss=1
        elif model=='gmix2':
            ngauss=2
        elif model=='gmix3':
            ngauss=3
        else:
            raise ValueError("bad psf model: '%s'" % model)

        ntry=1
        while ntry <= 4:
            res=gmix_image.gmix_fit.quick_fit_psf_coellip(im, self['skysig'], 
                                                          ngauss, ares=aresi)
            if res is None:
                # admom failed
                continue
            if res.flags==0:
                break
            ntry += 1

        out[model+'_flags'][ipsf] = res.flags
        out[model+'_numiter'][ipsf] = res.numiter
        out[model+'_ntry'][ipsf] = ntry

        out[model+'_s2n'][ipsf] = res.get_weighted_s2n()

        out[model+'_pars'][ipsf,:] = res.pars
        if res.flags==0:
            out[model+'_perr'][ipsf,:] = res.perr
            out[model+'_pcov'][ipsf,:,:] = res.pcov
            stats=res.get_stats()
            out[model+'_prob'][ipsf] = stats['fit_prob']
            out[model+'_aic'][ipsf] = stats['aic']
            out[model+'_bic'][ipsf] = stats['bic']


    def run_turb_psf(self,im, aresi, out, ipsf):

        fitter=MixMCPSF(im, self['ivar'], 'gturb',
                        nwalkers=self['nwalkers'],
                        nstep=self['nstep'], 
                        burnin=self['burnin'],
                        mca_a=self['mca_a'],
                        iter=self.get('iter',False),
                        ares=aresi)

        res=fitter.get_result()
        out['gturb_pars'][ipsf] = res['pars']
        out['gturb_pcov'][ipsf] = res['pcov']
        out['gturb_prob'][ipsf] = res['fit_prob']
        out['gturb_aic'][ipsf] = res['aic']
        out['gturb_bic'][ipsf] = res['bic']

    def run_admom(self):
        import admom

        print 'running admom',self.cat.size
        ares=self.get_admom_struct()
        ares['simid'][:] = self.cat['simid']
        ares['id'][:] = self.cat['id']

        for i in xrange(self.cat.size):
            #stdout.write(".")
            #stdout.flush()
            if (i % 100) == 0:
                print '%d/%d' % (i+1,self.cat.size)
            c=self.get_zerod_cutout(i)

            im=c.subimage
            cen=c.subcen

            nsub=self.get('psf_admom_nsub',1)
            res = admom.admom(im, cen[0], cen[1],
                              sigsky=self['skysig'],
                              guess=4.0,
                              nsub=nsub)

            if res['whyflag'] != 0:
                #print res['whyflag'],'%.16g %.16g' % (cen[0]-res['wrow'], cen[1]-res['wcol'])
                print 'index: %4d flag: %s' % (i,res['whystr'])

            ares['row_range'][i] = c.row_range
            ares['col_range'][i] = c.col_range

            ares['wrow'][i] = c.row_range[0] + res['wrow']
            ares['wcol'][i] = c.col_range[0] + res['wcol']
            
            #for n in ares.dtype.names:
            for n in res:
                if n in ares.dtype.names and n not in ['row','col','wrow','wcol']:
                    ares[n][i] = res[n]
        print

        w, = where(ares['whyflag'] != 0)
        print 'found %d/%d bad   %s' % (w.size, ares.size, w.size/float(ares.size))
        self.ares=ares

        files.write_fits_output(data=ares,
                                ftype='admom',
                                **self)


    def plot_admom_sizemag(self, show=False):
        """
        Plot the size-mag.  Show psf stars in a different color.
        """
        import biggles

        if not hasattr(self, 'ares'):
            raise ValueError("run_admom first")

        wstar=self.get_stars()
        wgal=self.get_gals()
        wpsf=self.get_psf_stars()


        plt=biggles.FramedPlot()
        sigma=sqrt( (self.ares['Irr']+self.ares['Icc'])/2 )
        mag=self.cat['mag_auto_r']
        psize=0.7
        pstar=biggles.Points(mag[wstar], sigma[wstar], 
                             type='filled circle', color='red',
                             size=psize)
        pgal=biggles.Points(mag[wgal], sigma[wgal], 
                            type='filled diamond', color='dark green',
                            size=psize)

        ppsf=biggles.Points(mag[wpsf], sigma[wpsf], 
                            type='filled circle', color='blue',
                            size=psize)

        pstar.label='star'
        pgal.label='gal'
        ppsf.label='psf star'

        key=biggles.PlotKey(0.95,0.92, [pstar,ppsf,pgal], halign='right')

        plt.add(pstar, ppsf, pgal, key)
        plt.xlabel='mag_auto_r'
        plt.ylabel=r'$\sigma_{AM}$ [pixels]'
        label='%s-p%s-s%s-%s' % (self['run'],self['psfnum'],self['shnum'],
                                 self['ccd'])
        plt.add(biggles.PlotLabel(0.05,0.92,label,halign='left'))

        plt.xrange=[16.5,25.5]
        plt.yrange=[0.75,8]
        if show:
            plt.show()

        epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)
        self._write_plot(plt, epsfile)



    def _write_plot(self, plt, epsfile):
        import converter
        d=os.path.dirname(epsfile)
        if not os.path.exists(d):
            try:
                os.makedirs(d)
            except:
                pass
        print 'writing eps file:',epsfile
        plt.write_eps(epsfile)

        converter.convert(epsfile, dpi=100)



    def show_many_cutouts(self, indices, **keys):
        for i in indices:
            self.show_cutout(i, **keys)
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                return


    def show_cutout(self, index, size=None, 
                    zero_seg=False):
        import biggles
        biggles.configure( 'default', 'fontsize_min', 2)
        if not hasattr(self,'ares'):
            self.run_admom()

        import biggles
        import images

        c=self.get_cutout(index, size=size)
        cz=self.get_zerod_cutout(index, size=size)

        minv=-2.5*self['skysig']

        subimage=c.get_subimage().copy()
        segsub=c.get_seg_subimage().copy()

        zsubimage=cz.get_subimage().copy()

        objtype='star'
        if self.cat['class'][index]==CLUSTERSTEP_GAL:
            objtype='gal'

        s2n=self.ares['s2n'][index]
        print 'index:',index
        print 's2n:',s2n
        print 'ranges:'
        print '\t',c.row_range
        print '\t',c.col_range
        print 'image shape:',subimage.shape
        print

        # make the other objects always show darker
        id=self.cat['id'][index]
        w=where((segsub != id) & (segsub != 0))
        if w[0].size > 0:
            segsub[w] = id*10

        implt=images.view(subimage,show=False, min=minv)
        segplt=images.view(segsub,show=False)

        title='%d %s S/N: %.1f [%d,%d]'
        title = title % (index, objtype,s2n,subimage.shape[0],subimage.shape[1])
        implt.title=title
        segplt.title='seg'

        zplt=images.view(zsubimage,show=False, min=minv)
        zplt.title='zerod'

        tab=biggles.Table(2,2)
        tab[0,0]=implt
        tab[0,1]=segplt
        tab[1,0]=zplt

        tab.show()

    def get_admom_struct(self):
        sxdt=self.cat.dtype.descr
        dt= [('wrow','f8'), # note simid,id in sxdt
             ('wcol','f8'),
             ('row_range','f8',2),
             ('col_range','f8',2),
             ('Irr','f8'),
             ('Irc','f8'),
             ('Icc','f8'),
             ('e1','f8'),
             ('e2','f8'),
             ('rho4','f8'),
             ('a4','f8'),
             ('s2','f8'),
             ('uncer','f8'),
             ('s2n','f8'),
             ('numiter','i4'),
             ('nsub','i4'),
             ('whyflag','i4'),
             ('whystr','S10'),
             ('shiftmax','f8')]

        dt=sxdt + dt
        data=zeros(self.cat.size, dtype=dt)
        for n in self.cat.dtype.names:
            data[n][:] = self.cat[n][:]

        return data

    def get_psf_struct(self, n):
        model_npars={'em1':6,'em2':2*6,'em2cocen':2*6,
                     'gmix1':6, 'gmix2':2*2+4,'gmix3':2*3+4,
                     'gturb':6,'admom':6}
        dt= [('simid','i4'),
             ('id','i4'),
             ('row_range','f8',2),
             ('col_range','f8',2)]

        for model in self['psf_models']:
            npars=model_npars[model]
            if model=='gturb':
                dt+=[('gturb_flags','i4'), # will always be zero
                     ('gturb_pars','f8',npars),
                     ('gturb_pcov','f8',(npars,npars)),
                     ('gturb_prob','f8'),
                     ('gturb_aic','f8'),
                     ('gturb_bic','f8')]

            elif model=='admom':
                dt+=[('admom_s2n','f8'),
                     ('admom_flags','i4'),
                     ('admom_pars','f8',npars)]
            else:
                for dti in [('flags','i4'),('numiter','i4'),
                            ('fdiff','f8'),('ntry','i4'),
                            ('pars','f8',npars),
                            ('prob','f8'),('aic','f8'),('bic','f8')]:
                    if len(dti) > 2:
                        dt += [(model+'_'+dti[0], dti[1], dti[2])]
                    else:
                        dt += [(model+'_'+dti[0], dti[1])]
                if 'em' in model:
                    dt += [(model+'_fdiff','f8')]
                if 'gmix' in model:
                    dt += [(model+'_s2n','f8'),
                           (model+'_perr','f8',npars),
                           (model+'_pcov','f8',(npars,npars))]

        data=zeros(n, dtype=dt)

        for model in self['psf_models']:
            if model in ['em1','em2','em2cocen']:
                data[model+'_fdiff']=9999.
                data[model+'_ntry']=9999
                data[model+'_pars']=-9999.
                data[model+'_aic'] = 1.e9
                data[model+'_bic'] = 1.e9
            elif 'gmix' in model:
                data[model+'_ntry']=9999
                data[model+'_pars']=-9999.
                data[model+'_perr']=9999.
                data[model+'_pcov']=9999.
                data[model+'_aic'] = 1.e9
                data[model+'_bic'] = 1.e9

            elif model=='gturb':
                data['gturb_aic'] = 1.e9
                data['gturb_bic'] = 1.e9

        return data

    def copy_shear_results(self, out, res, gmix_psf, igal, ares0):
        """
        Copy results into existing "out" structure
        """
        out['s2n_admom'][igal] = ares0['s2n']

        e1psf,e2psf,Tpsf=gmix_psf.get_e1e2T()
        Tobj=res['Tmean']
        s2 = Tpsf/Tobj
        sratio=sqrt(Tobj/Tpsf)

        out['Tpsf'][igal] = Tpsf
        out['e1psf'][igal] = e1psf
        out['e2psf'][igal] = e2psf
        out['s2'][igal] = s2
        out['sratio'][igal] = sratio

        #out['pars_psf'][igal] = gmix_psf.get_pars()

        for k in res:
            if k in out.dtype.names:
                out[k][igal] = res[k]


    def get_shear_struct(self, n):
        npars=6
        dt=[('simid','i4'),
            ('id','i4'),
            ('psfnum','i2'),
            ('shnum','i2'),
            ('ccd','i2'),
            ('mag_auto_r','f8'),
            ('row_range','f8',2),
            ('col_range','f8',2),
            ('model','S5'),
            ('flags','i4'),
            ('e1psf','f8'),
            ('e2psf','f8'),
            ('Tpsf','f8'),
            ('s2','f8'),
            ('sratio','f8'),
            ('g','f8',2),
            ('gsens','f8',2),
            ('gcov','f8',(2,2)),
            ('pars','f8',npars),
            ('pcov','f8',(npars,npars)),
            ('Tmean','f8','f8'),
            ('Terr','f8','f8'),
            ('Ts2n','f8','f8'),
            ('s2n_w','f8'),  # weighted s/n based on most likely point
            ('s2n_admom','f8'),
            ('arate','f8'),
            ('loglike','f8'),     # loglike of fit
            ('chi2per','f8'),     # chi^2/degree of freedom
            ('dof','i4'),         # degrees of freedom
            ('fit_prob','f8'),    # probability of the fit happening randomly
            ('aic','f8'),         # Akaike Information Criterion
            ('bic','f8')          # Bayesian Information Criterion
           ]

        data=zeros(n, dtype=dt)
        data['gcov']=9999
        data['pcov']=9999
        return data



class FFTStackPipe(dict):
    def __init__(self, **keys):
        """
        works similarly to the stack and shift method

        parameters
        ----------
        run:
            run id
        psfnum:
            psf number
        shnum:
            shear number
        ccd:
            ccd number

        version: optional
            version id

        ideas:
            shifting - non-linear make sense?  Different altogether?  
            
            Fourier space so centroid doesn't matter (no shifting)?

            does boosting really help?
        """

        self._check_keys(**keys)

        for k in keys:
            self[k] = keys[k]

        self._set_boost(**keys)

        conf=files.read_config(self['run'])

        for k in conf:
            self[k]=conf[k]

        self._load_data()

    def _check_keys(self, **keys):
        if ('run' not in keys
                or 'psfnum' not in keys
                or 'shnum' not in keys
                or 'ccd' not in keys):
            raise ValueError("send run=, psfnum=, "
                             "shnum=, ccd=")

    def _set_boost(self, **keys):
        boost=int(keys.get('stack_boost',1))
        if boost < 1:
            raise ValueError("stack_boost must be >= 1")
        self['stack_boost']=boost

    def _load_data(self):
        self.cat=files.read_cat(**self)
        image_orig, self.hdr=files.read_image(**self)

        shear_path=files.get_output_path(ftype='shear', **self)
        if os.path.exists(shear_path):
            self.shear_res=files.read_fits_output(ftype='shear',**self)

        self.seg, self.seg_hdr=files.read_image(ftype='seg', **self)

        skypersec=self.hdr['sky']
        exptime=self.hdr['exptime']

        # assume gain same for both amplifiers
        gain=self.hdr['gaina']
        read_noise=self.hdr['rdnoisea']

        # note unit mismatch? Both check out
        self['sky'] = skypersec*exptime/gain
        self['skyvar'] = self.hdr['sky'] + read_noise
        self['skysig'] = sqrt(self['skyvar'])
        self['ivar'] = 1.0/self['skyvar']

        self.image = image_orig.astype('f8') - self['sky']
        del image_orig

    def get_cutout(self, index, size=None):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]
        id=self.cat['id'][index]

        padding=self.get('seg_padding',0)
        include_all_seg=self.get('include_all_seg',True)

        cutout=CutoutWithSeg(self.image, self.seg, cen, id, size,
                             include_all_seg=include_all_seg,
                             padding=padding)

        return cutout

    def get_full_cutout(self, index, size=None):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]

        cutout=Cutout(self.image, cen, size)

        return cutout


    def get_zerod_cutout(self, index, **keys):
        try:
            c=self.get_cutout(index, **keys)
        except NoSegMatches as excpt:
            print >>stderr,str(excpt)
            # sometimes there are no associated
            # segment pixels!
            c=self.get_full_cutout(index, **keys)
            return c

        im=c.subimage
        seg=c.seg_subimage

        id=self.cat['id'][index]
        self.zero_seg(im, seg, id)

        return c

    def zero_seg(self, im, seg, id):
        w=where( (seg != id) & (seg != 0) )
        if w[0].size != 0:
            im[w] = self['skysig']*random.randn(w[0].size)


    def get_stars(self):
        """
        Get things labeled as stars.  if good, also trim
        to those with good adaptive moments
        """
        logic=self.cat['class']==CLUSTERSTEP_STAR

        am_logic = self.res['am_flags']==0
        logic = logic & am_logic

        w,=where(logic)
        return w

    def get_gals(self):
        """
        Get things labeled as galaxies.  if good, also trim
        to those with good adaptive moments
        """

        logic=self.cat['class']==CLUSTERSTEP_GAL

        am_logic = self.res['am_flags']==0
        logic = logic & am_logic

        w,=where(logic)
        return w

    def get_psf_stars(self):
        """
        Get things labeled as stars in the psf star mag range.  if good, also
        trim to those with good adaptive moments
        """
        star_logic = self.cat['class']==CLUSTERSTEP_STAR
        mag_logic  = ( (self.cat['mag_auto_r'] > self['star_minmag'])
                       & (self.cat['mag_auto_r'] < self['star_maxmag']))
        am_logic = self.res['am_flags']==0
        logic = star_logic & mag_logic & am_logic

        w,=where(logic)
        return w

    def _get_s2n_bins(self):
        s2n_bins=numpy.logspace(numpy.log10(self['stack_s2n_min']),
                                numpy.log10(self['stack_s2n_max']),
                                self['stack_nbin']+1)
        return s2n_bins

    def run(self):
        self._run_admom()

        self.stacks=self.get_stack_struct()

        print 'stacking psf images'
        self.psf_stack, self.npsf_stack, self.psf_skyvar\
                = self._stack_psf_stars()

        self._measure_stacked_psf_admom()

        print 'stacking gal images'

        s2n_bins=self._get_s2n_bins()
        nbin=len(s2n_bins)-1
        for i in xrange(nbin):
            s2n_min=s2n_bins[i]
            s2n_max=s2n_bins[i+1]
            print 's2n: [%s,%s]' % (s2n_min,s2n_max)
            gal_stack,nstack,skyvar=self._stack_galaxies(s2n_min, s2n_max)

            self.stacks['s2n_min'][i] = s2n_min
            self.stacks['s2n_max'][i] = s2n_max
            self.stacks['nstack'][i] = nstack
            self.stacks['skyvar'][i] = skyvar
            self.stacks['images_real'][i,:,:] = gal_stack.real
            self.stacks['images_imag'][i,:,:] = gal_stack.imag

        self._write_data()

    def _write_data(self):
        import fitsio
        path=files.get_output_path(ftype='shear', **self)
        print path

        dir=os.path.dirname(path)
        if not os.path.exists(dir):
            try:
                os.makedirs(dir)
            except:
                pass

        psf_header={'skyvar':self.psf_skyvar,
                    'nstack':self.npsf_stack,
                    'boost':self['stack_boost']}
        with fitsio.FITS(path, mode='rw', clobber=True) as fits:
            fits.write(self.res)
            fits.write(self.psf_stack.real, header=psf_header)
            fits.write(self.psf_stack.imag, header=psf_header)
            fits.write(self.stacks)
    
    def _measure_stacked_psf_admom(self):
        import admom

        psf=numpy.fft.fft(self.psf_stack).real

        cen=array(psf.shape)/2.
        guess=4.0*self['stack_boost']**2
        ares = admom.admom(psf, cen[0], cen[1],
                           sigsky=sqrt(self.psf_skyvar),
                           guess=guess,
                           nsub=1)
        self.psf_ares=ares


    def _stack_psf_stars(self):
        wpsf=self.get_psf_stars()
        imstack, nstack, skyvar = self._stack_images(wpsf)

        return imstack,nstack,skyvar

    def _stack_galaxies(self, s2n_min, s2n_max):
        wgal=self.get_gals()
        logic = ( (self.res['am_s2n'][wgal] > s2n_min) 
                 &(self.res['am_s2n'][wgal] < s2n_max) )

        T = self.res['am_irr'][wgal] + self.res['am_icc'][wgal]
        Tpsf = self.psf_ares['Irr'] + self.psf_ares['Icc']

        # Tpsf is in the boosted image
        Tpsf *= (1./self['stack_boost']**2)

        srat=float(self['stack_sratio'])
        logic = logic & (T > Tpsf* srat**2)

        wgal2=where(logic)
        wgal=wgal[wgal2]

        imstack,nstack,skyvar=self._stack_images(wgal)

        return imstack, nstack, skyvar

    def _stack_images(self, win):
        import images
        res=self.res
        nrow=1+res['row_range'][win,1]-res['row_range'][win,0]
        ncol=1+res['col_range'][win,1]-res['col_range'][win,0]

        w,=where((nrow==self['cutout_size']) & (ncol==self['cutout_size']))
        print '  using %s/%s in stack' % (w.size,win.size)

        w=win[w]
        
        size=self['cutout_size']*self['stack_boost']
        imstack=zeros( (size,size), dtype=numpy.complex64)
        
        skyvar=0.0

        for i in xrange(w.size):
            index=w[i]
            c=self.get_zerod_cutout(index)

            im=c.subimage
    
            if self['stack_boost'] != 1:
                im = images.boost(im, self['stack_boost'])

            if im.shape != imstack.shape:
                raise ValueError("shape incorrect; index error?")

            imstack += numpy.fft.fft(im)

            # sky variances add
            skyvar += self['skyvar']

        return imstack, w.size, skyvar

    def _run_admom(self, **keys):

        res=self.get_struct()
        self.res=res

        for i in xrange(self.cat.size):
            if (i % 100) == 0:
                print '%d/%d' % (i+1,self.cat.size)
            c=self.get_zerod_cutout(i)

            im=c.subimage
            cen=c.subcen
            self.res['row_range'][i] = c.row_range
            self.res['col_range'][i] = c.col_range

            ares=self._do_admom(im, i, cen)

            self.res['am_row'][i] = c.row_range[0] + ares['wrow']
            self.res['am_col'][i] = c.col_range[0] + ares['wcol']

        w, = where(self.res['am_flags'] != 0)
        print 'found %d/%d bad admom  %s' % (w.size,res.size,w.size/float(res.size))


    def _do_admom(self, im, i, cen):
        import admom
        ares = admom.admom(im, cen[0], cen[1],
                           sigsky=self['skysig'],
                           guess=4.0,
                           nsub=1)

        if ares['whyflag'] != 0:
            print 'index: %4d flag: %s' % (i,ares['whystr'])
        
        for n in ares:
            rn='am_'+n.lower()
            if (rn in self.res.dtype.names 
                    and n not in ['row','col','wrow','wcol']):
                self.res[rn][i] = ares[n]

        # two are differently named
        self.res['am_flags'][i] = ares['whyflag']
        self.res['am_flagstr'][i] = ares['whystr']

        return ares


    def get_stack_struct(self):
        size=self['cutout_size']*self['stack_boost']
        imshape=tuple([size]*2)
        nbin=self['stack_nbin']
        dt=[('s2n_min','f8'),
            ('s2n_max','f8'),
            ('nstack','i4'),
            ('skyvar','f8'),
            ('images_real', 'f8', imshape),
            ('images_imag', 'f8', imshape)]

        data=zeros(nbin, dtype=dt)
        return data

    def get_struct(self):
        sxdt=self.cat.dtype.descr

        dt= [('psfnum','i2'),
             ('shnum','i2'),
             ('ccd','i2'),
             ('row_range','f8',2),
             ('col_range','f8',2),
             ('am_row','f8'), # note simid,id in sxdt
             ('am_col','f8'),
             ('am_irr','f8'),
             ('am_irc','f8'),
             ('am_icc','f8'),
             ('am_e1','f8'),
             ('am_e2','f8'),
             ('am_rho4','f8'),
             ('am_a4','f8'),
             ('am_s2','f8'),
             ('am_uncer','f8'),
             ('am_s2n','f8'),
             ('am_numiter','i4'),
             ('am_flags','i4'),
             ('am_flagstr','S10'),
             ('am_shiftmax','f8')]
        dt=sxdt + dt
        data=zeros(self.cat.size, dtype=dt)
        for n in self.cat.dtype.names:
            data[n][:] = self.cat[n][:]

        data['psfnum'] = self['psfnum']
        data['shnum'] = self['shnum']
        data['ccd'] = self['ccd']
        return data



class StackPipe(dict):
    def __init__(self, **keys):
        """
        parameters
        ----------
        run:
            run id
        psfnum:
            psf number
        shnum:
            shear number
        ccd:
            ccd number

        version: optional
            version id

        """


        for k in keys:
            self[k] = keys[k]

        self._set_boost(**keys)

        conf=files.read_config(self['run'])

        for k in conf:
            self[k]=conf[k]

        self._check_inputs(**keys)

        self._load_data()

    def _check_inputs(self, **keys):
        if ('run' not in keys
                or 'psfnum' not in keys
                or 'shnum' not in keys
                or 'ccd' not in keys):
            raise ValueError("send run=, psfnum=, "
                             "shnum=, ccd=")


    def _set_boost(self, **keys):
        boost=int(keys.get('stack_boost',1))
        if boost < 1:
            raise ValueError("stack_boost must be >= 1")
        self['stack_boost']=boost

    def _load_data(self):
        self.cat=files.read_cat(**self)
        image_orig, self.hdr=files.read_image(**self)

        shear_path=files.get_output_path(ftype='shear', **self)
        if os.path.exists(shear_path):
            self.shear_res=files.read_fits_output(ftype='shear',**self)

        if 'mixmc_run' in self:
            self.mixmc_pipe=Pipe(run=self['mixmc_run'],
                                 psfnum=self['psfnum'], 
                                 shnum=self['shnum'], 
                                 ccd=self['ccd'])

        self.seg, self.seg_hdr=files.read_image(ftype='seg', **self)

        skypersec=self.hdr['sky']
        exptime=self.hdr['exptime']

        # assume gain same for both amplifiers
        gain=self.hdr['gaina']
        read_noise=self.hdr['rdnoisea']

        # note unit mismatch? Both check out
        self['sky'] = skypersec*exptime/gain
        self['skyvar'] = self.hdr['sky'] + read_noise
        self['skysig'] = sqrt(self['skyvar'])
        self['ivar'] = 1.0/self['skyvar']

        self.image = image_orig.astype('f8') - self['sky']
        del image_orig

    def get_cutout(self, index, size=None):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]
        id=self.cat['id'][index]

        padding=self.get('seg_padding',0)
        include_all_seg=self.get('include_all_seg',True)

        cutout=CutoutWithSeg(self.image, self.seg, cen, id, size,
                             include_all_seg=include_all_seg,
                             padding=padding)

        return cutout

    def get_full_cutout(self, index, size=None):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]

        cutout=Cutout(self.image, cen, size)

        return cutout


    def get_zerod_cutout(self, index, **keys):
        try:
            c=self.get_cutout(index, **keys)
        except NoSegMatches as excpt:
            print >>stderr,str(excpt)
            # sometimes there are no associated
            # segment pixels!
            c=self.get_full_cutout(index, **keys)
            return c

        im=c.subimage
        seg=c.seg_subimage

        id=self.cat['id'][index]
        self.zero_seg(im, seg, id)

        return c

    def zero_seg(self, im, seg, id):
        w=where( (seg != id) & (seg != 0) )
        if w[0].size != 0:
            im[w] = self['skysig']*random.randn(w[0].size)


    def get_stars(self):
        """
        Get things labeled as stars.  if good, also trim
        to those with good adaptive moments
        """
        logic=self.cat['class']==CLUSTERSTEP_STAR

        am_logic = self.res['am_flags']==0
        logic = logic & am_logic

        w,=where(logic)
        return w

    def get_gals(self):
        """
        Get things labeled as galaxies.  if good, also trim
        to those with good adaptive moments
        """

        logic=self.cat['class']==CLUSTERSTEP_GAL

        am_logic = self.res['am_flags']==0
        logic = logic & am_logic

        w,=where(logic)
        return w

    def get_psf_stars(self):
        """
        Get things labeled as stars in the psf star mag range.  if good, also
        trim to those with good adaptive moments
        """
        star_logic = self.cat['class']==CLUSTERSTEP_STAR
        mag_logic  = ( (self.cat['mag_auto_r'] > self['star_minmag'])
                       & (self.cat['mag_auto_r'] < self['star_maxmag']))
        am_logic = self.res['am_flags']==0
        logic = star_logic & mag_logic & am_logic

        w,=where(logic)
        return w

    def _get_s2n_bins(self):
        s2n_bins=numpy.logspace(numpy.log10(self['stack_s2n_min']),
                                numpy.log10(self['stack_s2n_max']),
                                self['stack_nbin']+1)
        return s2n_bins

    def run(self):
        self._run_admom()

        self.stacks=self.get_stack_struct()

        print 'stacking psf images'
        self.psf_stack, self.npsf_stack, self.psf_skyvar\
                = self._stack_psf_stars()

        self._measure_stacked_psf_admom()

        print 'stacking gal images'

        s2n_bins=self._get_s2n_bins()
        nbin=len(s2n_bins)-1
        for i in xrange(nbin):
            s2n_min=s2n_bins[i]
            s2n_max=s2n_bins[i+1]
            print 's2n: [%s,%s]' % (s2n_min,s2n_max)
            gal_stack,nstack,skyvar=self._stack_galaxies(s2n_min, s2n_max)

            self.stacks['s2n_min'][i] = s2n_min
            self.stacks['s2n_max'][i] = s2n_max
            self.stacks['nstack'][i] = nstack
            self.stacks['skyvar'][i] = skyvar
            self.stacks['images'][i,:,:] = gal_stack

        self._write_data()

    def _write_data(self):
        import fitsio
        path=files.get_output_path(ftype='shear', **self)
        print path

        dir=os.path.dirname(path)
        if not os.path.exists(dir):
            try:
                os.makedirs(dir)
            except:
                pass

        psf_header={'skyvar':self.psf_skyvar,
                    'nstack':self.npsf_stack,
                    'boost':self['stack_boost']}
        with fitsio.FITS(path, mode='rw', clobber=True) as fits:
            fits.write(self.res)
            fits.write(self.psf_stack, header=psf_header)
            fits.write(self.stacks)
    
    def _measure_stacked_psf_admom(self):
        import admom
        cen=array(self.psf_stack.shape)/2.
        guess=4.0*self['stack_boost']**2
        ares = admom.admom(self.psf_stack, cen[0], cen[1],
                           sigsky=sqrt(self.psf_skyvar),
                           guess=guess,
                           nsub=1)
        self.psf_ares=ares


    def _stack_psf_stars(self):
        wpsf=self.get_psf_stars()
        imstack, nstack, skyvar = self._stack_images(wpsf)

        if False:
            import images
            images.multiview(imstack/imstack.max(), nonlinear=1)

        return imstack,nstack,skyvar

    def _stack_galaxies(self, s2n_min, s2n_max):
        # adaptive good enough for psf since we know
        # it is gaussian
        Tpsf = self.psf_ares['Irr'] + self.psf_ares['Icc']
        if 'mixmc_run' in self:
            mpipe=self.mixmc_pipe
            wgal=mpipe.get_gals(s2n_min=mpipe['shear_s2n_min'])
            mres=mpipe.shear_res
            T = mres['Tmean']
            logic = (  (self.res['am_flags'][wgal]==0)
                     & (self.res['am_s2n'][wgal] > s2n_min) 
                     & (self.res['am_s2n'][wgal] < s2n_max) )
        else:
            wgal=self.get_gals()
            logic = ( (self.res['am_s2n'][wgal] > s2n_min) 
                     &(self.res['am_s2n'][wgal] < s2n_max) )

            T = self.res['am_irr'][wgal] + self.res['am_icc'][wgal]

        # Tpsf is in the boosted image
        Tpsf *= (1./self['stack_boost']**2)

        srat=float(self['stack_sratio'])
        logic = logic & (T > Tpsf* srat**2)

        wgal2,=where(logic)
        wgal=wgal[wgal2]

        imstack,nstack,skyvar=self._stack_images(wgal)

        if False:
            import images
            tit='s2n: [%.2f,%.2f]' % (s2n_min,s2n_max)
            images.multiview(imstack/imstack.max(), nonlinear=1,
                            title=tit)

        return imstack, nstack, skyvar

    def _stack_images(self, win):
        import images
        res=self.res
        nrow=1+res['row_range'][win,1]-res['row_range'][win,0]
        ncol=1+res['col_range'][win,1]-res['col_range'][win,0]

        w,=where((nrow==self['cutout_size']) & (ncol==self['cutout_size']))
        print '  using %s/%s in stack' % (w.size,win.size)

        w=win[w]
        
        row=res['am_row'][w]-res['row_range'][w,0]
        col=res['am_col'][w]-res['col_range'][w,0]

        row *= self['stack_boost']
        col *= self['stack_boost']

        row_mean=row.mean()
        col_mean=col.mean()


        print '    <row>:',row_mean,'<col>:',col_mean

        size=self['cutout_size']*self['stack_boost']
        imstack=zeros( (size,size) )
        
        skyvar=0.0

        use_seg=self.get('use_seg',True)
        for i in xrange(w.size):
            index=w[i]
            if use_seg:
                c=self.get_zerod_cutout(index)
            else:
                c=self.get_full_cutout(index)

            im=c.subimage
    
            if self['stack_boost'] != 1:
                im = images.boost(im, self['stack_boost'])

            if im.shape != imstack.shape:
                raise ValueError("shape incorrect; index error?")

            # rows are in boosted frame
            drow=row_mean-row[i]
            dcol=col_mean-col[i]

            imshift= self._shift_image(im, [drow, dcol])
            """
            guess=(self.res['am_irr'][index]+self.res['am_icc'][index])/2
            ares = admom.admom(imshift, row_mean, col_mean,
                               sigsky=self['skysig'],
                               guess=guess,
                               nsub=1)

            print row_mean-ares['wrow']
            print col_mean-ares['wcol']
            """

            imstack += imshift

            # sky variances add
            skyvar += self['skyvar']

        return imstack, w.size, skyvar

    def _shift_image(self, image, shift):
        import scipy.ndimage

        output = scipy.ndimage.interpolation.shift(image, shift, 
                                                   output='f8',
                                                   order=1,
                                                   mode='constant',
                                                   cval=0.0)
        return output

    def _run_admom(self, **keys):

        res=self.get_struct()
        self.res=res

        for i in xrange(self.cat.size):
            if (i % 100) == 0:
                print '%d/%d' % (i+1,self.cat.size)
            c=self.get_zerod_cutout(i)

            im=c.subimage
            cen=c.subcen
            self.res['row_range'][i] = c.row_range
            self.res['col_range'][i] = c.col_range

            ares=self._do_admom(im, i, cen)

            self.res['am_row'][i] = c.row_range[0] + ares['wrow']
            self.res['am_col'][i] = c.col_range[0] + ares['wcol']

        w, = where(self.res['am_flags'] != 0)
        print 'found %d/%d bad admom  %s' % (w.size,res.size,w.size/float(res.size))


    def _do_admom(self, im, i, cen):
        import admom
        ares = admom.admom(im, cen[0], cen[1],
                           sigsky=self['skysig'],
                           guess=4.0,
                           nsub=1)

        if ares['whyflag'] != 0:
            print 'index: %4d flag: %s' % (i,ares['whystr'])
        
        for n in ares:
            rn='am_'+n.lower()
            if (rn in self.res.dtype.names 
                    and n not in ['row','col','wrow','wcol']):
                self.res[rn][i] = ares[n]

        # two are differently named
        self.res['am_flags'][i] = ares['whyflag']
        self.res['am_flagstr'][i] = ares['whystr']

        return ares


    def get_stack_struct(self):
        size=self['cutout_size']*self['stack_boost']
        imshape=tuple([size]*2)
        nbin=self['stack_nbin']
        dt=[('s2n_min','f8'),
            ('s2n_max','f8'),
            ('nstack','i4'),
            ('skyvar','f8'),
            ('images', 'f8', imshape)]

        data=zeros(nbin, dtype=dt)
        return data

    def get_struct(self):
        sxdt=self.cat.dtype.descr

        dt= [('psfnum','i2'),
             ('shnum','i2'),
             ('ccd','i2'),
             ('row_range','f8',2),
             ('col_range','f8',2),
             ('am_row','f8'), # note simid,id in sxdt
             ('am_col','f8'),
             ('am_irr','f8'),
             ('am_irc','f8'),
             ('am_icc','f8'),
             ('am_e1','f8'),
             ('am_e2','f8'),
             ('am_rho4','f8'),
             ('am_a4','f8'),
             ('am_s2','f8'),
             ('am_uncer','f8'),
             ('am_s2n','f8'),
             ('am_numiter','i4'),
             ('am_flags','i4'),
             ('am_flagstr','S10'),
             ('am_shiftmax','f8')]
        dt=sxdt + dt
        data=zeros(self.cat.size, dtype=dt)
        for n in self.cat.dtype.names:
            data[n][:] = self.cat[n][:]

        data['psfnum'] = self['psfnum']
        data['shnum'] = self['shnum']
        data['ccd'] = self['ccd']
        return data



class MomentPipe(dict):
    def __init__(self, **keys):
        """
        parameters
        ----------
        run:
            run id
        psfnum:
            psf number
        shnum:
            shear number
        ccd:
            ccd number

        version: optional
            version id

        bugs
        ----

        the psf and shear row,col outputs are in the sub-image coordinates
        Currently need to add row_range[0],col_range[0] to get back to
        image coords

        """

        self._check_keys(**keys)

        for k in keys:
            self[k] = keys[k]

        conf=files.read_config(self['run'])

        for k in conf:
            self[k]=conf[k]

        self._load_data()

    def _check_keys(self, **keys):
        if ('run' not in keys
                or 'psfnum' not in keys
                or 'shnum' not in keys
                or 'ccd' not in keys):
            raise ValueError("send run=, psfnum=, "
                             "shnum=, ccd=")


    def _load_data(self):
        self.cat=files.read_cat(**self)
        image_orig, self.hdr=files.read_image(**self)

        shear_path=files.get_output_path(ftype='shear', **self)
        if os.path.exists(shear_path):
            self.shear_res=files.read_fits_output(ftype='shear',**self)

        self.seg, self.seg_hdr=files.read_image(ftype='seg', **self)

        skypersec=self.hdr['sky']
        exptime=self.hdr['exptime']

        # assume gain same for both amplifiers
        gain=self.hdr['gaina']
        read_noise=self.hdr['rdnoisea']

        # note unit mismatch? Both check out
        self['sky'] = skypersec*exptime/gain
        self['skyvar'] = self.hdr['sky'] + read_noise
        self['skysig'] = sqrt(self['skyvar'])
        self['ivar'] = 1.0/self['skyvar']

        self.image = image_orig.astype('f8') - self['sky']
        del image_orig

    def _set_psfres(self, psfres):
        self.psfres=psfres

    def get_cutout(self, index, size=None):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]
        id=self.cat['id'][index]

        padding=self.get('seg_padding',0)
        include_all_seg=self.get('include_all_seg',True)

        cutout=CutoutWithSeg(self.image, self.seg, cen, id, size,
                             include_all_seg=include_all_seg,
                             padding=padding)

        return cutout

    def get_full_cutout(self, index, size=None):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]

        cutout=Cutout(self.image, cen, size)

        return cutout


    def get_zerod_cutout(self, index, **keys):
        try:
            c=self.get_cutout(index, **keys)
        except NoSegMatches as excpt:
            print >>stderr,str(excpt)
            # sometimes there are no associated
            # segment pixels!
            c=self.get_full_cutout(index, **keys)
            return c

        im=c.subimage
        seg=c.seg_subimage

        id=self.cat['id'][index]
        self.zero_seg(im, seg, id)

        return c

    def zero_seg(self, im, seg, id):
        w=where( (seg != id) & (seg != 0) )
        if w[0].size != 0:
            im[w] = self['skysig']*random.randn(w[0].size)


    def get_stars(self):
        """
        Get things labeled as stars.  if good, also trim
        to those with good adaptive moments
        """
        logic=self.cat['class']==CLUSTERSTEP_STAR

        w,=where(logic)
        return w

    def get_gals(self):
        """
        Get things labeled as galaxies.  if good, also trim
        to those with good adaptive moments
        """

        logic=self.cat['class']==CLUSTERSTEP_GAL

        w,=where(logic)
        return w

    def get_psf_stars(self):
        """
        Get things labeled as stars in the psf star mag range.  if good, also
        trim to those with good adaptive moments
        """
        star_logic = self.cat['class']==CLUSTERSTEP_STAR
        mag_logic  = ( (self.cat['mag_auto_r'] > self['star_minmag'])
                       & (self.cat['mag_auto_r'] < self['star_maxmag']))

        logic = star_logic & mag_logic

        w,=where(logic)
        return w

    
    def run(self, **keys):

        use_seg=self.get('use_seg',True)

        res=self.get_struct()
        self.res=res

        do_isolated=self.get('isolated',False)
        for i in xrange(self.cat.size):
            if (i % 100) == 0:
                print '%d/%d' % (i+1,self.cat.size)

            # we might use full below
            c=self.get_zerod_cutout(i)

            im=c.subimage
            cen=c.subcen
            self.res['row_range'][i] = c.row_range
            self.res['col_range'][i] = c.col_range

            # this modifies self.res with the values, flags, etc.
            ares=self._do_admom(im, i, cen)

            if ares['whyflag'] != 0:
                # we need a center!
                self.res['uw_flags'][i] = UW_BAD_ADMOM
                continue

            self.res['am_row'][i] = c.row_range[0] + ares['wrow']
            self.res['am_col'][i] = c.col_range[0] + ares['wcol']

            acen=[ares['wrow'], ares['wcol']]

            if do_isolated:
                if not self._is_isolated(c, i):
                    self.res['uw_flags'][i] += UW_CROWDED
                    continue

            if not use_seg:
                cfull=self.get_full_cutout(i)
                self._do_uw_moments(cfull.subimage, i, acen)
            else:
                self._do_uw_moments(im, i, acen)

        w, = where(self.res['uw_flags'] != 0)
        print 'found %d/%d bad   %s' % (w.size,res.size,w.size/float(res.size))

        files.write_fits_output(data=self.res, ftype='shear', **self)

    def _is_isolated(self, cutout, index):
        id=self.cat['id'][index]
        w=where(  (cutout.seg_subimage != id) 
                & (cutout.seg_subimage != 0))
        if w[0].size==0:
            return True
        else:
            #print '%s/%s' % (w[0].size, cutout.seg_subimage.size)
            return False

    def _do_admom(self, im, i, cen):
        import admom
        ares = admom.admom(im, cen[0], cen[1],
                          sigsky=self['skysig'],
                          guess=4.0,
                          nsub=1)

        if ares['whyflag'] != 0:
            print 'index: %4d flag: %s' % (i,ares['whystr'])
        
        for n in ares:
            rn='am_'+n.lower()
            if (rn in self.res.dtype.names 
                    and n not in ['row','col','wrow','wcol']):
                self.res[rn][i] = ares[n]

        # two are differently named
        self.res['am_flags'][i] = ares['whyflag']
        self.res['am_flagstr'][i] = ares['whystr']

        return ares


    def _do_uw_moments(self, im, i, cen):
        """
        not currently used
        """

        row,col=numpy.mgrid[0:im.shape[0], 0:im.shape[1]]
        rm=row.astype('f8')-cen[0]
        cm=col.astype('f8')-cen[1]

        isum=im.sum()

        irrvals = im*rm**2
        ircvals = im*rm*cm
        iccvals = im*cm**2

        irrsum = irrvals.sum()
        ircsum = ircvals.sum()
        iccsum = iccvals.sum()

        self.res['uw_rmax'][i] = 9999
        self.res['uw_isum'][i] = isum
        self.res['uw_irrsum'][i] = irrsum
        self.res['uw_ircsum'][i] = ircsum
        self.res['uw_iccsum'][i] = iccsum



    def _do_uw_moments_trim(self, im, i, cen):
        """
        not currently used
        """
        if self.res['am_flags'][i] != 0:
            self.res['uw_flags'][i] = UW_BAD_ADMOM
            return

        nsig=self['uw_nsig']
        T=self.res['am_irr'][i] + self.res['am_icc'][i]

        arad=sqrt(T/2)
        rmax = nsig*arad
        rmax2=rmax**2

        row,col=numpy.mgrid[0:im.shape[0], 0:im.shape[1]]
        rm=row.astype('f8')-cen[0]
        cm=col.astype('f8')-cen[1]
        rad2=rm**2 + cm**2

        w=where(rad2 <= rmax2)
        if w[0].size < self['uw_min_pixels']:
            self.res['uw_flags'][i] += UW_TOO_FEW_PIXELS
            return

        isum=im[w].sum()

        rr = im*rm**2
        rc = im*rm*cm
        cc = im*cm**2

        irrsum = rr[w].sum()
        ircsum = rc[w].sum()
        iccsum = cc[w].sum()

        self.res['uw_rmax'][i] = rmax
        self.res['uw_npix'][i] = w[0].size
        self.res['uw_isum'][i] = isum
        self.res['uw_irrsum'][i] = irrsum
        self.res['uw_ircsum'][i] = ircsum
        self.res['uw_iccsum'][i] = iccsum


    def get_struct(self):
        sxdt=self.cat.dtype.descr
        dt= [('psfnum','i2'),
             ('shnum','i2'),
             ('ccd','i2'),
             ('row_range','f8',2),
             ('col_range','f8',2),
             ('am_row','f8'), # note simid,id in sxdt
             ('am_col','f8'),
             ('am_irr','f8'),
             ('am_irc','f8'),
             ('am_icc','f8'),
             ('am_e1','f8'),
             ('am_e2','f8'),
             ('am_rho4','f8'),
             ('am_a4','f8'),
             ('am_s2','f8'),
             ('am_uncer','f8'),
             ('am_s2n','f8'),
             ('am_numiter','i4'),
             ('am_flags','i4'),
             ('am_flagstr','S10'),
             ('am_shiftmax','f8'),
             ('uw_flags','i4'),
             ('uw_type','i1'), # CLUSTERSTEP_GAL,CLUSTERSTEP_STAR,PSF_STAR
             ('uw_nsig','f8'),
             ('uw_rmax','f8'), # nsig*radius from adaptive moments
             ('uw_npix','i4'),
             ('uw_isum','f8'),
             ('uw_irrsum','f8'),
             ('uw_ircsum','f8'),
             ('uw_iccsum','f8')]

        dt=sxdt + dt
        data=zeros(self.cat.size, dtype=dt)
        for n in self.cat.dtype.names:
            data[n][:] = self.cat[n][:]

        wpsf=self.get_psf_stars() 
        data['uw_nsig'] = self.get('uw_nsig',9999)
        data['uw_type'] = data['class']
        data['uw_type'][wpsf] += CLUSTERSTEP_PSF_STAR

        data['psfnum'] = self['psfnum']
        data['shnum'] = self['shnum']
        data['ccd'] = self['ccd']
        return data


    def plot_admom_sizemag(self, show=False):
        """
        Plot the size-mag.  Show psf stars in a different color.
        """
        import biggles

        if not hasattr(self, 'ares'):
            raise ValueError("run_admom first")

        wstar=self.get_stars()
        wgal=self.get_gals()
        wpsf=self.get_psf_stars()


        plt=biggles.FramedPlot()
        sigma=sqrt( (self.ares['Irr']+self.ares['Icc'])/2 )
        mag=self.cat['mag_auto_r']
        psize=0.7
        pstar=biggles.Points(mag[wstar], sigma[wstar], 
                             type='filled circle', color='red',
                             size=psize)
        pgal=biggles.Points(mag[wgal], sigma[wgal], 
                            type='filled diamond', color='dark green',
                            size=psize)

        ppsf=biggles.Points(mag[wpsf], sigma[wpsf], 
                            type='filled circle', color='blue',
                            size=psize)

        pstar.label='star'
        pgal.label='gal'
        ppsf.label='psf star'

        key=biggles.PlotKey(0.95,0.92, [pstar,ppsf,pgal], halign='right')

        plt.add(pstar, ppsf, pgal, key)
        plt.xlabel='mag_auto_r'
        plt.ylabel=r'$\sigma_{AM}$ [pixels]'
        label='%s-p%s-s%s-%s' % (self['run'],self['psfnum'],self['shnum'],
                                 self['ccd'])
        plt.add(biggles.PlotLabel(0.05,0.92,label,halign='left'))

        plt.xrange=[16.5,25.5]
        plt.yrange=[0.75,8]
        if show:
            plt.show()

        epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)
        self._write_plot(plt, epsfile)



    def _write_plot(self, plt, epsfile):
        import converter
        d=os.path.dirname(epsfile)
        if not os.path.exists(d):
            try:
                os.makedirs(d)
            except:
                pass
        print 'writing eps file:',epsfile
        plt.write_eps(epsfile)

        converter.convert(epsfile, dpi=100)



    def show_many_cutouts(self, indices, **keys):
        for i in indices:
            self.show_cutout(i, **keys)
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                return


    def show_cutout(self, index, size=None, 
                    zero_seg=False):
        import biggles
        biggles.configure( 'default', 'fontsize_min', 2)
        if not hasattr(self,'ares'):
            self.run_admom()

        import biggles
        import images

        c=self.get_cutout(index, size=size)
        cz=self.get_zerod_cutout(index, size=size)

        minv=-2.5*self['skysig']

        subimage=c.get_subimage().copy()
        segsub=c.get_seg_subimage().copy()

        zsubimage=cz.get_subimage().copy()

        objtype='star'
        if self.cat['class'][index]==CLUSTERSTEP_GAL:
            objtype='gal'

        s2n=self.ares['s2n'][index]
        print 'index:',index
        print 's2n:',s2n
        print 'ranges:'
        print '\t',c.row_range
        print '\t',c.col_range
        print 'image shape:',subimage.shape
        print

        # make the other objects always show darker
        id=self.cat['id'][index]
        w=where((segsub != id) & (segsub != 0))
        if w[0].size > 0:
            segsub[w] = id*10

        implt=images.view(subimage,show=False, min=minv)
        segplt=images.view(segsub,show=False)

        title='%d %s S/N: %.1f [%d,%d]'
        title = title % (index, objtype,s2n,subimage.shape[0],subimage.shape[1])
        implt.title=title
        segplt.title='seg'

        zplt=images.view(zsubimage,show=False, min=minv)
        zplt.title='zerod'

        tab=biggles.Table(2,2)
        tab[0,0]=implt
        tab[0,1]=segplt
        tab[1,0]=zplt

        tab.show()


class NoSegMatches(Exception):
    def __init__(self, value):
        self.value = str(value)
    def __str__(self):
        return self.value

class CutoutWithSeg:
    def __init__(self, image, seg, cen, id, minsize, 
                 include_all_seg=True, padding=0):
        self.image=image
        self.padding=padding
        self.seg=seg
        self.cen=cen
        self.minsize=minsize
        self.id=id

        self.include_all_seg=include_all_seg

        self._set_box()
        self._make_cutout()

    def get_subimage(self):
        return self.subimage
    def get_seg_subimage(self):
        return self.seg_subimage

    def get_subcen(self):
        return self.subcen

    def _get_seg_box(self):
        """
        Region containing the seg pixels
        """
        w=where(self.seg == self.id)

        if w[0].size == 0:
            mess="no seg pixels with id %s" % self.id
            raise NoSegMatches(mess)

        minrow = w[0].min() - self.padding
        maxrow = w[0].max() + self.padding
        mincol = w[1].min() - self.padding
        maxcol = w[1].max() + self.padding

        return minrow,maxrow,mincol,maxcol

    def _get_minimal_box(self):
        sh=self.image.shape
        cen=self.cen
        size=self.minsize

        if (cen[0] < 0 or cen[1] < 0
                or cen[0] > (sh[0]-1)
                or cen[1] > (sh[1]-1) ):
            mess=("center [%s,%s] is out of "
                  "bounds of image: [%s,%s] ")
            mess=mess % (cen[0],cen[1],sh[0],sh[1])
            raise ValueError(mess)

        sz2 = (size-1)/2.
        minrow=int(cen[0]-sz2 )
        maxrow=int(cen[0]+sz2)
        mincol=int(cen[1]-sz2)
        maxcol=int(cen[1]+sz2)
        
        return minrow,maxrow,mincol,maxcol

    def _set_box(self):

        sh=self.image.shape
        minrow0,maxrow0,mincol0,maxcol0 = self._get_minimal_box()
        if self.include_all_seg:
            minrow,maxrow,mincol,maxcol     = self._get_seg_box()

            if minrow0 < minrow:
                minrow = minrow0
            if maxrow0 > maxrow:
                maxrow=maxrow0
            if mincol0 < mincol:
                mincol = mincol0
            if maxcol0 > maxcol:
                maxcol=maxcol0
        else:
            minrow,maxrow,mincol,maxcol = minrow0,maxrow0,mincol0,maxcol0

        if minrow < 0:
            minrow=0
        if maxrow > (sh[0]-1):
            maxrow=sh[0]-1

        if mincol < 0:
            mincol=0
        if maxcol > (sh[1]-1):
            maxcol=sh[1]-1


        self.row_range=[minrow,maxrow]
        self.col_range=[mincol,maxcol]


    def _make_cutout(self):
        cen=self.cen

        minrow,maxrow=self.row_range
        mincol,maxcol=self.col_range

        # note +1 for python slices
        self.subimage=self.image[minrow:maxrow+1, mincol:maxcol+1].copy()
        self.seg_subimage=self.seg[minrow:maxrow+1, mincol:maxcol+1].copy()
        self.subcen=[cen[0]-minrow, cen[1]-mincol]


class Cutout:
    def __init__(self, image, cen, size):
        self.image=image
        self.cen=cen
        self.size=size
        self._make_cutout()

    def get_subimage(self):
        return self.subimage

    def get_subcen(self):
        return self.subcen

    def _make_cutout(self):
        sh=self.image.shape
        cen=self.cen
        size=self.size

        if (cen[0] < 0 or cen[1] < 0
                or cen[0] > (sh[0]-1)
                or cen[1] > (sh[1]-1) ):
            mess=("center [%s,%s] is out of "
                  "bounds of image: [%s,%s] ")
            mess=mess % (cen[0],cen[1],sh[0],sh[1])
            raise ValueError(mess)

        #minrow=int(     cen[0]-size/2.-0.5)
        #maxrow=int(ceil(cen[0]+size/2.+0.5))
        #mincol=int(     cen[1]-size/2.-0.5)
        #maxcol=int(ceil(cen[1]+size/2.+0.5))
        sz2 = (size-1)/2.
        minrow=int(cen[0]-sz2 )
        maxrow=int(cen[0]+sz2)
        mincol=int(cen[1]-sz2)
        maxcol=int(cen[1]+sz2)
 
        if minrow < 0:
            minrow=0
        if maxrow > (sh[0]-1):
            maxrow=sh[0]-1

        if mincol < 0:
            mincol=0
        if maxcol > (sh[1]-1):
            maxcol=sh[1]-1
        
        # note +1 for python slices
        self.subimage=self.image[minrow:maxrow+1, mincol:maxcol+1].copy()
        self.subcen=[cen[0]-minrow, cen[1]-mincol]

        self.row_range=[minrow,maxrow]
        self.col_range=[mincol,maxcol]

