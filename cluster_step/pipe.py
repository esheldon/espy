import os
from sys import stdout
import pprint
from math import ceil
from numpy import where, sqrt, random, zeros, arange, median

from . import files
from . import prior

import gmix_image
from gmix_image.gmix import GMix
from gmix_image.gmix_mcmc import MixMCStandAlone, MixMCPSF
from gmix_image.gmix_em import GMixEMPSF

CLUSTERSTEP_GAL=1
CLUSTERSTEP_STAR=2

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
        """

        self._check_keys(**keys)

        for k in keys:
            self[k] = keys[k]

        conf=files.read_config(self['run'])

        for k in conf:
            self[k]=conf[k]
        
        self._set_priors()

        random.seed(self['seed'])
        self._load_data()

    def _check_keys(self, **keys):
        if ('run' not in keys
                or 'psfnum' not in keys
                or 'shnum' not in keys
                or 'ccd' not in keys):
            raise ValueError("send run=, psfnum=, "
                             "shnum=, ccd=")

    def _set_priors(self):
        priors={}
        if self['gprior_type'] == 'old':
            priors['gexp']=gmix_image.priors.GPrior(A=12.25,
                                                    B=0.2,
                                                    C=1.05,
                                                    D=13.)
            priors['gdev'] = priors['gexp']
        else:
            priors['gexp']=prior.GPriorExp(self['gprior_pars_exp'])
            priors['gdev']=prior.GPriorDev(self['gprior_pars_dev'])
        self.gpriors=priors


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

    def get_cutout(self, index, size=None, with_seg=False):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]

        cutout=Cutout(self.image, cen, size)

        if with_seg:
            segcut=Cutout(self.seg, cen, size)
            return cutout, segcut

        return cutout

    def get_zerod_cutout(self, index, **keys):
        c,cseg=self.get_cutout(index, with_seg=True, **keys)

        im=c.subimage
        seg=cseg.subimage

        id=self.cat['id'][index]
        self.zero_seg(c.subimage, cseg.subimage, id)

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
        Get things labeled as stars.  if good, also trim
        to those with good adaptive moments
        """

        logic=self.cat['class']==CLUSTERSTEP_GAL

        if good:
            logic = logic & (self.ares['whyflag']==0)

        if s2n_min is not None:
            logic=logic & (self.ares['s2n'] > s2n_min)

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
            self.plot_admom_sizemag()
        if not hasattr(self,'psfres') or run_psf:
            self.run_psf()

        ares=self.ares

        #wpsf=self.get_psf_stars()
        wgal=self.get_gals(s2n_min=self['shear_s2n_min'])
        #wgal=wgal[0:10]

        out=self.get_shear_struct(wgal.size)

        out['id'][:] = self.cat['id'][wgal]
        out['simid'][:] = self.cat['simid'][wgal]

        for igal in xrange(wgal.size):
            from esutil.numpy_util import aprint
            #if ((igal % 10) == 0):
            #    print "  %s/%s done" % (igal+1,wgal.size)

            index=wgal[igal]

            print "%s/%s s2n: %s" % (igal+1,wgal.size,ares['s2n'][index])

            c=self.get_zerod_cutout(index)

            im=c.subimage

            gmix_psf=self.get_gmix_psf()

            res=self.fit_shear_models(im, ares[index].copy(), gmix_psf)
            self.copy_shear_results(out, res, gmix_psf, igal)

        self.shear_res=out
        files.write_fits_output(data=out,
                                ftype='shear',
                                **self)


    def fit_shear_models(self, im, ares, gmix_psf):
        """
        Fit all listed models, return the best fitting
        """
        probrand=-9999.
        fitmodels=self.get_fitmodels()
        for fitmodel in fitmodels:
            fitter=self.run_shear_model(im, ares, gmix_psf, fitmodel)

            res0 = fitter.get_result()
            if len(fitmodels) > 1:
                print '  model: %s prob: %.6f aic: %.6f bic: %.6f Ts/n: %.6f ' % \
                    (fitmodel,res0['fit_prob'],res0['aic'],res0['bic'],res0['Ts2n'])
            if res0['fit_prob'] > probrand:
                res=res0
                probrand=res0['fit_prob']

        if len(fitmodels) > 1:
            print '    best model:',res['model']
        return res


    def run_shear_model(self, im, ares0, gmix_psf, fitmodel):
        """
        Run the shear model though the mixmc code
        """

        ares={'wrow':ares0['wrow']-ares0['row_range'][0],
              'wcol':ares0['wcol']-ares0['col_range'][0],
              'Irr':ares0['Irr'],
              'Irc':ares0['Irc'],
              'Icc':ares0['Icc'],
              'whyflag':ares0['whyflag']}
 
        gprior=self.gpriors[fitmodel]

        fitter=MixMCStandAlone(im, self['ivar'],
                               gmix_psf, gprior, fitmodel,
                               nwalkers=self['nwalkers'],
                               nstep=self['nstep'], 
                               burnin=self['burnin'],
                               mca_a=self['mca_a'],
                               iter=self.get('iter',False),
                               draw_gprior=self['draw_gprior'],
                               ares=ares)
        return fitter

    def get_fitmodels(self):
        fitmodels=self['fitmodel']
        if not isinstance(fitmodels,list):
            fitmodels=[fitmodels]
        return fitmodels

    def get_gmix_psf(self):
        """
        Get a psf model
        """
        if self['psf_interp']=='random':
            return self.get_random_gmix_psf()
        elif self['psf_interp']in ['mean','median']:
            return self.get_mean_gmix_psf()

    def _set_best_psf_model(self):
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
                if self['psf_interp']=='median':
                    pars[i]=median(allpars[:,i])
                else:
                    pars[i]=allpars[:,i].mean()

        gmix=self._get_model_psf_gmix(pars,model)
        return gmix

    def _get_model_psf_gmix(self, pars, model):
        if model=='gturb':
            # note the MixMCPSF uses e1,e2 not g1,g2
            gmix=GMix(pars, model='gturb')
        else:
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

        ares=self.ares

        wpsf=self.get_psf_stars()

        out=self.get_psf_struct(wpsf.size)
        out['id'][:] = self.cat['id'][wpsf]
        out['simid'][:] = self.cat['simid'][wpsf]

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
                         'whyflag':whyflag}

                c=self.get_zerod_cutout(index)

                im=c.subimage

                if model in ['em1','em2','em2cocen']:
                    self.run_em_psf(im, aresi, model, out, ipsf)
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

        files.write_fits_output(data=out,
                                ftype='psf',
                                **self)


    def run_em_psf(self, im, aresi, model, out, ipsf):
        cocenter=False
        if model=='em1':
            ngauss=1
        elif model=='em2':
            ngauss=2
        elif model=='em2cocen':
            ngauss=2
            cocenter=True
        else:
            raise ValueError("bad psf model: '%s'" % model)
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

        print 'running admom'
        ares=self.get_admom_struct()

        for i in xrange(self.cat.size):
            c=self.get_zerod_cutout(i)

            im=c.subimage
            cen=c.subcen

            res = admom.admom(im, cen[0], cen[1],
                              sigsky=self['skysig'],
                              guess=4.0,
                              nsub=1)

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


    def show_cutout(self, index, size=None, with_seg=False,
                    zero_seg=False):
        import biggles
        import images


        minv=-2.5*self['skysig']
        if not with_seg:
            c=self.get_cutout(index, size=size)
            subimage=c.get_subimage()


            images.multiview(subimage, min=minv)
        else:
            c,segc=self.get_cutout(index, size=size, with_seg=True)

            subimage=c.get_subimage()
            segsub=segc.get_subimage()


            implt=images.view(subimage,show=False, min=minv)
            segplt=images.view(segsub,show=False)

            implt.title='image cutout'
            segplt.title='seg cutout'

            if zero_seg:
                self.zero_seg(subimage, segsub, self.cat['id'][index])
                zplt=images.view(subimage,show=False, min=minv)
                zplt.title='zerod'

                tab=biggles.Table(2,2)
                tab[0,0]=implt
                tab[0,1]=segplt
                tab[1,0]=zplt
            else:
                tab=biggles.Table(1,2)
                tab[0,0]=implt
                tab[0,1]=segplt

            tab.show()

    def get_admom_struct(self):
        sxdt=self.cat.dtype.descr
        dt= [('wrow','f8'),
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
        npars1=6
        npars2=2*6
        npars_turb=6
        dt= [('id','i4'),
             ('simid','i4'),
             ('em1_flags','i4'),
             ('em1_numiter','i4'),
             ('em1_fdiff','f8'),
             ('em1_ntry','i4'),
             ('em1_pars','f8',npars1),
             ('em1_prob','f8'),
             ('em1_aic','f8'),
             ('em1_bic','f8'),
             ('em2_flags','i4'),
             ('em2_numiter','i4'),
             ('em2_fdiff','f8'),
             ('em2_ntry','i4'),
             ('em2_pars','f8',npars2),
             ('em2_prob','f8'),
             ('em2_aic','f8'),
             ('em2_bic','f8'),
             ('em2cocen_flags','i4'),
             ('em2cocen_numiter','i4'),
             ('em2cocen_fdiff','f8'),
             ('em2cocen_ntry','i4'),
             ('em2cocen_pars','f8',npars2),
             ('em2cocen_prob','f8'),
             ('em2cocen_aic','f8'),
             ('em2cocen_bic','f8'),
             ('gturb_flags','i4'), # will always be zero
             ('gturb_pars','f8',npars_turb),
             ('gturb_pcov','f8',(npars_turb,npars_turb)),
             ('gturb_prob','f8'),
             ('gturb_aic','f8'),
             ('gturb_bic','f8')
             ]


        data=zeros(n, dtype=dt)

        for model in ['em1','em2','em2cocen']:
            data[model+'_fdiff']=9999.
            data[model+'_ntry']=9999
            data[model+'_pars']=-9999.
            data[model+'_aic'] = 1.e9
            data[model+'_bic'] = 1.e9
        data['gturb_aic'] = 1.e9
        data['gturb_bic'] = 1.e9

        return data

    def copy_shear_results(self, out, res, gmix_psf, igal):
        """
        Copy results into existing "out" structure
        """
        out['s2n_admom'][igal] = self.ares['s2n'][igal]

        e1psf,e2psf,Tpsf=gmix_psf.get_e1e2T()
        Tobj=res['Tmean']
        s2 = Tpsf/Tobj

        out['Tpsf'][igal] = Tpsf
        out['e1psf'][igal] = e1psf
        out['e2psf'][igal] = e2psf
        out['s2'][igal] = s2

        #out['pars_psf'][igal] = gmix_psf.get_pars()

        for k in res:
            if k in out.dtype.names:
                out[k][igal] = res[k]


    def get_shear_struct(self, n):
        npars=6
        dt=[('id','i4'),
            ('simid','i4'),
            ('model','S20'),
            ('flags','i4'),
            ('e1psf','f8'),
            ('e2psf','f8'),
            ('Tpsf','f8'),
            ('s2','f8'),
            ('g','f8',2),
            ('gsens','f8',2),
            ('gcov','f8',(2,2)),
            ('pars','f8',npars),
            ('pcov','f8',(npars,npars)),
            #('pars_psf','f8',npars_psf),
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
        return data

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

        minrow=int(     cen[0]-size/2.-0.5)
        maxrow=int(ceil(cen[0]+size/2.+0.5))
        mincol=int(     cen[1]-size/2.-0.5)
        maxcol=int(ceil(cen[1]+size/2.+0.5))

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

