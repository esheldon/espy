from sys import stderr

import numpy
import sdsspy
import es_sdsspy

from esutil.random import srandu

import gmix_image
import admom


from . import files

class GMixSweep(dict):
    def __init__(self, gmix_run, run, camcol, **keys):
        
        conf=files.read_config(gmix_run)
        self.update(conf)

        self.run=run
        self.camcol=camcol

        self._load_data()

        if self['make_plots'] or self['make_psf_plots']:
            import biggles
            biggles.configure("screen","width",1100)
            biggles.configure("screen","height",1100)


    def process(self):
        
        min_field=self.objs['field'].min()
        max_field=self.objs['field'].max()
        for field in xrange(min_field, max_field+1):
            self._process_field(field)

    def _process_field(self, field):
        w,=numpy.where(self.objs['field']==field)
        if w.size==0:
            return

        psfield,psf_kl = self._read_psfield(field)

        print 'field:',field,'objects',w.size
        for i in xrange(w.size):
            stderr.write('%s/%s\n' % ((i+1), w.size))
            index=w[i]
            res=self._process_object(index, psfield, psf_kl)

    def _process_object(self, index, psfield, psf_kl):

        obj=self.objs[index]

        result={}
        for filt in self['filters']:
            fnum=sdsspy.FILTERNUM[filt]

            psf_res=self._measure_psf(obj,psf_kl,fnum)
            if psf_res is None:
                continue

            psf_gmix=psf_res.get_gmix()
            res=self._measure_obj(obj, psfield, fnum, psf_gmix)

            result[filt] = res
        
        return result

    def _get_sigma_guess(self, psf_gmix, m_rr_cc):

        psf_sigma=numpy.sqrt(psf_gmix.get_T()/2.)

        if m_rr_cc < 0:
            sigma_guess=psf_sigma
        else:
            sigma_guess=numpy.sqrt(m_rr_cc/2.)
            if sigma_guess < psf_sigma:
                sigma_guess=psf_sigma

        if sigma_guess > 20:
            sigma_guess=20

        return sigma_guess
            
    def _prepare_atlas(self, im, background, row, col):
        """
        Trim empty space, convert to 'f8', and background
        subtract

        can do better than this
        """
        w1,w2 = numpy.where(im != background)
        if w1.size > 0:
            minrow = w1.min()
            maxrow = w1.max()
            mincol = w2.min()
            maxcol = w2.max()
            
            im_out=im[minrow:maxrow+1, mincol:maxcol+1].astype('f8')
            row_out = row-minrow
            col_out = col-mincol
        else:
            # no signal
            im_out=im.astype('f8')

        im_out -= background
        return im_out,row_out,col_out

    def _measure_psf(self, obj, psf_kl, fnum):
        rowc=obj['rowc'][fnum]
        colc=obj['colc'][fnum]

        psf_im = psf_kl[fnum].rec(rowc, colc, trim=True)
        guess_psf=admom.fwhm2mom( obj['psf_fwhm'][fnum], pixscale=0.4)/2.

        psf_res=self._fit_psf(psf_im, guess_psf, 1.0)

        return psf_res


    def _fit_psf(self, im, sigma_guess, skysig):

        row=im.shape[0]/2.
        col=im.shape[1]/2.

        ares=self._run_admom(im, row, col, sigma_guess, 1.0)
        if ares is None:
            print >>stderr,'failed to get psf admom'
            return None

        res = self._fit_psf_model(im, ares, self['psf_model'], skysig)
        if res is None:
            print >>stderr,'failed to fit psf model'

        return res

    def _fit_psf_model(self, im, ares, model, skysig):
        if model=='gmix1':
            ngauss=1
        elif model=='gmix2':
            ngauss=2
        elif model=='gmix3':
            ngauss=3
        else:
            raise ValueError("bad psf model: '%s'" % model)

        ntry=1
        res=None
        while ntry <= self['psf_lm_ntry']:
            tmp=gmix_image.gmix_fit.quick_fit_psf_coellip(im, 
                                                          skysig, 
                                                          ngauss, 
                                                          ares=ares)
            if tmp is None:
                continue
            if tmp.flags==0:
                res=tmp
                break
            ntry += 1

        if self['make_psf_plots'] and res is not None:
            import images
            model=gmix_image.render.gmix2image(res.get_gmix(),im.shape)
            images.compare_images(im,model,label1='image',label2='model')
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                stop

        return res


    def _measure_obj(self, obj, psfield, fnum, psf_gmix):

        atlas=self._read_atlas(obj)
        background=atlas['SOFT_BIAS']

        rowc=obj['rowc'][fnum]
        colc=obj['colc'][fnum]
        row_notrim = rowc - atlas['row0'][fnum] - 0.5
        col_notrim = colc - atlas['col0'][fnum] - 0.5
        
        im_notrim=atlas['images'][fnum]
        im,row,col=self._prepare_atlas(im_notrim, background, 
                                       row_notrim, col_notrim)
        
        sigma_guess=self._get_sigma_guess(psf_gmix,obj['m_rr_cc'][fnum])

        skysig = psfield['skysig'][0,fnum]

        res=self._fit_object(obj, im, skysig, row, col, sigma_guess, psf_gmix)

        return res

    def _fit_object(self, obj, im, skysig, row, col, sigma_guess, psf_gmix):

        ares=self._run_admom(im, row, col, sigma_guess, skysig)
        if ares is None:
            print >>stderr,'failed to run object admom'
            return None

        mag=None

        results={'ares':ares}
        for fitter_type in self['obj_fitters']:
            results[fitter_type] = {}
            for model in self['obj_models']:
                if fitter_type=='lm':
                    print 'lm:',model,
                    fitter=self._fit_object_model_lm(
                            im, mag, skysig, ares, model, psf_gmix)
                elif fitter_type=='mcmc':

                    fitter=self._fit_object_model_mcmc(
                            im, mag, skysig, ares, model, psf_gmix)
                    res=fitter.get_result()
                    print '  model: %s prob: %.6f aic: %.6f bic: %.6f Ts/n: %.6f ' % \
                            (model,res['fit_prob'],res['aic'],res['bic'],res['Ts2n'])
                else:
                    raise ValueError("bad fitter type '%s'" % fitter_type)
                results[fitter_type][model] = fitter

            if fitter_type=='lm':
                print
        return results

    def _fit_object_model_lm(self, im, mag, skysig, ares, model, psf_gmix):
        from gmix_image.gmix_fit import GMixFitSimple
        ivar=1./skysig**2
        cen_width=self.get('cen_width',1.0)

        # can be None
        gprior=self._get_gprior(mag, model)

        aic=9.999e9
        for i in xrange(self['obj_lm_ntry']):
            fitter0=GMixFitSimple(im, 
                                  ivar,
                                  psf_gmix, 
                                  model,
                                  ares,
                                  cen_width=cen_width,
                                  gprior=gprior)
            res=fitter0.get_result()
            if res['aic'] < aic:
                aic=res['aic']
                fitter=fitter0
        return fitter


    def _fit_object_model_mcmc(self, im, mag, skysig, ares, model, psf_gmix):
        from gmix_image.gmix_mcmc import MixMCStandAlone

        ivar=1./skysig**2
        cen_width=self.get('cen_width',1.0)

        # cannot be None
        gprior=self._get_gprior(mag, model)

        nwalkers=self['nwalkers']
        burnin=self['burnin']
        nstep=self['nstep']

        fitter=MixMCStandAlone(im,
                               ivar,
                               psf_gmix, 
                               gprior, 
                               model,
                               nwalkers=nwalkers,
                               nstep=nstep,
                               burnin=burnin,
                               mca_a=self['mca_a'],
                               iter=self['iter'],
                               draw_gprior=True,
                               ares=ares,
                               cen_width=cen_width,
                               make_plots=self['make_plots'])

        return fitter
    def _run_admom(self, im, row, col, guess, skysig):

        i=0
        Tguess=4.

        res=None
        while i < self['admom_ntry']: 
            tmp = admom.admom(im, row, col,
                              sigsky=skysig,
                              guess=guess)
            if tmp['whyflag']==0:
                res=tmp
                break
            row = row*(1.+0.01*srandu())
            col = col*(1.+0.01*srandu())
            guess = guess*(1.+0.01*srandu())

            i += 1


        return res

    def _read_psfield(self, field):
        psfield=sdsspy.read('psfield',
                            run=self.run,
                            camcol=self.camcol,
                            field=field)
        psf_kl = []
        for filter in sdsspy.FILTERCHARS:
            kl = sdsspy.read('psField',
                             run=self.run,
                             camcol=self.camcol,
                             field=field,
                             filter=filter)
            psf_kl.append(kl)

        return psfield, psf_kl


    def _read_atlas(self, obj):
        """
        Read atlas images as sky-subtracted, 'f8'
        """

        field=obj['field']
        id=obj['id']
        atlas = sdsspy.read('fpAtlas',
                            run=self.run,
                            camcol=self.camcol,
                            field=field,
                            id=id)
        #for i in xrange(5):
        #    im=atlas['images'][i].astype('f8')-atlas['SOFT_BIAS']
        #    atlas['images'][i] = im
        return atlas

    def _select(self, objs):
        print "Selecting objects"
        s = es_sdsspy.select.Selector(objs)

        print "  getting resolve logic"
        resolve_logic = s.resolve_logic()

        print "  getting flag logic"
        flag_logic = s.flag_logic()

        print "  getting rmag logic"
        rmag_logic = s.cmodelmag_logic("r", self['max_cmodelmag_r'])

        #rmag_logic=rmag_logic & (s.mags['cmodel_dered'][:,2] <18)

        logic = resolve_logic & flag_logic & rmag_logic

        keep, = numpy.where(logic)
        print "  keeping %i/%i" % (keep.size, objs.size)

        return keep

    def _load_data(self):
        objs=sdsspy.read('calibObj.gal', 
                         run=self.run, 
                         camcol=self.camcol, 
                         lower=True,
                         verbose=True)
        w=self._select(objs)
        objs=objs[w]

        self.objs=objs

    def _get_gprior(self, mag, model):
        if not hasattr(self, 'gpriors'):
            self._set_gpriors_vs_mag()

        gprior=None
        for pdict in self.gpriors[model]:
            if mag >= pdict['minmag'] and mag <= pdict['maxmag']:
                gprior=pdict['gprior']
        if gprior is None:
            if mag < self.gpriors[model][0]['minmag']:
                gprior=self.gpriors[model][0]['gprior']
            elif mag > self.gpriors[model][-1]['maxmag']:
                gprior=self.gpriors[model][-1]['gprior']
            else:
                raise ValueError("not possible error finding mag: nan?")

        return gprior


    def _set_gpriors_vs_mag(self):
        """
        Note we use the GPriorExp for both dev and exp, but
        different galaxies were used to train it
        """
        import cluster_step
        exp_prior_pars=cluster_step.files.read_gprior(type='gexp')
        dev_prior_pars=cluster_step.files.read_gprior(type='gdev')

        gpriors={}

        exp_plist=[]
        for i in xrange(exp_prior_pars.size):
            exp_pdict={}

            pexp=cluster_step.prior.GPriorExp(exp_prior_pars['pars'][i])

            exp_pdict['gprior'] = pexp
            exp_pdict['minmag'] = exp_prior_pars['minmag'][i]
            exp_pdict['maxmag'] = exp_prior_pars['maxmag'][i]

            exp_plist.append(exp_pdict)


        dev_plist=[]
        for i in xrange(dev_prior_pars.size):
            dev_pdict={}

            pdev=cluster_step.prior.GPriorExp(dev_prior_pars['pars'][i])

            dev_pdict['gprior'] = pdev
            dev_pdict['minmag'] = dev_prior_pars['minmag'][i]
            dev_pdict['maxmag'] = dev_prior_pars['maxmag'][i]

            dev_plist.append(dev_pdict)

        gpriors['gexp']=exp_plist
        gpriors['gdev']=dev_plist

        self.gpriors=gpriors


