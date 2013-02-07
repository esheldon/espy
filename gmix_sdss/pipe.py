"""
TODO
    limit to BOSS footprint
        - its ok to write empty files
    limit galaxy measure to am s/n > 10 
        - mag limit stricter as well?

    Produce wq for each field of the minscore runs
        - add fieldlist method to window in sdsspy

    output structure
        - mark the best model!

    multiple psf models?


    Collation
"""
from sys import stderr

import numpy
import sdsspy
import es_sdsspy

import esutil as eu
from esutil.random import srandu

import gmix_image
import admom
import time

from . import files

NO_ATLAS=2**0
AM_PSF_FAILED=2**1
PSF_FAILED=2**2
AM_OBJ_FAILED=2**3
AM_FAINT=2**4  # am s/n < min_s2n


class GMixField(dict):
    def __init__(self, gmix_run, run, camcol, field):
        
        conf=files.read_config(gmix_run)
        self.update(conf)

        self['gmix_run'] = self['run']
        self['run']=run
        self['camcol']=camcol
        self['field']=field
        self['fnum']=sdsspy.FILTERNUM[self['filter']]

        self._load_data()

        if self['make_plots'] or self['make_psf_plots']:
            self._setup_plotting()

    def go(self):
        if self.objs is not None:

            nobj=self.objs.size
            st=self._get_struct()

            self._set_start_time()
            for i in xrange(self.objs.size):
                stderr.write('%s/%s\n' % ((i+1), nobj))

                obj=self.objs[i]

                # setting this means the image gets the same random numbers
                # added each time we run the code; but note the emcee sampler
                # uses its own internal random number generator

                self._set_object_seed(obj)

                res=self._process_object(obj)

                self._copy_to_output(st, i, res)
            self._print_time_stats()
        else:
            st=None

        files.write_output(data=st,**self)


    def _process_object(self, obj):

        result = {'flags':0}

        psf_res=self._measure_psf(obj)
        result['psf_res'] = psf_res
        result['flags'] += psf_res['flags']

        if psf_res['flags'] == 0:
            psf_gmix=psf_res['fitter'].get_gmix()
            res=self._measure_obj(obj, psf_gmix)
            result['obj_res'] = res
            result['flags'] += res['flags']
        else:
            result['obj_res'] = {'flags':PSF_FAILED}

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
            
    def _prepare_atlas(self, im, background, skysig, row, col):
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
            
            im_out=im[minrow:maxrow+1, mincol:maxcol+1]
            row_out = row-minrow
            col_out = col-mincol
        else:
            print 'no signal found'
            im_out=im
            row_out=row
            col_out=col

        wb=numpy.where(im_out==background)

        im_out=im_out.astype('f8')
        if wb[0].size != 0:
            im_out[wb] += skysig*numpy.random.randn(wb[0].size)

        im_out -= background
        return im_out,row_out,col_out

    def _set_object_seed(self, obj):
        pid=sdsspy.get_photoid(obj)
        numpy.random.seed(int(pid))

    def _measure_psf(self, obj):
        fnum=self['fnum']
        rowc=obj['rowc'][fnum]
        colc=obj['colc'][fnum]

        psf_im = self.psf_kl[fnum].rec(rowc, colc, trim=True)
        psf_im=psf_im.astype('f8')

        guess_psf=admom.fwhm2mom( obj['psf_fwhm'][fnum], pixscale=0.4)/2.

        psf_res=self._fit_psf(psf_im, guess_psf, 1.0)

        return psf_res


    def _fit_psf(self, im, sigma_guess, skysig):

        row=im.shape[0]/2.
        col=im.shape[1]/2.

        ares=self._run_admom(im, row, col, sigma_guess, 1.0)
        if ares is None:
            print >>stderr,'failed to get psf admom'
            return {'flags':AM_PSF_FAILED}

        res = self._fit_psf_model(im, ares, self['psf_model'], skysig)
        if res is None:
            print >>stderr,'failed to fit psf model'
            return {'flags':PSF_FAILED}

        return {'flags':0, 'fitter':res}

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


    def _measure_obj(self, obj, psf_gmix):

        fnum=self['fnum']

        atlas=self._read_atlas(obj)
        if atlas is None:
            return {'flags':NO_ATLAS}

        background=atlas['SOFT_BIAS']
        skysig = self.psfield['skysig'][0,fnum]

        rowc=obj['rowc'][fnum]
        colc=obj['colc'][fnum]
        row_notrim = rowc - atlas['row0'][fnum] - 0.5
        col_notrim = colc - atlas['col0'][fnum] - 0.5
        
        im_notrim=atlas['images'][fnum]

        im,row,col=self._prepare_atlas(im_notrim, background, skysig,
                                       row_notrim, col_notrim)
        
        sigma_guess=self._get_sigma_guess(psf_gmix,obj['m_rr_cc'][fnum])


        res=self._fit_object(obj, im, skysig, row, col, sigma_guess, psf_gmix)

        return res

    def _fit_object(self, obj, im, skysig, row, col, sigma_guess, psf_gmix):

        ares=self._run_admom(im, row, col, sigma_guess, skysig)
        if ares is None:
            print >>stderr,'failed to run object admom'
            return {'flags':AM_OBJ_FAILED}

        if ares['s2n'] < self['min_s2n']:
            mess='s/n %s is less than minimum %s' %(ares['s2n'],self['min_s2n'])
            print >>stderr,mess
            return {'flags':AM_FAINT,'ares':ares}

        flags=0
        mag=None

        results={'flags':0,'ares':ares}
        for fitter_type in self['obj_fitters']:
            results[fitter_type] = {}

            for model in self['obj_models']:

                if fitter_type=='lm':
                    fitter=self._fit_object_model_lm(
                            im, mag, skysig, ares, model, psf_gmix)
                elif fitter_type=='mcmc':

                    fitter=self._fit_object_model_mcmc(
                            im, mag, skysig, ares, model, psf_gmix)
                    res=fitter.get_result()
                    print '  model: %s prob: %.6f aic: %.6f bic: %.6f s/n: %.6f Ts/n: %.6f ' % \
                            (model,res['fit_prob'],res['aic'],res['bic'],res['s2n_w'],res['Ts2n'])
                else:
                    raise ValueError("bad fitter type '%s'" % fitter_type)
                results[fitter_type][model] = fitter

        return results

    def _fit_object_model_lm(self, im, mag, skysig, ares, model, psf_gmix):
        from gmix_image.gmix_fit import GMixFitSimple
        ivar=1./skysig**2
        cen_width=self.get('cen_width',1.0)

        # can be None
        gprior=self._get_gprior(mag, model)

        aic=9.999e9
        fitter=None
        for i in xrange(self['obj_lm_ntry']):
            fitter0=GMixFitSimple(im, 
                                  ivar,
                                  psf_gmix, 
                                  model,
                                  ares,
                                  cen_width=cen_width,
                                  gprior=gprior)
            res=fitter0.get_result()
            if res['flags'] != 0:
                continue
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
                            run=self['run'],
                            camcol=self['camcol'],
                            field=field)
        psf_kl = []
        for filter in sdsspy.FILTERCHARS:
            kl = sdsspy.read('psField',
                             run=self['run'],
                             camcol=self['camcol'],
                             field=field,
                             filter=filter)
            psf_kl.append(kl)

        return psfield, psf_kl


    def _read_atlas(self, obj):
        """
        Read atlas images as sky-subtracted, 'f8'
        """
        from sdsspy.atlas import NoAtlasImageError

        field=obj['field']
        id=obj['id']
        try:
            atlas = sdsspy.read('fpAtlas',
                                run=self['run'],
                                camcol=self['camcol'],
                                field=field,
                                id=id)
        except NoAtlasImageError:
            print >>stderr,"object has not atlas image"
            atlas=None
        #for i in xrange(5):
        #    im=atlas['images'][i].astype('f8')-atlas['SOFT_BIAS']
        #    atlas['images'][i] = im
        return atlas

    def _select(self, objs):

        # first trim the photoObj down to what appears in a data sweep
        print 'Trimming to sweep objects'
        sweep_selector=es_sdsspy.select.SweepSelector(objs, 'gal')
        sweep_ind=sweep_selector.get_indices()
        if sweep_ind.size == 0:
            return sweep_ind

        print "Selecting objects"
        s = es_sdsspy.select.Selector(objs[sweep_ind])

        print "  getting resolve logic"
        resolve_logic = s.resolve_logic()

        print "  getting flag logic"
        flag_logic = s.flag_logic()

        print "  getting rmag logic"
        rmag_logic = s.cmodelmag_logic("r", self['max_cmodelmag_r'])

        logic = resolve_logic & flag_logic & rmag_logic

        keep, = numpy.where(logic)

        if keep.size > 0:
            keep=sweep_ind[keep]

        print "  keeping %i/%i" % (keep.size, objs.size)

        return keep

    def _load_data(self):
        objs=sdsspy.read('photoObj', 
                         run=self['run'], 
                         camcol=self['camcol'], 
                         field=self['field'],
                         lower=True,
                         verbose=True)
        if objs is None:
            print >>stderr,"no data read"
            self.objs=objs
            return

        w=self._select(objs)

        if w.size == 0:
            print >>stderr,"no objects passed cuts"
            self.objs=None
            return

        objs=objs[w]
        self.objs=objs

        self.psfield,self.psf_kl = self._read_psfield(self['field'])

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

    def _set_start_time(self):
        self._t0=time.time()

    def _print_time_stats(self):
        tm=time.time()-self._t0
        print 'total time:',tm/60.,'minutes'
        print 'time per object:',tm/self.objs.size

    def _setup_plotting(self):
        import biggles
        biggles.configure("screen","width",1100)
        biggles.configure("screen","height",1100)

    def _get_psf_npars(self):
        if self['psf_model'] == 'gmix3':
            ngauss=3
        elif self['psf_model'] == 'gmix2':
            ngauss=2
        elif self['psf_model'] == 'gmix1':
            ngauss=1

        return 2*ngauss+4

    def _copy_to_output(self, st, i, res):

        st['flags'][i] = res['flags']

        if res['psf_res']['flags'] == 0:
            self._copy_psf_result(st, i, res)

        if 'ares' in res['obj_res']:
            ares=res['obj_res']['ares']
            st['am_s2n'][i] = ares['s2n']

        if res['obj_res']['flags'] == 0:
            self._copy_obj_result(st, i, res)

    def _copy_psf_result(self, st, i, res):
        pres=res['psf_res']['fitter'].get_result()
        st['psf_pars'][i,:] = pres['pars']
        st['psf_perr'][i,:] = pres['perr']
        st['psf_pcov'][i,:,:] = pres['pcov']

        for n in ['flags', 'numiter','loglike','chi2per',
                  'dof','fit_prob','aic','bic']:
            nn='psf_%s' % n
            st[nn][i] = pres[n]

    def _copy_obj_result(self, st, i, res):
        """
        We get here as long as adaptive moments didn't fail
        """
        ores=res['obj_res']

        for fitter_type in self['obj_fitters']:
            for model in self['obj_models']:

                fitter=ores[fitter_type][model]
                if fitter is None:
                    continue

                r=fitter.get_result()

                if fitter_type=='lm':
                    flags=r['flags']
                else:
                    flags=0

                if flags == 0:
                    front='%s_%s' % (fitter_type, model)


                    st[front+'_pars'][i,:] = r['pars']
                    st[front+'_perr'][i,:] = r['perr']
                    st[front+'_pcov'][i,:] = r['pcov']
                    for n in r:
                        if n in ['pars','perr','pcov']:
                            continue

                        if n=='s2n_w':
                            nn='s2n'
                        else:
                            nn=n

                        nn = '%s_%s' % (front,nn)

                        if nn in st.dtype.names:
                            st[nn][i] = r[n]


    def _get_struct(self):

        nobj=self.objs.size

        npars_psf=self._get_psf_npars()
        npars=6

        dt=[('photoid','i8'),
            ('run','i2'),
            ('rerun','i2'),
            ('camcol','i2'),
            ('field','i2'),
            ('id','i2'),
            ('filter','S1'),
            ('fnum','i1'),

            ('flags','i4'), # overall flag indicators
           
            ('psf_pars','f8',npars_psf),
            ('psf_perr','f8',npars_psf),
            ('psf_pcov','f8',(npars_psf,npars_psf)),
            ('psf_flags','i4'),
            ('psf_numiter','i4'),
            ('psf_loglike','f8'),
            ('psf_chi2per','f8'),
            ('psf_dof','i4'),
            ('psf_fit_prob','f8'),
            ('psf_aic','f8'),
            ('psf_bic','f8'),
            
            ('am_s2n','f8')

           ]

        for fitter_type in self['obj_fitters']:
            for model in self['obj_models']:

                front='%s_%s' % (fitter_type, model)

                dt += [(front+'_flags','i4'),
                       (front+'_g','f8',2),
                       (front+'_gcov','f8',(2,2)),
                       (front+'_pars','f8',npars),
                       (front+'_perr','f8',npars),
                       (front+'_pcov','f8',(npars,npars))]

                dt += [(front+'_s2n','f8'),
                       (front+'_loglike','f8'),
                       (front+'_chi2per','f8'),
                       (front+'_dof','i4'),
                       (front+'_fit_prob','f8'),
                       (front+'_aic','f8'),
                       (front+'_bic','f8'),
                       (front+'_Tmean','f8'),
                       (front+'_Terr','f8'),
                       (front+'_Ts2n','f8')
                      ]

                if fitter_type=='lm':
                    dt += [(front+'_numiter','i4')]

                if fitter_type=='mcmc':
                    dt += [(front+'_gsens','f8',2),
                           (front+'_arate','f8')]

        st=numpy.zeros(nobj, dtype=dt)

        for n in st.dtype.names:
            st[n] = -9999
            
        st['photoid']=sdsspy.get_photoid(self.objs)
        st['run'] = self.objs['run']
        st['rerun'] = self.objs['rerun']
        st['camcol'] = self.objs['camcol']
        st['field'] = self.objs['field']
        st['id'] = self.objs['id']

        st['filter'] = self['filter']
        st['fnum'] = self['fnum']

        return st
