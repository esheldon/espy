import os
from sys import stderr,stdout
import time
import numpy

import meds

# starting new values for these
DEFVAL=-9999
PDEFVAL=9999
BIG_DEFVAL=-9.999e9
BIG_PDEFVAL=9.999e9


NO_CUTOUTS=2**0
PSF_FIT_FAILURE=2**1
PSF_LARGE_OFFSETS=2**2
EXP_FIT_FAILURE=2**3
DEV_FIT_FAILURE=2**4

BOX_SIZE_TOO_BIG=2**5

NO_ATTEMPT=2**30

#PSF_S2N=1.e6
PSF_OFFSET_MAX=0.25
PSF_TOL=1.0e-5
EM_MAX_TRY=3
EM_MAX_ITER=100

SIMPLE_MODELS_DEFAULT = ['exp','dev']

_CHECKPOINTS_DEFAULT_MINUTES=[10,30,60,90]

class MedsFit(object):
    def __init__(self, meds_file, truth_file, psf_file, **keys):
        """
        This differs from nfit in gmix_meds in that there is no coadd

        parameters
        ----------
        meds_file: string
            meds path
        truth_file: string
            truth info to get the PSF file
        psf_file: string
            file holding the psf images

        The following are sent through keywords

        fit_types: list of strings
            ['simple']
        obj_range: optional
            a 2-element sequence or None.  If not None, only the objects in the
            specified range are processed and get_data() will only return those
            objects. The range is inclusive unlike slices.
        psf_model: string, int
            e.g. "em2"
        psf_offset_max: optional
            max offset between multi-component gaussians in psf models

        checkpoints: number, optional
            Times after which to checkpoint, seconds
        checkpoint_file: string, optional
            File which will hold a checkpoint.
        checkpoint_data: dict, optional
            The data representing a previous checkpoint, object and
            psf fits
        """

        self.meds_file=meds_file
        self.truth_file=truth_file
        self.psf_file=psf_file

        self._load_meds()
        self._load_truth()
        self._load_psf_fobj()

        self.conf={}
        self.conf.update(keys)

        self.fit_types=self.conf['fit_types']
        self.simple_models=keys.get('simple_models',SIMPLE_MODELS_DEFAULT )

        self.nwalkers=keys.get('nwalkers',80)
        self.burnin=keys.get('burnin',400)
        self.nstep=keys.get('nstep',800)
        self.mca_a=keys.get('mca_a',2.0)

        self._unpack_priors()

        self._setup_checkpoints()


        self.obj_range=keys.get('obj_range',None)
        self._set_index_list()

        self.psf_model=keys.get('psf_model','em2')
        self.psf_offset_max=keys.get('psf_offset_max',PSF_OFFSET_MAX)
        self.psf_ngauss=get_psf_ngauss(self.psf_model)

        self.debug=keys.get('debug',0)

        self.psf_ntry=keys.get('psf_ntry', EM_MAX_TRY)
        self.psf_maxiter=keys.get('psf_maxiter', EM_MAX_ITER)
        self.psf_tol=keys.get('psf_tol', PSF_TOL)

        self.region=keys.get('region','seg_and_sky')
        self.max_box_size=keys.get('max_box_size',2048)

        self.make_plots=keys.get('make_plots',False)
        self.prompt=keys.get('prompt',True)

        if self._checkpoint_data is None:
            self._make_struct()
            self._make_epoch_struct()

    def _unpack_priors(self):
        conf=self.conf

        nmod=len(self.simple_models)
        T_priors=conf['T_priors']
        counts_priors=conf['counts_priors']
        g_priors=conf['g_priors']

        if (len(T_priors) != nmod or len(g_priors) != nmod ):
            raise ValueError("models and T,g priors must be same length")

        priors={}
        models=self.simple_models
        for i in xrange(nmod):
            model=models[i]
            T_prior=T_priors[i]

            # note it is a list
            counts_prior=[ counts_priors[i] ]

            g_prior=g_priors[i]
            
            modlist={'T':T_prior, 'counts':counts_prior,'g':g_prior}
            priors[model] = modlist

        self.priors=priors
        self.draw_g_prior=conf.get('draw_g_prior',True)

        # in arcsec (or units of jacobian)
        self.cen_prior=conf.get("cen_prior",None)


    def get_data(self):
        """
        Get the data structure.  If a subset was requested, only those rows are
        returned.
        """
        return self.data

    def get_epoch_data(self):
        """
        Get the epoch data structure, including psf fitting
        """
        return self.epoch_data

    def get_meds_meta(self):
        """
        get copies of the meta data
        """
        return self.meds_meta.copy()

    def get_magzp(self):
        """
        Get the magnitude zero point.
        """
        return self.meds_meta['magzp_ref'][0]

    def do_fits(self):
        """
        Fit all objects in our list
        """

        t0=time.time()

        last=self.index_list[-1]
        num=len(self.index_list)

        for dindex in xrange(num):
            if self.data['processed'][dindex]==1:
                # was checkpointed
                continue

            mindex = self.index_list[dindex]
            print >>stderr,'index: %d:%d' % (mindex,last),
            self.fit_obj(dindex)

            tm=time.time()-t0

            self._try_checkpoint(tm) # only at certain intervals

        tm=time.time()-t0
        print >>stderr,"time:",tm
        print >>stderr,"time per:",tm/num


    def fit_obj(self, dindex):
        """
        Process the indicated object through the requested fits
        """

        t0=time.time()

        # for checkpointing
        self.data['processed'][dindex]=1

        mindex = self.index_list[dindex]

        # need to do this because we work on subset files
        if 'number' not in self.meds._cat.dtype.names:
            self.data['number'][dindex] = mindex+1
        else:
            self.data['number'][dindex] = self.meds['number'][mindex]

        self.data['flags'][dindex] = self._obj_check(mindex)
        if self.data['flags'][dindex] != 0:
            return 0

        # lists of lists.
        im_list,wt_list= self._get_imlist_wtlist(dindex,mindex)
        jacob_list=self._get_jacobian_list(mindex)

        print >>stderr,im_list[0].shape
    

        print >>stderr,'    fitting: psf models'
        keep_list,psf_gmix_list,flags=self._fit_psfs(dindex,jacob_list)
        if flags != 0:
            self.data['flags'][dindex] = PSF_FIT_FAILURE 
            return


        im_list, wt_list, jacob_list = \
            self._extract_sub_lists(keep_list,im_list,wt_list,jacob_list)

        nim = len(im_list)
        if nim==0:
            print >>stderr,'wierd lists zero length!'
            self.data['flags'][dindex] = PSF_FIT_FAILURE 
            return

        self.data['nimage_use'][dindex] = nim

        sdata={'keep_list':keep_list,
               'im_list':im_list,
               'wt_list':wt_list,
               'jacob_list':jacob_list,
               'psf_gmix_list':psf_gmix_list}

        self._fit_all_models(dindex, sdata)

        self.data['time'][dindex] = time.time()-t0



    def _obj_check(self, mindex):
        """
        Check box sizes, number of cutouts
        """
        flags=0

        box_size=self.meds['box_size'][mindex]
        if box_size > self.max_box_size:
            print >>stderr,'Box size too big:',box_size
            flags |= BOX_SIZE_TOO_BIG

        if self.meds['ncutout'][mindex] < 1:
            print >>stderr,'No cutouts'
            flags |= NO_CUTOUTS
        return flags


    def _get_imlist_wtlist(self, dindex, mindex):
        """
        get lists of lists
        """

        # inherited functions
        imlist=self._get_imlist(mindex)
        wtlist=self._get_wtlist(mindex)

        self.data['nimage_tot'][dindex] = len(imlist)
        
        return imlist,wtlist

    def _get_imlist(self, mindex, type='image'):
        """
        get the image list
        """
        imlist0 = self.meds.get_cutout_list(mindex,type=type)

        imlist = [im.astype('f8') for im in imlist0]
        return imlist


    def _get_wtlist(self, mindex):
        """
        get the weight list.
        """
        if self.region=='seg_and_sky':
            wtlist_all=self.meds.get_cweight_cutout_list(mindex)

            wtlist=[wt.astype('f8') for wt in wtlist_all]
        else:
            raise ValueError("support other region types")

        return wtlist


    def _get_jacobian_list(self, mindex):
        """
        Get a list of the jocobians for this object
        """
        import ngmix
        jlist0=self.meds.get_jacobian_list(mindex)

        jlist=[]
        for jdict in jlist0:
            #print jdict
            j=ngmix.Jacobian(jdict['row0'],
                             jdict['col0'],
                             jdict['dudrow'],
                             jdict['dudcol'],
                             jdict['dvdrow'],
                             jdict['dvdcol'])
            jlist.append(j)

        return jlist


    def _fit_psfs(self,dindex,jacob_list):
        """
        Get images for all SE images and fit
        them to gaussian mixture models
        """

        from ngmix import GMixMaxIterEM

        # guess fwhm=0.9''
        sigma=0.9/2.35

        imlist = self._get_psf_imlist(dindex)

        keep_list=[]
        gmix_list=[]

        flags=0

        gmix_psf=None
        mindex = self.index_list[dindex]

        for i in xrange(len(imlist)):

            im=imlist[i]
            jacob0=jacob_list[i]

            cen0=[ (im.shape[0]-1.0)/2.,
                   (im.shape[1]-1.0)/2. ]

            # the dimensions of the psfs are different, need
            # new center
            jacob=jacob0.copy()
            jacob._data['row0'] = cen0[0]
            jacob._data['col0'] = cen0[1]

            tflags=0
            try:
                # previous result gets re-used
                fitter=self._do_fit_psf(im,jacob,sigma,first_guess=gmix_psf)

                gmix_psf=fitter.get_gmix()
                self._set_psf_result(gmix_psf)

                keep,offset_arcsec=self._should_keep_psf(gmix_psf)
                if keep:
                    gmix_list.append( gmix_psf )
                    keep_list.append(i)
                else:
                    print >>stderr,'large psf offset: '
                    tflags |= PSF_LARGE_OFFSETS 
                    gmix_psf=None

                
            except GMixMaxIterEM:
                print >>stderr,'psf fail',flist[i]

                tflags = PSF_FIT_FAILURE

                gmix_psf=None

            flags |= tflags

            self._set_psf_data(mindex, i, tflags)
            self.psf_index += 1

        return keep_list, gmix_list, flags

    def _set_psf_result(self, gm):
        """
        Set psf fit data. Index can be got from the main model
        fits struct
        """

        psf_index=self.psf_index

        pars=gm.get_full_pars()
        g1,g2,T=gm.get_g1g2T()

        ed=self.epoch_data
        ed['psf_fit_g'][psf_index,0]    = g1
        ed['psf_fit_g'][psf_index,1]    = g2
        ed['psf_fit_T'][psf_index]      = T
        ed['psf_fit_pars'][psf_index,:] = pars

    def _set_psf_data(self, mindex, icut, flags):
        """
        Set all the meta data for the psf result
        """
        psf_index=self.psf_index
        ed=self.epoch_data

        # mindex can be an index into a sub-range meds
        if 'number' not in self.meds._cat.dtype.names:
            ed['number'][psf_index] = mindex+1
        else:
            ed['number'][psf_index] = self.meds['number'][mindex]

        ed['cutout_index'][psf_index] = icut
        ed['file_id'][psf_index]  = self.meds['file_id'][mindex,icut].astype('i4')
        ed['orig_row'][psf_index] = self.meds['orig_row'][mindex,icut]
        ed['orig_col'][psf_index] = self.meds['orig_col'][mindex,icut]
        ed['psf_fit_flags'][psf_index] = flags

    def _get_psf_imlist(self, dindex):
        """
        Get psf images for the SE images
        associated with the cutouts
        """

        mindex = self.index_list[dindex]
        ncut=self.meds['ncutout'][mindex]
        imlist=[]

        for icut in xrange(ncut):
            
            psf_id = self.truth['id_psf'][mindex]
            im = self.psf_fobj[psf_id].read()
            imlist.append( im )

        return imlist

    def _should_keep_psf(self, gm):
        """
        For double gauss we limit the separation
        """
        keep=True
        if self.psf_ngauss == 2:
            offset_arcsec = calc_offset_arcsec(gm)
            if offset_arcsec > self.psf_offset_max:
                keep=False

        return keep, offset_arcsec

    def _do_fit_psf(self, im, jacob, sigma_guess, first_guess=None):
        """
        Fit a single psf
        """
        import ngmix
        from ngmix import GMixMaxIterEM

        s2=sigma_guess**2
        im_with_sky, sky = ngmix.em.prep_image(im)

        fitter=ngmix.em.GMixEM(im_with_sky, jacobian=jacob)

        for i in xrange(self.psf_ntry):

            if i == 0 and first_guess is not None:
                gm_guess=first_guess.copy()
            else:
                s2guess=s2*jacob._data['det'][0]
                gm_guess=self._get_em_guess(s2guess)
            try:
                fitter.go(gm_guess, sky,
                          maxiter=self.psf_maxiter,
                          tol=self.psf_tol)
                break
            except GMixMaxIterEM:
                res=fitter.get_result()
                print >>stderr,'last fit:'
                print >>stderr,fitter.get_gmix()
                print >>stderr,'try:',i+1,'fdiff:',res['fdiff'],'numiter:',res['numiter']
                if i == (self.psf_ntry-1):
                    raise

        return fitter

    def _get_em_guess(self, sigma2):
        """
        Guess for the EM algorithm
        """
        import ngmix
        from ngmix import srandu

        if self.psf_ngauss==1:
            pars=numpy.array( [1.0, 0.0, 0.0, 
                               sigma2*(1.0 + 0.1*srandu()),
                               0.0,
                               sigma2*(1.0 + 0.1*srandu())] )
        else:

            pars=numpy.array( [_em2_pguess[0],
                               0.1*srandu(),
                               0.1*srandu(),
                               _em2_fguess[0]*sigma2*(1.0 + 0.1*srandu()),
                               0.0,
                               _em2_fguess[0]*sigma2*(1.0 + 0.1*srandu()),

                               _em2_pguess[1],
                               0.1*srandu(),
                               0.1*srandu(),
                               _em2_fguess[1]*sigma2*(1.0 + 0.1*srandu()),
                               0.0,
                               _em2_fguess[1]*sigma2*(1.0 + 0.1*srandu())] )

        return ngmix.gmix.GMix(pars=pars)




    def _extract_sub_lists(self,
                           keep_list,
                           im_list0,
                           wt_list0,
                           jacob_list0):
        """
        extract those that passed some previous cuts
        """

        im_list = [im_list0[i] for i in keep_list]
        wt_list = [wt_list0[i] for i in keep_list]
        jacob_list = [jacob_list0[i] for i in keep_list]

        return im_list, wt_list, jacob_list


    def _fit_all_models(self, dindex, sdata):
        """
        Fit psf flux and other models
        """
        self._fit_psf_flux(dindex, sdata)
 
        s2n=self.data['psf_flux'][dindex]/self.data['psf_flux_err'][dindex]
        max_psf_s2n=numpy.nanmax(s2n)
         
        if max_psf_s2n >= self.conf['min_psf_s2n']:
            if 'simple' in self.fit_types:
                self._fit_simple_models(dindex, sdata)
        else:
            mess="    psf s/n too low: %s (%s)"
            mess=mess % (max_psf_s2n,self.conf['min_psf_s2n'])
            print >>stderr,mess

    def _fit_psf_flux(self, dindex, sdata):
        """
        Fit the PSF flux
        """
        import ngmix

        print >>stderr,'    fitting: psf'

        fitter=ngmix.fitting.PSFFluxFitter(sdata['im_list'],
                                           sdata['wt_list'],
                                           sdata['jacob_list'],
                                           sdata['psf_gmix_list'])
        fitter.go()
        res=fitter.get_result()
        self.data['psf_flags'][dindex] = res['flags']
        self.data['psf_flux'][dindex] = res['flux']
        self.data['psf_flux_err'][dindex] = res['flux_err']
        self.data['psf_chi2per'][dindex] = res['chi2per']
        self.data['psf_dof'][dindex] = res['dof']

        mess='         %s +/- %s'
        mess=mess % (self.data['psf_flux'][dindex],
                     self.data['psf_flux_err'][dindex])
        print >>stderr,mess


    def _fit_simple_models(self, dindex, sdata):
        """
        Fit all the simple models
        """

        for model in self.simple_models:
            print >>stderr,'    fitting:',model

            gm=self._fit_simple(dindex, model, sdata)

            res=gm.get_result()

            self._copy_simple_pars(dindex, res)
            self._print_simple_res(res)

            if self.make_plots:
                mindex = self.index_list[dindex]
                ptrials,wptrials,resplots=\
                        gm.make_plots(title='%s multi-epoch' % model,
                                      do_residual=True)
                ptrials.write_img(800,800,
                                  'trials-%06d-%s.png' % (mindex,model))
                wptrials.write_img(800,800,
                                   'trials-%06d-%s-weighted.png' % (mindex,model))
                for i,p in enumerate(resplots):
                    p.write_img(1400,800,
                                'resid-%06d-%s-%02i.png' % (mindex,model,i))

    def _fit_simple(self, dindex, model, sdata):
        """
        Fit one of the "simple" models, e.g. exp or dev
        """

        full_guess=None
        counts_guess=self._get_counts_guess(dindex,sdata)
        T_guess=self._get_T_guess(dindex,sdata)
        
        gm=self._do_fit_simple(model, 
                               sdata['im_list'],
                               sdata['wt_list'],
                               sdata['jacob_list'],
                               sdata['psf_gmix_list'],
                               self.burnin,
                               self.nstep,
                               T_guess=T_guess,
                               counts_guess=counts_guess,
                               full_guess=full_guess)
        return gm


    def _get_counts_guess(self, dindex, sdata):
        """
        Based on the psf flux guess
        """
        psf_flux=self.data['psf_flux'][dindex].clip(min=0.1, max=None)
        return psf_flux

    def _get_T_guess(self, dindex, sdata):
        """
        Guess at T in arcsec**2

        Guess corresponds to FWHM=2.0 arcsec

        Assuming scale is 0.27''/pixel
        """
        return 1.44

    def _do_fit_simple(self, model, im_list, wt_list, jacob_list, psf_gmix_list,
                       burnin,nstep,
                       T_guess=None,
                       counts_guess=None,
                       full_guess=None):
        import ngmix

        priors=self.priors[model]
        g_prior=priors['g']
        T_prior=priors['T']
        counts_prior=priors['counts']

        gm=ngmix.fitting.MCMCSimple(im_list,
                                    wt_list,
                                    jacob_list,
                                    model,
                                    psf=psf_gmix_list,

                                    nwalkers=self.nwalkers,
                                    burnin=burnin,
                                    nstep=nstep,
                                    mca_a=self.mca_a,

                                    iter=True,

                                    T_guess=T_guess,
                                    counts_guess=counts_guess,

                                    full_guess=full_guess,

                                    cen_prior=self.cen_prior,
                                    T_prior=T_prior,
                                    counts_prior=counts_prior,
                                    g_prior=g_prior,
                                    draw_g_prior=self.draw_g_prior,
                                    do_pqr=True)
        gm.go()
        return gm


    def _copy_simple_pars(self, dindex, res):
        """
        Copy from the result dict to the output array
        """
        model=res['model']
        n=get_model_names(model)

        self.data[n['flags']][dindex] = res['flags']

        if res['flags'] == 0:
            pars=res['pars']
            pars_cov=res['pars_cov']

            flux=pars[5:]
            flux_cov=pars_cov[5:, 5:]

            self.data[n['pars']][dindex,:] = pars
            self.data[n['pars_cov']][dindex,:,:] = pars_cov

            self.data[n['flux']][dindex] = flux
            self.data[n['flux_cov']][dindex] = flux_cov

            self.data[n['g']][dindex,:] = res['g']
            self.data[n['g_cov']][dindex,:,:] = res['g_cov']

            self.data[n['arate']][dindex] = res['arate']
            if res['tau'] is not None:
                self.data[n['tau']][dindex] = res['tau']

            for sn in _stat_names:
                self.data[n[sn]][dindex] = res[sn]

            self.data[n['P']][dindex] = res['P']
            self.data[n['Q']][dindex,:] = res['Q']
            self.data[n['R']][dindex,:,:] = res['R']
                


    def _load_meds(self):
        """
        Load all listed meds files
        """

        print >>stderr,self.meds_file
        self.meds = meds.MEDS(self.meds_file)
        self.meds_meta=self.meds.get_meta()
        self.nobj_tot = self.meds.size

    def _load_truth(self):
        """
        load the truth file for getting the psf index
        """
        import fitsio
        print >>stderr,self.truth_file
        self.truth = fitsio.read(self.truth_file,lower=True)
        
    def _load_psf_fobj(self):
        """
        Load the psf file as a FITS object
        """
        import fitsio
        print >>stderr,self.psf_file
        self.psf_fobj = fitsio.FITS(self.psf_file)

    def _set_index_list(self):
        """
        set the list of indices to be processed
        """
        if self.obj_range is None:
            start=0
            end=self.nobj_tot-1
        else:
            start=self.obj_range[0]
            end=self.obj_range[1]

        self.index_list = numpy.arange(start,end+1)


    def _print_simple_res(self, res):
        if res['flags']==0:
            self._print_simple_fluxes(res)
            self._print_simple_T(res)
            self._print_simple_shape(res)
            print >>stderr,'        arate:',res['arate']

    def _print_simple_shape(self, res):
        g1=res['pars'][2]
        g1err=numpy.sqrt(res['pars_cov'][2,2])
        g2=res['pars'][3]
        g2err=numpy.sqrt(res['pars_cov'][3,3])

        print >>stderr,'        g1: %.4g +/- %.4g g2: %.4g +/- %.4g' % (g1,g1err,g2,g2err)

    def _print_simple_fluxes(self, res):
        """
        print in a nice format
        """
        from numpy import sqrt
        flux=res['pars'][5]
        flux_err=sqrt( res['pars_cov'][5, 5] )

        print >>stderr,'        flux: %s +/- %s' % (flux,flux_err)

    def _print_simple_T(self, res):
        """
        print T, Terr, Ts2n and sigma
        """
        T = res['pars'][4]
        Terr = numpy.sqrt( res['pars_cov'][4,4] )

        if Terr > 0:
            Ts2n=T/Terr
        else:
            Ts2n=-9999.0
        if T > 0:
            sigma=numpy.sqrt(T/2.)
        else:
            sigma=-9999.0

        tup=(T,Terr,Ts2n,sigma)
        print >>stderr, '        T: %s +/- %s Ts2n: %s sigma: %s' % tup

    def _setup_checkpoints(self):
        """
        Set up the checkpoint times in minutes and data
        """
        self.checkpoints = self.conf.get('checkpoints',_CHECKPOINTS_DEFAULT_MINUTES)
        self.n_checkpoint    = len(self.checkpoints)
        self.checkpointed    = [0]*self.n_checkpoint
        self.checkpoint_file = self.conf.get('checkpoint_file',None)

        self._set_checkpoint_data()

        if self.checkpoint_file is not None:
            self.do_checkpoint=True
        else:
            self.do_checkpoint=False

    def _set_checkpoint_data(self):
        """
        See if checkpoint data was sent
        """
        self._checkpoint_data=self.conf.get('checkpoint_data',None)
        if self._checkpoint_data is not None:
            self.data=self._checkpoint_data['data']
            self.epoch_data=self._checkpoint_data['epoch_data']

            if self.epoch_data.dtype.names is not None:
                # start where we left off
                w,=numpy.where( self.epoch_data['number'] == -1)
                self.psf_index = w.min()

    def _try_checkpoint(self, tm):
        """
        Checkpoint at certain intervals.  
        Potentially modified self.checkpointed
        """

        should_checkpoint, icheck = self._should_checkpoint(tm)

        if should_checkpoint:
            self._write_checkpoint(tm)
            self.checkpointed[icheck]=1

    def _should_checkpoint(self, tm):
        """
        Should we write a checkpoint file?
        """

        should_checkpoint=False
        icheck=-1

        if self.do_checkpoint:
            tm_minutes=tm/60

            for i in xrange(self.n_checkpoint):

                checkpoint=self.checkpoints[i]
                checkpointed=self.checkpointed[i]

                if tm_minutes > checkpoint and not checkpointed:
                    should_checkpoint=True
                    icheck=i

        return should_checkpoint, icheck

    def _write_checkpoint(self, tm):
        """
        Write out the current data structure to a temporary
        checkpoint file.
        """
        import fitsio

        print >>stderr,'checkpointing at',tm/60,'minutes'
        print >>stderr,self.checkpoint_file

        with fitsio.FITS(self.checkpoint_file,'rw',clobber=True) as fobj:
            fobj.write(self.data, extname="model_fits")
            fobj.write(self.epoch_data, extname="epoch_data")

    def _count_all_cutouts(self):
        """
        Count the cutouts for the objects
        """
        ncutout = self.meds['ncutout'][self.index_list].sum()
        return ncutout


    def _make_epoch_struct(self):
        """
        We will make the maximum number of possible psfs according
        to the cutout count
        """

        npars=self.psf_ngauss*6
        dt=[('number','i4'), # 1-n as in sextractor
            ('cutout_index','i4'), # this is the index into e.g. m['orig_row'][3,index]
            ('orig_row','f8'),
            ('orig_col','f8'),
            ('file_id','i4'),
            ('psf_fit_flags','i4'),
            ('psf_fit_g','f8',2),
            ('psf_fit_T','f8'),
            ('psf_fit_pars','f8',npars)]

        ncutout=self._count_all_cutouts()
        if ncutout > 0:
            epoch_data = numpy.zeros(ncutout, dtype=dt)

            epoch_data['number'] = -1
            epoch_data['cutout_index'] = -1
            epoch_data['file_id'] = -1
            epoch_data['psf_fit_g'] = PDEFVAL
            epoch_data['psf_fit_T'] = PDEFVAL
            epoch_data['psf_fit_pars'] = PDEFVAL
            epoch_data['psf_fit_flags'] = NO_ATTEMPT

            self.epoch_data=epoch_data
        else:
            self.epoch_data=numpy.zeros(1)

        # where the next psf data will be written
        self.psf_index = 0


    def _make_struct(self):
        """
        make the output structure
        """

        dt=[('number','i4'),
            ('processed','i1'),
            ('flags','i4'),
            ('nimage_tot','i4'),
            ('nimage_use','i4'),
            ('time','f8')]

        n=get_model_names('psf')
        dt += [(n['flags'],   'i4'),
               (n['flux'],    'f8'),
               (n['flux_err'],'f8'),
               (n['chi2per'],'f8'),
               (n['dof'],'f8')]
       
        if 'simple' in self.fit_types:
    
            simple_npars=6
            simple_models=self.simple_models

            for model in simple_models:
                n=get_model_names(model)

                np=simple_npars


                dt+=[(n['flags'],'i4'),
                     (n['pars'],'f8',np),
                     (n['pars_cov'],'f8',(np,np)),
                     (n['flux'],'f8'),
                     (n['flux_cov'],'f8'),
                     (n['g'],'f8',2),
                     (n['g_cov'],'f8',(2,2)),
                    
                     (n['s2n_w'],'f8'),
                     (n['chi2per'],'f8'),
                     (n['dof'],'f8'),
                     (n['aic'],'f8'),
                     (n['bic'],'f8'),
                     (n['arate'],'f8'),
                     (n['tau'],'f8'),
                    ]

                dt += [(n['P'], 'f8'),
                       (n['Q'], 'f8', 2),
                       (n['R'], 'f8', (2,2))]


        num=self.index_list.size
        data=numpy.zeros(num, dtype=dt)

        data['psf_flags'] = NO_ATTEMPT
        data['psf_flux'] = DEFVAL
        data['psf_flux_err'] = PDEFVAL

        if 'simple' in self.fit_types:
            for model in simple_models:
                n=get_model_names(model)

                data[n['flags']] = NO_ATTEMPT

                data[n['pars']] = DEFVAL
                data[n['pars_cov']] = PDEFVAL
                data[n['flux']] = DEFVAL
                data[n['flux_cov']] = PDEFVAL
                data[n['g']] = DEFVAL
                data[n['g_cov']] = PDEFVAL

                data[n['s2n_w']] = DEFVAL
                data[n['chi2per']] = PDEFVAL
                data[n['aic']] = BIG_PDEFVAL
                data[n['bic']] = BIG_PDEFVAL

                data[n['tau']] = PDEFVAL

                data[n['P']] = DEFVAL
                data[n['Q']] = DEFVAL
                data[n['R']] = DEFVAL

     
        self.data=data


_psf_ngauss_map={'em1':1, 'em2':2}
def get_psf_ngauss(psf_model):
    if psf_model not in _psf_ngauss_map:
        raise ValueError("bad psf model: '%s'" % psf_model)
    return _psf_ngauss_map[psf_model]



def calc_offset_arcsec(gm, scale=1.0):
    data=gm.get_data()

    offset=numpy.sqrt( (data['row'][0]-data['row'][1])**2 + 
                       (data['col'][0]-data['col'][1])**2 )
    offset_arcsec=offset*scale
    return offset_arcsec


_em2_fguess=numpy.array([0.5793612389470884,1.621860687127999])
_em2_pguess=numpy.array([0.596510042804182,0.4034898268889178])
#_em2_fguess=numpy.array([12.6,3.8])
#_em2_fguess[:] /= _em2_fguess.sum()
#_em2_pguess=numpy.array([0.30, 0.70])


_stat_names=['s2n_w',
             'chi2per',
             'dof',
             'aic',
             'bic']


def get_model_names(model):
    names=['rfc_flags',
           'rfc_tries',
           'rfc_iter',
           'rfc_pars',
           'rfc_pars_cov',
           'flags',
           'pars',
           'pars_cov',
           'flux',
           'flux_err',
           'flux_cov',
           'g',
           'g_cov',
           'g_sens',
           'e',
           'e_cov',
           'e_sens',
           'P',
           'Q',
           'R',
           'iter',
           'tries',
           'arate',
           'tau']
    names += _stat_names

    ndict={}
    for n in names:
        ndict[n] = '%s_%s' % (model,n)

    return ndict


