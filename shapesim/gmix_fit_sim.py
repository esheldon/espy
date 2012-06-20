"""
Generate image simulations and process them with the
gmix fitting pipeline
"""

import os
import numpy
from numpy import random, zeros, sqrt, array, ceil, isfinite, \
        where, diag, arctan2, median
from numpy.random import random as randu
from numpy.linalg import eig, LinAlgError
import sys
from sys import stderr
from lensing.util import e2gamma, e1e2_to_g1g2
from . import shapesim
from fimage import mom2sigma, cov2sigma
from pprint import pprint 
import copy
import images
import esutil as eu
from esutil.misc import wlog

try:
    import gmix_image
    from gmix_image import pars2gmix_coellip
    from gmix_image import GMIXFIT_SINGULAR_MATRIX
    from gmix_image import GMIXFIT_NEG_COV_EIG
    from gmix_image import GMIXFIT_NEG_COV_DIAG

    from gmix_image import ellip2eta, eta2ellip, print_pars
except ImportError:
    wlog("could not import gmix_image")

try:
    import admom
except ImportError:
    wlog("could not import admom")

class GMixFitSim(shapesim.BaseSim):
    """
    We only override

        .run()
        .out_dtype()
        .copy_output()

    """
    def __init__(self, run):
        super(GMixFitSim,self).__init__(run)
        if 'verbose' not in self:
            self['verbose'] = False

    def run(self, ci):
        """
        Process the input convolved image

        Output will be a dict with
        --------------------------
        flags:
            Flags of the last processing
        psf_res:
            Result of psf processing.
        res:
            Result of image processing, if psf processing succeeded.
        """
        show=False
        dostop=True
        out={}

        #images.multiview(ci.psf)
        #stop
        ptype=self.get('ptype','cov')

        coellip_psf=self['coellip_psf']
        coellip_obj=self['coellip_obj']

        out['psf_res'] = self.process_image(ci.psf, 
                                            self['ngauss_psf'],
                                            ci['cen_psf_admom'],
                                            ci['cov_psf_admom'],
                                            skysig=ci['skysig_psf'],
                                            coellip=coellip_psf)
        out['flags'] = out['psf_res']['flags']
        if out['flags'] == 0:
            cov_admom = ci['cov_admom']
            cov_psf_admom = ci['cov_psf_admom']
            Tadmom=cov_admom[0]+cov_admom[2]
            Tpsf_admom = cov_psf_admom[0]+cov_psf_admom[2]
            T0admom = Tadmom-Tpsf_admom
            Radmom = T0admom/Tadmom
            out['res'] = self.process_image(ci.image, 
                                            self['ngauss_obj'],
                                            ci['cen_admom'],
                                            cov_admom,
                                            psf=out['psf_res']['gmix'],
                                            skysig=ci['skysig'],
                                            Radmom=Radmom,
                                            coellip=coellip_obj)
            out['flags'] = out['res']['flags']
            if show and out['flags'] == 0:
                self.show_residual(ci, out['psf_res']['gmix'], 
                                   objmix=out['res']['gmix'],
                                   dostop=dostop)
            elif show:
                self.show_residual(ci, out['psf_res']['gmix'],dostop=dostop)

            wlog("\ne1true:",ci['e1true'],"e2true:",ci['e2true'])
        else:
            if self['verbose']:
                wlog('failed PSF flags:')
                gmix_image.printflags(ptype,out['flags'])

        if out['flags'] != 0 and self['verbose']:
            wlog('flags:')
            gmix_image.printflags(ptype,out['flags'])
        return out

    def process_image(self, image, ngauss, cen, cov, psf=None,
                      Radmom=None,
                      skysig=None, coellip=True):
        if not coellip:
            raise ValueError("must use coellip for now")

        counts = image.sum()
        method=self.get('method','lm')

        # will only apply when psf is present
        randomize=self.get('randomize',False)
        if psf:
            verbose=False
            maxtry=self.get('maxtry',1)
        else:
            verbose=False
            maxtry=1
        ntry=0
        chi2arr=zeros(maxtry) + 1.e9
        chi2perarr=zeros(maxtry) + 1.e9
        gmlist=[]

        ptype = self.get('ptype','cov')
        use_jacob=self.get('use_jacob',True)
        wlog("ptype:",ptype)
        wlog("use_jacob:",use_jacob)
        while ntry < maxtry:
            if ptype == 'cov':
                guess = self.get_guess_coellip_cov(counts, ngauss, cen, cov, 
                                                   randomize=randomize,
                                                   psf=psf)
            elif ptype == 'e1e2':
                if Radmom is not None:
                    Tfac=None
                    if ntry==0:
                        eguess=None
                    if ntry == 1:
                        eguess=[0,0]
                else:
                    #Tfac=1.
                    Tfac=None
                    eguess=None

                    # first try we guess the admom ellip
                    if ntry == 1:
                        eguess=[0,0]
                    elif ntry == 2:
                        Tfac = 5.
                    elif ntry == 3:
                        eguess=[0.,0.]
                        Tfac = 5.
                    elif ntry == 4:
                        eguess=[0.3,0.3]
                        Tfac = 5.
                    elif ntry == 5:
                        eguess=[-0.3,-0.3]
                        Tfac = 5.

                    elif ntry == 6:
                        Tfac = .2
                    elif ntry == 7:
                        eguess=[0.,0.]
                        Tfac = .2
                    elif ntry == 8:
                        eguess=[0.3,0.3]
                        Tfac = .2
                    elif ntry == 9:
                        eguess=[-0.3,-0.3]
                        Tfac = .2

                    elif ntry == 10:
                        eguess=[0.0,0.0]
                        Tfac = .05  
                    elif ntry == 11:
                        eguess=[0.3,0.3]
                        Tfac = .05  
                    elif ntry == 12:
                        eguess=[-0.3,-0.3]
                        Tfac = .05
                    else:
                        eguess=None
                guess = self.get_guess_coellip_e1e2(counts, ngauss, cen, cov, 
                                                    randomize=randomize,
                                                    psf=psf,
                                                    Radmom=Radmom,
                                                    eguess=eguess,
                                                    Tfac=Tfac)
                #if psf:
                #    print_pars(guess,front="guess: ")
                #    stop
            else:
                raise ValueError("ptype should be 'cov','e1e2'")

            gm = gmix_image.GMixFitCoellip(image,guess,
                                           psf=psf,
                                           method=method,
                                           ptype=ptype,
                                           use_jacob=use_jacob,
                                           verbose=verbose)

            print_pars(gm.popt,front="pars:  ")
            if skysig is not None:
                chi2arr[ntry] = gm.chi2(gm.popt)
                chi2perarr[ntry] = gm.chi2per(gm.popt,skysig)
                wlog("chi2/pdeg:",chi2perarr[ntry])
            else:
                chi2arr[ntry] = gm.chi2(gm.popt)
                wlog("chi2:",chi2arr[ntry])
            gmlist.append(gm)

            ntry += 1
                
        #if psf:
        #    stop
        w=chi2arr.argmin()
        gm = gmlist[w]

        print_pars(gm.popt,front='popt: ')
        print_pars(gm.perr,front='perr: ')
        #wlog('chi2arr:',chi2arr)
        if skysig is not None:
            #wlog('chi2arr/perdeg:',chi2perarr)
            print_pars(chi2perarr,front='chi2/deg: ')
        else:
            print_pars(chi2arr,front='chi2: ')
        wlog("numiter gmix:",gm.numiter)
        wlog("poptT/Tguess:",gm.popt[ngauss+4:]/guess[ngauss+4:])

        out={'gmix':    gm.gmix,
             'pars':    gm.popt,
             'perr':    gm.perr,
             'pcov':    gm.pcov,
             'flags':   gm.flags,
             'ier':     gm.ier,
             'numiter': gm.numiter,
             'coellip': coellip,
             'chi2per': chi2perarr[w]}
        return out

    def get_guess_coellip_e1e2(self, counts, ngauss, cen, cov, 
                               randomize=False, 
                               eguess=None,
                               Radmom=None,
                               Tfac=1.,
                               psf=None):
        wlog("\nusing coellip e1e2")
        npars = 2*ngauss+4
        guess=zeros(npars)
        guess[0] = cen[0]
        guess[1] = cen[1]

        T = cov[2]+cov[0]
        e1=(cov[2]-cov[0])/T
        e2=2*cov[1]/T

        guess[2] = e1
        guess[3] = e2

        wlog("Tadmom:",T)
        if ngauss==1:
            wlog("    using ngauss==1")

            guess[4] = 1.0
            guess[5] = T

            if psf is not None:
                wlog("======> with psf")
                psfmoms = gmix_image.total_moms(psf)
                tcov=cov.copy()
                tcov[0] -= psfmoms['irr']
                tcov[1] -= psfmoms['irc']
                tcov[2] -= psfmoms['icc']
                tdet=tcov[0]*tcov[2] - tcov[1]**2
                if tdet > 1.e-5:
                    wlog("using special guesses")
                    tT = tcov[0]+tcov[2]
                    te1=(tcov[2]-tcov[0])/tT
                    te2=2*tcov[1]/tT

                    guess[2] = te1
                    guess[3] = te2
                    guess[4] = 1.0
                    guess[5] = tT
                else:
                    # use defaults
                    wlog("NOT USING special guesses")
                    pass
            if randomize:
                guess[0] += 1*(randu()-0.5)  # cen0
                guess[1] += 1*(randu()-0.5)  # cen1
                guess[2] += 0.2*(randu()-0.5)  # e1
                guess[3] += 0.2*(randu()-0.5)  # e2
                guess[4] += 0.1*(randu()-0.5)  # p
                guess[5] += 1*(randu()-0.5)   # T

        elif ngauss==3:
                wlog("    using ngauss=3")
            
                if eguess is not None:
                    guess[2],guess[3] = eguess
                else:
                    guess[2] = e1# + 0.05*(randu()-0.5)
                    guess[3] = e2# + 0.05*(randu()-0.5)
                    """
                    while abs(guess[2]) > 0.95:
                      guess[2] = e1 + 0.05*(randu()-0.5)
                    while abs(guess[3]) > 0.95:
                      guess[3] = e2 + 0.05*(randu()-0.5)
                    """

                wlog("    starting e1,e2:",guess[2],guess[3])
                guess[4] = 0.62
                guess[5] = 0.34
                guess[6] = 0.04
                guess[7] = T*2.7
                guess[8] = T*0.62
                guess[9] = T*0.09

        elif ngauss==4:
                wlog("    using ngauss=4")

                guess[4] = 1./ngauss
                guess[5] = 1./ngauss
                guess[6] = 1./ngauss
                guess[7] = 1./ngauss

                if eguess is not None:
                    guess[2],guess[3] = eguess
                else:
                    guess[2] = e1# + 0.05*(randu()-0.5)
                    guess[3] = e2# + 0.05*(randu()-0.5)

                wlog("    starting e1,e2:",guess[2],guess[3])

                if Tfac is None and Radmom > 0.01:
                    # from gmix_fit_sim.plot_admom_max_tratio
                    wlog("Using Radmom fit")
                    ply=numpy.poly1d([ 72.51258096,   4.14321833])
                    tratio = ply(Radmom)

                    # this T is Tadmom
                    Tmax = tratio*T

                    # 1.000000 0.182911 0.035527 0.002705
                    guess[8] = Tmax
                    guess[9] = Tmax*0.18
                    guess[10] = Tmax*0.035
                    guess[11] = Tmax*0.0027
                else:
                    if Tfac is None:
                        Tfac=1
                    # Tfac is to help cover the cases where
                    # the object is much smaller or much larger
                    # than the PSF

                    """
                    guess[4] = 0.21*Tfac
                    guess[5] = 0.27*Tfac
                    guess[6] = 0.25*Tfac
                    guess[7] = 0.24*Tfac
                    """
                    """
                    guess[8] = T*4.5*Tfac
                    guess[9] = T*.8*Tfac
                    guess[10] = T*.18*Tfac
                    guess[11] = T*.04*Tfac
                    """
                    guess[8] = T*2.6*Tfac
                    guess[9] = T*.38*Tfac
                    guess[10] = T*.076*Tfac
                    guess[11] = T*.015*Tfac
        else:
            raise RuntimeError("implement other guesses!")
 
        return guess

    def get_guess_coellip_cov(self, counts, ngauss, cen, cov, 
                              randomize=False, psf=None):
        npars = 2*ngauss+4
        guess=zeros(npars)
        guess[0] = cen[0]
        guess[1] = cen[1]

        guess[2] = cov[0]
        guess[3] = cov[1]
        guess[4] = cov[2]

        """
        if psf is not None:
            psfmoms = gmix_image.total_moms(psf)
            tcov=cov.copy()
            tcov[0] -= psfmoms['irr']
            tcov[1] -= psfmoms['irc']
            tcov[2] -= psfmoms['icc']
            tdet=tcov[0]*tcov[2] - tcov[1]**2
            T = tcov[0]+tcov[2]
            if tdet > 0:
                print 'tcov:',tcov
                guess[2:2+3] = tcov
        """ 
        # If psf sent, this is an object. If ngauss==3, 
        # make guesses good for an exp disk
        if psf is not None and ngauss == 3:
            if len(psf) == 3:
                #pvals=array([0.419696,0.0725887,0.499471])
                pvals = array([1./ngauss]*3,dtype='f8')
                if randomize:
                    pvals += 0.1*(randu(3)-0.5)
                pvals *= counts/pvals.sum()
                guess[5:5+3] = pvals
                       
                guess[8] = 0.3 + 0.01*(randu()-0.5)
                guess[9] = 3.5 + 1.*(randu()-0.5)
                guess[2:2+3] += 0.1*(randu(3)-0.5)
                if randomize:
                    guess[8] += 0.01*(randu()-0.5)
                    guess[9] += 1.*(randu()-0.5)
                    guess[2:2+3] += 0.1*(randu(3)-0.5)
            else:
                pvals = array([1./ngauss]*3,dtype='f8')
                if randomize:
                    pvals += 0.4*(randu(3)-0.5)
                pvals *= counts/pvals.sum()
                guess[5:5+3] = pvals
                       
                # try: - more trials and broader sampleing
                #      - a cutoff on fi at the high and low end?
                #      - priors?
                #      - no randomizing pvals guess?
                #      - randomize cov with single number, or not at all?
                guess[8] = 0.7# + 0.4*(random.random()-0.5)
                guess[9] = 3.0# + 3.*(random.random()-0.5)

                if randomize:
                    guess[8] += .6*(randu()-0.5)
                    guess[9] += 3.*(randu()-0.5)

                # these (good) guesses seem to force us to go through
                # some local (wrong) minimum
                #guess[5] = 0.4
                #guess[6] = 0.065
                #guess[7] = 0.54
                #guess[8] = 0.2
                #guess[9] = 3.8

        elif ngauss == 4:
            """
            guess[5] = .25 + 0.1*(randu()-0.5)
            guess[6] = .35 + 0.1*(randu()-0.5)
            guess[7] = .25 + 0.1*(randu()-0.5)
            guess[8] = .15 + 0.1*(randu()-0.5)

            guess[9] = .10 + 0.1*(randu()-0.5) # start on "inside"
            guess[10] = .22 + 0.1*(randu()-0.5)
            guess[11] = 6. + 1*(randu()-0.5)
            """

            # for low res case; this works, so the guess really
            # does matter, even in no noise case
            #guess[5] = .16
            #guess[6] = .5
            #guess[7] = .245
            #guess[8] = .08
            guess[5] = 1./ngauss
            guess[6] = 1./ngauss
            guess[7] = 1./ngauss
            guess[8] = 1./ngauss

            guess[9] = .05
            guess[10] = .25
            guess[11] = 5.2

            if randomize:
                # these move together
                guess[2:2+3] += 1.*(randu()-0.5)

                guess[5] += .1*(randu()-0.5)
                guess[6] += .1*(randu()-0.5)
                guess[7] += .1*(randu()-0.5)
                guess[8] += .1*(randu()-0.5)

                guess[9]  += .1*randu() # only upward 
                guess[10] += .1*(randu()-0.5)
                guess[11] += 1.*(randu()-0.5)

        else:
            # generic guesses
            guess[5:5+ngauss] = counts/ngauss

            if ngauss > 1:
                if ngauss == 2:
                    guess[5+ngauss] = 3.0
                elif ngauss == 3:
                    guess[5+ngauss] = 0.5
                    guess[5+ngauss+1] = 3.0
                else:
                    # 4 or mor
                    guess[5+ngauss] = 0.5
                    guess[5+ngauss+1] = 3.0
                    guess[5+ngauss+2:] = 4.0

        return guess

    def show_residual(self, ci, psfmix, dostop=True, objmix=None):
        """
        Show plots of the input compared with the fit gaussian mixtures.
        """
        
        psfmodel = gmix_image.gmix2image(psfmix,ci.psf.shape,
                                         renorm=False) 
        images.compare_images(ci.psf,psfmodel,
                              label1='psf',label2='gmix')

        if objmix is not None:
            skysig=None
            if ci['skysig'] > 0:
                skysig=ci['skysig']
            #model0 = gmix_image.gmix2image(objmix,ci.image0.shape,
            #                               renorm=False) 
            model = gmix_image.gmix2image(objmix,ci.image.shape,
                                          psf=psfmix,
                                          renorm=False) 

            gmix_image.gmix_print(objmix)
            for i,g in enumerate(objmix):
                wlog('g',i)
                tim=gmix_image.gmix2image([g],ci.image0.shape)
                w=where(isfinite(tim) == False)
                if w[0].size > 0:
                    wlog('found NaN')
                wlog(tim)
            #images.compare_images(ci.image0,model0,
            #                      label1='object0',label2='gmix',
            #                      skysig=skysig)
            images.compare_images(ci.image,model,
                                  label1='object+psf',label2='gmix',
                                  skysig=skysig)
        if dostop:
            stop
        else:
            key=raw_input('hit a key (q to quit): ')
            if key == 'q':
                stop

    def copy_output(self, s2, ellip, s2n, ci, res):
        if not self['coellip_psf'] or not self['coellip_obj']:
            raise ValueError("must use coellip for psf/obj for now")

        st = zeros(1, dtype=self.out_dtype())

        # first copy inputs and data from the CI
        st['s2'] = s2
        st['s2n'] = s2n
        st['ellip'] = ellip
        st['e1true'] = ci['e1true']
        st['e2true'] = ci['e2true']
        st['etrue']  = ci['etrue']
        st['gamma'] = e2gamma(st['etrue'])
        st['gamma1'],st['gamma2'] = e1e2_to_g1g2(st['e1true'],st['e2true'])

        st['irr_uw'] = ci['cov_uw'][0]
        st['irc_uw'] = ci['cov_uw'][1]
        st['icc_uw'] = ci['cov_uw'][2]

        st['irr_psf_uw'] = ci['cov_psf_uw'][0]
        st['irc_psf_uw'] = ci['cov_psf_uw'][1]
        st['icc_psf_uw'] = ci['cov_psf_uw'][2]

        st['e1_uw'] = ci['e1_image0_uw']
        st['e2_uw'] = ci['e2_image0_uw']
        st['e_uw']  = ci['e_image0_uw']

        size2psf = ci['cov_psf_uw'][0]+ci['cov_psf_uw'][2]
        size2obj = ci['cov_image0_uw'][0]+ci['cov_image0_uw'][2]
        st['s2_uw'] = size2psf/size2obj

        s2psf_am = ci['cov_psf_admom'][0]+ci['cov_psf_admom'][2]
        s2obj_am = ci['cov_image0_admom'][0]+ci['cov_image0_admom'][2]
        st['s2admom'] = s2psf_am/s2obj_am
        st['T_psf_admom'] = ci['cov_psf_admom'][0]+ci['cov_psf_admom'][2]
        st['T0_admom'] = ci['cov_image0_admom'][0]+ci['cov_image0_admom'][2]
        st['T_admom'] = ci['cov_admom'][0]+ci['cov_admom'][2]

        if 'psf_res' in res:
            st['pars_psf']     = res['psf_res']['pars']
            st['pars_psf_err'] = res['psf_res']['perr']
            st['pars_psf_cov'] = res['psf_res']['pcov']

            ng = len(res['psf_res']['gmix'])
            for i,r in enumerate(res['psf_res']['gmix']):
                for k in ['p','row','col','irr','irc','icc']:
                    if ng == 1:
                        st['gmix_psf'][k][0] = r[k] 
                    else:
                        st['gmix_psf'][k][0,i] = r[k] 

            psf_moms = gmix_image.total_moms(res['psf_res']['gmix'])
            st['irr_psf_meas'] = psf_moms['irr']
            st['irc_psf_meas'] = psf_moms['irc']
            st['icc_psf_meas'] = psf_moms['icc']
            st['sigma_psf_meas'] = 0.5*(psf_moms['irr']+psf_moms['icc'])

            st['numiter_psf'] = res['psf_res']['numiter']

            st['chi2per_psf'] = res['psf_res']['chi2per']
        if 'res' in res:
            st['pars']     = res['res']['pars']
            st['pars_err'] = res['res']['perr']
            st['pars_cov'] = res['res']['pcov']

            ng = len(res['res']['gmix'])
            for i,r in enumerate(res['res']['gmix']):
                for k in ['p','row','col','irr','irc','icc']:
                    if ng == 1:
                        st['gmix'][k][0] = r[k] 
                    else:
                        st['gmix'][k][0,i] = r[k] 

            moms = gmix_image.total_moms(res['res']['gmix'])
            st['irr_meas'] = moms['irr']
            st['irc_meas'] = moms['irc']
            st['icc_meas'] = moms['icc']
            st['s2_meas'] = \
                (psf_moms['irr']+psf_moms['icc'])/(moms['irr']+moms['icc'])
            st['sigma_meas'] = 0.5*(moms['irr']+moms['icc'])

            st['e1_meas'] = (moms['icc']-moms['irr'])/(moms['icc']+moms['irr']) 
            st['e2_meas'] = 2*moms['irc']/(moms['icc']+moms['irr']) 
            st['e_meas'] = sqrt(st['e1_meas']**2 + st['e2_meas']**2)


            st['gamma_meas'] = e2gamma(st['e_meas'])
            st['gamma1_meas'],st['gamma2_meas'] = \
                    e1e2_to_g1g2(st['e1_meas'],st['e2_meas'])

            st['flags'] = res['res']['flags']
            st['numiter'] = res['res']['numiter']

            mcov = res['res']['pcov'][2:2+3,2:2+3]
            st['e1_chol_err'] = numpy.inf
            st['e2_chol_err'] = numpy.inf
            st['e_chol_err'] = numpy.inf
            """
            try:
                e1, e1_err, e2, e2_err, e, e_err = \
                    get_ellip_cholesky(res['res']['pars'][2:2+3], 
                                       res['res']['pcov'][2:2+3,2:2+3])
                st['e1_chol'] = e1
                st['e1_chol_err'] = e1_err
                st['e2_chol'] = e2
                st['e2_chol_err'] = e2_err
                st['e_chol'] = e
                st['e_chol_err'] = e_err
            except LinAlgError as e:
                wlog(e)
                st['e1_chol_err'] = numpy.inf
                st['e2_chol_err'] = numpy.inf
                st['e_chol_err'] = numpy.inf
                stop
            """
            st['chi2per'] = res['res']['chi2per']
        else:
            st['s2_meas'] = -9999



        # figure out how to measure this
        st['s2n_meas'] = st['s2n']


        return st


    def out_dtype(self):
        npars_psf = 2*self['ngauss_psf']+4
        npars_obj = 2*self['ngauss_obj']+4

        gmix_dt = [('p','f8'),('row','f8'),('col','f8'),
                   ('irr','f8'),('irc','f8'),('icc','f8')]
        dt=[('s2n','f8'),
            ('ellip','f8'),

            ('s2','f8'),         # requested (spsf/sobj)**2
            ('s2_uw','f8'), # unweighted s2 of object before noise
            ('T_psf_admom','f8'),
            ('T_admom','f8'),
            ('T0_admom','f8'),
            ('s2admom','f8'),    # s2 from admom, generally different

            ('irr_uw','f8'),
            ('irc_uw','f8'),
            ('icc_uw','f8'),
            ('irr_psf_uw','f8'),
            ('irc_psf_uw','f8'),
            ('icc_psf_uw','f8'),
            ('e_uw','f8'),
            ('e1_uw','f8'),
            ('e2_uw','f8'),

            ('etrue','f8'),
            ('e1true','f8'),
            ('e2true','f8'),
            ('gamma','f8'),
            ('gamma1','f8'),
            ('gamma2','f8'),

            ('numiter','i8'),
            ('numiter_psf','i8'),

            ('flags','i8'),

            ('s2n_meas','f8'),    # use admom s2n

            ('Tadmom','f8'),      # of convolved object

            ('s2_meas','f8'),
            ('irr_psf_meas','f8'),
            ('irc_psf_meas','f8'),
            ('icc_psf_meas','f8'),
            ('irr_meas','f8'),
            ('irc_meas','f8'),
            ('icc_meas','f8'),
            ('sigma_meas','f8'),
            ('sigma_psf_meas','f8'),
            ('e_meas','f8'),
            ('e1_meas','f8'),
            ('e2_meas','f8'),

            ('e_chol','f8'),
            ('e_chol_err','f8'),
            ('e1_chol','f8'),
            ('e1_chol_err','f8'),
            ('e2_chol','f8'),
            ('e2_chol_err','f8'),

            ('gamma_meas','f8'),
            ('gamma1_meas','f8'),
            ('gamma2_meas','f8'),

            ('chi2per_psf','f8'),
            ('gmix_psf',gmix_dt,self['ngauss_psf']),
            ('pars_psf','f8',npars_psf),
            ('pars_psf_err','f8',npars_psf),
            ('pars_psf_cov','f8',(npars_psf,npars_psf)),

            ('chi2per','f8'),
            ('gmix',gmix_dt,self['ngauss_obj']),
            ('pars','f8',npars_obj),
            ('pars_err','f8',npars_obj),
            ('pars_cov','f8',(npars_obj,npars_obj)) ]

        return dt


def get_ellip_cholesky(means, cov, n=100000):
    r = eu.stat.cholesky_sample(cov, n, means=means)

    Tr = r[2,:]+r[0,:]
    e1r = (r[2,:]-r[0,:])/Tr
    e2r = 2*r[1,:]/Tr
    er = sqrt(e1r**2 + e2r**2)

    e1 = e1r.mean()
    e1_err = e1r.std()

    e2 = e2r.mean()
    e2_err = e2r.std()

    e = er.mean()
    e_err = er.std()

    return e1, e1_err, e2, e2_err, e, e_err

def plot_admom_max_tratio(run, ei=0):
    """
    Plot the ratio of the maximum moment to the adaptive moment
    """
    import glob

    d=shapesim.get_output_dir(run)
    pattern='%s-*-%03i.rec' % (run,ei)
    pattern=os.path.join(d,pattern)

    wlog(pattern)
    flist = glob.glob(pattern)
    nf=len(flist)

    Radmoms = zeros(nf)
    Tratios = zeros(nf)

    for j,f in enumerate(flist):
        wlog(f)
        t=eu.io.read(f)

        # had a typo!  Is fixed now, and saving T instead of sigma
        Tam = 2*t['sigma0_admom']**2
        Tam_psf = 2*t['sigma_psf_admom']**2

        Tmax = zeros(len(t))
        for i in xrange(len(t)):
            Tmax[i] = t['pars'][i,8:].max()
        
        Radmoms[j] = median((Tam-Tam_psf)/Tam)
        Tratios[j] = median(Tmax/Tam)

    w,=where(Tratios > 10)
    pfit = numpy.polyfit(Radmoms[w], Tratios[w], 1)
    p = numpy.poly1d(pfit)
    print pfit

    plt=eu.plotting.bscatter(Radmoms,Tratios,show=False,
                             xlabel=r'$R_{AM}$',
                             ylabel=r'$T/T_{AM}$')
    eu.plotting.bscatter(Radmoms[w],p(Radmoms[w]),type='solid',color='red',
                         plt=plt)

