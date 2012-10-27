"""
Generate image simulations and process them with the
gmix fitting pipeline
"""

import os
from math import copysign
import numpy
from numpy import random, zeros, sqrt, array, ceil, isfinite, \
        where, diag, arctan2, median, poly1d
from numpy.random import random as randu
from numpy.random import randn
from numpy.linalg import eig, LinAlgError
import sys
from sys import stderr
from lensing.util import e2gamma, e1e2_to_g1g2
from . import shapesim
from fimage import mom2sigma, cov2sigma, mom2ellip
from pprint import pprint 
import copy
import images
import esutil as eu
from esutil.misc import wlog
import scipy.stats

import fimage

try:
    import gmix_image
    from gmix_image import GMIXFIT_SINGULAR_MATRIX
    from gmix_image import GMIXFIT_NEG_COV_EIG
    from gmix_image import GMIXFIT_NEG_COV_DIAG

    from gmix_image import ellip2eta, eta2ellip, print_pars, pars2gmix, gmix2pars
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

        coellip_psf=self['coellip_psf']
        coellip_obj=self['coellip_obj']

        out['psf_res'] = self.process_image(ci.psf, ci['skysig_psf'],
                                            self['ngauss_psf'],
                                            ci['cen_psf_admom'],
                                            ci['cov_psf_admom'],
                                            coellip=coellip_psf)


        out['flags'] = out['psf_res']['flags']
        if out['flags'] == 0:
            cov_admom = ci['cov_admom']
            cov_psf_admom = ci['cov_psf_admom']
            Tadmom=cov_admom[0]+cov_admom[2]
            Tpsf_admom = cov_psf_admom[0]+cov_psf_admom[2]
            pfrac_am = Tpsf_admom/Tadmom

            pfrac_am2 = get_admom_pfrac(ci)

            out['res'] = self.process_image(ci.image, ci['skysig'],
                                            self['ngauss_obj'],
                                            ci['cen_admom'],
                                            cov_admom,
                                            psf=out['psf_res']['pars'],
                                            psf_gmix=out['psf_res']['gmix'],
                                            pfrac_am=pfrac_am,
                                            coellip=coellip_obj)


            out['flags'] = out['res']['flags']
            if show and out['flags'] == 0:
                self.show_residual(ci, out['psf_res']['gmix'], 
                                   objmix=out['res']['gmix'],
                                   dostop=dostop)
            elif show:
                self.show_residual(ci, out['psf_res']['gmix'],dostop=dostop)

            if out['flags'] == 0:
                mess=("e1true: %.6g e1: %.6g +/- %.6g diff: %.6g\n"
                      "e2true: %.6g e2: %.6g +/- %.6g diff: %.6g")
                mess = mess % (ci['e1true'], 
                               out['res']['pars'][2],
                               out['res']['perr'][2],
                               out['res']['pars'][2]-ci['e1true'],
                               ci['e2true'], 
                               out['res']['pars'][3],
                               out['res']['perr'][3],
                               out['res']['pars'][3]-ci['e2true'])
                self.wlog(mess)
        else:
            self.wlog('failed PSF flags:')
            if self['verbose']:
                gmix_image.printflags("flags:",out['flags'])

        if out['flags'] != 0 and self['verbose']:
            self.wlog('flags:')
            if self['verbose']:
                gmix_image.printflags("flags:",out['flags'])
        return out


    def process_image(self, image, skysig, ngauss, cen, cov, 
                      psf=None,psf_gmix=None,
                      pfrac_am=None, coellip=True):
        if not coellip:
            raise ValueError("must use coellip for now")

        counts = image.sum()

        randomize=self['randomize']

        verbose=False

        if psf is not None:
            maxtry=self['maxtry']
        else:
            maxtry=self['maxtry_psf']

        ntry=0
        chi2perarr=zeros(maxtry) + 1.e9
        guess_chi2perarr=zeros(maxtry) + 1.e9

        gmlist=[]


        generic = self['generic_prior']
        use_nlsolver = self.get("use_nlsolver",False)
        while ntry < maxtry:
            if psf is not None:
                if generic:
                    uniform_p=False
                    eguess=None
                else:
                    eguess=[0,0]
                    uniform_p=False
                    if ntry > 0:
                        uniform_p=True
            else:
                uniform_p=False
                if ntry == 0:
                    eguess=[0,0]
                else:
                    eguess=None

                # overriding, can't remember why
                if ngauss==1:
                    randomize=True
                #else:
                #    randomize=False

                if ngauss==3:
                    eguess=[0,0]


            if generic:
                admom_mult = 4.
                guess,width = self.get_prior_generic(ngauss,
                                                     counts, cen, cov, admom_mult,
                                                     eguess=eguess,
                                                     psf=psf,
                                                     uniform_p=uniform_p,
                                                     randomize=randomize)
            else:
                guess,width = self.get_prior(ngauss,
                                             counts, cen, cov, pfrac_am,
                                             eguess=eguess,
                                             psf=psf,
                                             uniform_p=uniform_p,
                                             randomize=randomize)

            if self['verbose']:
                print_pars(guess,front="guess: ")

            if use_nlsolver:
                maxiter=2000
                if psf is not None:
                    # full representation 6*ngauss
                    psf_pars=gmix2pars(psf_gmix)
                else:
                    psf_pars=None
                gm = gmix_image.gmix_nlsolve.GMixCoellipSolver(image, 
                                                               guess,
                                                               1.,
                                                               maxiter,
                                                               psf_pars,
                                                               False)


            else:
                Tpositive=True
                if psf is not None and ngauss > 1:
                    Tpositive=False
                gm = gmix_image.GMixFitCoellip(image, skysig,
                                               guess,width,
                                               psf=psf_gmix,
                                               verbose=verbose,
                                               Tpositive=Tpositive)


            pars=gm.get_pars()
            perr=gm.get_perr()
            chi2per = gm.get_chi2per()
            if use_nlsolver:
                perr *= skysig
                chi2per /= skysig**2

            if self['verbose']:
                print_pars(pars,front="pars:  ")
                print_pars(perr,front="perr:  ")

            chi2perarr[ntry] = chi2per
            self.wlog("chi2/pdeg:",chi2perarr[ntry])
            gmlist.append(gm)

            ntry += 1
                
        w=chi2perarr.argmin()

        gm = gmlist[w]
        if self['verbose']:
            print_pars(chi2perarr,front='chi2/deg: ')
        s2n = gm.get_s2n()

        dof=gm.get_dof()
        chi2per=chi2perarr[w]
        prob = scipy.stats.chisqprob(chi2per*dof, dof)


        popt=gm.get_pars()
        pcov=gm.get_pcov()
        perr=gm.get_perr()
        gmix=pars2gmix(popt, coellip=True)

        if use_nlsolver:
            perr *= skysig
            pcov *= skysig**2
            s2n /= skysig

        if self['verbose']:
            wlog("\n")

            print_pars(popt,front='popt: ')
            print_pars(perr,front='perr: ')
            if psf is not None:
                wlog("chi^2/deg:",chi2per,"prob:",prob)
            wlog("s2n:",s2n)
            #wlog("numiter gmix:",gm.numiter)
            wlog("Topt/Tguess:",popt[ngauss+4:]/guess[ngauss+4:])

        if not hasattr(gm,'ier'):
            ier=0
            numiter=0
        else:
            ier=gm.ier
            numiter=gm.get_numiter()
            if self['verbose']:
                wlog("numiter:",numiter)
        out={'gmix':    gmix,
             'pars':    popt,
             'perr':    perr,
             'pcov':    pcov,
             'flags':   gm.get_flags(),
             'ier':     ier,
             'numiter': numiter,
             'coellip': coellip,
             's2n':     s2n,
             'chi2per': chi2per,
             'dof':dof,
             'prob': prob}
        return out


    def get_prior(self, ngauss, counts, cen, cov, pfrac_am,
                  psf=None,
                  randomize=False,
                  uniform_p=False,
                  eguess=None):

        npsf=0
        if psf is not None:
            npsf = gmix_image.get_ngauss_coellip(psf)

        tight_priors=self.get('tight_priors',False)
        npars=2*ngauss+4
        prior=zeros(npars)
        width=zeros(npars)
        prior[0] = cen[0]
        prior[1] = cen[1]
        width[0] = 1
        width[1] = 1

        T = cov[2]+cov[0]
        e1=(cov[2]-cov[0])/T
        e2=2*cov[1]/T

        self.wlog("pfrac_am:",pfrac_am)

        if ngauss==4:

            if eguess:
                prior[2],prior[3]=eguess
            else:
                prior[2],prior[3]=e1,e2

            width[2] = 10
            width[3] = 10

            if psf is not None:
                if npsf == 3:
                    #Tmax = T*55
                    #p0,p1,p2,p3=(0.0251877,0.0381688,0.0979805,0.835861)
 
                    #p0_poly = poly1d([0.3177065 -0.65823369 0.36093455])
                    #p1_poly = poly1d([0.72778915 -1.5281274 0.82178094])
                    #p2_poly = poly1d([-1.1024044 1.0738875 0.016592872])
                    #p3_poly = poly1d([-0.17724095 1.5283721 -0.38947674])

                    Trat_poly = poly1d([-2032.5693, 3281.6972, -1269.4376])
                    p0_poly = poly1d([0.3177065, -0.65823369, 0.36093455])
                    p1_poly = poly1d([0.72778915, -1.5281274, 0.82178094])
                    p2_poly = poly1d([-1.1024044, 1.0738875, 0.016592872])
                    p3_poly = poly1d([-0.17724095, 1.5283721, -0.38947674])
                    T1_poly = poly1d([2.0163789, -6.9140866, 6.7484268, 
                                      -1.7787715])
                    T2_poly = poly1d([-0.41175713, 0.85397488, -0.56542694, 
                                      0.153802])
                    T3_poly = poly1d([0.0093353597, -0.022581187, 0.015612458, 
                                      -0.0020975941])

                    p0 = p0_poly(pfrac_am)
                    p1 = p1_poly(pfrac_am)
                    p2 = p2_poly(pfrac_am)
                    p3 = p3_poly(pfrac_am)

                    Tfrac1 = T1_poly(pfrac_am)
                    Tfrac2 = T2_poly(pfrac_am)
                    Tfrac3 = T3_poly(pfrac_am)

                    # test forcing it
                    force=55
                    if False:
                        # try the single gauss psf case
                        #Tply=poly1d([-72.51258096,  76.65579929])
                        #tratio = Tply(pfrac_am)
                        #p0_poly = poly1d([-0.24772161,  0.22211086])
                        #p1_poly = poly1d([-0.50867256,  0.48633848])
                        #p2_poly = poly1d([-0.50787093,  0.6208675 ])
                        #p3_poly = poly1d([ 1.28199182, -0.34435804])
                        tratio=55
                        p0,p1,p2,p3=0.0838067,0.146986,0.201162,0.558154
                        Tmax = tratio*T
                        Tfrac1 = .18
                        Tfrac2 = .035
                        Tfrac3 = .0027
                    else:
                        self.wlog("forcing",force)
                        Tmax = T*force
                elif npsf == 1:
                    # this works quite well for a single gaussian psf
                    Tply=poly1d([-72.51258096,  76.65579929])
                    # from low ellip
                    p0_poly = poly1d([-0.24772161,  0.22211086])
                    p1_poly = poly1d([-0.50867256,  0.48633848])
                    p2_poly = poly1d([-0.50787093,  0.6208675 ])
                    p3_poly = poly1d([ 1.28199182, -0.34435804])

                    tratio = Tply(pfrac_am)
                    Tmax = tratio*T
                    p0 = p0_poly(pfrac_am)
                    p1 = p1_poly(pfrac_am)
                    p2 = p2_poly(pfrac_am)
                    p3 = p3_poly(pfrac_am)
                    Tfrac1 = .18
                    Tfrac2 = .035
                    Tfrac3 = .0027

            else:
                raise ValueError("need psf for now")
            #prior[4:4+4] = counts/ngauss
            prior[4] = Tmax
            prior[5] = Tfrac1
            prior[6] = Tfrac2
            prior[7] = Tfrac3

            if uniform_p:
                self.wlog("    uniform p")
                prior[8] = counts/ngauss
                prior[9] = counts/ngauss
                prior[10] = counts/ngauss
                prior[11] = counts/ngauss
            else:
                prior[8] = p0
                prior[9] = p1
                prior[10] = p2
                prior[11] = p3

            if npsf==3:
                if tight_priors:
                    width[4] = .01                # informative
                    #width[4] = 1.e-5                # informative
                else:
                    width[4] = 5                # informative
            else:
                width[4] = 100              # uninformative
            if tight_priors:
                wfac=.01
                #wfac=1.e-5
            else:
                wfac=0.2
            width[5] = prior[5]*wfac
            width[6] = prior[6]*wfac
            width[7] = prior[7]*wfac

            width[8] = prior[8]*wfac
            width[9] = prior[9]*wfac
            width[10] = prior[10]*wfac
            width[11] = prior[11]*wfac

            if randomize:
                e1start=prior[2]
                e2start=prior[3]
                prior[2],prior[3] = randomize_e1e2(e1start,e2start)

                prior[4] += prior[4]*0.05*(randu()-0.5)
                prior[5] += prior[5]*0.05*(randu()-0.5)
                prior[6] += prior[6]*0.05*(randu()-0.5)
                prior[7] += prior[7]*0.05*(randu()-0.5)
                prior[8] += prior[8]*0.05*(randu()-0.5)
                prior[9] += prior[9]*0.05*(randu()-0.5)
                prior[10] += prior[10]*0.05*(randu()-0.5)
                prior[11] += prior[11]*0.05*(randu()-0.5)

        elif ngauss==3 and psf is None:
            # these guesses are for turbulent
            self.wlog("    using ngauss=3")

            if eguess is not None:
                prior[2],prior[3] = eguess
            else:
                prior[2] = e1# + 0.05*(randu()-0.5)
                prior[3] = e2# + 0.05*(randu()-0.5)

            Tmax = T*8.3
            Tfrac1 = 1.7/8.3
            Tfrac2 = 0.8/8.3
            prior[4] = Tmax
            prior[5] = Tfrac1
            prior[6] = Tfrac2

            if uniform_p:
                self.wlog("    uniform p")
                prior[7] = counts/ngauss
                prior[8] = counts/ngauss
                prior[9] = counts/ngauss
            else:
                prior[7] = counts*0.08
                prior[8] = counts*0.38
                prior[9] = counts*0.53

            # uninformative priors for PSF, might want to revisit for real stars
            width[2] = 10
            width[3] = 10
            width[4] = 100
            width[5:] = 10

            if randomize:
                self.wlog("    randomizing")
                e1start=prior[2]
                e2start=prior[3]
                prior[2],prior[3] = randomize_e1e2(e1start,e2start)

                prior[4] += prior[4]*0.05*(randu()-0.5)
                prior[5] += prior[5]*0.05*(randu()-0.5)
                prior[6] += prior[6]*0.05*(randu()-0.5)
                prior[7] += prior[7]*0.05*(randu()-0.5)
                prior[8] += prior[8]*0.05*(randu()-0.5)
                prior[9] += prior[9]*0.05*(randu()-0.5)

        elif ngauss==3 and psf is not None:
            self.wlog("    using psf ngauss=3")
        
            if eguess is not None:
                prior[2],prior[3] = eguess 
            else:
                prior[2] = e1
                prior[3] = e2

            self.wlog("    starting e1,e2:",prior[2],prior[3])


            #self.wlog("Using pfrac_am fit:",pfrac_am)
            # need to do this at higher S/N
            #ply=poly1d([-3.20824373,  3.40727954])
            #tratio = ply(pfrac_am)

            # this T is Tadmom
            #Tmax = tratio*T
            #Tfrac1 = 0.35
            #Tfrac2 = 0.07
            
            Tmax = T*5
            #Tmax = T*8
            Tfrac1 = .3
            Tfrac2 = .06
            prior[4] = Tmax
            prior[5] = Tfrac1
            prior[6] = Tfrac2

            # these ratios are important when there is noise
            #prior[7] = counts*0.14
            #prior[8] = counts*0.53
            #prior[9] = counts*0.33

            if uniform_p:
                self.wlog("    uniform p")
                prior[7] = counts/ngauss
                prior[8] = counts/ngauss
                prior[9] = counts/ngauss
            else:
                prior[7] = counts*0.26
                prior[8] = counts*0.55
                prior[9] = counts*0.18


            # uninformative
            width[2] = 10
            width[3] = 10
            width[4] = 100 # Tmax
            width[5:] = 10 # Ti/pi

            if randomize:
                self.wlog("    randomizing")
                e1start=prior[2]
                e2start=prior[3]
                prior[2],prior[3] = randomize_e1e2(e1start,e2start)

                prior[4] += prior[4]*0.05*(randu()-0.5)
                prior[5] += prior[5]*0.05*(randu()-0.5)
                prior[6] += prior[6]*0.05*(randu()-0.5)
                prior[7] += prior[7]*0.05*(randu()-0.5)
                prior[8] += prior[8]*0.05*(randu()-0.5)
                prior[9] += prior[9]*0.05*(randu()-0.5)


        elif ngauss==1:
            self.wlog("    using ngauss==1")

            prior[4] = T
            prior[5] = counts

            # uninformative
            width[2] = 10
            width[3] = 10
            width[4] = 100 # Tmax
            width[5] = 10  # p


            if psf is not None:
                self.wlog("======> with psf")
                psf_gmix = gmix_image.pars2gmix(psf,coellip=True)
                psfmoms = gmix_image.total_moms(psf_gmix)
                tcov=cov.copy()
                tcov[0] -= psfmoms['irr']
                tcov[1] -= psfmoms['irc']
                tcov[2] -= psfmoms['icc']
                tdet=tcov[0]*tcov[2] - tcov[1]**2
                if tdet > 1.e-5:
                    self.wlog("using special guesses")
                    tT = tcov[0]+tcov[2]
                    te1=(tcov[2]-tcov[0])/tT
                    te2=2*tcov[1]/tT

                    prior[2] = te1
                    prior[3] = te2
                    prior[4] = tT
                    prior[5] = counts
                else:
                    # use defaults
                    self.wlog("NOT USING special guesses")
                    pass
            if randomize:
                prior[0] += 1*(randu()-0.5)  # cen0
                prior[1] += 1*(randu()-0.5)  # cen1
                prior[2] += 0.2*(randu()-0.5)  # e1
                prior[3] += 0.2*(randu()-0.5)  # e2
                prior[4] += 0.1*(randu()-0.5)  # p
                prior[5] += 1*(randu()-0.5)   # T

        elif ngauss==2 and psf is None:
            self.wlog("    ngauss:",ngauss)

            if eguess is not None:
                prior[2],prior[3] = eguess
            else:
                prior[2] = e1# + 0.05*(randu()-0.5)
                prior[3] = e2# + 0.05*(randu()-0.5)

            e1start,e2start = prior[2],prior[3]
            prior[2],prior[3] = randomize_e1e2(e1start,e2start)
            prior[4] = T*5 # Tmax
            prior[5] = 1./5.08
            prior[6]= counts*0.1
            prior[7]= counts*0.9

            prior[4] += 0.1*prior[4]*(randu()-0.5)
            prior[5] += 0.1*prior[5]*(randu()-0.5)
            prior[6] += 0.1*prior[6]*(randu()-0.5)
            prior[7] += 0.1*prior[7]*(randu()-0.5)

            # uninformative
            width[2] = 10  # e1
            width[3] = 10  # e2
            width[4] = 100 # Tmax
            width[5] = 10  # Tfrac2
            width[6] = 10  # p1
            width[7] = 10  # p2
 
        else:
            raise ValueError("implement other guesses")

        return prior, width


    def get_prior_generic(self, ngauss, counts, cen, cov, admom_mult,
                          psf=None,
                          randomize=False,
                          uniform_p=False,
                          eguess=None):
        """
        admom_mult is used for the prepsf fit ngauss==3
        """

        npsf=0
        if psf is not None:
            npsf = gmix_image.get_ngauss_coellip(psf)

        npars=2*ngauss+4

        prior=zeros(npars)
        width=zeros(npars) + 1.e20

        prior[0] = cen[0]
        prior[1] = cen[1]

        T = cov[2]+cov[0]
        e1=(cov[2]-cov[0])/T
        e2=2*cov[1]/T

        if psf is not None:
            if ngauss == 3:
                prior[2],prior[3] = randomize_e1e2(0., 0.)

                # should really be 0.08, 0.0 for dev
                Tfracs = array([0.3, 0.0])

                prior[4] = T*admom_mult
                prior[5] = .3 # should be ~0.08 for dev galaxies
                prior[6] = 0.
                
                # should be ~.4, .45, .05 or something for dev
                prior[7] = 0.26
                prior[8] = 0.55
                prior[9] = 0.18

                width[6] = 1.e-8
                prior[6] = width[6]*randn()

                rand_fac=0.2
                prior[4] += prior[4]*rand_fac*(randu()-0.5)
                prior[5] += prior[5]*rand_fac*(randu()-0.5)
                prior[7] += prior[7]*rand_fac*(randu()-0.5)
                prior[8] += prior[8]*rand_fac*(randu()-0.5)
                prior[9] += prior[9]*rand_fac*(randu()-0.5)


            elif ngauss==1:
                psf_gmix = gmix_image.pars2gmix(psf,coellip=True)
                psfmoms = gmix_image.total_moms(psf_gmix)
                tcov=cov.copy()
                tcov[0] -= psfmoms['irr']
                tcov[1] -= psfmoms['irc']
                tcov[2] -= psfmoms['icc']
                tdet=tcov[0]*tcov[2] - tcov[1]**2
                if tdet > 1.e-5:
                    self.wlog("using special guesses")
                    tT = tcov[0]+tcov[2]
                    te1=(tcov[2]-tcov[0])/tT
                    te2=2*tcov[1]/tT

                    prior[2] = te1
                    prior[3] = te2
                    prior[4] = 1.0
                    prior[5] = tT
                else:
                    # use defaults
                    self.wlog("NOT USING special guesses")
                    pass

                if randomize:
                    prior[0] += 1*(randu()-0.5)  # cen0
                    prior[1] += 1*(randu()-0.5)  # cen1
                    prior[2] += 0.2*(randu()-0.5)  # e1
                    prior[3] += 0.2*(randu()-0.5)  # e2
                    prior[4] += 0.1*(randu()-0.5)  # p
                    prior[5] += 1*(randu()-0.5)   # T

            else:
                raise ValueError("only 1 or 3 gauss prepsf for now")
        else:
            if ngauss==1:
                prior[0] += 0.02*(randu()-0.5)  # cen0
                prior[1] += 0.02*(randu()-0.5)  # cen1
                prior[2],prior[3] = randomize_e1e2(e1,e2)
                prior[4] = T + T*0.1*(randu()-0.5)
                prior[5] = counts + 0.1*counts*(randu()-0.5)
            elif ngauss==2:
                self.wlog("    ngauss:",ngauss)

                if eguess is not None:
                    prior[2],prior[3] = eguess
                else:
                    prior[2] = e1# + 0.05*(randu()-0.5)
                    prior[3] = e2# + 0.05*(randu()-0.5)

                e1start,e2start = prior[2],prior[3]
                prior[2],prior[3] = randomize_e1e2(e1start,e2start)
                prior[4] = T*5 # Tmax
                prior[5] = 1./5.08
                prior[6]= counts*0.1
                prior[7]= counts*0.9

                prior[4] += 0.1*prior[4]*(randu()-0.5)
                prior[5] += 0.1*prior[5]*(randu()-0.5)
                prior[6] += 0.1*prior[6]*(randu()-0.5)
                prior[7] += 0.1*prior[7]*(randu()-0.5)

            elif ngauss==3:

                if eguess is not None:
                    prior[2],prior[3] = eguess
                else:
                    prior[2] = e1# + 0.05*(randu()-0.5)
                    prior[3] = e2# + 0.05*(randu()-0.5)

                Tmax = T*8.3
                Tfrac1 = 1.7/8.3
                Tfrac2 = 0.8/8.3
                prior[4] = Tmax
                prior[5] = Tfrac1
                prior[6] = Tfrac2

                if uniform_p:
                    self.wlog("    uniform p")
                    prior[7] = counts/ngauss
                    prior[8] = counts/ngauss
                    prior[9] = counts/ngauss
                else:
                    prior[7] = counts*0.08
                    prior[8] = counts*0.38
                    prior[9] = counts*0.53

                if randomize:
                    self.wlog("    randomizing")
                    e1start=prior[2]
                    e2start=prior[3]
                    prior[2],prior[3] = randomize_e1e2(e1start,e2start)

                    prior[4] += prior[4]*0.05*(randu()-0.5)
                    prior[5] += prior[5]*0.05*(randu()-0.5)
                    prior[6] += prior[6]*0.05*(randu()-0.5)
                    prior[7] += prior[7]*0.05*(randu()-0.5)
                    prior[8] += prior[8]*0.05*(randu()-0.5)
                    prior[9] += prior[9]*0.05*(randu()-0.5)
 
            else:
                raise ValueError("implement other psf guesses")

        return prior, width


    def show_residual(self, ci, psfmix, dostop=True, objmix=None):
        """
        Show plots of the input compared with the fit gaussian mixtures.
        """
        
        psfmodel = gmix_image.gmix2image(psfmix,ci.psf.shape) 
        images.compare_images(ci.psf,psfmodel,
                              label1='psf',label2='gmix')

        if objmix is not None:
            skysig=None
            if ci['skysig'] > 0:
                skysig=ci['skysig']
            model = gmix_image.gmix2image(objmix,ci.image.shape, psf=psfmix) 

            gmix_image.gmix_print(objmix)
            for i,g in enumerate(objmix):
                self.wlog('g',i)
                tim=gmix_image.gmix2image([g],ci.image0.shape)
                w=where(isfinite(tim) == False)
                if w[0].size > 0:
                    self.wlog('found NaN')
                self.wlog(tim)
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
        if not self['coellip_psf']:
            raise ValueError("must use coellip for psf for now")

        st = zeros(1, dtype=self.out_dtype())

        # first copy inputs and data from the CI
        st['s2'] = s2
        st['s2n_uw'] = ci['s2n_uw']
        st['s2n_matched'] = ci['s2n_matched']
        st['s2n_admom'] = ci['s2n_admom']
        st['s2n_uw_psf'] = ci['s2n_uw_psf']
        st['s2n_matched_psf'] = ci['s2n_matched_psf']
        st['s2n_admom_psf'] = ci['s2n_admom_psf']
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
        #wlog(" -- s2_uw:",st['s2_uw'],"goal:",st['s2'])

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
            st['sigma_psf_meas'] = sqrt(0.5*(psf_moms['irr']+psf_moms['icc']))

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

            self.wlog(" -- s2_meas:",st['s2_meas'],"goal:",st['s2'])
            self.wlog(" -- psf sigma:",sqrt((psf_moms['irr']+psf_moms['icc'])/2))
            self.wlog(" -- obj sigma:",sqrt((moms['irr']+moms['icc'])/2))

            st['sigma_meas'] = sqrt(0.5*(moms['irr']+moms['icc']))

            st['e1_meas'] = (moms['icc']-moms['irr'])/(moms['icc']+moms['irr']) 
            st['e2_meas'] = 2*moms['irc']/(moms['icc']+moms['irr']) 
            st['e_meas'] = sqrt(st['e1_meas']**2 + st['e2_meas']**2)


            st['gamma_meas'] = e2gamma(st['e_meas'])
            st['gamma1_meas'],st['gamma2_meas'] = \
                    e1e2_to_g1g2(st['e1_meas'],st['e2_meas'])

            st['s2n_meas'] = res['res']['s2n']
            st['flags'] = res['res']['flags']
            st['numiter'] = res['res']['numiter']

            """
            st['e1_chol_err'] = numpy.inf
            st['e2_chol_err'] = numpy.inf
            st['e_chol_err'] = numpy.inf
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
            st['prob'] = res['res']['prob']
            st['dof'] = res['res']['dof']
        else:
            st['s2_meas'] = -9999





        return st


    def out_dtype(self):

        ngauss_psf=self['ngauss_psf']
        ngauss_obj=self['ngauss_obj']
        npars_psf = 2*ngauss_psf+4
        npars_obj = 2*ngauss_obj+4

        gmix_dt = [('p','f8'),('row','f8'),('col','f8'),
                   ('irr','f8'),('irc','f8'),('icc','f8')]
        dt=[('s2n_uw','f8'),
            ('s2n_uw_psf','f8'),
            ('s2n_matched','f8'),
            ('s2n_admom','f8'),
            ('s2n_matched_psf','f8'),
            ('s2n_admom_psf','f8'),
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

            #('e_chol','f8'),
            #('e_chol_err','f8'),
            #('e1_chol','f8'),
            #('e1_chol_err','f8'),
            #('e2_chol','f8'),
            #('e2_chol_err','f8'),

            ('gamma_meas','f8'),
            ('gamma1_meas','f8'),
            ('gamma2_meas','f8'),

            ('chi2per_psf','f8'),
            ('gmix_psf',gmix_dt,ngauss_psf),
            ('pars_psf','f8',npars_psf),
            ('pars_psf_err','f8',npars_psf),
            ('pars_psf_cov','f8',(npars_psf,npars_psf)),

            ('chi2per','f8'),
            ('dof','f8'),
            ('prob','f8'),
            ('gmix',gmix_dt,ngauss_obj),
            ('pars','f8',npars_obj),
            ('pars_err','f8',npars_obj),
            ('pars_cov','f8',(npars_obj,npars_obj)) ]

        return dt


class GMixGalSim(dict):
    """
    Only works for high s/n right now
    """
    def __init__(self, run, profile, ellip, res):
        self['run']=run
        self['profile']=profile
        self['ellip']=int(ellip)
        self['resolution']=int(res)

        self.load_config()


    def load_config(self):
        f=self.get_config_file()
        wlog(f)
        c=eu.io.read(f)

        sim_name='%s%s' % (self['run'],self['profile'])
        if c['name'] != sim_name:
            raise ValueError("name in config does not match "
                             "itself: '%s' instead of '%s'" % (c['name'],sim_name))

        for k,v in c.iteritems():
            self[k] = v

    def get_config_dir(self):
        dir=os.environ['ESPY_DIR']
        dir=os.path.join(dir,'shapesim/galsim-config')
        return dir
    def get_config_file(self):
        import yaml
        dir=self.get_config_dir()
        return os.path.join(dir,'gmix-%s%s.yaml' % (self['run'],self['profile']))

    def get_run_dir(self):
        return os.path.expanduser('~/galsim/%s' % self['run'])

    def get_image_dir(self):
        dir=self.get_run_dir()
        profile=self['profile'].replace('p','')
        dir=os.path.join(dir, 'profile%s' % profile)
        return dir
    def get_image_basename(self):
        return 'ellip%02ihlr%02i.fits' % (self['ellip'],self['resolution'])

    def get_image_file(self):
        dir=self.get_image_dir()
        bname=self.get_image_basename()
        return os.path.join(dir,bname)

    def get_psf_file(self):
        dir=self.get_image_dir()
        return os.path.join(dir,'psf.fits')

    def get_output_dir(self):
        dir=os.path.expanduser('~/galsim-outputs')
        profile=self['profile'].replace('p','')
        dir=os.path.join(dir,'%s/profile%s' % (self['run'],profile))
        return dir


    def get_output_file(self):
        dir=self.get_output_dir()
        bname=self.get_image_basename()
        bname=bname.replace('.fits','-gmix.fits')
        return os.path.join(dir,bname)
    def read_output(self):
        f=self.get_output_file()
        print f
        return eu.io.read(f)

    def load_data(self):
        im_file=self.get_image_file()
        psf_file=self.get_psf_file()

        wlog(im_file)
        image=eu.io.read(im_file)
        self.image = array(image,dtype='f8',copy=False)
        wlog(str(self.image.shape))

        wlog(psf_file)
        psf=eu.io.read(psf_file)
        self.psf=array(psf,dtype='f8',copy=False)
        wlog(str(self.psf.shape))

    def run(self):
        self.load_data()
        ntot=self['nobj_row']*self['nobj_col']
        output=zeros(ntot,dtype=self.get_dtype())

        ii=0
        for orow in xrange(self['nobj_row']):
            for ocol in xrange(self['nobj_col']):
                if self['verbose']:
                    wlog('-'*70)
                wlog(orow,ocol)
                tres=self.run_obj(orow, ocol)
                
                output[ii] = self.get_output(orow,ocol,tres)
                
                ii+=1
                #self.write_output(output)
                #stop
        self.write_output(output)



    def run_obj_old(self, orow, ocol):
        """
        Process all objects in the image and psf
        """

        rows_per=self.image.shape[0]/self['nobj_row']
        cols_per=self.image.shape[1]/self['nobj_col']

        row1=orow*rows_per
        row2=(orow+1)*rows_per

        col1=ocol*cols_per
        col2=(ocol+1)*cols_per


        image = self.image[row1:row2, col1:col2].copy()
        psf0 = self.psf[row1:row2, col1:col2].copy()
        psf,skysig_psf = fimage.add_noise_admom(psf0,self['s2n_psf'])
        self['skysig_psf']=skysig_psf

        cenrow=image.shape[0]/2.
        cencol=image.shape[1]/2.

        psf_admom=self.run_admom(psf,cenrow,cencol,skysig_psf,4.0)
        obj_admom=self.run_admom(image,cenrow,cencol,self['skysig_obj'],6.0)

        cen_psf_admom=array([psf_admom['row'],psf_admom['col']])
        cov_psf_admom = array([psf_admom['Irr'],psf_admom['Irc'],psf_admom['Icc']])
        cen_admom = array([obj_admom['row'], obj_admom['col']])
        cov_admom = array([obj_admom['Irr'],obj_admom['Irc'],obj_admom['Icc']])

        out={}
        out['psf_res'] = self.process_image(psf, self['skysig_psf'],
                                            self['ngauss_psf'],
                                            cen_psf_admom, cov_psf_admom)
        out['psf_res']['s2n_admom'] = psf_admom['s2n']


        out['flags'] = out['psf_res']['flags']
        if out['flags'] == 0:
            Tadmom=cov_admom[0]+cov_admom[2]
            Tpsf_admom = cov_psf_admom[0]+cov_psf_admom[2]

            psf_gmix=pars2gmix(out['psf_res']['pars'], coellip=True)
            out['res'] = self.process_image(image, self['skysig_obj'],
                                            self['ngauss_obj'],
                                            cen_admom,
                                            cov_admom,
                                            psf=psf_gmix)

            out['res']['s2n_admom'] = obj_admom['s2n']

            out['flags'] = out['res']['flags']

            if out['flags'] != 0:
                wlog('failed PSF flags:')

        else:
            wlog('failed PSF flags:')

        # this is a ring test, we cannot have failures
        if out['flags'] != 0:
            gmix_image.printflags("flags:",out['flags'])
            raise RuntimeError("error, halting")
        return out

    def run_admom(self, im, rowguess, colguess, skysig, Tguess):
        res=admom.admom(im, rowguess, colguess, 
                        sigsky=skysig, guess=Tguess/2, nsub=1)  
        return res

    def run_obj(self, orow, ocol):
        """
        Process all objects in the image and psf
        """

        rows_per=self.image.shape[0]/self['nobj_row']
        cols_per=self.image.shape[1]/self['nobj_col']

        row1=orow*rows_per
        row2=(orow+1)*rows_per

        col1=ocol*cols_per
        col2=(ocol+1)*cols_per


        image = self.image[row1:row2, col1:col2].copy()
        psf0 = self.psf[row1:row2, col1:col2].copy()
        psf,skysig_psf = fimage.add_noise_admom(psf0,self['s2n_psf'])
        self['skysig_psf']=skysig_psf

        cenrow=image.shape[0]/2.
        cencol=image.shape[1]/2.

        psf_admom=self.run_admom(psf,cenrow,cencol,skysig_psf,4.0)
        obj_admom=self.run_admom(image,cenrow,cencol,self['skysig_obj'],6.0)

        cen_psf_admom=array([psf_admom['row'],psf_admom['col']])
        cov_psf_admom = array([psf_admom['Irr'],psf_admom['Irc'],psf_admom['Icc']])
        cen_admom = array([obj_admom['row'], obj_admom['col']])
        cov_admom = array([obj_admom['Irr'],obj_admom['Irc'],obj_admom['Icc']])

        out={}
        out['flags'] = 9999
        i=1
        while out['flags'] != 0:
            wlog("  psf try:",i)
            out['psf_res'] = self.process_image(psf, self['skysig_psf'],
                                                self['ngauss_psf'],
                                                cen_psf_admom, cov_psf_admom)
            out['psf_res']['s2n_admom'] = psf_admom['s2n']
            out['flags'] = out['psf_res']['flags']
            i+=1

        Tadmom=cov_admom[0]+cov_admom[2]
        Tpsf_admom = cov_psf_admom[0]+cov_psf_admom[2]
        psf_gmix=pars2gmix(out['psf_res']['pars'], coellip=True)

        i=1
        out['flags'] = 9999
        while out['flags'] != 0:
            wlog("  obj try:",i)
            out['res'] = self.process_image(image, self['skysig_obj'],
                                            self['ngauss_obj'],
                                            cen_admom,
                                            cov_admom,
                                            psf=psf_gmix)

            out['res']['s2n_admom'] = obj_admom['s2n']

            out['flags'] = out['res']['flags']
            i+=1

        return out


    def process_image(self, image, skysig, ngauss, cen, cov, psf=None):
        counts = image.sum()

        if psf is not None:
            maxtry=self['maxtry']
        else:
            maxtry=self['maxtry_psf']

        ntry=0
        chi2perarr=zeros(maxtry) + 1.e9
        guess_chi2perarr=zeros(maxtry) + 1.e9

        gmlist=[]


        while ntry < maxtry:
            if ntry == 0:
                eguess=[0,0]
            else:
                eguess=None

            admom_mult = 4.
            guess,width = self.get_prior_generic(ngauss,
                                                 counts, cen, cov, admom_mult,
                                                 eguess=eguess,
                                                 psf=psf)


            Tpositive=True
            if psf is not None and ngauss > 1:
                Tpositive=False
            gm = gmix_image.GMixFitCoellip(image, skysig,
                                           guess,width,
                                           psf=psf,
                                           verbose=self['verbose'],
                                           Tpositive=Tpositive)


            pars=gm.get_pars()
            perr=gm.get_perr()
            chi2per = gm.get_chi2per()

            if self['verbose']:
                print_pars(pars,front="pars:  ")
                print_pars(perr,front="perr:  ")

            chi2perarr[ntry] = chi2per
            if self['verbose']:
                wlog("chi2/pdeg:",chi2perarr[ntry])
            gmlist.append(gm)

            ntry += 1
                
        w=chi2perarr.argmin()

        gm = gmlist[w]
        if self['verbose']:
            print_pars(chi2perarr,front='chi2/deg: ')
        s2n_uw = gm.get_s2n()
        s2n_w = gm.get_weighted_s2n()

        dof=gm.get_dof()
        chi2per=chi2perarr[w]
        prob = scipy.stats.chisqprob(chi2per*dof, dof)


        popt=gm.get_pars()
        pcov=gm.get_pcov()
        perr=gm.get_perr()
        #if pcov is None:
        #    raise ValueError("best fit has failed cov calculation")
        gmix=pars2gmix(popt, coellip=True)


        if self['verbose']:
            wlog("\n")

            print_pars(popt,front='popt: ')
            print_pars(perr,front='perr: ')
            if psf is not None:
                wlog("chi^2/deg:",chi2per,"prob:",prob)
            wlog("s2n_w:",s2n_w)
            #wlog("numiter gmix:",gm.numiter)
            wlog("Topt/Tguess:",popt[ngauss+4:]/guess[ngauss+4:])

        if not hasattr(gm,'ier'):
            ier=0
            numiter=0
        else:
            ier=gm.ier
            numiter=gm.get_numiter()
            if self['verbose']:
                wlog("numiter:",numiter)
        out={'gmix':    gmix,
             'pars':    popt,
             'perr':    perr,
             'pcov':    pcov,
             'flags':   gm.get_flags(),
             'ier':     ier,
             'numiter': numiter,
             's2n_uw':     s2n_uw,
             's2n_w':     s2n_w,
             'chi2per': chi2per,
             'dof':dof,
             'prob': prob}
        return out




    def get_prior_generic(self, ngauss, counts, cen, cov, admom_mult,
                          psf=None,
                          eguess=None):
        """
        admom_mult is used for the prepsf fit ngauss==3
        """

        npars=2*ngauss+4

        prior=zeros(npars)
        width=zeros(npars) + 1.e20

        prior[0] = cen[0]
        prior[1] = cen[1]

        T = cov[2]+cov[0]
        e1=(cov[2]-cov[0])/T
        e2=2*cov[1]/T

        if psf is not None:

            if ngauss == 4:
                prior[2] = 0.4*(randu()-0.5)
                prior[3] = 0.4*(randu()-0.5)

                prior[4] = T*admom_mult*2
                prior[5] = .5
                prior[6] = 0.1
                prior[7] = 0.002
                
                # should be ~.4, .45, .05 or something for dev
                prior[8]  = 0.05*counts
                prior[9]  = 0.26*counts
                prior[10]  = 0.55*counts
                prior[11] = 0.18*counts

                rand_fac=0.2
                prior[4]  += prior[4]*rand_fac*(randu()-0.5)
                prior[5]  += prior[5]*rand_fac*(randu()-0.5)
                prior[6]  += prior[6]*rand_fac*(randu()-0.5)
                prior[7]  += prior[7]*rand_fac*(randu()-0.5)
                prior[8]  += prior[8]*rand_fac*(randu()-0.5)
                prior[9]  += prior[9]*rand_fac*(randu()-0.5)
                prior[10] += prior[10]*rand_fac*(randu()-0.5)
                prior[11] += prior[11]*rand_fac*(randu()-0.5)

            elif ngauss == 3:
                prior[2] = 0.4*(randu()-0.5)
                prior[3] = 0.4*(randu()-0.5)

                prior[4] = T*admom_mult
                prior[5] = .3
                prior[6] = 0.02
                
                prior[7] = 0.26*counts
                prior[8] = 0.55*counts
                prior[9] = 0.18*counts

                rand_fac=0.2
                prior[4] += prior[4]*rand_fac*(randu()-0.5)
                prior[5] += prior[5]*rand_fac*(randu()-0.5)
                prior[6] += prior[6]*rand_fac*(randu()-0.5)
                prior[7] += prior[7]*rand_fac*(randu()-0.5)
                prior[8] += prior[8]*rand_fac*(randu()-0.5)
                prior[9] += prior[9]*rand_fac*(randu()-0.5)


            elif ngauss==1:
                prior[0] += 1*(randu()-0.5)
                prior[1] += 1*(randu()-0.5)
                prior[2],prior[3] = randomize_e1e2(0., 0.)
                prior[4] =      T*(1. + .1*(randu()-0.5))
                prior[5] = counts*(1. + .1*(randu()-0.5))

            else:
                raise ValueError("only 1 or 3 gauss prepsf for now")
        else:
            if ngauss==1:
                prior[0] += 0.02*(randu()-0.5)  # cen0
                prior[1] += 0.02*(randu()-0.5)  # cen1
                prior[2],prior[3] = randomize_e1e2(e1,e2)
                prior[4] = T + T*0.1*(randu()-0.5)
                prior[5] = counts + 0.1*counts*(randu()-0.5)
            elif ngauss==2:
                if self['verbose']:
                    wlog("    ngauss:",ngauss)

                if eguess is not None:
                    prior[2],prior[3] = eguess
                else:
                    prior[2] = e1# + 0.05*(randu()-0.5)
                    prior[3] = e2# + 0.05*(randu()-0.5)

                e1start,e2start = prior[2],prior[3]
                prior[2],prior[3] = randomize_e1e2(e1start,e2start)
                prior[4] = T*5 # Tmax
                prior[5] = 1./5.08
                prior[6]= counts*0.1
                prior[7]= counts*0.9

                prior[4] += 0.1*prior[4]*(randu()-0.5)
                prior[5] += 0.1*prior[5]*(randu()-0.5)
                prior[6] += 0.1*prior[6]*(randu()-0.5)
                prior[7] += 0.1*prior[7]*(randu()-0.5)

            elif ngauss==3:

                if eguess is not None:
                    prior[2],prior[3] = eguess
                else:
                    prior[2] = e1# + 0.05*(randu()-0.5)
                    prior[3] = e2# + 0.05*(randu()-0.5)

                Tmax = T*8.3
                Tfrac1 = 1.7/8.3
                Tfrac2 = 0.8/8.3
                prior[4] = Tmax
                prior[5] = Tfrac1
                prior[6] = Tfrac2

                prior[7] = counts*0.08
                prior[8] = counts*0.38
                prior[9] = counts*0.53

                if self['verbose']:
                    wlog("    randomizing")
                e1start=prior[2]
                e2start=prior[3]
                prior[2],prior[3] = randomize_e1e2(e1start,e2start)

                prior[4] += prior[4]*0.05*(randu()-0.5)
                prior[5] += prior[5]*0.05*(randu()-0.5)
                prior[6] += prior[6]*0.05*(randu()-0.5)
                prior[7] += prior[7]*0.05*(randu()-0.5)
                prior[8] += prior[8]*0.05*(randu()-0.5)
                prior[9] += prior[9]*0.05*(randu()-0.5)

            else:
                raise ValueError("implement other psf guesses")

        return prior, width



    def get_output(self, orow, ocol, res):

        st = zeros(1, dtype=self.get_dtype())
        st['orow']=orow
        st['ocol']=ocol
        st['flags']=res['flags']

        st['pars_psf']     = res['psf_res']['pars']
        st['pars_psf_err'] = res['psf_res']['perr']
        st['pars_psf_cov'] = res['psf_res']['pcov']

        psf_moms = gmix_image.total_moms(res['psf_res']['gmix'])

        st['numiter_psf'] = res['psf_res']['numiter']
        st['chi2per_psf'] = res['psf_res']['chi2per']
        st['dof_psf'] = res['psf_res']['dof']

        st['s2n_admom_psf'] = res['psf_res']['s2n_admom']
        st['s2n_uw_psf'] = res['psf_res']['s2n_uw']
        st['s2n_w_psf'] = res['psf_res']['s2n_w']

        st['pars']=res['res']['pars']
        st['pars_err']=res['res']['perr']
        st['pars_cov']=res['res']['pcov']

        st['e1']    = res['res']['pars'][2]
        st['e1err'] = res['res']['perr'][2]
        st['e2']    = res['res']['pars'][3]
        st['e2err'] = res['res']['perr'][3]
        st['ecov']  = res['res']['pcov'][2:4,2:4]

        st['numiter'] = res['res']['numiter']
        st['chi2per'] = res['res']['chi2per']
        st['dof'] = res['res']['dof']

        st['s2n_admom'] = res['res']['s2n_admom']
        st['s2n_uw'] = res['res']['s2n_uw']
        st['s2n_w'] = res['res']['s2n_w']

        moms = gmix_image.total_moms(res['res']['gmix'])
        st['s2']=(psf_moms['irr']+psf_moms['icc'])/(moms['irr']+moms['icc'])

        return st


    def get_dtype(self):

        ngauss_psf=self['ngauss_psf']
        ngauss_obj=self['ngauss_obj']
        npars_psf = 2*ngauss_psf+4
        npars_obj = 2*ngauss_obj+4

        dt=[('orow','i4'),
            ('ocol','i4'),
            ('flags','i8'),
            ('e1','f8'),
            ('e1err','f8'),
            ('e2','f8'),
            ('e2err','f8'),
            ('ecov','f8',(2,2)),
            ('pars','f8',npars_obj),
            ('pars_err','f8',npars_obj),
            ('pars_cov','f8',(npars_obj,npars_obj)),
            ('chi2per','f8'),
            ('dof','f8'),
            ('prob','f8'),
            ('s2n_admom','f8'),
            ('s2n_uw','f8'), # unweighted s/n
            ('s2n_w','f8'),  # weighted by best-fit model
            ('s2','f8'),
            ('numiter','i8'),
            ('s2n_admom_psf','f8'),
            ('s2n_uw_psf','f8'),
            ('s2n_w_psf','f8'),
            ('pars_psf','f8',npars_psf),
            ('pars_psf_err','f8',npars_psf),
            ('pars_psf_cov','f8',(npars_psf,npars_psf)),
            ('chi2per_psf','f8'),
            ('dof_psf','f8'),
            ('numiter_psf','i8')]

        return dt



    def write_output(self, output):
        fname=self.get_output_file()
        if 'hdfs' not in fname:
            dname=self.get_output_dir()
            if not os.path.exists(dname):
                os.makedirs(dname)

        wlog(fname)
        eu.io.write(fname,output,clobber=True)

def randomize_e1e2(e1start,e2start):
    if e1start == 0 and e2start==0:
        e1rand = 0.05*(randu()-0.5)
        e2rand = 0.05*(randu()-0.5)
    else:
        nmax=100
        ii=0
        while True:
            e1rand = e1start*(1 + 0.2*(randu()-0.5))
            e2rand = e2start*(1 + 0.2*(randu()-0.5))
            etot = sqrt(e1rand**2 + e2rand**2)
            if etot < 0.95:
                break
            ii += 1
            if ii==nmax:
                wlog("---- hit max try on randomize e1e2, setting zero and restart")
                return randomize_e1e2(0.0,0.0)

    return e1rand, e2rand

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

def get_admom_pfrac(ci):
    cov = ci['cov_admom']
    a4o = ci['a4']
    covp = ci['cov_psf_admom']
    a4p = ci['a4_psf']

    To = cov[0]+cov[2]
    e1o = (cov[2]-cov[1])/To
    e2o = 2*cov[1]/To

    Tp = covp[0]+covp[2]
    e1p = (covp[2]-covp[1])/Tp
    e2p = 2*covp[1]/Tp


    e1,e2,R,flags = admom.correct(To, e1o, e2o, a4o,
                                  Tp, e1p, e2p, a4p)
    pfrac = 1-R
    return pfrac

def plot_admom_max_tratio(run, ei=0, tfrac=False, pfrac_am_max=None):
    """
    Plot the ratio of the maximum moment to the adaptive moment

    if tfrac==True, assume this is a Tfrac parametrization
    """
    import glob
    import biggles

    title='ei: %d' % ei
    if tfrac:
        title += ' Tfrac'

    d=shapesim.get_output_dir(run)
    pattern='%s-*-%03i.rec' % (run,ei)
    pattern=os.path.join(d,pattern)

    wlog(pattern)
    flist = glob.glob(pattern)
    nf=len(flist)


    t=eu.io.read(flist[0])
    ngauss = (t['pars'][0,:].size-4)/2

    pfrac_admoms = zeros(nf) - 9999
    Tratios = zeros(nf) -9999
    TbyTmaxs = zeros((nf,ngauss))
    Ps = zeros((nf,ngauss))

    for j,f in enumerate(flist):

        wlog(f)
        t=eu.io.read(f)


        if 'sigma0_admom' in t.dtype.names:
            # had a typo!  Is fixed now, and saving T instead of sigma
            Tam = 2*t['sigma0_admom']**2
            Tam_psf = 2*t['sigma_psf_admom']**2
        else:
            Tam = t['T_admom']
            Tam_psf = t['T_psf_admom']


        Tmax = zeros(len(t))
        TbyTmax = TbyTmax = zeros((len(t),ngauss))

        Pmax = zeros(len(t))
        P = zeros((len(t),ngauss))
        for i in xrange(len(t)):
            if tfrac:
                Tmax[i] = t['pars'][i,4].max()
                pvals = t['pars'][i,4+ngauss:].copy()
                Tfracs = t['pars'][i,4:4+ngauss].copy()
                Tfracs[0] = 1
            else:
                Tmax[i] = t['pars'][i,4+ngauss:].max()
                tvals = t['pars'][i,4+ngauss:].copy()
                s = tvals.argsort()
                tvals=array(list(reversed(tvals[s])))
                Tfracs =  tvals/Tmax[i]

                pvals = t['pars'][i,4:4+ngauss].copy()
                pvals=array(list(reversed(pvals[s])))

            Pmax[i] = pvals.max()

            #wlog("T/Tmax:",)
            print_pars(Tfracs, front="T/Tmax: ")

            TbyTmax[i,:] = Tfracs


            P[i,:] = pvals
        
        pfrac_admoms[j] = median(Tam_psf/Tam)
        Tratios[j] = median(Tmax/Tam)
        #Tratios[j] = (Tmax/Tam).max()
        for i in xrange(ngauss):
            TbyTmaxs[j,i] = median(TbyTmax[:,i])
            Ps[j,i] = median(P[:,i])

        print Tratios[j]

    s=pfrac_admoms.argsort()

    if ngauss == 4:
        w,=where(Tratios[s] > 10)
        #w,=where(Tratios > 0)
    else:
        w,=where(Tratios[s] > 0)
    w=s[w]

    if pfrac_am_max is not None:
        ww,=where(pfrac_admoms[w] < pfrac_am_max)
        if ww.size == 0:
            raise ValueError("none with pfrac_admoms < %s" % pfrac_am_max)
        w=w[ww]

    wlog("pfrac_admoms:",pfrac_admoms)
    Trat_order=2
    if run == 'gmix-fit-dt01r03':
        wtp=where(pfrac_admoms[w] > 0.8)
        Tpfit = numpy.polyfit(pfrac_admoms[w[wtp]], Tratios[w[wtp]], Trat_order)
    else:
        wtp=numpy.arange(w.size)
        Tpfit = numpy.polyfit(pfrac_admoms[w], Tratios[w], Trat_order)
    Tp = numpy.poly1d(Tpfit)

    plt=eu.plotting.bscatter(pfrac_admoms[w],Tratios[w],show=False,
                             xlabel=r'$pfrac_{AM}$',
                             ylabel=r'$T_{MAX}^{guess}/T_{AM}$')
    #tmp=numpy.linspace(pfrac_admoms[w].min(),
    #                   pfrac_admoms[w].max(),
    #                   1000)
    tmp = pfrac_admoms[w]
    plt=eu.plotting.bscatter(tmp[wtp],Tp(tmp[wtp]),
                             type='solid',color='red',
                             show=False,
                             plt=plt)


    ppfit_list = []
    pply_list = []

    Trelfit_list = []
    Trelply_list = []
    colors=['magenta','blue','red','orange']
    pplt=biggles.FramedPlot()

    pplt.xlabel=r'$pfrac_{AM}$'
    pplt.ylabel='P'
    pplts=[]

    porder=2
    Torder=3

    tab=biggles.Table(ngauss-1,1)
    for i in xrange(ngauss):
        pfit = numpy.polyfit(pfrac_admoms[w], Ps[w,i], porder)
        ply = numpy.poly1d(pfit)

        ppfit_list.append(pfit)
        pply_list.append(ply)


        pfit = numpy.polyfit(pfrac_admoms[w], TbyTmaxs[w,i], Torder)
        ply = numpy.poly1d(pfit)

        Trelfit_list.append(pfit)
        Trelply_list.append(ply)

        pply=pply_list[i]
        Trelply=Trelply_list[i]

        p=biggles.Points(pfrac_admoms[w], Ps[w,i], 
                         type='filled circle',color=colors[i])
        p.label=r'$P_{%d}$' % i
        c=biggles.Curve(tmp,pply(tmp),color=colors[i])

        pplt.add(c,p)
        pplts.append(p)
        """
        pplt=eu.plotting.bscatter(pfrac_admoms[w],
                                 Ps[w,i],show=False,plt=pplt,
                                 color=colors[i],
                                 xlabel=r'$pfrac_{AM}$',
                                 ylabel=r'$P_%d$' % i)
        pplt=eu.plotting.bscatter(tmp,
                                  pply(tmp),
                                 type='solid',color=colors[i],
                                 show=False, plt=pplt)
        """
        if i > 0:

            Trelplt=biggles.FramedPlot()
            Trelplt.xlabel=r'$pfrac_{AM}$'
            Trelplt.ylabel=r'$T/T_{max}$'
            Trelplt.title=r'$T_{%i}$' % i

            p=biggles.Points(pfrac_admoms[w], TbyTmaxs[w,i], 
                             type='filled circle',color=colors[i])
            c=biggles.Curve(tmp,Trelply(tmp),color=colors[i])
            Trelplt.add(p,c)
            tab[i-1,0] = Trelplt

    key=biggles.PlotKey(0.1,0.2,pplts,halign='left')
    pplt.add(key)

   
    wlog("P")
    wlog(Ps[w])
    wlog("T/Tmax")
    wlog(TbyTmaxs[w])

    Trat_fmt=', '.join(['%.8g']*(Trat_order+1))
    Trat_fmt= 'Trat_poly = poly1d(['+Trat_fmt+'])'
    print Trat_fmt % tuple(Tpfit)
    
    pfmt=', '.join(['%.8g']*(porder+1))
    pfmt = 'p%d_poly = poly1d(['+pfmt+'])'

    print pfmt % ((0,) + tuple(ppfit_list[0]))
    print pfmt % ((1,) + tuple(ppfit_list[1]))
    print pfmt % ((2,) + tuple(ppfit_list[2]))
    if ngauss == 4:
        print pfmt % ((3,) + tuple(ppfit_list[3]))

    Tfmt=', '.join(['%.8g']*(Torder+1))
    Tfmt = 'T%d_poly = poly1d(['+Tfmt+'])'
    print Tfmt % ( (1,) + tuple(Trelfit_list[1]) )
    print Tfmt % ( (2,) + tuple(Trelfit_list[2]) )
    if ngauss == 4:
        print Tfmt % ( (3,) + tuple(Trelfit_list[3]) )

    plt.title=title
    pplt.title=title
    tab.title=title
    plt.show()
    pplt.show()
    tab.show()
    #Trelplt.show()
