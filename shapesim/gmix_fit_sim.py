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

try:
    import gmix_image
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

        ptype=self.get('ptype','e1e2')

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
            pfrac_am = Tpsf_admom/Tadmom

            pfrac_am2 = get_admom_pfrac(ci)

            out['res'] = self.process_image(ci.image, 
                                            self['ngauss_obj'],
                                            ci['cen_admom'],
                                            cov_admom,
                                            psf=out['psf_res']['gmix'],
                                            skysig=ci['skysig'],
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
                gmix_image.printflags(ptype,out['flags'])

        if out['flags'] != 0 and self['verbose']:
            self.wlog('flags:')
            if self['verbose']:
                gmix_image.printflags(ptype,out['flags'])
        return out

    def process_image(self, image, ngauss, cen, cov, psf=None,
                      pfrac_am=None,
                      skysig=None, coellip=True,
                      e1true=None, e2true=None):
        if not coellip:
            raise ValueError("must use coellip for now")

        counts = image.sum()
        method=self.get('method','lm')

        randomize=self.get('randomize',False)

        if psf:
            verbose=False
            maxtry=self.get('maxtry',1)
        else:
            maxtry=self.get('maxtry_psf',1)
            # this is a psf
            verbose=False
        ntry=0
        chi2arr=zeros(maxtry) + 1.e9
        chi2perarr=zeros(maxtry) + 1.e9
        guess_chi2perarr=zeros(maxtry) + 1.e9

        gmlist=[]

        ptype = self.get('ptype','e1e2')
        use_jacob=self.get('use_jacob',False)
        Tmin = self.get('Tmin',0.0)
        self.wlog("ptype:",ptype,"use_jacob:",use_jacob,"Tmin:",Tmin)
        while ntry < maxtry:
            if psf and ptype == 'dev':
                if pfrac_am is None:
                    raise ValueError("Must have pfrac_am for dev fit")
                guess = self.get_guess_dev(counts, cen, cov, pfrac_am,
                                           randomize=randomize)
            elif psf and ptype=='Tfrac':
            #elif ptype=='Tfrac':
                #eguess=None
                eguess=[0,0]
                uniform_p=False
                if ntry > 0:
                    uniform_p=True
                guess,width = self.get_prior_Tfrac(ngauss,
                                                   counts, cen, cov, pfrac_am,
                                                   eguess=eguess,
                                                   psf=psf,
                                                   uniform_p=uniform_p,
                                                   randomize=randomize)
            elif (psf is None) or (ptype == 'e1e2'):
                # we always go here for psf measurement
                Tfac=None
                eguess=None

                if ngauss==1 or psf is None:
                    randomize=True
                else:
                    randomize=False

                if pfrac_am is not None:
                    Tfac=None
                    if ntry==0:
                        eguess=None
                    elif ntry == 1:
                        eguess=[0,0]
                    else:
                        if (ntry % 2) == 0:
                            # needed this for eg, e=0.8, very high S/N
                            # over-ride randomize
                            randomize=True
                            eguess=None
                        else:
                            # needed this for eg, e=0.8, very high S/N
                            # over-ride randomize
                            randomize=True
                            eguess=[0,0]
                else:
                    if ngauss==3 and psf is None:
                        eguess=[0,0]

                guess = self.get_guess_coellip_e1e2(counts, 
                                                    ngauss, cen, cov, 
                                                    randomize=randomize,
                                                    psf=psf,
                                                    pfrac_am=pfrac_am,
                                                    eguess=eguess,
                                                    Tfac=Tfac)
            else:
                raise ValueError("ptype should be 'Tfrac' or 'e1e2' or 'dev'")

            if self['verbose']:
                print_pars(guess,front="guess: ")
            if psf and ptype=='dev':
                gm = gmix_image.GMixFitDev(image,guess,
                                           psf=psf,
                                           verbose=verbose)

            if psf and ptype=='Tfrac':
                #verbose=True
                if skysig is None:
                    raise ValueError("skysig must not be None")
                gm = gmix_image.GMixFitCoellipTfrac(image,
                                                    skysig,
                                                    guess,width,
                                                    psf=psf,
                                                    verbose=verbose)

            else:
                gm = gmix_image.GMixFitCoellip(image,guess,
                                               psf=psf,
                                               method=method,
                                               ptype='e1e2',
                                               Tmin=Tmin,
                                               use_jacob=use_jacob,
                                               verbose=verbose)

            if self['verbose']:
                print_pars(gm.popt,front="pars:  ")
                print_pars(gm.perr,front="perr:  ")
            if skysig is not None:
                chi2arr[ntry] = gm.get_chi2(gm.popt)
                if ptype=='Tfrac' and psf is not None:
                    chi2perarr[ntry] = gm.get_chi2per(gm.popt)
                else:
                    chi2perarr[ntry] = gm.get_chi2per(gm.popt,skysig)
                self.wlog("chi2/pdeg:",chi2perarr[ntry])
            else:
                chi2arr[ntry] = gm.get_chi2(gm.popt)
                self.wlog("chi2:",chi2arr[ntry])
            gmlist.append(gm)

            ntry += 1
                
        #if psf:
        #    stop
        w=chi2arr.argmin()
        gm = gmlist[w]
        if skysig is not None:
            if self['verbose']:
                print_pars(chi2perarr,front='chi2/deg: ')
            if ptype=='Tfrac' and psf is not None:
                s2n = gm.get_s2n(gm.popt)
            else:
                s2n = gm.get_s2n(gm.popt, skysig)
        else:
            if self['verbose']:
                print_pars(chi2arr,front='chi2: ')
            s2n = -9999
        if self['verbose']:
            wlog("\n")

            print_pars(gm.popt,front='popt: ')
            print_pars(gm.perr,front='perr: ')
            wlog("s2n:",s2n)
            wlog("numiter gmix:",gm.numiter)
            wlog("Topt/Tguess:",gm.popt[ngauss+4:]/guess[ngauss+4:])

        out={'gmix':    gm.gmix,
             'pars':    gm.popt,
             'perr':    gm.perr,
             'pcov':    gm.pcov,
             'flags':   gm.flags,
             'ier':     gm.ier,
             'numiter': gm.numiter,
             'coellip': coellip,
             's2n':     s2n,
             'chi2per': chi2perarr[w]}
        return out

    def get_prior_Tfrac(self, ngauss, counts, cen, cov, pfrac_am,
                        psf=None,
                        randomize=False,
                        uniform_p=False,
                        eguess=None):
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

            if psf:
                if len(psf)==3:
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
                elif len(psf) == 1:
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

            if len(psf)==3:
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

        elif ngauss==3 and psf:
            self.wlog("    using psf ngauss=3")
        
            if eguess is not None:
                prior[2],prior[3] = eguess 
            else:
                prior[2] = e1
                prior[3] = e2

            self.wlog("    starting e1,e2:",prior[2],prior[3])


            self.wlog("Using pfrac_am fit:",pfrac_am)
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

            prior[4] = counts
            prior[5] = T

            # uninformative
            width[2] = 10
            width[3] = 10
            width[4] = 100 # Tmax
            width[5] = 10  # p


            if psf is not None:
                self.wlog("======> with psf")
                psfmoms = gmix_image.total_moms(psf)
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
            raise ValueError("implement other guesses")

        return prior, width

    def get_guess_dev(self, counts, cen, cov, pfrac_am,
                      eguess=None,
                      randomize=False):

        if pfrac_am > 0.8:
            pfrac_am = 0.8
        ngauss=4
                                
        T = cov[2]+cov[0]
        if eguess:
            e1,e2 = eguess
        else:
            e1=(cov[2]-cov[0])/T
            e2=2*cov[1]/T

        guess=zeros(9)
        guess[0] = cen[0]
        guess[1] = cen[1]
        guess[2] = e1
        guess[3] = e2

        # we know this is wrong for turbulent!
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

        #guess[4:4+4] = counts/ngauss
        guess[4] = counts*p0
        guess[5] = counts*p1
        guess[6] = counts*p2
        guess[7] = counts*p3
        guess[8] = Tmax

        if randomize:
            self.wlog("    randomizing")
            e1start=guess[2]
            e2start=guess[3]
            if e1start == 0:
                doabs=True
            else:
                doabs=False
            while True:
                if not doabs:
                    guess[2] += 0.2*e1start*(randu()-0.5)
                    guess[3] += 0.2*e2start*(randu()-0.5)
                else:
                    guess[2] = 0.05*(randu()-0.5)
                    guess[3] = 0.05*(randu()-0.5)
                etot = sqrt(guess[2]**2 + guess[3]**2)
                if etot < 0.95:
                    break
            guess[4] += 0.2*guess[4]*(randu()-0.5)
            guess[5] += 0.2*guess[5]*(randu()-0.5)
            guess[6] += 0.2*guess[6]*(randu()-0.5)
            guess[7] += 0.2*guess[7]*(randu()-0.5)
            guess[8] += 0.2*guess[8]*(randu()-0.5)

        return guess

    def get_guess_coellip_e1e2(self, counts, ngauss, cen, cov, 
                               randomize=False, 
                               eguess=None,
                               pfrac_am=None,
                               Tfac=1.,
                               psf=None):
        self.wlog("\nusing coellip e1e2")
        npars = 2*ngauss+4
        guess=zeros(npars)
        guess[0] = cen[0]
        guess[1] = cen[1]

        T = cov[2]+cov[0]
        e1=(cov[2]-cov[0])/T
        e2=2*cov[1]/T

        guess[2] = e1
        guess[3] = e2

        self.wlog("Tadmom:",T)
        if ngauss==1:
            self.wlog("    using ngauss==1")

            guess[4] = counts
            guess[5] = T

            if psf is not None:
                self.wlog("======> with psf")
                psfmoms = gmix_image.total_moms(psf)
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

                    guess[2] = te1
                    guess[3] = te2
                    guess[4] = 1.0
                    guess[5] = tT
                else:
                    # use defaults
                    self.wlog("NOT USING special guesses")
                    pass
            if randomize:
                guess[0] += 1*(randu()-0.5)  # cen0
                guess[1] += 1*(randu()-0.5)  # cen1
                guess[2] += 0.2*(randu()-0.5)  # e1
                guess[3] += 0.2*(randu()-0.5)  # e2
                guess[4] += 0.1*(randu()-0.5)  # p
                guess[5] += 1*(randu()-0.5)   # T

        elif ngauss==3 and psf is None:
            # these are for turbulent
            self.wlog("    using ngauss=3")

            if eguess is not None:
                guess[2],guess[3] = eguess
            else:
                guess[2] = e1# + 0.05*(randu()-0.5)
                guess[3] = e2# + 0.05*(randu()-0.5)

            guess[4] = counts*0.08
            guess[5] = counts*0.38
            guess[6] = counts*0.53
            guess[7] = T*8.3
            guess[8] = T*1.7
            guess[9] = T*0.8
            if randomize:
                self.wlog("    randomizing")
                e1start=guess[2]
                e2start=guess[3]
                guess[2],guess[3] = randomize_e1e2(e1start,e2start)

                guess[4] += 0.05*guess[4]*(randu()-0.5)
                guess[5] += 0.05*guess[5]*(randu()-0.5)
                guess[6] += 0.05*guess[6]*(randu()-0.5)
                guess[7] += 0.05*guess[7]*(randu()-0.5)
                guess[8] += 0.05*guess[8]*(randu()-0.5)
                guess[9] += 0.05*guess[9]*(randu()-0.5)

        elif ngauss==3 and psf:
            self.wlog("    using psf ngauss=3")
        
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

            self.wlog("    starting e1,e2:",guess[2],guess[3])

            if psf:
                guess[4] = counts*0.62
                guess[5] = counts*0.34
                guess[6] = counts*0.04
                if pfrac_am is not None:
                    self.wlog("Using pfrac_am fit:",pfrac_am)
                    # need to do this at higher S/N
                    ply=poly1d([-3.20824373,  3.40727954])
                    tratio = ply(pfrac_am)

                    # this T is Tadmom
                    Tmax = tratio*T

                    guess[7] = Tmax
                    guess[8] = Tmax*0.35
                    guess[9] = Tmax*0.07

                    # these ratios are important when there is noise
                    guess[4] = counts*0.14
                    guess[5] = counts*0.53
                    guess[6] = counts*0.33
                    #0.13999945  0.52742499  0.33407491

                    if randomize:
                        e1start=guess[2]
                        e2start=guess[3]
                        if e1start == 0:
                            doabs=True
                        else:
                            doabs=False
                        while True:
                            if not doabs:
                                guess[2] += 0.2*e1start*(randu()-0.5)
                                guess[3] += 0.2*e2start*(randu()-0.5)
                            else:
                                guess[2] = 0.05*(randu()-0.5)
                                guess[3] = 0.05*(randu()-0.5)
                            if sqrt(guess[2]**2 + guess[3]**2) < 0.95:
                                break
                        while True:
                            guess[4] += 0.2*guess[4]*(randu()-0.5)
                            guess[5] += 0.2*guess[5]*(randu()-0.5)
                            guess[6] += 0.2*guess[6]*(randu()-0.5)
                            if guess[4] > 0 and guess[5] > 0 and guess[6] > 0:
                                break

                        guess[7] += 0.2*guess[7]*(randu()-0.5)
                        guess[8] += 0.2*guess[8]*(randu()-0.5)
                        guess[9] += 0.2*guess[9]*(randu()-0.5)

                else:
                    guess[7] = T*2.7
                    guess[8] = T*0.62
                    guess[9] = T*0.09
            else:
                guess[4:4+ngauss] = counts/ngauss    
                guess[7] = T*3.0
                guess[8] = T
                guess[9] = T*.5
        elif ngauss==4:
                self.wlog("    using ngauss=4")

                guess[4] = counts/ngauss
                guess[5] = counts/ngauss
                guess[6] = counts/ngauss
                guess[7] = counts/ngauss

                if eguess is not None:
                    guess[2],guess[3] = eguess
                else:
                    guess[2] = e1# + 0.05*(randu()-0.5)
                    guess[3] = e2# + 0.05*(randu()-0.5)

                self.wlog("    starting e1,e2:",guess[2],guess[3])

                if Tfac is None and pfrac_am is not None:
                    # from gmix_fit_sim.plot_admom_max_tratio
                    self.wlog("Using pfrac_am fit:",pfrac_am)
                    Tply=poly1d([-72.51258096,  76.65579929])

                    if pfrac_am > 0.8:
                        pfrac_am = 0.8
                    """
                    if pfrac_am > 0.8:
                        pfrac_am_p = 0.8
                    else:
                        pfrac_am_p = pfrac_am
                    """
                    pfrac_am_p = pfrac_am

                    # from low ellip
                    p0_poly = poly1d([-0.24772161,  0.22211086])
                    p1_poly = poly1d([-0.50867256,  0.48633848])
                    p2_poly = poly1d([-0.50787093,  0.6208675 ])
                    p3_poly = poly1d([ 1.28199182, -0.34435804])

                    
                    # from high ellip
                    """
                    p0_poly = poly1d([-0.25917109,  0.23568394])
                    p1_poly = poly1d([-0.52025075,  0.51105923])
                    p2_poly = poly1d([-0.44028671,  0.5890214])
                    p3_poly = poly1d([1.2392543,  -0.35254647])
                    T0_poly = poly1d([1.1559449e-16,  1])
                    T1_poly = poly1d([0.15519577,  0.067887437])
                    T2_poly = poly1d([0.042111311,  0.0032187106])
                    T3_poly = poly1d([0.00083959131,  0.0018759731])
                    """


                    tratio = Tply(pfrac_am)
                    p0 = p0_poly(pfrac_am_p)
                    p1 = p1_poly(pfrac_am_p)
                    p2 = p2_poly(pfrac_am_p)
                    p3 = p3_poly(pfrac_am_p)
  
                    #  p0 poly: [-0.24772161  0.22211086]
                    #  p1 poly: [-0.50867256  0.48633848]
                    #  p2 poly: [-0.50787093  0.6208675 ]
                    #  p3 poly: [ 1.28199182 -0.34435804]

                    # this T is Tadmom
                    Tmax = tratio*T

                    # 1.000000 0.182911 0.035527 0.002705
                    guess[8] = Tmax
                    guess[9] = Tmax*0.18
                    guess[10] = Tmax*0.035
                    guess[11] = Tmax*0.0027

                    # these ratios are important when there is noise
                    #0.03982061  0.11440515  0.25216625  0.59183215
                    guess[4] = counts*p0
                    guess[5] = counts*p1
                    guess[6] = counts*p2
                    guess[7] = counts*p3

                    # old
                    #guess[4] = counts*0.04
                    #guess[5] = counts*0.11
                    #guess[6] = counts*0.25
                    #guess[7] = counts*0.60

                    if randomize:
                        self.wlog("    randomizing")
                        e1start=guess[2]
                        e2start=guess[3]
                        if e1start == 0:
                            doabs=True
                        else:
                            doabs=False
                        while True:
                            if not doabs:
                                guess[2] += 0.2*e1start*(randu()-0.5)
                                guess[3] += 0.2*e2start*(randu()-0.5)
                            else:
                                guess[2] = 0.05*(randu()-0.5)
                                guess[3] = 0.05*(randu()-0.5)
                            etot = sqrt(guess[2]**2 + guess[3]**2)
                            if etot < 0.95:
                                break
                        guess[4] += 0.2*guess[4]*(randu()-0.5)
                        guess[5] += 0.2*guess[5]*(randu()-0.5)
                        guess[6] += 0.2*guess[6]*(randu()-0.5)
                        guess[7] += 0.2*guess[7]*(randu()-0.5)
                        guess[8] += 0.2*guess[8]*(randu()-0.5)
                        guess[9] += 0.2*guess[9]*(randu()-0.5)
                        guess[10] += 0.2*guess[10]*(randu()-0.5)
                        guess[11] += 0.2*guess[11]*(randu()-0.5)

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
        elif ngauss==2:
            self.wlog("    using ngauss==1")
            # generic guesses
            guess[4]= counts/ngauss
            guess[5]= counts/ngauss
            guess[6] = T*2
            guess[7] = T*0.5

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
        ptype=self['ptype']
        if not self['coellip_psf']:
            raise ValueError("must use coellip for psf for now")
        if (ptype != 'dev') and (not self['coellip_obj']):
            raise ValueError("must use coellip or dev for obj for now")

        st = zeros(1, dtype=self.out_dtype())

        # first copy inputs and data from the CI
        st['s2'] = s2
        st['s2n_uw'] = ci['s2n_uw']
        st['s2n_matched'] = ci['s2n_matched']
        st['s2n_uw_psf'] = ci['s2n_uw_psf']
        st['s2n_matched_psf'] = ci['s2n_matched_psf']
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
        else:
            st['s2_meas'] = -9999





        return st


    def out_dtype(self):
        ptype=self['ptype']

        ngauss_psf=self['ngauss_psf']
        ngauss_obj=self['ngauss_obj']
        npars_psf = 2*ngauss_psf+4
        if ptype == 'dev':
            npars_obj = 9
        else:
            npars_obj = 2*ngauss_obj+4

        gmix_dt = [('p','f8'),('row','f8'),('col','f8'),
                   ('irr','f8'),('irc','f8'),('icc','f8')]
        dt=[('s2n_uw','f8'),
            ('s2n_uw_psf','f8'),
            ('s2n_matched','f8'),
            ('s2n_matched_psf','f8'),
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
            ('gmix',gmix_dt,ngauss_obj),
            ('pars','f8',npars_obj),
            ('pars_err','f8',npars_obj),
            ('pars_cov','f8',(npars_obj,npars_obj)) ]

        return dt

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
                e1start=0
                e2start=0
                ii=0

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
