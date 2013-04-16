import os, sys
import sdsspy
from sdsspy.atlas.atlas import NoAtlasImageError
import es_sdsspy

import numpy
from numpy import where,sqrt,array,zeros,diag
from numpy.random import random as randu

import esutil as eu
from esutil.numpy_util import where1
from esutil.misc import wlog


from .util import randomize_e1e2

import gmix_image
from gmix_image import gmix2image, gmix2image_em, print_pars, gmix_print, printflags
import admom

class PSFCache:
    def __init__(self):
        self.psf_key=''
        self.psfield=None

    def cache(self, run, camcol, field, verbose=True):
        key = '%06i-%i-%04i' % (run,camcol,field)
        if key != self.psf_key:
            self.psfield = sdsspy.read('psfield',run,camcol,field, verbose=verbose)
            self.psf_kls = []
            for filter in sdsspy.FILTERCHARS:
                kl = sdsspy.read('psField',run,camcol,field,filter=filter, verbose=verbose)
                self.psf_kls.append(kl)

            self.psf_key = key

    def get_psf(self, run, camcol, field, filter, rowc, colc, verbose=True):
        """
        This will cache the psf for the input field as a list for each band.
        """
        self.cache(run, camcol, field)
        filternum=sdsspy.FILTERNUM[filter]
        kl = self.psf_kls[filternum]
        psf = kl.rec(rowc, colc, trim=True)
        return psf

    def get_psfield(self, run, camcol, field):
        self.cache(run, camcol, field)
        return self.psfield

_psf_cache = PSFCache()
def get_psf(run, camcol, field, filter, rowc, colc):
    return _psf_cache.get_psf(run, camcol, field, filter, rowc, colc)

def get_psfield(run, camcol, field):
    return _psf_cache.get_psfield(run, camcol, field)


class GMixPSF(dict):
    """
    Process SDSS psf kl decompositions.

    A cache of the psfield info for the field is kept, so process a whole field
    at once for efficiency.
    """
    def __init__(self, **keys):
        self.check_args(**keys)
        for k,v in keys.iteritems():
            self[k] = v

        self['verbose'] = self.get('verbose',False)
        conf=read_config(self['sample'])
        for k,v in conf.iteritems():
            self[k] = v

        self['filternum'] = sdsspy.FILTERNUM[self['filter']]

        psf = get_psf(self['run'],self['camcol'],self['field'],
                           self['filter'],self['rowc'],self['colc'])
        psf = psf.astype('f8')
        psf /= psf.sum()
        self.psf = psf

        self['counts'] = self.psf.sum()

        # This is typical little bit of noise in outskirts of kl decomp
        self['skysig'] = 5.6e-5

        self.psfield = get_psfield(self['run'],self['camcol'],self['field'])

        self['cen'] = (array(self.psf.shape,dtype='f8')-1.)/2.

        use_nlsolver=keys.get('use_nlsolver',False)
        if use_nlsolver:
            self.process_nlsolver()
        else:
            self.process()

    def check_args(self, **keys):
        anames=['sample','run','camcol','field','filter','rowc','colc']
        for n in anames:
            if n not in keys:
                raise ValueError("send '%s' in constructor" % n)


    def process(self):
        c=self['filternum']
        self.set_admom()

        ntrials = len(self['psf']['trials'])

        chi2arr=zeros(ntrials) + 1.e9
        gmlist=[]

        for i,trial in enumerate(self['psf']['trials']):
            ngauss=trial['ngauss']
            if ngauss==3:
                prior,width = self.get_prior_turb(trial)
            elif ngauss==2:
                prior,width = self.get_prior_2generic(trial)
            else:
                raise ValueError("only have ngauss in [2,3] now")
            #prior,width = self.get_prior_2psfield(trial)
            #prior,width = self.get_prior_test(trial)

            if self['verbose']:
                print_pars(prior,front="guess: ")
            gm = gmix_image.GMixFitCoellip(self.psf, self['skysig'],
                                           prior,width, verbose=False)
            #gm = gmix_image.gmix_fit.GMixFitCoellipNoPrior(self.psf, self['skysig'],
            #                                               prior, verbose=False)

            if gm.flags != 0:
                gmix_image.printflags("flags:",gm.flags)
                raise ValueError("error")
            chi2arr[i] = gm.get_chi2per(gm.popt)
            if self['verbose']:
                print_pars(gm.popt,front="pars:  ")
                print_pars(gm.perr,front="perr:  ")
                wlog("chi2/pdeg:",chi2arr[i])

            gmlist.append(gm)

        w=chi2arr.argmin()
        self.gm = gmlist[w]
        wlog('w:',w)

        if self['verbose']:
            print_pars(chi2arr,front='chi2/deg: ')

            wlog("\n")

            print_pars(gm.popt,front='popt: ')
            print_pars(gm.perr,front='perr: ')
            #wlog("s2n:",s2n)
            wlog("numiter gmix:",gm.numiter)
            ngauss=(len(gm.popt)-4)/2

    def process_nlsolver(self):
        wlog("Using jarvis solver")
        c=self['filternum']
        self.set_admom()

        ntrials = len(self['psf']['trials'])

        chi2arr=zeros(ntrials) + 1.e9
        gmlist=[]

        for i,trial in enumerate(self['psf']['trials']):
            ngauss=trial['ngauss']
            if ngauss==3:
                prior,width = self.get_prior_turb(trial)
            elif ngauss==2:
                prior,width = self.get_prior_2generic(trial)
            else:
                raise ValueError("only have ngauss in [2,3] now")

            #prior,width = self.get_prior_2psfield(trial)
            #prior,width = self.get_prior_test(trial)
            #prior,width = self.get_prior_2generic(trial)

            if self['verbose']:
                print_pars(prior,front="guess: ")
            # kludge using skysig=1 for now
            maxiter=2000
            psf=None
            gm=gmix_image.gmix_nlsolve.GMixCoellipSolver(self.psf, prior, 1., maxiter, psf, False)

            success=gm.get_success()
            if not success:
                raise ValueError("error")

            chi2per=gm.get_chi2per()

            # kludge
            chi2per /= self['skysig']**2

            chi2arr[i] = chi2per
            if self['verbose']:
                popt = gm.get_pars()
                cov=gm.get_pcov()
                perr=sqrt(diag(cov))

                print_pars(popt,front="pars:  ")
                # kludge
                perr *= self['skysig']
                print_pars(perr,front="perr:  ")
                wlog("chi2/pdeg:",chi2arr[i])

            gmlist.append(gm)

        w=chi2arr.argmin()
        self.gm = gmlist[w]
        wlog('w:',w)

        if self['verbose']:
            print_pars(chi2arr,front='chi2/deg: ')

            wlog("\n")

            print_pars(gm.get_pars(),front='popt: ')
            cov=gm.get_pcov()
            perr=sqrt(diag(cov))
            # kludge
            perr *= self['skysig']
            print_pars(perr,front='perr: ')



    def set_admom(self):
        c=self['filternum']
        Tguess = admom.sigma2mom(self.psfield['psf_sigma1'][0,c])
        self['amguess'] = Tguess/2.
        row,col=self['cen']
        self.amres = admom.admom(self.psf, row, col, sky=0., guess=self['amguess'], sigsky=self['skysig'])

    def get_prior_test(self, trial):
        import images
        eguess=trial['eguess']
        randomize=trial['randomize']
        uniform_p=trial['uniform_p']

        # good answer from gmix
        dguess=[{'p':0.58,'row':15.,'col':15.,'irr':1.3,'irc':-0.05,'icc':1.16},
                {'p':0.19,'row':15.,'col':15.,'irr':7.7,'irc':-0.23,'icc':6.8}]
        #psf = gmix_image.gmix2image(dguess,self.psf.shape)
        #psf += self['skysig']*numpy.random.random(self.psf.size).reshape(self.psf.shape)
        #psf *= self['counts']/psf.sum()
        #self.psf=psf
        #stop

        e1 = (dguess[0]['icc']-dguess[0]['irr'])/(dguess[0]['icc']+dguess[0]['irr'])
        e2 = 2.*dguess[0]['irc']/(dguess[0]['icc']+dguess[0]['irr'])

        Tmax = dguess[1]['icc'] + dguess[1]['irr']
        Tfrac1 = (dguess[0]['icc'] + dguess[0]['irr'])/Tmax

        pvals = array([d['p'] for d in dguess])
        #pvals *= self['counts']/pvals.sum()

        ngauss=len(dguess)
        npars = ngauss*2+4
        prior=zeros(npars)
        width=zeros(npars) + 1.e20

        prior[0] = dguess[0]['row']
        prior[1] = dguess[0]['col']
        if eguess:
            prior[2],prior[3] = eguess
        else:
            prior[2] = e1
            prior[3] = e2

        prior[4] = Tmax
        prior[5] = Tfrac1

        if uniform_p:
            prior[6] = self['counts']/ngauss
            prior[7] = self['counts']/ngauss
        else:
            prior[6] = pvals[1]
            prior[7] = pvals[0]

        if randomize:
            e1start,e2start = prior[2],prior[3]
            prior[2],prior[3] = randomize_e1e2(e1start,e2start)

            prior[4] += prior[4]*0.05*(randu()-0.5)
            prior[5] += prior[5]*0.05*(randu()-0.5)
            prior[6] += prior[6]*0.05*(randu()-0.5)
            prior[7] += prior[7]*0.05*(randu()-0.5)

            #pvals = prior[ [6,7] ].copy()
            #pvals *= self['counts']/pvals.sum()
            #prior[6] = pvals[0]
            #prior[7] = pvals[1]

        """
        width[0] = 0.001
        width[1] = 0.001
        width[2] = 0.001
        width[3] = 0.001
        width[4] = 0.001
        width[5] = 0.001
        width[6] = 0.001
        width[7] = 0.001

        width[:] = 1.e-7
        """

        return prior, width



    def get_prior_2psfield(self, trial):
        """
        Take two of the guesses from the psfield sigma1,sigma2
        """
        ngauss=trial['ngauss']
        eguess=trial['eguess']
        uniform_p=trial['uniform_p']
        randomize=trial['randomize']

        if ngauss != 2:
            raise ValueError("ngauss==2 only for now")

        npars=2*ngauss+4
        prior=zeros(npars)
        width=zeros(npars) + 1.e20
        
        prior[0] = self.amres['row']
        prior[1] = self.amres['col']

        if eguess is not None:
            prior[2],prior[3] = eguess
        else:
            prior[2] = self.amres['e1']
            prior[3] = self.amres['e2']

        T = self.amres['Irr'] + self.amres['Icc']

        c = self['filternum']
        T1 = 2*self.psfield['psf_sigma1'][0,c]**2
        T2 = 2*self.psfield['psf_sigma2'][0,c]**2

        Tmax = T2
        Tfrac1 = T1/T2
        prior[4] = Tmax
        prior[5] = Tfrac1

        if uniform_p:
            wlog("    uniform p")
            prior[6] = self['counts']/ngauss
            prior[7] = self['counts']/ngauss
        else:
            # psf_b is p2/p1
            pvals = array([self.psfield['psf_b'][0,c], 1.])
            pvals /= self['counts']*pvals.sum()
            prior[6] = pvals[0]
            prior[7] = pvals[1]


        if randomize:
            wlog("    randomizing")
            e1start=prior[2]
            e2start=prior[3]
            prior[2],prior[3] = randomize_e1e2(e1start,e2start)

            prior[4] += prior[4]*0.05*(randu()-0.5)
            prior[5] += prior[5]*0.05*(randu()-0.5)
            prior[6] += prior[6]*0.05*(randu()-0.5)
            prior[7] += prior[7]*0.05*(randu()-0.5)

            pvals = prior[ [6,7] ].copy()
            pvals *= self['counts']/pvals.sum()
            prior[6] = pvals[0]
            prior[7] = pvals[1]


        return prior, width

    def get_prior_2generic(self, trial):
        """
        generic guesses
        """

        wlog("using generic guesses")
        ngauss=trial['ngauss']
        eguess=trial['eguess']
        uniform_p=trial['uniform_p']
        randomize=trial['randomize']

        if ngauss != 2:
            raise ValueError("ngauss==2 only for now")

        npars=2*ngauss+4
        prior=zeros(npars)
        width=zeros(npars) + 1.e20
        
        prior[0] = self.amres['row']
        prior[1] = self.amres['col']

        if eguess is not None:
            prior[2],prior[3] = eguess
        else:
            prior[2] = self.amres['e1']
            prior[3] = self.amres['e2']

        T = self.amres['Irr'] + self.amres['Icc']

        Tmax = T*3
        T1 = T*0.5
        Tfrac1 = T1/Tmax

        prior[4] = Tmax
        prior[5] = Tfrac1

        if uniform_p:
            wlog("    uniform p")
            prior[6] = self['counts']/ngauss
            prior[7] = self['counts']/ngauss
        else:
            prior[6] = self['counts']*0.2
            prior[7] = self['counts']*0.8


        if randomize:
            wlog("    randomizing")
            e1start=prior[2]
            e2start=prior[3]
            prior[2],prior[3] = randomize_e1e2(e1start,e2start)

            prior[4] += prior[4]*0.05*(randu()-0.5)
            prior[5] += prior[5]*0.05*(randu()-0.5)
            prior[6] += prior[6]*0.05*(randu()-0.5)
            prior[7] += prior[7]*0.05*(randu()-0.5)

            pvals = prior[ [6,7] ].copy()
            pvals *= self['counts']/pvals.sum()
            prior[6] = pvals[0]
            prior[7] = pvals[1]


        return prior, width



    def get_prior_turb(self, trial):
        ngauss=trial['ngauss']
        eguess=trial['eguess']
        uniform_p=trial['uniform_p']
        randomize=trial['randomize']

        if ngauss != 3:
            raise ValueError("ngauss==3 only for now")

        npars=2*ngauss+4
        prior=zeros(npars)
        width=zeros(npars)
        
        prior[0] = self.amres['row']
        prior[1] = self.amres['col']
        width[0] = 1000
        width[1] = 1000

        if eguess is not None:
            prior[2],prior[3] = eguess
        else:
            prior[2] = self.amres['e1']
            prior[3] = self.amres['e2']

        T = self.amres['Irr'] + self.amres['Icc']

        # turbulent psf guess
        Tmax = T*8.3
        Tfrac1 = 1.7/8.3
        Tfrac2 = 0.8/8.3
        prior[4] = Tmax
        prior[5] = Tfrac1
        prior[6] = Tfrac2

        if uniform_p:
            wlog("    uniform p")
            prior[7] = self['counts']/ngauss
            prior[8] = self['counts']/ngauss
            prior[9] = self['counts']/ngauss
        else:
            prior[7] = self['counts']*0.08
            prior[8] = self['counts']*0.38
            prior[9] = self['counts']*0.53

        # uninformative priors
        width[2] = 1000
        width[3] = 1000
        width[4] = 1000
        width[5:] = 1000

        if randomize:
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

        return prior, width

    def show(self):
        import images
        if hasattr(self.gm,'popt'):
            pars=self.gm.popt
        else:
            pars=self.gm.get_pars()
        model=gmix2image(pars,self.psf.shape, coellip=True)
        #print((model-self.psf).max())
        #images.multiview(self.psf)
        #images.multiview(model)
        #model *= self['counts']/model.sum()
        images.compare_images(self.psf, model)


class GMixEMPSF(dict):
    """
    Process SDSS psf kl decompositions.

    A cache of the psfield info for the field is kept, so process a whole field
    at once for efficiency.
    """
    def __init__(self, **keys):
        self.check_args(**keys)
        for k,v in keys.iteritems():
            self[k] = v

        self['verbose'] = self.get('verbose',False)
        conf=read_config(self['sample'])
        for k,v in conf.iteritems():
            self[k] = v

        self['filternum'] = sdsspy.FILTERNUM[self['filter']]

        self.psf = get_psf(self['run'],self['camcol'],self['field'],
                           self['filter'],self['rowc'],self['colc'])
        self['counts'] = self.psf.sum()

        self.psfield = get_psfield(self['run'],self['camcol'],self['field'])

        self['cen'] = (array(self.psf.shape,dtype='f8')-1.)/2.

        self.process()

    def check_args(self, **keys):
        anames=['sample','run','camcol','field','filter','rowc','colc']
        for n in anames:
            if n not in keys:
                raise ValueError("send '%s' in constructor" % n)


    def process(self):
        c=self['filternum']
        self.set_admom()

        ntrials = len(self['psf']['trials'])

        chi2arr=zeros(ntrials) + 1.e9
        gmlist=[]

        im=self.psf.copy()
        im_min = im.min()
        if im_min <= 0:
            im -= im_min
            sky=0.001*im.max()
            im += sky
        else:
            sky = im_min


        for i,trial in enumerate(self['psf']['trials']):
            #prior,width = self.get_prior_turb(trial)
            guess = self.get_em_guess(trial)

            if self['verbose']:
                wlog('guess')
                gmix_print(guess,title='guess:')

            gm = gmix_image.GMixEM(im,guess,
                                   sky=sky,
                                   maxiter=4000,
                                   tol=1.e-6,
                                   coellip=False,
                                   cocenter=False) # true required for deconv



            chi2arr[i] = gm.get_fdiff()
            if self['verbose']:
                gmix_print(gm.pars,title='pars:')
                wlog("chi2/pdeg:",chi2arr[i])

            gmlist.append(gm)

        w=chi2arr.argmin()
        self.gm = gmlist[w]

        flags = gm.get_flags()
        if flags != 0:
            printflags('em',flags)
            raise ValueError("error")
        if self['verbose']:
            gmix_print(gm.pars,title='popt:')

            wlog("\n")

            wlog("numiter gmix:",gm.numiter)

    def set_admom(self):
        c=self['filternum']
        self['skysig'] = self.psfield['skysig'][0,c]
        Tguess = admom.sigma2mom(self.psfield['psf_sigma1'][0,c])
        self['amguess'] = Tguess/2.
        row,col=self['cen']
        self.amres = admom.admom(self.psf, row, col, sky=0., guess=self['amguess'], sigsky=self['skysig'])


    def get_em_guess(self, trial):
        """
        Take two of the guesses from the psfield sigma1,sigma2
        """
        ngauss=trial['ngauss']
        if ngauss > 3:
            raise ValueError("implement ngauss > 3")
        eguess=trial['eguess']
        uniform_p=trial['uniform_p']
        randomize=trial['randomize']

        c = self['filternum']

        guess=[]
        for i in xrange(ngauss):
            g={'p':1./ngauss,
               'row':self.amres['row'],
               'col':self.amres['col']}
            if i == 0:
                T = 2*self.psfield['psf_sigma1'][0,c]**2
                g['irr'] = T/2.
                g['irc'] = 0.0
                g['icc'] = T/2.
            elif i == 1:
                T = 2*self.psfield['psf_sigma2'][0,c]**2
                g['irr'] = T/2.
                g['irc'] = 0.0
                g['icc'] = T/2.
            elif i == 2:
                T1 = 2*self.psfield['psf_sigma1'][0,c]**2
                T = T1/4.
                g['irr'] = T/2.
                g['irc'] = 0.0
                g['icc'] = T/2.
            guess.append(g)

        return guess

    def show(self):
        import images
        model=gmix2image_em(self.gm.pars, self.psf.shape, counts=self['counts'])
        images.compare_images(self.psf, model)




def get_config_file(sample):
    dir=os.environ['ESPY_DIR']
    dir=os.path.join(dir,'lensing','gmix_sdss_config')
    f='gmix-sdss-%s.yaml' % sample
    f=os.path.join(dir,f)
    return f
def read_config(sample):
    f=get_config_file(sample)
    return eu.io.read(f)


def get_psf_skysig(psf):
    p1 = psf[0:5,:]
    p2 = psf[:,0:5]
    p3 = psf[30-5:,:]
    p4 = psf[:,30-5:]

    ntot = p1.size + p2.size + p3.size + p4.size
    data=zeros(ntot)

    i1 = 0
    i2 = p1.size
    data[i1:i2] = p1.ravel()
    i1 += p1.size
    i2 += p2.size
    data[i1:i2] = p2.ravel()
    i1 += p2.size
    i2 += p3.size
    data[i1:i2] = p3.ravel()
    i1 += p3.size
    i2 += p4.size
    data[i1:i2] = p4.ravel()

    print(data.std())

def test():
    import images
    dguess=[{'p':0.58,'row':15.,'col':15.,'irr':1.3,'irc':-0.05,'icc':1.16},
            {'p':0.19,'row':15.,'col':15.,'irr':7.7,'irc':-0.23,'icc':6.8}]
    im=gmix_image.gmix2image(dguess, [31,31])
    print 'counts:',im.sum()
    images.multiview(im)
