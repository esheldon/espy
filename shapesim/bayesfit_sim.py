from numpy import sqrt, cos, sin, exp, pi, zeros, ones, empty, \
        random, where, array, linspace
from . import shapesim
import lensing
from fimage.convolved import NoisyConvolvedImage
import gmix_image
import images
from esutil.misc import wlog

from sys import stderr


class BayesFitSim(shapesim.BaseSim):
    def __init__(self, run):
        """
        use config files 
        """
        super(BayesFitSim,self).__init__(run)
        if 'verbose' not in self:
            self['verbose'] = False

        self.gprior = GPrior()
        self.fitter = BayesFitter(None)


    def process_trials_by_s2n(self, is2, is2n):
        """
        ring test

        Generates a random total ellipticity from the 
        assumed prior

        Run this many times to sample the prior distribution
        """

        s2n = shapesim.get_s2n(self, is2n)
        s2n_fac = self['s2n_fac']

        nring = self.simc['nring']
        nrepeat = shapesim.get_s2n_nrepeat(s2n, fac=s2n_fac)

        ntot = nring*nrepeat
        out = zeros(ntot, dtype=self.out_dtype())

        # random g value
        g = self.gprior.sample1d(1)[0]

        ii = 0
        for i in xrange(nring):
            itheta=i

            if self['verbose']:
                stderr.write('-'*70)
                stderr.write('\n')
            # we always write this
            stderr.write("%d/%d %d%% done\n" % ((i+1),nring,
                                                100.*(i+1)/float(nring)))

            dolog=False
            if i==0:
                dolog=True
            st = self.process_trial_by_s2n(is2, is2n, itheta, g, dolog=dolog)
            out[ii:ii+nrepeat] = st
            ii += nrepeat

        #write_output(self['run'], is2, is2n, out, fs=self.fs)
        return out



    def process_trial_by_s2n(self, is2, is2n, itheta, g,
                             dolog=False):
        """
        Process a singe element in the ring, with nrepeat
        possible noise realizations
        """
        ellip=lensing.util.g2e(g)

        s2 = linspace(self.simc['mins2'],
                      self.simc['maxs2'], 
                      self.simc['nums2'])[is2]
        s2n_psf = self['s2n_psf']
        s2n = shapesim.get_s2n(self, is2n)
        theta = shapesim.get_theta(self.simc, itheta=itheta)

        print 'true g:',g,'theta:',theta,'s/n:',s2n,'s2:',s2

        s2n_fac = self['s2n_fac']
        s2n_method = self['s2n_method']
        s2ncalc_fluxfrac =self['s2ncalc_fluxfrac']

        nrepeat = shapesim.get_s2n_nrepeat(s2n, fac=s2n_fac)

        # do this once and get noise realizations.
        ci_nonoise = self.shapesim.get_trial(s2, ellip, theta)

        e1true=ci_nonoise['e1true']
        e2true=ci_nonoise['e2true']
        g1true,g2true=lensing.util.e1e2_to_g1g2(e1true,e2true)
        print 'true sheared g1,g2:',g1true,g2true

        if dolog:
            wlog("ring theta: %s/%s" % (itheta+1,self.simc['nring']))
            wlog('ellip:',ellip,'s2n:',s2n,
                 's2n_psf:',s2n_psf,'s2n_method:',s2n_method)

        out = zeros(nrepeat, dtype=self.out_dtype())
        out['gtrue'][:,0] = g1true
        out['gtrue'][:,1] = g2true
        out['shear_true'][:,0] = self.simc['shear'][0]
        out['shear_true'][:,1] = self.simc['shear'][1]
        for irepeat in xrange(nrepeat):
            
            if self['verbose']:
                stderr.write('-'*70 + '\n')
            # we always write this, although slower when not verbose
            if (nrepeat > 1) and (( (irepeat+1) % 10) == 0 or irepeat == 0):
                stderr.write("  %s/%s repeat done\n" % ((irepeat+1),nrepeat))

            ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf,
                                     s2n_method=s2n_method,
                                     fluxfrac=s2ncalc_fluxfrac)
            if self['verbose']:
                wlog("s2n_admom:",ci['s2n_admom'],"s2n_uw:",ci['s2n_uw'],"s2n_matched:",ci['s2n_matched'])

            # self.fitter is used in here
            self._run_fitter(ci)

            # just likelihood right now
            like = self.fitter.get_like()
            like = like/like.sum()

            res = self.fitter.get_result()

            out['g'][irepeat,:] = res['g']
            out['gcov'][irepeat,:,:] = res['gcov']
            """
            print 'g1: %.4g +/- %.4g' % (res['g'][0],sqrt(res['cov'][0,0]))
            print 'g2: %.4g +/- %.4g' % (res['g'][1],sqrt(res['cov'][1,1]))
            print 'covar:',res['cov'][0,1]
            """

            #images.multiview(like)
            #stop
            """
            if res['flags'] == 0:
                st = self.copy_output(s2, ellip, s2n, ci, res)
                out[irepeat] = st
                break
            else:
                raise ValueError("error encountered")
            """
        if self['verbose']:
            stderr.write("niter: %d\n" % (iter+1))

        return out

    def _run_fitter(self, ci):
        """
        cheat on psf, T and cen for now
        """
        cov=ci['covtrue']
        T = cov[0] + cov[2]
        cen = ci['cen']

        # p/row/col unimportant
        psf_cov = ci.psfpars['cov']
        psf=[{'p':1,'row':1,'col':1,
              'irr':psf_cov[0],'irc':psf_cov[1],'icc':psf_cov[2]}]
        self.fitter.process_image(image=ci.image,
                                  pixerr=ci['skysig'],
                                  cen=ci['cen'],
                                  T=T,
                                  psf=psf)

    def out_dtype(self):
        dt=[('shear_true','f8',2),
            ('gtrue','f8',2),
            ('g','f8',2),
            ('gcov','f8',(2,2))]
        return dt


class BayesFitter:
    """
    Desiged to be re-used since the models are all
    stored in a cache.

    likelihood is exp(0.5*A*B^2)

    A is 
        sum( (model/err)^2 ) which we will fix at 1
    B is
        sum( model*data/err^2 )/A

    We must fix A as a constant for every model we generate, so that
    the relative height of the likelihood is valid give the other
    simplifications we have made.  We can take A=1 for example

    So we should normalize the model according to the A condition.
    create
        ymod = model/err
    which will have a normalization S. Then the renormalization is given by
        N = sqrt( S*A/sum(ymod^2) )
    (i.e. set ymod = ymod*N/S) 
    """
    def __init__(self, prior):
        self.prior=prior

        # range for search
        self.gmin=-.98
        self.gmax= .98
        #self.ngrid=20 # in both g1 and g2
        self.ngrid=100 # in both g1 and g2

        self.g1vals = linspace(self.gmin, self.gmax, self.ngrid)
        self.g2vals = self.g1vals.copy()

        self.g1matrix = self.g1vals.reshape(self.ngrid,1)*ones(self.ngrid)
        self.g2matrix = self.g2vals*ones(self.ngrid).reshape(self.ngrid,1)

        self.delta_g = self.g1vals[1]-self.g1vals[0]

        self.models={}

    def process_image(self, 
                      image=None, 
                      pixerr=None, 
                      psf=None,
                      cen=None, # send to fix cen
                      T=None):  # send to fix T
        self._set_data(image=image,
                       pixerr=pixerr,
                       psf=psf,
                       cen=cen,
                       T=T)

        self._go()

    def _go(self):
        T     = self.T
        cen   = self.cen
        dims = self.image.shape

        loglike=zeros((self.ngrid, self.ngrid))

        for i1 in xrange(self.ngrid):
            g1 = self.g1vals[i1]
            for i2 in xrange(self.ngrid):
                g2 = self.g2vals[i2]

                g=sqrt(g1**2 + g2**2)
                if g >= 1:
                    # leave like[i1,i2]==0
                    continue

                # like is exp(0.5*A*B^2) where
                # a is sum((model/err)^2)
                # and B is sum(model*image/err^2)/A
                # since A is 1, all we need is sum (model*image)/err^2

                # this is model/err^2/sqrt(2)
                mod_over_err2=self._get_normalized_model(i1, i2, T, dims, cen)

                # this creates a temporary, can speed up
                chi2=( mod_over_err2*self.image ).sum()
                #like[i1,i2] = exp(chi2)
                loglike[i1,i2] = chi2

        loglike -= loglike.max()
        like = exp(loglike)

        like /= like.sum()
        self._like=like
        # get the expectation values, sensitiviey and errors
        self._calc_expected()

    def get_like(self):
        return self._like
    def get_result(self):
        return self._result

    def _calc_expected(self):
        if self.prior is None:
            #g1 = self.delta_g**2 * (self._like*self.g1matrix).sum()
            #g2 = self.delta_g**2 * (self._like*self.g2matrix).sum()
            # assuming likelihood is normalized
            g1 = (self._like*self.g1matrix).sum()
            g2 = (self._like*self.g2matrix).sum()


            g1diff = g1-self.g1matrix
            g2diff = g2-self.g2matrix
            g11var = (self._like*g1diff**2).sum()
            g22var = (self._like*g2diff**2).sum()
            g12var = (self._like*g1diff*g2diff).sum()

        else:
            raise ValueError("implement working with prior")
        
        g=zeros(2)
        cov=zeros((2,2))
        g[0] = g1
        g[1] = g2
        cov[0,0] = g11var
        cov[0,1] = g12var
        cov[1,0] = g12var
        cov[1,1] = g22var

        self._result={'g':g,'gcov':cov}



    def _get_normalized_model(self, i1, i2, T, dims, cen):
        #import images
        model = self._get_model(i1, i2, T, dims, cen)

        #images.multiview(model)
        #key=raw_input('hit a key: ')
        ymod = model/self.pixerr

        ymodsum = ymod.sum()
        ymod2sum = (ymod**2).sum()

        norm = sqrt(ymodsum**2*self.A/ymod2sum)

        # also divide again by err since we actually need model/err^2
        ymod *= (norm/ymodsum)

        #print 'A:',self.A,'sum(ymod^2):',(ymod**2).sum()
        #print ymod.shape
        #print ymod
        # also divide again by err since we actually need model/err^2
        # put this above after we do the check
        # also put 1/sqrt(2) since we need 0.5*B^2 in the likelihood
        ymod *= (1./self.pixerr/sqrt(2))

        return ymod


    def _get_model(self, i1, i2, T, dims, cen):
        """
        for now don't bother with cache
        """

        sh = self.get_shape(i1,i2)
        #print 'shape:',sh

        model=empty(dims,dtype='f8')
        #model=zeros(dims,dtype='f8')

        # for now pars are only for single gaussian
        # note setting p=1 since we must renormalize anyway
        pars=array([cen[0],cen[1],sh.e1,sh.e2,T,1.0],dtype='f8')
        #print 'pars:',pars
        #print 'psf_pars:',self.psf_pars
        gmix_image.render._render.fill_model_coellip(model, pars, self.psf_pars, None)

        return model

    def get_shape(self, i1, i2):
        g1 = self.g1vals[i1]
        g2 = self.g2vals[i2]
        #print 'g1,g2:',g1,g2
        shape = lensing.Shear(g1=g1,g2=g2)
        return shape

    def _create_model_container(self):
        if not hasattr(self,'model'):
            pass
    def _set_data(self, 
                  image=None, 
                  pixerr=None, 
                  psf=None,
                  cen=None,
                  T=None):
        self.image=image
        self.pixerr=pixerr
        self.T=T
        self.cen=cen

        self._set_psf(psf)

        if image is None or pixerr is None:
            raise ValueError("send at least image and pixerr")
        if T is None or cen is None:
            raise ValueError("for now send T and cen to fix them")

        self.A = ( (image/pixerr)**2 ).sum()

    def _set_psf(self, psf):
        self.psf_gmix = psf
        self.psf_pars = None

        if psf is not None:
            if not isinstance(psf[0],dict):
                raise ValueError("psf must be list of dicts")
            self.psf_pars = gmix_image.gmix2pars(psf)


class GPrior:
    """
    This is in g1,g2 space

    2D
    Prob = A cos(|g| pi/2) exp( - [ 2 |g| / B / (1 + |g|^D) ]^C )

    For 1D prob, you need to multiply by 2*pi*|g|
    """
    def __init__(self):
        # A actually depends on norm when doing full thing
        self.A = 12.25
        self.B = 0.03
        self.C = 0.45
        self.D = 13.

        self.maxval = self.prior(0., 0.)

    def prior(self, g1, g2):
        """
        Get the prior for the input g1,g2 value(s)

        This is the 2 dimensional prior.  If you want to just
        generate the |g| values, use prior1d
        """
        g = sqrt(g1**2 + g2**2)
        return self.prior_gabs(g)

    def prior_gabs(self, g):
        """
        Get the 2d prior for the input |g| value(s)
        """
        g = array(g, ndmin=1, copy=False)
        prior = zeros(g.size)

        w,=where(g < 1)
        if w.size > 0:
            prior[w] = self.A * cos(g[w]*pi/2)*exp( - ( 2*g[w] / self.B / (1 + g[w]**self.D) )**self.C )
        return prior



    def sample(self, nrand, as_shear=False):
        """
        Get random g1,g2 values

        parameters
        ----------
        nrand: int
            Number to generate
        as_shear: bool, optional
            If True, get a list of Shear objects
        """
        g1 = zeros(nrand)
        g2 = zeros(nrand)

        ngood=0
        nleft=nrand
        while ngood < nrand:

            # generate total g**2 in [0,1)
            grand2 = random.random(nleft)
            grand = sqrt(grand2)
            # now uniform angles
            rangle = random.random(nleft)*2*pi

            # now get cartesion locations in g1,g2 plane
            g1rand = grand*cos(rangle)
            g2rand = grand*sin(rangle)

            # now finally the height from [0,maxval)
            h = self.maxval*random.random(nleft)

            pvals = self.prior(g1rand, g2rand)

            w,=where(h < pvals)
            if w.size > 0:
                g1[ngood:ngood+w.size] = g1rand[w]
                g2[ngood:ngood+w.size] = g2rand[w]
                ngood += w.size
                nleft -= w.size

        if as_shear:
            from lensing.shear import Shear
            shlist=[]
            for g1i,g2i in zip(g1,g2):
                shlist.append(Shear(g1=g1i,g2=g2i))
            return shlist
        else:
            return g1, g2


    def prior1d(self, g):
        """
        Get the 1d prior for an input |g| value(s).

        To generate 2-d g1,g2, use prior()
        """
        return 2*pi*g*self.prior_gabs(g)

    def sample1d(self, nrand):
        """
        Get random |g| from the 1d distribution

        parameters
        ----------
        nrand: int
            Number to generate
        """

        if not hasattr(self,'maxval1d'):
            self.set_maxval1d()

        g = zeros(nrand)

        ngood=0
        nleft=nrand
        while ngood < nrand:

            # generate total g**2 in [0,1)
            grand = random.random(nleft)

            # now finally the height from [0,maxval)
            h = self.maxval1d*random.random(nleft)

            pvals = self.prior1d(grand)

            w,=where(h < pvals)
            if w.size > 0:
                g[ngood:ngood+w.size] = grand[w]
                ngood += w.size
                nleft -= w.size
   
        return g

    def set_maxval1d(self):
        import scipy.optimize
        
        (minvalx, fval, iterations, fcalls, warnflag) \
                = scipy.optimize.fmin(self.prior1dneg, 0.1, full_output=True, 
                                      disp=False)
        if warnflag != 0:
            raise ValueError("failed to find min: warnflag %d" % warnflag)
        self.maxval1d = -fval

    def prior1dneg(self, g, *args):
        """
        So we can use the minimizer
        """
        return -self.prior1d(g)

