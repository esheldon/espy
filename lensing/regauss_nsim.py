import os
from sys import stderr
from pprint import pprint
import numpy

import nsim
from nsim.sim import TryAgainError
import admom

MAX_ELLIP=4.0
MAX_R = 1.0
MIN_R = 0.01

class RegaussNSim(nsim.sim.NGMixSim):
    """
    Use ngmix for the sim, regauss for fitting
    """

    def process_pair(self):
        """
        Create a simulated image pair and perform the fit
        """

        imdicts = self.get_noisy_image_pair()
        psf_image = imdicts['psf']['im']

        reslist=[]
        for key in ['im1','im2']:

            imd = imdicts[key]

            res=self.fit_galaxy(imd, psf_image)
            if res['flags'] != 0:
                raise TryAgainError("failed at %s" % key)

            self.print_res(res)

            reslist.append(res)

        return reslist


    def fit_galaxy(self, imdict, psf_image):

        image=imdict['image']

        cen_guess=image.shape[0]/2.0

        pars=imdict['pars']
        
        T_guess_psf = self.simc['psf_T']
        T_guess_obj = T_guess_psf + pars[4]

        rg = admom.ReGauss(imdict['image'],
                           cen_guess,
                           cen_guess,
                           psf_image,
                           guess_psf=T_guess_psf/2.,
                           guess=T_guess_obj/2.,
                           sigsky=self.skysig)
        rg.do_all()
        res=rg['rgcorrstats']
        if res== None:
            #pprint(rg)
            #stop
            res={'flags':1}
        else:
            e1,e2,R=res['e1'],res['e2'],res['R']
            if (abs(e1) > MAX_ELLIP
                    or abs(e2) > MAX_ELLIP
                    or R <= MIN_R
                    or R >  MAX_R):
                res['flags']=1
            else:
                # err is per component.  Boost by 1/R
                res['err']=rg['rgstats']['uncer']/R
                weight,ssh=self._get_weight_and_ssh(res['e1'],
                                                    res['e2'],
                                                    res['err'])

                res['weight'] = weight
                res['ssh'] = ssh

        return res

    def _get_weight_and_ssh(self, e1, e2, err):
        """
        err is per component, as is the shape noise
        """
        esq = e1**2 + e2**2
        err2=err**2

        # for responsivity. 
        #   Shear = 0.5 * sum(w*e1)/sum(w)/R
        # where
        #   ssh = sum(w*ssh)/sum(w)
        # this shape noise is per component

        sn2 = self.get_sn2() 

        # coefficients (eq 5-35 Bern02) 

        f = sn2/(sn2 + err2)

        ssh = 1.0 - (1-f)*sn2 - 0.5*f**2 * esq

        weight = 1.0/(sn2 + err2)

        return weight, ssh


    def get_sn2(self):
        """
        get the shape noise per component
        """
        if not hasattr(self,'SN2'):
            import ngmix
            n=10000
            g1,g2 = self.g_prior.sample2d(n)
            e1=g1.copy()
            e2=g2.copy()
            for i in xrange(n):
                e1[i], e2[i] = ngmix.shape.g1g2_to_e1e2(g1[i],g2[i])
            self.SN2 = e1.var()

            print >>stderr,'SN2:',self.SN2

        return self.SN2

    def print_res(self,res):
        print >>stderr,'    e1: %.6g e2: %.6g err: %.6g' % (res['e1'],res['e2'],res['err'])
    def get_dims_cen(self, T):
        """
        For regauss we want odd for some reason, can't remember
        """
        sigma=numpy.sqrt(T/2.)

        dim=int( 2.*sigma*nsim.sim.NSIGMA_RENDER )

        if (dim % 2) == 0:
            dim += 1

        dims = [dim]*2
        cen = [(dims[0]-1.)/2.]*2

        return dims, cen


    def copy_to_output(self, res, i):
        """
        Copy results into the output
        """
        d=self.data
        d['processed'][i] = 1
        d['e1'][i] = res['e1']
        d['e2'][i] = res['e2']
        d['err'][i] = res['err']
        d['R'][i] = res['R']

        d['weight'][i] = res['weight']
        d['ssh'][i] = res['ssh']

    def make_struct(self):
        """
        Make the output array
        """

        dt=[('processed','i2'),
            ('e1','f8'),
            ('e2','f8'),
            ('err','f8'),
            ('R','f8'),
            ('ssh','f8'),
            ('weight','f8')]

        self.data=numpy.zeros(self.npairs*2, dtype=dt)

def get_shear(data, weighted=True, shape_noise=True):
    """
    If shape_noise==False then it is assumed there is
    no shape noise (e.g. a ring) so only the errors are used
    """

    if weighted:
        return get_shear_weighted(data, shape_noise=shape_noise)
    else:
        return get_shear_unweighted(data, shape_noise=shape_noise)

def get_shear_unweighted(data, shape_noise=True):
    """
    If shape_noise==False then it is assumed there is
    no shape noise (e.g. a ring) so only the errors are used
    """

    e1mean = data['e1'].mean()
    e2mean = data['e2'].mean()

    if not shape_noise:
        err2=data['err']**2
        e1ivar = ( 1.0/err2 ).sum()
        e1err = numpy.sqrt( 1.0/e1ivar )
        e2err = e1err
    else:
        e1err = data['e1'].std()/sqrt(data.size)
        e2err = data['e2'].std()/sqrt(data.size)

    ssh = data['ssh'].mean()
    R = data['R'].mean()

    sh1=0.5*e1mean/ssh
    sh2=0.5*e2mean/ssh

    sh1err = 0.5*e1err/ssh
    sh2err = 0.5*e2err/ssh

    return sh1, sh1err, sh2, sh2err, ssh, R


def get_shear_weighted(data, shape_noise=True):
    """
    If shape_noise==False then it is assumed there is
    no shape noise (e.g. a ring) so only the errors are used
    """
    from esutil.stat import wmom

    wt=data['weight']

    e1mean,e1err = wmom(data['e1'], wt, calcerr=True)
    e2mean,e2err = wmom(data['e2'], wt, calcerr=True)
    R,Rerr = wmom(data['R'], wt, calcerr=True)
    ssh,ssherr = wmom(data['ssh'], wt, calcerr=True)

    if not shape_noise:
        err2=data['err']**2
        e1ivar = ( 1.0/err2 ).sum()
        e1err = numpy.sqrt( 1.0/e1ivar )
        e2err = e1err

    sh1=0.5*e1mean/ssh
    sh2=0.5*e2mean/ssh

    sh1err = 0.5*e1err/ssh
    sh2err = 0.5*e2err/ssh

    return sh1, sh1err, sh2, sh2err, ssh, R

def get_config_dir():
    d=os.environ['ESPY_DIR']
    return os.path.join(d,'lensing','regauss_sim_config')

def get_config_file(run):
    d=get_config_dir()
    name='%s.yaml' % run
    return os.path.join(d, name)

def read_config(run):
    """
    run could be 'name' in sim
    """
    import yaml
    f=get_config_file(run)
    with open(f) as fobj:
        c=yaml.load(fobj)

    if 'name' in c:
        n=c['name']
    else:
        n=c['run']
    if n != run:
        raise ValueError("%s in config does not match "
                         "itself: '%s' instead of '%s'" % (n,c[n],run))
    return c


def test_err(n=10000, skysig=0.1):
    """
    this shows that the uncer is per component
    """
    import ngmix
    import admom
    from numpy.random import randn

    T=4.0
    sigma=numpy.sqrt(T/2.)
    dim = int( round(2.0*5.0*sigma) )
    dims = [dim]*2

    cen = [(dim-1.)/2.]*2

    pars=[cen[0], cen[1], 0.0, 0.0, T, 100.0]

    gm = ngmix.gmix.GMixModel(pars, 'gauss')

    im0=gm.make_image(dims, nsub=16.0)

    e1s=numpy.zeros(n)
    e2s=numpy.zeros(n)
    errs = numpy.zeros(n)
    for i in xrange(n):
        
        if ( (i % 100) == 0):
            print '%d/%d' % (i+1,n)

        im = im0 + skysig*randn(im0.size).reshape(im0.shape)

        res=admom.admom(im, cen[0], cen[1], guess=T/2.,
                        sigsky=skysig)

        e1s[i] = res['e1']
        e2s[i] = res['e2']
        errs[i] = res['uncer']

    print 'mean err:    ',errs.mean()
    print 'e1 shape std:',e1s.std()
    print 'e2 shape std:',e2s.std()
