import os
from sys import stderr
from pprint import pprint
import numpy

import nsim
from nsim.sim import TryAgainError
import admom

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
        T_guess=pars[4]/2.

        rg = admom.ReGauss(imdict['image'],
                           cen_guess,
                           cen_guess,
                           psf_image,
                           guess_psf=self.simc['psf_T']/2.,
                           guess=T_guess,
                           sigsky=self.skysig)
        rg.do_all()
        res=rg['rgcorrstats']
        if res== None:
            #pprint(rg)
            #stop
            res={'flags':1}
        else:
            # err is per component
            res['err']=rg['rgstats']['uncer']
            weight,wssh=self._get_weight_and_ssh(res['e1'],res['err'])

            res['weight'] = weight

            # ssh will be sum(wssh)/sum(weight)
            res['wssh'] = wssh
        return res

    def _get_weight_and_ssh(self, e1, err):
        err2=err**2

        # for responsivity. 
        #   Shear = 0.5 * sum(w*e1)/sum(w)/R
        # where
        #   R = sum(w*F)/sum(w)
        # this shape noise is per component

        sn2 = self.get_sn2() 
        weight = 1.0/(sn2 + err2)


        # coefficients (p 596 Bern02) 
        # there is a k1*e^2/2 in Bern02 because
        # its the total ellipticity he is using
        # 
        # factor of (1/2)^2 does not cancel here, 4 converts to shapenoise

        f = sn2/(sn2 + err2)
        k0 = (1-f)*sn2
        k1 = f**2

        F = 1. - k0 - k1*e1**2

        wssh = weight*F
        return weight, wssh




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
        d['wssh'][i] = res['wssh']

    def make_struct(self):
        """
        Make the output array
        """

        dt=[('processed','i2'),
            ('e1','f8'),
            ('e2','f8'),
            ('err','f8'),
            ('R','f8'),
            ('weight','f8'),
            ('wssh','f8')]

        self.data=numpy.zeros(self.npairs*2, dtype=dt)

def get_shear(data, shape_noise=True):
    """
    If shape_noise==False then it is assumed there is
    no shape noise (e.g. a ring) so only the errors are used
    """
    from esutil.stat import wmom

    w=data['weight']
    wssh = data['wssh']

    e1=data['e1']
    e2=data['e2']

    e1mean,e1err = wmom(e1, w, calcerr=True)
    e2mean,e2err = wmom(e2, w, calcerr=True)

    if not shape_noise:
        e1ivar = ( 1.0/data['err']**2 ).sum()
        e1err = numpy.sqrt( 1.0/e1ivar )
        e2err = e1err

    wsum=w.sum()

    ssh = data['wssh'].sum()/wsum
    R = (data['R']*data['weight']).sum()/wsum

    sh1=0.5*e1mean/ssh
    sh2=0.5*e2mean/ssh

    sh1err = 0.5*e1err/ssh
    sh2err = 0.5*e2err/ssh

    return sh1, sh1err, sh2, sh2err, R

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


