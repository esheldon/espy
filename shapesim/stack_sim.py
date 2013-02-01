"""

For fits and shear=[0.04,,0.00] Getting <e> of 0.0785 at all s/n instead of
0.08 For unweighted moments corrected by 1/Rshear I get the right answer.  The
residuals of model and image look very good.  See similar bias for "et" model
as gauss model.

    - max like bias some how? no, tried mcmc
    - working in e instead of g?  no
    - pixel problem

        * made both psf and object *much* bigger but got same bias!  But making
        the object bigger relative to the PSF did change the bias.

        * made objects round and I see no bias (and no need for Rshear)

        * for gauss, there is no bias if I make the shear pure e2.  If there is
        any component that is e1 I see a bias in both shear1 and shear2.  For "et"
        I still saw the bias for pure e2.

        * moved the center around but got same bias

        * maybe nsub=16 isn't good enough in sim?  Trying 24.  Nope.

        - maybe it is fortran itself? nope
            * used gmix_image for sims.
            * added more padding around image


"""

import os
import numpy
from numpy import sqrt, diag, log10, exp
import gmix_image
from gmix_image.util import get_estyle_pars, calculate_some_stats, print_pars
import lensing
from lensing.shear import Shear
from lensing.util import ShapeRangeError
from .shapesim import read_config, ShapeSim
from . import shapesim

import copy

import esutil as eu
from esutil.random import srandu

LOWVAL=-9999.9e9
RSHEAR_DEFAULT = 1.-0.126228

def combine_stacks(run, is2=None, is2n=None):
    conf=read_config(run)
    simc = read_config(conf['sim'])

    num_s2 = simc['nums2']
    num_s2n=len(conf['s2nvals'])
    nsplit=conf['nsplit']

    if is2 is None:
        is2_vals=range(num_s2)
    else:
        is2_vals=numpy.array(is2,ndmin=1)
    if is2n is None:
        is2n_vals=range(num_s2n)
    else:
        is2n_vals=numpy.array(is2n,ndmin=1)

    for is2 in is2_vals:
        for is2n in is2n_vals:
            print 'combining is2:',is2,'is2n:',is2n
            combiner=StackCombiner(run, is2=is2, is2n=is2n)
            combiner.load_splits()
            combiner.fit_stacks()
            combiner.write()


def plot_is2_stacks(run, is2=None, type='fit', true_shear=None, 
                    yrange=None):
    import biggles
    data=load_is2_stacks(run,is2)
    
    if is2 is None or type is None:
        raise ValueError("send is2= and type=")
    if type=='uw':
        sh1=0.5*data['e1_uw']/data['Rshear']
        sh2=0.5*data['e2_uw']/data['Rshear']
    elif type=='fit':
        #sh1=0.5*data['pars'][:,2]/RSHEAR_DEFAULT
        #sh2=0.5*data['pars'][:,3]/RSHEAR_DEFAULT
        # wierd, this gives a closer answer for ngauss=2,3!
        sh1=0.5*data['pars'][:,2]
        sh2=0.5*data['pars'][:,3]
    elif type=='admom':
        sh1=0.5*data['e1_admom']/data['Rshear']
        sh2=0.5*data['e2_admom']/data['Rshear']
    else:
        raise ValueError("bad type: '%s'" % type)

    sherr=0.16/sqrt(data['nimage'])
    err1=sqrt(sherr**2 + data['perr'][:,2]**2)
    err2=sqrt(sherr**2 + data['perr'][:,3]**2)

    plt=biggles.FramedPlot()
    plt.add( biggles.Points(data['s2n'], sh1, color='blue',type='filled circle'))
    plt.add( biggles.SymmetricErrorBarsY(data['s2n'], sh1, err1, color='blue'))
    plt.add( biggles.Points(data['s2n'], sh2, color='red',type='filled circle'))
    plt.add( biggles.SymmetricErrorBarsY(data['s2n'], sh2, err2, color='red'))

    plt.title="run: %s  shear type: %s" % (run,type)

    if true_shear is not None:
        t1=data['s2n']*0 + true_shear[0]
        t2=data['s2n']*0 + true_shear[1]
        plt.add(biggles.Curve(data['s2n'], t1))
        plt.add(biggles.Curve(data['s2n'], t2))

    if yrange is not None:
        plt.yrange=yrange
    plt.show()

def load_is2_stacks(run, is2):
    """
    load stacks for the input is2
    """
    conf=read_config(run)
    simc = read_config(conf['sim'])

    num_s2 = simc['nums2']
    if is2 < 0 or is2 > (num_s2-1):
        raise ValueError("s2 %d out of range: [%d,%d)" % \
                         (is2, 0, num_s2))
    num_s2n=len(conf['s2nvals'])

    data=[]
    for is2n in xrange(num_s2n):
        data0 = shapesim.read_output(run, is2, is2n, verbose=True)

        data.append(data0)

    data=eu.numpy_util.combine_arrlist(data)
    return data

class StackSimBase(dict):
    """
    Generate or load and combine stacks.  Fit for shear.

    Send itrial= on construction if you plan to generate the images.  This is a
    "split", an organizational tool for splitting the data across multiple
    jobs.

    But you can also just run load_splits() to load all the "trials" or
    splits.  and combine them before running.  In that case you don't need to
    send itrial

    """
    def __init__(self, run, **keys):
        conf=read_config(run)
        self.update(conf)

        self._extract_keys(**keys)

        self.simc = read_config(self['sim'])

        self.s2n = shapesim.get_s2n(self, self.is2n)
        self.s2 = numpy.linspace(self.simc['mins2'],
                                 self.simc['maxs2'], 
                                 self.simc['nums2'])[self.is2]
        self.sratio=numpy.sqrt(1./self.s2)
        self.nimage=self._get_nimage()

        self._set_gprior()

        simpars=self.get('simpars',{})
        self.shapesim = ShapeSim(self['sim'], **simpars)

        self.nsub=1
        self.ngauss_psf=3
        self.ngauss_obj=3

    def fit_stacks(self):
        self._fit_stacks()
        self._set_result()

    def write(self):
        shapesim.write_output(self['run'], self.is2, self.is2n, self._data, 
                              itrial=self.itrial)

    def get_data(self):
        return self._data.copy()

    def _extract_keys(self, **keys):
        self.is2=keys.get('is2',None)
        if self.is2 is None:
            raise ValueError("send is2=")

        self.is2n=keys.get('is2n',None)
        if self.is2n is None:
            raise ValueError("send is2n=")

        self.itrial=keys.get('itrial',None)

        self.verbose=keys.get('verbose',False)


    def _get_nimage(self):
        s2n_fac = self['s2n_fac']
        nellip = shapesim.get_s2n_nrepeat(self.s2n, fac=s2n_fac)

        if nellip < self['min_gcount']:
            nellip=self['min_gcount']
        return nellip


    def _set_result(self):

        res=self._fitter.get_result()
        pres=self._psf_fitter.get_result()

        pars     = res['pars']
        perr     = res['perr']
        pcov     = res['pcov']

        psf_pars = pres['pars']
        psf_perr = pres['perr']
        psf_pcov = pres['perr']


        s2n=res['s2n_w']
        psf_s2n=pres['s2n_w']

        data=self._get_struct(pars.size, psf_pars.size)

        data['model'] = self['fitmodel']
        data['nimage'] = self.nimage
        data['s2n'] = self.s2n
        data['s2'] = self.s2
        data['sratio'] = self.sratio

        data['shear_true'] = self.simc['shear']


        data['e'] = pars[2:2+2]
        data['err'] = perr[2:2+2]
        data['ecov'] = pcov[2:2+2, 2:2+2]
        data['pars'] = pars
        data['perr'] = perr
        data['pcov'] = pcov
        data['Tmean'] = pars[4]
        data['Terr'] = numpy.sqrt(pcov[4,4])
        data['Ts2n'] = data['Tmean']/data['Terr']

        data['psf_pars'] = psf_pars
        data['psf_perr'] = psf_perr
        data['psf_pcov'] = psf_pcov

        data['loglike'] = res['loglike']
        data['chi2per'] = res['chi2per']
        data['dof'] = res['dof']
        data['fit_prob'] = res['fit_prob']

        data['e1_uw'] = self.res_uw['e1']
        data['e2_uw'] = self.res_uw['e2']
        data['e1_admom'] = self._ares['e1corr']
        data['e2_admom'] = self._ares['e2corr']

        data['s2n_stack'] = s2n
        data['image_stack'][0,:,:] = self._image_stack
        data['image_skyvar'] = self._image_skyvar

        data['psf_s2n_stack'] = psf_s2n
        data['psf_stack'][0,:,:] = self._psf_stack
        data['psf_skyvar'] = self._psf_skyvar

        data['Rshear'] = self.Rshear
        self._data=data


    def _get_struct(self, npars, psf_npars):
        shape=self._image_stack.shape
        psf_shape=self._psf_stack.shape

        dt=[('model','S5'),
            ('nimage','i4'),
            ('s2n','f8'),  # per object
            ('s2','f8'),
            ('sratio','f8'),
            ('shear_true','f8',2),

            ('e','f8',2),
            ('err','f8',2),
            ('ecov','f8',(2,2)),
            ('Tmean','f8'),
            ('Terr','f8'),
            ('Ts2n','f8'),

            ('pars','f8',npars),
            ('perr','f8',npars),
            ('pcov','f8',(npars,npars)),
            ('psf_pars','f8',psf_npars),
            ('psf_perr','f8',psf_npars),
            ('psf_pcov','f8',(psf_npars,psf_npars)),

            ('loglike','f8'),     # loglike of fit
            ('chi2per','f8'),     # chi^2/degree of freedom
            ('dof','i4'),         # degrees of freedom
            ('fit_prob','f8'),    # probability of the fit happening randomly

            ('e1_uw','f8'),
            ('e2_uw','f8'),
            ('e1_admom','f8'),
            ('e2_admom','f8'),

            ('s2n_stack','f8'),
            ('image_stack','f8',shape),
            ('image_skyvar','f8'),

            ('psf_s2n_stack','f8'),
            ('psf_stack','f8',psf_shape),
            ('psf_skyvar','f8'),
            ('Rshear','f8')
           ]

        data=numpy.zeros(1, dtype=dt)
        return data


    def _fit_stacks(self):
        import admom

        self.res_uw=self._run_uw(self._image_stack, self._psf_stack)

        ares = self._run_admom(self._image_stack, self._image_skyvar)
        psf_ares = self._run_admom(self._psf_stack, self._psf_skyvar)
        
        T=ares['Irr']+ares['Icc']
        Tpsf=psf_ares['Irr'] + psf_ares['Icc']
        e1corr,e2corr,R,flags=admom.correct(T, ares['e1'], ares['e2'],ares['a4'],
                                            Tpsf, psf_ares['e1'], psf_ares['e2'], psf_ares['a4'])
        ares['e1corr'] = e1corr
        ares['e2corr'] = e2corr
        ares['R'] = R
        ares['compea4_flags']=flags

        self._ares=ares
        self._psf_ares=psf_ares
        print 'admom sh1:',0.5*e1corr/self.Rshear,'e2:',0.5*e2corr/self.Rshear


        psf_fitter = self._fit_stack(self._psf_stack, self._psf_skyvar, psf_ares, 
                                     ngauss=self.ngauss_psf, name='PSF')

        psf_gmix=psf_fitter.get_gmix()

        fitter=self._fit_stack(self._image_stack, self._image_skyvar, ares, psf=psf_gmix, 
                               ngauss=self.ngauss_obj, name='Images')

        self._fitter=fitter
        self._psf_fitter=psf_fitter


        res=self._fitter.get_result()
        sh1,sh2=0.5*res['pars'][2], 0.5*res['pars'][3]
        err1=0.5*res['perr'][2]
        err2=0.5*res['perr'][3]
        print 'fit:  sh1: %.6g +/- %.6g  sh2: %.6g +/- %.6g' % (sh1,err1,sh2,err2)
        print 'fit/Rshear: ',sh1/self.Rshear,sh2/self.Rshear

    def _run_uw(self, image, psf):
        import fimage
        res=fimage.fmom(image)
        pres=fimage.fmom(psf)

        irr = res['cov'][0]-pres['cov'][0]
        irc = res['cov'][1]-pres['cov'][1]
        icc = res['cov'][2]-pres['cov'][2]

        e1=(icc-irr)/(irr+icc)
        e2=2.*irc/(irr+icc)

        sh1=0.5*e1
        sh2=0.5*e2
        print 'uw:',sh1,sh2
        print 'uw/Rshear:',sh1/self.Rshear,sh2/self.Rshear
        return {'e1':e1, 'e2':e2}

    def _fit_stack(self, image_stack, skyvar, ares, psf=None, 
                        ngauss=1, name=''):
        return self._fit_stack_lm(image_stack, skyvar, ares, psf=psf, 
                                  ngauss=ngauss, name=name)
        #return self._fit_stack_mcmc(image_stack, skyvar, ares, psf=psf, 
        #                            ngauss=ngauss, name=name)

    def _fit_stack_mcmc(self, image_stack, skyvar, ares, psf=None, 
                        ngauss=1, name=''):


        if psf is not None:
            # MCM is a pure likelihood code, no gprior
            from .mcm_sim import MCM
            fitter=MCM(image_stack, 1./skyvar, psf, self.gprior, 'gauss', 
                       200, 200, 200, ares=ares)
        else:
            from gmix_image.gmix_mcmc import MixMCPSF
            fitter=MixMCPSF(image_stack, 1./skyvar, 'gauss',ares=ares)

        if self.verbose:
            res=fitter.get_result()
            print_pars(res['pars'], front='pars: ')
            print_pars(res['perr'], front='perr: ')

        return fitter


    def _fit_stack_lm(self, image_stack, skyvar, ares, psf=None, 
                      ngauss=1, name=''):

        row_guess=ares['wrow'] 
        col_guess=ares['wcol'] 
        e1guess=ares['e1']
        e2guess=ares['e2']
        Tguess=ares['Irr']+ares['Icc']
        counts=image_stack.sum()

        while True:
            if ngauss==3:
                prior=numpy.array( [row_guess+0.01*srandu(),
                                    col_guess+0.01*srandu(),
                                    e1guess+0.05*srandu(),
                                    e2guess+0.05*srandu(),
                                    Tguess*8.0*(1.+0.1*srandu()),
                                    Tguess*2.42*(1.+0.1*srandu()),
                                    Tguess*0.20*(1.+0.1*srandu()),
                                    counts*0.28*(1.+0.1*srandu()), 
                                    counts*0.34*(1.+0.1*srandu()), 
                                    counts*0.38*(1.+0.1*srandu())])
            elif ngauss==2:
                prior=numpy.array( [row_guess+0.01*srandu(),
                                    col_guess+0.01*srandu(),
                                    e1guess+0.05*srandu(),
                                    e2guess+0.05*srandu(),
                                    Tguess*3.0*(1.+0.1*srandu()),
                                    Tguess*0.9*(1.+0.1*srandu()),
                                    counts*0.6*(1.+0.1*srandu()), 
                                    counts*0.4*(1.+0.1*srandu())])

            elif ngauss==1:
                prior=numpy.array( [row_guess+0.01*srandu(),
                                    col_guess+0.01*srandu(),
                                    e1guess+0.05*srandu(),
                                    e2guess+0.05*srandu(),
                                    Tguess*(1.+0.1*srandu()),
                                    counts*(1.+0.1*srandu())] )
            else:
                raise ValueError("send ngauss 1,2,3")
            width=numpy.abs(prior)*1.e6

            fitter=gmix_image.gmix_fit.GMixFitCoellip(image_stack, sqrt(skyvar), prior, width, 
                                                      model='coellip',
                                                      nsub=self.nsub, psf=psf, verbose=False)
            flags=fitter.get_flags()

            if flags==0:
                break
            else:
                print flags

        res=fitter.get_result()
        print 'numiter:',res['numiter']
        if self.verbose:
            print_pars(prior,front='prior: ')
            print_pars(res['pars'], front='pars: ')
            print_pars(res['perr'], front='perr: ')

        if False:
            import images
            model=fitter.get_model()
            #images.multiview(im/im.max(), nonlinear=1)
            images.compare_images(image_stack/image_stack.max(), 
                                  model/model.max(),
                                  nonlinear=1,
                                  title=name)

            key=raw_input('hit a key: ')
            if key=='q':
                stop
        return fitter


    def _run_admom(self, image, skyvar):
        import admom

        cen_guess=numpy.array(image.shape,dtype='f8')/2.
        n=0
        while n < 100:
            ares = admom.admom(image, 
                               cen_guess[0]+0.1*srandu(), 
                               cen_guess[1]+0.1*srandu(),
                               sigsky=sqrt(skyvar),
                               guess=4.0*(1+0.1*srandu()),
                               nsub=self.nsub)

            if ares['whyflag'] == 0:
                break
            n += 1
        return ares 

    def _make_pair(self, g):

        # make sure these are totally random
        numpy.random.seed(None)
        theta1 = 360.0*numpy.random.random()
        theta2 = theta1 + 90.0

        imd1=self._make_image(g, theta1)
        imd2=self._make_image(g, theta2)

        return imd1,imd2

    def _make_image(self, g, theta):
        s2n_psf = self['s2n_psf']

        ellip=lensing.util.g2e(g)

        ci,ares = self._get_ci_ares_psf(self.s2, ellip, theta, self.s2n, s2n_psf)

        imd={'image':ci.image,
              'psf_image':ci.psf,
              'image_skyvar':ci['skysig']**2,
              'ivar':1./ci['skysig']**2,
              'psf_skyvar':ci['skysig_psf']**2,
              'psf_ivar':1./ci['skysig_psf']**2,
              'ares':ares,
              'model':self.simc['psfmodel']}

        return imd

    def _get_ci_ares_psf(self, s2, ellip, theta, s2n, s2n_psf):
        from fimage.convolved import NoisyConvolvedImage

        ci_nonoise = self.shapesim.get_trial(s2, ellip, theta)

        ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf,
                                 s2n_method=self['s2n_method'])
        # build up a fake ares 
        ares={'wrow':ci['cen_admom'][0],
              'wcol':ci['cen_admom'][1],
              'Irr':ci['cov_admom'][0],
              'Irc':ci['cov_admom'][1],
              'Icc':ci['cov_admom'][2],
              'e1':ci['e1_admom'],
              'e2':ci['e2_admom'],
              'whyflag':0}

        return ci, ares

    def _set_gprior(self):
        if self.simc['gprior_type']=='fits-vs-mag':
            import cluster_step
            exp_prior_pars=cluster_step.files.read_gprior(type='gexp')

            index=3
            self.gprior=cluster_step.prior.GPriorExp(exp_prior_pars['pars'][3])
            self.Rshear = RSHEAR_DEFAULT
            
            if index != 3:
                raise ValueError("re-calculate Rshear")
        elif self.simc['gprior_type']=='round':
            self.Rshear = 1.0
        else:
            raise ValueError("expected gprior type 'fits-vs-mag' or 'round'")

class StackSim(StackSimBase):
    """
    Send itrial= on construction

    Generate stacks for a single "split".  Inherits all methods from
    StackSimBase.  You will want to call 

        s.fit_stacks() and
        s.write()

    """
    def __init__(self, run, **keys):
        super(StackSim,self).__init__(run, **keys)

        if self.itrial is None:
            raise ValueError("send itrial= to genrate stacks")

    def generate_stacks(self):
        if self.simc['orient']=='ring':
            self._generate_ring_stacks()
        else:
            self._generate_random_stacks()

    def _generate_random_stacks(self):

        gvals=self._get_gvals()
        for i,g in enumerate(gvals):
            if (((i+1) % 100) == 0):
                print '%d/%d' % (i+1,self.nimage)

            theta=360*numpy.random.random()
            imd=self._make_image(g,theta)

            if i==0:
                image_stack  = numpy.zeros(imd['image'].shape)
                psf_stack    = numpy.zeros(imd['psf_image'].shape)
                image_skyvar = 0.0
                psf_skyvar   = 0.0

            image_stack += imd['image']
            image_skyvar += imd['image_skyvar']
            psf_stack += imd['psf_image']
            psf_skyvar += imd['psf_skyvar']
            

        self._image_stack=image_stack
        self._image_skyvar=image_skyvar
        self._psf_stack=psf_stack
        self._psf_skyvar=psf_skyvar

        if False:
            import images
            images.multiview(image_stack)


    def _generate_ring_stacks(self):

        gvals=self._get_gvals()
        for i,g in enumerate(gvals):
            if (((i+1) % 100) == 0):
                print '%d/%d' % (i+1,self.nimage)

            imd1,imd2=self._make_pair(g)

            if i==0:
                image_stack  = numpy.zeros(imd1['image'].shape)
                psf_stack    = numpy.zeros(imd1['psf_image'].shape)
                image_skyvar = 0.0
                psf_skyvar   = 0.0

            image_stack += imd1['image']
            image_stack += imd2['image']

            image_skyvar += imd1['image_skyvar']
            image_skyvar += imd2['image_skyvar']

            psf_stack += imd1['psf_image']
            psf_stack += imd2['psf_image']

            psf_skyvar += imd1['psf_skyvar']
            psf_skyvar += imd2['psf_skyvar']
            

        self._image_stack=image_stack
        self._image_skyvar=image_skyvar
        self._psf_stack=psf_stack
        self._psf_skyvar=psf_skyvar

        if False:
            import images
            images.multiview(image_stack)

    def _get_gvals(self):
        n=self.nimage

        if self.simc['gprior_type']=='round':
            print 'using round galaxies'
            gvals=numpy.zeros(n)
        else:
            print 'sampling gprior'
            gvals=self.gprior.sample1d(n)

        return gvals

class StackCombiner(StackSimBase):
    """

    Combine stacks for all splits. Inherits all methods from StackSimBase.  You
    will want to call 

        s.fit_stacks() and
        s.write()

    """
    def __init__(self, run, **keys):
        super(StackCombiner,self).__init__(run, **keys)

        if self.itrial is not None:
            raise ValueError("itrial %s was sent to a StackCombiner" % self.itrial)

    def load_splits(self):
        """
        Load and combine all splits for the input run,is2,is2n
        """
        nsplit=self['nsplit']

        for isplit in xrange(nsplit):
            data0=shapesim.read_output(self['run'], 
                                       self.is2, 
                                       self.is2n, 
                                       itrial=isplit,
                                       verbose=self.verbose)
            if isplit==0:
                image_stack=numpy.zeros(data0['image_stack'][0,:,:].shape)
                psf_stack=numpy.zeros(data0['psf_stack'][0,:,:].shape)
                image_skyvar=0.0
                psf_skyvar=0.0
                nimage=0
                nimage=0

            image_stack += data0['image_stack'][0,:,:]
            psf_stack += data0['psf_stack'][0,:,:]

            image_skyvar += data0['image_skyvar'][0]
            psf_skyvar += data0['psf_skyvar'][0]

            nimage += data0['nimage'][0]

        self._image_stack=image_stack
        self._image_skyvar=image_skyvar
        self._psf_stack=psf_stack
        self._psf_skyvar=psf_skyvar
        self.nimage=nimage



def load_admom_is2_stacks(run, is2):
    """
    load stacks for the input is2
    """
    conf=read_config(run)
    simc = read_config(conf['sim'])

    num_s2 = simc['nums2']
    if is2 < 0 or is2 > (num_s2-1):
        raise ValueError("s2 %d out of range: [%d,%d)" % \
                         (is2, 0, num_s2))
    num_s2n=len(conf['s2nvals'])

    data=[]
    for is2n in xrange(num_s2n):
        data0 = shapesim.read_output(run, is2, is2n, verbose=True)

        data.append(data0)

    data=eu.numpy_util.combine_arrlist(data)
    return data



class AdmomSim(StackSimBase):
    """
    Send itrial= on construction

    currently assumes a ring test
    """
    def __init__(self, run, **keys):
        super(AdmomSim,self).__init__(run, **keys)

        if self.itrial is None:
            raise ValueError("send itrial= to genrate stacks")

    def run(self):
        import lensing
        import admom
        i=0
        j=0

        e1vals=numpy.zeros(self.nimage)
        e2vals=numpy.zeros(self.nimage)
        gvals=numpy.zeros(self.nimage/2)

        while i < self.nimage:
            if (((i+1) % 100) == 0):
                print '%d/%d' % (i+1,self.nimage)

            g=self.gprior.sample1d(1)

            imd1,imd2=self._make_pair(g[0])

            ares1 = self._run_admom(imd1['image'], imd1['image_skyvar'])
            if ares1['whyflag'] != 0:
                print 'bad1'
                continue

            psfres1 = self._run_admom(imd1['psf_image'], imd1['psf_skyvar'])
            if psfres1['whyflag'] != 0:
                print 'badpsf1'
                continue

            ares2 = self._run_admom(imd2['image'], imd2['image_skyvar'])
            if ares2['whyflag'] != 0:
                print 'bad2'
                continue

            psfres2 = self._run_admom(imd2['psf_image'], imd2['psf_skyvar'])
            if psfres2['whyflag'] != 0:
                print 'badpsf2'
                continue

            gvals[i] = g

            for ares,pres in zip([ares1,ares2],[psfres1,psfres2]):
                T=ares['Irr']+ares['Icc']
                Tpsf=pres['Irr'] + pres['Icc']
                e1,e2,R,flags=admom.correct(T, ares['e1'], ares['e2'],ares['a4'],
                                            Tpsf, pres['e1'], pres['e2'], pres['a4'])
                if flags != 0:
                    raise ValueError("compea4 failed")

                e1vals[j] = e1
                e2vals[j] = e2
                j += 1
            

            i += 1

        data=numpy.zeros(e1vals.size, dtype=[('e1','f8'),('e2','f8')])
        data['e1'] = e1vals
        data['e2'] = e2vals
        self._data=data

        self.calc_shear()

    def calc_shear(self):
        e1mean=self._data['e1'].mean()
        e2mean=self._data['e2'].mean()

        g1samp,g2samp=self.gprior.sample2d(100000)
        e1samp,e2samp = lensing.util.g1g2_to_e1e2(g1samp,g2samp)

        e1var=e1samp.var()
        Rsh = 1.0 - e1var

        sh1 = 0.5*e1mean/Rsh
        sh2 = 0.5*e2mean/Rsh

        sh1err=0.5*self._data['e1'].std()/sqrt(self._data.size)
        sh2err=0.5*self._data['e2'].std()/sqrt(self._data.size)
        print 'shear1: %s +/- %s  shear2: %s +/- %s' % (sh1,sh1err,sh2,sh2err)

        self.sh1=sh1
        self.sh2=sh2
        self.sh1err=sh1err
        self.sh2err=sh2err

class AdmomCombiner(AdmomSim):
    """

    Combine stacks for all splits. Inherits all methods from StackSimBase.  You
    will want to call 

        s.load_splits()
        s.write()

    """
    def __init__(self, run, **keys):
        # fake itrial here
        super(AdmomCombiner,self).__init__(run, itrial=1, **keys)

        self.itrial=None

    def load_splits(self):
        """
        Load and combine all splits for the input run,is2,is2n
        """
        nsplit=self['nsplit']

        data_list=[]
        for isplit in xrange(nsplit):
            data0=shapesim.read_output(self['run'], 
                                       self.is2, 
                                       self.is2n, 
                                       itrial=isplit,
                                       verbose=self.verbose)
            data_list.append(data0)

        self._data=eu.numpy_util.combine_arrlist(data_list)
