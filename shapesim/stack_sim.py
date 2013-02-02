"""

unweighted seems to be unbiased when stacked!  Maybe we don't care
that they are noisy when we may be dominated by other effects.  It *is*
much noiser, don't get me wrong... radically noisier
    
    - run with more trials to see if it is actually biased

    - try trimming to 4-sigma to reduce noise.  Might be more biased since this
    is effectively a weight function.  Looks biased low, perhaps because of
    the gaussian weighting underestimating the needed region.

    WOW! much much less noisy, but way biased :)

    - try same aperture for both object and psf
        - circular 4 sigma - looks biased
        - circular 5 sigma - looks decent
            get05r02
        - might be circularity: 4 sigma box
        - 5 sigma box
    - try fixed circular/box aperture
        - try 12.5 (image radius) 
        - try 10.0
    - try fixed square aperture

    - try with same exact fixed aperture for both galaxies and stars.  Most
    galaxies are using radius of ~9, so let's set it to 10

    - found bug in cluster step version of this, fixed part of it; see how the
    result looks.  Don't expect it to be great: there are much bigger problems
    there.


For fits and shear=[0.04,,0.00] Getting <e> of 0.0785 at all s/n instead of
0.08 For unweighted moments corrected by 1/Rshear I think it does better
(although noisier).  The residuals of model and image look very good.  See
similar bias for "et" model as gauss model.

    - is it straight percentage?  If I increase the shear is the ratio
    off by the same?  Why is no sensitivity correction needed, because
    we do a fit instead of moments?

    * max like bias some how? no, tried mcmc
    * working in e instead of g?  no
    * pixel problem

        * made both psf and object *much* bigger but got same bias!  But making
        the object bigger relative to the PSF did change the bias.

        * made objects round and I see no bias (and no need for Rshear)

        * for gauss, there is no bias if I make the shear pure e2.  If there is
        any component that is e1 I see a bias in both shear1 and shear2.  For "et"
        I still saw the bias for pure e2.

        * moved the center around but got same bias

        * maybe nsub=16 isn't good enough in sim?  Trying 24.  Nope.

        * maybe it is fortran itself? nope
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
                    xlog=False,yrange=None, ignore_missing=False):
    import biggles
    data=load_is2_stacks(run,is2)
    
    if is2 is None or type is None:
        raise ValueError("send is2= and type=")
    if type=='sums':
        sh1=0.5*data['e1_bysum']/data['Rshear']
        sh2=0.5*data['e2_bysum']/data['Rshear']
        err1=0.16/sqrt(data['nimage'])
        err2=err1
    elif type=='fmom':
        if 'e1_fmom' in data.dtype.names:
            sh1=0.5*data['e1_fmom']/data['Rshear']
            sh2=0.5*data['e2_fmom']/data['Rshear']
        else:
            sh1=0.5*data['e1_uw']/data['Rshear']
            sh2=0.5*data['e2_uw']/data['Rshear']
        err1=0.16/sqrt(data['nimage'])
        err2=err1
    elif type=='fit':
        #sh1=0.5*data['pars'][:,2]/RSHEAR_DEFAULT
        #sh2=0.5*data['pars'][:,3]/RSHEAR_DEFAULT
        # wierd, this gives a closer answer for ngauss=2,3!
        sh1=0.5*data['pars'][:,2]
        sh2=0.5*data['pars'][:,3]
        err1=0.5*data['perr'][:,2]
        err2=0.5*data['perr'][:,3]
    elif type=='admom':
        sh1=0.5*data['e1_admom']/data['Rshear']
        sh2=0.5*data['e2_admom']/data['Rshear']
        err1=0.5*data['perr'][:,2]
        err2=0.5*data['perr'][:,3]
    else:
        raise ValueError("bad type: '%s'" % type)


    plt=biggles.FramedPlot()
    plt.xlog=xlog

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
    plt.xrange=[0.9*data['s2n'].min(), 1.1*data['s2n'].max()]
    plt.show()

def load_is2_stacks(run, is2, ignore_missing=False):
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
        fname=shapesim.get_output_url(run, is2, is2n)
        if not os.path.exists(fname):
            if not ignore_missing:
                raise ValueError("file not found: %s" % fname)
            else:
                continue

        data0 = shapesim.read_output(run, is2, is2n, verbose=True)
        data.append(data0)

    if len(data) == 0:
        raise ValueError("no data read")

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

        sdata=self._stack_data

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

        data['e1_fmom'] = self.res_fmom['e1']
        data['e2_fmom'] = self.res_fmom['e2']
        data['e1_admom'] = self._ares['e1corr']
        data['e2_admom'] = self._ares['e2corr']

        if 'imsum' in sdata:
            data['imsum'] = sdata['imsum']
            """
            data['rowcen2_sum'] = sdata['rowcen2_sum']
            data['rowcolcen_sum'] = sdata['rowcolcen_sum']
            data['colcen2_sum'] = sdata['colcen2_sum']
            """
            data['row2_sum'] = sdata['row2_sum']
            data['rowcol_sum'] = sdata['rowcol_sum']
            data['col2_sum'] = sdata['col2_sum']

            data['psf_imsum'] = sdata['psf_imsum']
            """
            data['psf_rowcen2_sum'] = sdata['psf_rowcen2_sum']
            data['psf_rowcolcen_sum'] = sdata['psf_rowcolcen_sum']
            data['psf_colcen2_sum'] = sdata['psf_colcen2_sum']
            """
            data['psf_row2_sum'] = sdata['psf_row2_sum']
            data['psf_rowcol_sum'] = sdata['psf_rowcol_sum']
            data['psf_col2_sum'] = sdata['psf_col2_sum']

            data['e1_bysum'] = self.res_sums['e1']
            data['e2_bysum'] = self.res_sums['e2']

        data['s2n_stack'] = s2n
        data['image_stack'][0,:,:] = sdata['image_stack']
        data['image_skyvar'] = sdata['image_skyvar']

        data['psf_s2n_stack'] = psf_s2n
        data['psf_stack'][0,:,:] = sdata['psf_stack']
        data['psf_skyvar'] = sdata['psf_skyvar']

        data['Rshear'] = self.Rshear
        self._data=data


    def _get_struct(self, npars, psf_npars):
        shape=self._stack_data['image_stack'].shape
        psf_shape=self._stack_data['psf_stack'].shape

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

            ('imsum','f8'),
            #('rowcen2_sum','f8'),
            #('rowcolcen_sum','f8'),
            #('colcen2_sum','f8'),
            ('row2_sum','f8'),
            ('rowcol_sum','f8'),
            ('col2_sum','f8'),

            ('psf_imsum','f8'),
            #('psf_rowcen2_sum','f8'),
            #('psf_rowcolcen_sum','f8'),
            #('psf_colcen2_sum','f8'),
            ('psf_row2_sum','f8'),
            ('psf_rowcol_sum','f8'),
            ('psf_col2_sum','f8'),

            ('e1_bysum','f8'),
            ('e2_bysum','f8'),


            ('e1_fmom','f8'),
            ('e2_fmom','f8'),
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

        sdata=self._stack_data

        self.res_fmom=self._run_fmom(sdata['image_stack'], sdata['psf_stack'])
        if 'imsum' in sdata:
            self.res_sums=self._run_sums_shear()

        ares = self._run_admom(sdata['image_stack'], sdata['image_skyvar'])
        psf_ares = self._run_admom(sdata['psf_stack'], sdata['psf_skyvar'])
        
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


        psf_fitter = self._fit_stack(sdata['psf_stack'], sdata['psf_skyvar'], psf_ares, 
                                     ngauss=self.ngauss_psf, name='PSF')

        psf_gmix=psf_fitter.get_gmix()

        fitter=self._fit_stack(sdata['image_stack'], sdata['image_skyvar'], ares, psf=psf_gmix, 
                               ngauss=self.ngauss_obj, name='Images')

        self._fitter=fitter
        self._psf_fitter=psf_fitter


        res=self._fitter.get_result()
        sh1,sh2=0.5*res['pars'][2], 0.5*res['pars'][3]
        err1=0.5*res['perr'][2]
        err2=0.5*res['perr'][3]
        fmt='%15.6f'
        print_pars(res['pars'], front='pars: ', fmt='%15.6f')
        print_pars(res['perr'], front='perr: ', fmt='%15.6f')
        print 'fit:  sh1: %.6g +/- %.6g  sh2: %.6g +/- %.6g' % (sh1,err1,sh2,err2)
        print 'fit/Rshear: ',sh1/self.Rshear,sh2/self.Rshear

    def _run_fmom(self, image, psf):
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
        print 'fmom:       ',sh1,sh2
        print 'fmom/Rshear:',sh1/self.Rshear,sh2/self.Rshear
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
            if ngauss==4:
                prior=numpy.array( [row_guess+0.01*srandu(),
                                    col_guess+0.01*srandu(),
                                    e1guess+0.05*srandu(),
                                    e2guess+0.05*srandu(),
                                    Tguess*16.0*(1.+0.1*srandu()),
                                    Tguess*8.0*(1.+0.1*srandu()),
                                    Tguess*2.42*(1.+0.1*srandu()),
                                    Tguess*0.20*(1.+0.1*srandu()),
                                    counts*0.10*(1.+0.1*srandu()), 
                                    counts*0.20*(1.+0.1*srandu()), 
                                    counts*0.25*(1.+0.1*srandu()), 
                                    counts*0.31*(1.+0.1*srandu())])

            elif ngauss==3:
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
             'model':self.simc['psfmodel'],
             'flags':0}

        ares = self._run_admom(ci.image, ci['skysig']**2)

        if ares['whyflag'] != 0:
            imd['flags'] = ares['whyflag']
            return imd

        psf_ares = self._run_admom(ci.psf, ci['skysig_psf']**2)
        if psf_ares['whyflag'] != 0:
            imd['flags'] = psf_ares['whyflag']
            return imd

        imd['ares'] = ares
        imd['psf_ares'] = psf_ares

        rmax=self._get_sums_aperture(ares, psf_ares)
        sums=self._get_sums(ci.image, ares, rmax)
        imd.update(sums)

        psf_sums=self._get_sums(ci.psf, psf_ares,rmax)

        for k in psf_sums:
            imd['psf_'+k] = psf_sums[k]

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

    def _get_sums_aperture(self, ares, psf_ares):
        if 'trim_nsig' in self:
            nsig=self['trim_nsig']
            rad=sqrt( max(ares['Irr'],ares['Icc']) )
            rad_psf=sqrt( max(psf_ares['Irr'],psf_ares['Icc']) )
            rad=max(rad, rad_psf)

            nsig=float(nsig)

            rmax = nsig*rad
            if rmax < 4:
                rmax=4

        elif 'trim_rad' in self:
            rmax=float(self['trim_rad'])
        else:
            rmax=1.e9

        minrad=self.get('trim_min_rad',4)
        if rmax < minrad:
            rmax=minrad
        return rmax

    def _get_sums(self, image, ares, rmax):

        row0=ares['wrow']
        col0=ares['wcol']

        dims=image.shape
        row,col=numpy.mgrid[0:dims[0], 0:dims[1]]

        row = row.astype('f8') - row0
        col = col.astype('f8') - col0

        rad2=row**2 + col**2

        rmax2=rmax**2

        w=numpy.where(rad2 <= rmax2)
        if w[0].size == 0:
            raise ValueError("rmax %s too small, no pixels" % rmax)

        #print rmax,image.shape[0]/2.

        image=image[w]
        row=row[w]
        col=col[w]

        imsum=image.sum()
        row2_sum = (image*row**2).sum()
        rowcol_sum = (image*row*col).sum()
        col2_sum = (image*col**2).sum()

        return {'imsum':imsum,
                'row2_sum':row2_sum,
                'rowcol_sum':rowcol_sum,
                'col2_sum':col2_sum}

    def _run_sums_shear(self):
        
        sd=self._stack_data
        psf_imsum=sd['psf_imsum']
        psf_irr = sd['psf_row2_sum']/psf_imsum
        psf_irc = sd['psf_rowcol_sum']/psf_imsum
        psf_icc = sd['psf_col2_sum']/psf_imsum

        imsum=sd['imsum']
        irrO = sd['row2_sum']/imsum
        ircO = sd['rowcol_sum']/imsum
        iccO = sd['col2_sum']/imsum

        irr = irrO - psf_irr
        irc = ircO - psf_irc
        icc = iccO - psf_icc

        e1=(icc-irr)/(irr+icc)
        e2=2.*irc/(irr+icc)

        sh1=0.5*e1
        sh2=0.5*e2
        print 'sums:       ',sh1,sh2
        print 'sums/Rshear:',sh1/self.Rshear,sh2/self.Rshear
        return {'e1':e1, 'e2':e2}


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

        sdata={}

        gvals=self._get_gvals()

        i=0
        while i < self.nimage:
            if (((i+1) % 100) == 0):
                print '%d/%d' % (i+1,self.nimage)

            g=gvals[i]
            theta=360*numpy.random.random()
            imd=self._make_image(g,theta)

            if imd['flags'] != 0:
                print 'bad'
                continue

            if i==0:
                sdata['image_stack']  = numpy.zeros(imd['image'].shape)
                sdata['psf_stack']    = numpy.zeros(imd['psf_image'].shape)
                sdata['image_skyvar'] = 0.0
                sdata['psf_skyvar']   = 0.0

                sdata['imsum']=0.0
                """
                sdata['rowcen2_sum'] = 0.0
                sdata['rowcolcen_sum'] = 0.0
                sdata['colcen2_sum'] = 0.0
                """

                sdata['row2_sum'] = 0.0
                sdata['rowcol_sum'] = 0.0
                sdata['col2_sum'] = 0.0

                sdata['psf_imsum']=0.0
                """
                sdata['psf_rowcen2_sum'] = 0.0
                sdata['psf_rowcolcen_sum'] = 0.0
                sdata['psf_colcen2_sum'] = 0.0
                """

                sdata['psf_row2_sum'] = 0.0
                sdata['psf_rowcol_sum'] = 0.0
                sdata['psf_col2_sum'] = 0.0




            sdata['image_stack'] += imd['image']
            sdata['image_skyvar'] += imd['image_skyvar']
            sdata['psf_stack'] += imd['psf_image']
            sdata['psf_skyvar'] += imd['psf_skyvar']

            sdata['imsum'] += imd['imsum']
            """
            sdata['rowcen2_sum'] += imd['rowcen2_sum']
            sdata['rowcolcen_sum'] += imd['rowcolcen_sum']
            sdata['colcen2_sum'] += imd['colcen2_sum']
            """

            sdata['row2_sum'] += imd['row2_sum']
            sdata['rowcol_sum'] += imd['rowcol_sum']
            sdata['col2_sum'] += imd['col2_sum']

            sdata['psf_imsum'] += imd['psf_imsum']
            """
            sdata['psf_rowcen2_sum'] += imd['psf_rowcen2_sum']
            sdata['psf_rowcolcen_sum'] += imd['psf_rowcolcen_sum']
            sdata['psf_colcen2_sum'] += imd['psf_colcen2_sum']
            """

            sdata['psf_row2_sum'] += imd['psf_row2_sum']
            sdata['psf_rowcol_sum'] += imd['psf_rowcol_sum']
            sdata['psf_col2_sum'] += imd['psf_col2_sum']


            i+= 1
           

        self._stack_data=sdata

        if False:
            import images
            images.multiview(image_stack)


    def _generate_ring_stacks(self):

        sdata={}
        gvals=self._get_gvals()
        
        i=0
        while i < self.nimage:
            if (((i+1) % 100) == 0):
                print '%d/%d' % (i+1,self.nimage)

            g=gvals[i]
            imd1,imd2=self._make_pair(g)
            if imd1['flags'] != 0 or imd2['flags'] != 0:
                print 'bad'
                continue

            if i==0:
                sdata['image_stack']  = numpy.zeros(imd1['image'].shape)
                sdata['psf_stack']    = numpy.zeros(imd1['psf_image'].shape)
                sdata['image_skyvar'] = 0.0
                sdata['psf_skyvar']   = 0.0

                sdata['imsum']=0.0
                """
                sdata['rowcen2_sum'] = 0.0
                sdata['rowcolcen_sum'] = 0.0
                sdata['colcen2_sum'] = 0.0
                """

                sdata['row2_sum'] = 0.0
                sdata['rowcol_sum'] = 0.0
                sdata['col2_sum'] = 0.0

                sdata['psf_imsum']=0.0
                """
                sdata['psf_rowcen2_sum'] = 0.0
                sdata['psf_rowcolcen_sum'] = 0.0
                sdata['psf_colcen2_sum'] = 0.0
                """

                sdata['psf_row2_sum'] = 0.0
                sdata['psf_rowcol_sum'] = 0.0
                sdata['psf_col2_sum'] = 0.0



            sdata['image_stack'] += imd1['image']
            sdata['image_skyvar'] += imd1['image_skyvar']
            sdata['psf_stack'] += imd1['psf_image']
            sdata['psf_skyvar'] += imd1['psf_skyvar']


            sdata['imsum'] += imd1['imsum']
            """
            sdata['rowcen2_sum'] += imd1['rowcen2_sum']
            sdata['rowcolcen_sum'] += imd1['rowcolcen_sum']
            sdata['colcen2_sum'] += imd1['colcen2_sum']
            """

            sdata['row2_sum'] += imd1['row2_sum']
            sdata['rowcol_sum'] += imd1['rowcol_sum']
            sdata['col2_sum'] += imd1['col2_sum']

            sdata['psf_imsum'] += imd1['psf_imsum']
            """
            sdata['psf_rowcen2_sum'] += imd1['psf_rowcen2_sum']
            sdata['psf_rowcolcen_sum'] += imd1['psf_rowcolcen_sum']
            sdata['psf_colcen2_sum'] += imd1['psf_colcen2_sum']
            """

            sdata['psf_row2_sum'] += imd1['psf_row2_sum']
            sdata['psf_rowcol_sum'] += imd1['psf_rowcol_sum']
            sdata['psf_col2_sum'] += imd1['psf_col2_sum']



            sdata['image_stack'] += imd2['image']
            sdata['image_skyvar'] += imd2['image_skyvar']
            sdata['psf_stack'] += imd2['psf_image']
            sdata['psf_skyvar'] += imd2['psf_skyvar']

            sdata['imsum'] += imd2['imsum']
            """
            sdata['rowcen2_sum'] += imd2['rowcen2_sum']
            sdata['rowcolcen_sum'] += imd2['rowcolcen_sum']
            sdata['colcen2_sum'] += imd2['colcen2_sum']
            """

            sdata['row2_sum'] += imd2['row2_sum']
            sdata['rowcol_sum'] += imd2['rowcol_sum']
            sdata['col2_sum'] += imd2['col2_sum']
 
            sdata['psf_imsum'] += imd2['psf_imsum']
            """
            sdata['psf_rowcen2_sum'] += imd2['psf_rowcen2_sum']
            sdata['psf_rowcolcen_sum'] += imd2['psf_rowcolcen_sum']
            sdata['psf_colcen2_sum'] += imd2['psf_colcen2_sum']
            """

            sdata['psf_row2_sum'] += imd2['psf_row2_sum']
            sdata['psf_rowcol_sum'] += imd2['psf_rowcol_sum']
            sdata['psf_col2_sum'] += imd2['psf_col2_sum']

            i += 1
 
        
        self._stack_data=sdata

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

        sdata={}
        for isplit in xrange(nsplit):
            data0=shapesim.read_output(self['run'], 
                                       self.is2, 
                                       self.is2n, 
                                       itrial=isplit,
                                       verbose=self.verbose)
            if isplit==0:
                sdata['image_stack']  = numpy.zeros(data0['image_stack'][0,:,:].shape)
                sdata['psf_stack']    = numpy.zeros(data0['psf_stack'][0,:,:].shape)

                sdata['image_skyvar'] = 0.0
                sdata['psf_skyvar']   = 0.0

                nimage = 0
                if 'imsum' in data0.dtype.names:
                    sdata['imsum']=0.0
                    """
                    sdata['rowcen2_sum'] = 0.0
                    sdata['rowcolcen_sum'] = 0.0
                    sdata['colcen2_sum'] = 0.0
                    """

                    sdata['row2_sum'] = 0.0
                    sdata['rowcol_sum'] = 0.0
                    sdata['col2_sum'] = 0.0

                    sdata['psf_imsum']=0.0
                    """
                    sdata['psf_rowcen2_sum'] = 0.0
                    sdata['psf_rowcolcen_sum'] = 0.0
                    sdata['psf_colcen2_sum'] = 0.0
                    """

                    sdata['psf_row2_sum'] = 0.0
                    sdata['psf_rowcol_sum'] = 0.0
                    sdata['psf_col2_sum'] = 0.0


            sdata['image_stack'] += data0['image_stack'][0,:,:]
            sdata['psf_stack'] += data0['psf_stack'][0,:,:]
            sdata['image_skyvar'] += data0['image_skyvar'][0]
            sdata['psf_skyvar'] += data0['psf_skyvar'][0]

            nimage += data0['nimage'][0]

            if 'imsum' in data0.dtype.names:
                sdata['imsum'] += data0['imsum'][0]
                """
                sdata['rowcen2_sum'] += data0['rowcen2_sum'][0]
                sdata['rowcolcen_sum'] += data0['rowcolcen_sum'][0]
                sdata['colcen2_sum'] += data0['colcen2_sum'][0]
                """

                sdata['row2_sum'] += data0['row2_sum'][0]
                sdata['rowcol_sum'] += data0['rowcol_sum'][0]
                sdata['col2_sum'] += data0['col2_sum'][0]

                sdata['psf_imsum'] += data0['psf_imsum'][0]
                """
                sdata['psf_rowcen2_sum'] += data0['psf_rowcen2_sum'][0]
                sdata['psf_rowcolcen_sum'] += data0['psf_rowcolcen_sum'][0]
                sdata['psf_colcen2_sum'] += data0['psf_colcen2_sum'][0]
                """

                sdata['psf_row2_sum'] += data0['psf_row2_sum'][0]
                sdata['psf_rowcol_sum'] += data0['psf_rowcol_sum'][0]
                sdata['psf_col2_sum'] += data0['psf_col2_sum'][0]


        self.nimage=nimage
        self._stack_data=sdata



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
