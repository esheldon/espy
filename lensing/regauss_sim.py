from __future__ import print_function
import os
import sys
from sys import stderr
import numpy
from numpy import sqrt, array, ogrid, where, median
from numpy.random import standard_normal
import images

import admom
import fimage
import glob

import esutil as eu
from esutil.ostools import path_join
from esutil.misc import wlog
from esutil import hdfs
import biggles

import pcolors

from pprint import pprint
import converter

try:
    import scipy.signal
except:
    print("Could not import scipy.signal, cannot do convolutions")


def shear_fracdiff(e, em, deriv=1.0):
    """
    e=etrue
    em=emeasured

    Hirata & Seljak eq 27
    putting in 1 for derivative d(emeas)/d(etrue)

    e=etrue
    em=emeas
    deriv deriviative of measured e with true e

    """
    return ((1-e**2)*deriv + em/e)/(2-em**2) - 1.0

def shear_fracdiff_stupid(e, em):
    return em/e-1.0


def plot_many():
    runs=['%02i' % r for r in xrange(12)]
    for run in runs:
        rsp = RegaussSimPlotter(run)

        for Rmin in [0.0, 0.33]:
            for yrange in [None, [-0.1,0.1]]:
                for dotitle in [False,True]:
                    rsp.plot('regauss',show=False, Rmin=Rmin,yrange=yrange,dotitle=dotitle)
                    rsp.plot('noweight',show=False, Rmin=Rmin,yrange=yrange,dotitle=dotitle)
                    rsp.plot('am+',show=False, Rmin=Rmin,yrange=yrange,dotitle=dotitle)

                    rsp.plot_with_alt('am+',show=False, Rmin=Rmin,yrange=yrange,dotitle=dotitle)
                    rsp.plot_with_alt('noweight',show=False, Rmin=Rmin,yrange=yrange,dotitle=dotitle)




class RegaussSimPlotter(dict):
    def __init__(self, run):

        c = read_config(run)
        self['run'] = run
        self.config = c
        self['objmodel'] = c['objmodel']
        self['psfmodel'] = c['psfmodel']
        self['psf_sigma'] = c['psf_sigma']

        self['psf_ellip'] = None
        if self['psfmodel'] != 'sdss':
            self['psf_ellip'] = c['psf_ellip']

        self.read_data()

        pdir = get_plotdir(self['run'])
        if not os.path.exists(pdir):
            print("Making plot dir:",pdir)
            os.makedirs(pdir)

    def plotfile(self, run, objmodel, psfmodel, type, alt=None, 
                 s2=None, Rmin=None, dotitle=False, yrange=None):
        dir = get_plotdir(run)
        f = 'rgsim-%sobj-%spsf-%s' % (objmodel, psfmodel, type)
        if alt is not None:
            f += '-alt'+alt
        if s2 is not None:
            f += '-s2-%0.3f' % s2
        if Rmin is not None:
            f += '-Rmin%0.3f' % Rmin
        if dotitle:
            f += '-tit'

        if yrange is not None:
            f += '-yr%0.3f-%0.3f' % tuple(yrange)
        f += '.eps'
        f = path_join(dir, f)
        return f

    def plot(self, type, Rmin=0.0, show=True, yrange=None, dotitle=False, 
             reduce_plot_key=True):
        """
        Send Rmin to limit the plotted R values
        if too many s2 use reduce_plot_key=False
        """
        pfile = self.plotfile(self['run'],
                              self['objmodel'],
                              self['psfmodel'], 
                              type,
                              Rmin=Rmin,
                              yrange=yrange, 
                              dotitle=dotitle)

        if 'R_rg' in self.alldata[0].dtype.names:
            print("using R_rg")
            rname = 'R_rg'
        else:
            rname = 'R'
        #rname='R'
        keepdata = []
        for st in self.alldata:
            if numpy.median(st[rname]) > Rmin:
                keepdata.append(st)


        #keepdata = self.alldata
        ndata = len(keepdata)
        colors=pcolors.rainbow(ndata, 'hex')

        biggles.configure('PlotKey','key_vsep',1.0)
        plt = biggles.FramedPlot()
        plt.aspect_ratio=1
        plt.xlabel='object ellipticity'
        plt.ylabel=r'$\Delta \gamma/\gamma$'

        allplots=[]
        i=0
        for st in keepdata:
            # this s2 is the value we were aiming for, could be pretty
            # far off for some models
            if 's2noweight' in st.dtype.names:
                s2 = st['s2noweight'][0]
            else:
                s2 = st['s2'][0]

            # this "etrue" is adaptive moments of pre-psf image
            s = st['etrue'].argsort()
            etrue = st['etrue'][s]

            if type == 'regauss':
                emeas = st['ecorr_rg'][s]
            elif type == 'am+':
                emeas = st['ecorr'][s]
            elif type == 'noweight':
                emeas = st['ecorr_uw'][s]
            else:
                raise ValueError("type should be 'regauss','am+', or 'noweight'")

            gamma_frac_rg = shear_fracdiff(etrue, emeas, deriv=1.0)

            Rmean = numpy.median( st[rname] )

            label = 's2: %0.3f R: %0.3f' % (s2,Rmean)

            crg = biggles.Curve(etrue, gamma_frac_rg, color=colors[i])
            crg.label=label

            plt.add(crg)
            
            allplots.append(crg)
            i += 1

        if dotitle:
            title='obj: %s psf: %s run: %s' \
                    % (self['objmodel'],self['psfmodel'],self['run'])

            if 'forcegauss' in self.config:
                if self.config['forcegauss']:
                    title += ' forcegauss'
            plt.title=title


        fsize=1.5
        if not reduce_plot_key:
            key = biggles.PlotKey(0.9,0.9, allplots, halign='right', fontsize=fsize)
        else:
            # pick a few
            nplot=len(allplots)
            tplots = [allplots[0], 
                      allplots[nplot*1/4], 
                      allplots[nplot/2], 
                      allplots[nplot*3/4], 
                      allplots[-1]]
            key = biggles.PlotKey(0.9,0.9, tplots, halign='right', fontsize=fsize)

        plt.add(key)

        plab='%s %s %s' % (type,self['objmodel'],self['psfmodel'])
        l = biggles.PlotLabel(0.1,0.9, plab, halign='left')
        plt.add(l)

        siglab=r'$\sigma_{PSF}: %.1f$ pix' % self['psf_sigma']
        if 's2n' in self.config:
            siglab+=r'$ S/N: %(s2n)d N_{trial}: %(ntrial)d$' % self.config
        elif 'ntrial' in self.config:
            siglab+=r'$  N_{trial}: %(ntrial)d$' % self.config
        sl = biggles.PlotLabel(0.1,0.1, siglab, halign='left')
        plt.add(sl)

        if not reduce_plot_key:
            plt.xrange = [0,1.4]
        if yrange is not None:
            plt.yrange = yrange
        print("Writing plot file:",pfile)
        if show:
            plt.show()
        plt.write_eps(pfile)
        converter.convert(pfile,dpi=100)




    def plot_with_alt(self, alt, Rmin=0.0, show=True, yrange=None, dotitle=False):
        '''
        Plot the regauss with another for comparison, e.g. noweight or AM+
        '''
        pfile = self.plotfile(self['run'],
                              self['objmodel'],
                              self['psfmodel'], 
                              'regauss-'+alt,
                              alt=alt, Rmin=Rmin, yrange=yrange,
                              dotitle=dotitle)

        keepdata = []
        for st in self.alldata:
            if numpy.median(st['R']) > Rmin:
                keepdata.append(st)


        #keepdata = self.alldata
        ndata = len(keepdata)
        colors=pcolors.rainbow(ndata, 'hex')

        allplots=[]
        i=0
        for st in keepdata:
            # this s2 is the value we were aiming for, could be pretty
            # far off for some models
            if 's2noweight' in st.dtype.names:
                s2 = st['s2noweight'][0]
            else:
                s2 = st['s2'][0]

            # this "etrue" is adaptive moments of pre-psf image
            s = st['etrue'].argsort()
            etrue = st['etrue'][s]

            # the top plot will be the regular analytic correction or the
            # unweighted correction. 
            if alt == 'noweight':
                gamma_frac_alt = shear_fracdiff(etrue, st['ecorr_uw'][s])
            elif alt == 'am+':
                gamma_frac_alt = shear_fracdiff(etrue, st['ecorr'][s])
            else:
                raise ValueError("alt should be am+ or noweight")

            gamma_frac_rg = shear_fracdiff(etrue, st['ecorr_rg'][s])

            Rmean = numpy.median( st['R'] )

            if Rmean >= Rmin:
                pdict={}
                label = 's2: %0.3f R: %0.3f' % (s2,Rmean)

                calt = biggles.Curve(etrue, gamma_frac_alt, color=colors[i])
                calt.label = label 
                pdict['calt'] = calt

                crg = biggles.Curve(etrue, gamma_frac_rg, color=colors[i])
                crg.label=label
                pdict['crg'] = crg
                
                allplots.append(pdict)
                i += 1

        conv=self.config.get('sim_conv', 'analytic')
        if dotitle:
            title='obj: %s psf: %s run: %s' \
                    % (self['objmodel'],self['psfmodel'],self['run'])

            if 'forcegauss' in self.config:
                if self.config['forcegauss']:
                    title += ' forcegauss'
        else:
            title=None

        arr=biggles.FramedArray(2,1, title=title)
        arr.xlabel='object ellipticity'
        arr.ylabel=r'$\Delta \gamma/\gamma$'


        allc1 = []
        allc2 = []
        i=0
        for p in allplots:
            arr[0,0].add( p['calt'] )
            arr[1,0].add( p['crg'] )
            if i < (ndata//2):
                allc1.append(p['calt'])
            else:
                allc2.append(p['calt'])
            i+=1

        fsize=2
        key1 = biggles.PlotKey(0.95,0.9, allc1, halign='right', fontsize=fsize)
        arr[0,0].add(key1)
        key2 = biggles.PlotKey(0.95,0.9, allc2, halign='right', fontsize=fsize)
        arr[1,0].add(key2)

        if alt == 'noweight':
            l1 = biggles.PlotLabel(0.1,0.9,'UnWeighted', halign='left')
        else:
            l1 = biggles.PlotLabel(0.1,0.9,'AM+', halign='left')
        l2 = biggles.PlotLabel(0.1,0.9,'ReGauss',  halign='left')
        arr[0,0].add(l1)
        arr[1,0].add(l2)

        arr.uniform_limits = 1

        arr.xrange = [0,1.2]
        if yrange is not None:
            arr.yrange = yrange
        print("Writing plot file:",pfile)
        arr.write_eps(pfile)
        if show:
            arr.show()


    def read_data(self):
        c=self.config
        s2vals=sample_s2(c['mins2'],c['maxs2'],c['ns2'])
        #s2vals=s2vals[1:]
        self.alldata = []
        for s2 in s2vals:
            #print('-'*72)
            #print("s2:",s2)

            pattern=get_simfile(self['run'], self['objmodel'], self['psfmodel'],
                            s2, '*', self['psf_ellip'])
            print(pattern)
            if hdfs.is_in_hdfs(pattern):
                flist = hdfs.ls(pattern,full=True)
            else:
                flist = glob.glob(pattern)
            nf=len(flist)
            if nf != 0:
                st = self.struct(nf)
                for i in xrange(nf):
                    f=flist[i]
                    #print("    Found %s" % f)
                    t=eu.io.read(f)

                    for n in st.dtype.names:
                        if n in t.dtype.names:
                            st[n][i] = t[n].mean()
                self.alldata.append(st)
            else:
                print("FOUND NO FILES")
 

    def struct(self, n):
        st=numpy.zeros(n, dtype=[('s2','f8'),
                                 ('s2noweight','f8'),
                                 ('etrue','f8'),
                                 ('econv','f8'),
                                 ('etrue_uw','f8'),
                                 ('ecorr_uw','f8'),
                                 ('ecorr','f8'),
                                 ('ecorr_rg','f8'),
                                 ('R','f8'),
                                 ('R_rg','f8')])
        return st

def run_many_s2(run, verbose=False):
    """

    Run lots of s2

    Not yet adjusting double gauss parameters, just using default

    """

    c=read_config(run)
    s2list=sample_s2(c['mins2'],c['maxs2'],c['ns2'])
    print('s2list\n    ',s2list)
    for s2 in s2list:
        print("-"*70)
        print("s2:",s2)
        rs = RegaussSimulatorRescontrol(run, s2, verbose=verbose)
        rs.run_many_ellip()

    print("Done")

def sample_s2(mins2,maxs2,ns2):
    s2vals = numpy.linspace(mins2, maxs2, ns2)
    return s2vals


class RegaussSimulatorRescontrol(dict):
    """
    Simulate objects convolved with the PSF. 

    This is distinguished from the SDSS simulator because you have precise
    control over the resolution effects in the PSF.

    Random realizations of the orientation are used based on the ntrial
    parameter in the config.  
    
    If s2n is also present, each of these realizations will get a different
    noise realization.  Note s2n is the *unweighted* signal to noise,
    the weighted s2n will be different!

    For high resolution and no noise objects, only a single trial is necessary.
    
    """
    def __init__(self, run, s2, debug=False, verbose=False):

        c = read_config(run)
        self['run']=run
        self['s2']=s2
        self['objmodel']=c['objmodel']
        self['psfmodel']=c['psfmodel']
        self['psf_ellip']=c.get('psf_ellip',0.0)
        self['psf_sigma']=c.get('psf_sigma',3.0)

        # for dgauss PSF
        self['psf_sigrat']=c.get('psf_sigrat',2.3)
        self['psf_cenrat']=c.get('psf_cenrat',0.09)
        
        if 's2n' in c:
            self['s2n'] = c['s2n']

        if 'ntrial' not in c or 'itmax' not in c:
            raise ValueError("you must send ntrial/itmax")
        self['ntrial'] = c['ntrial']
        self['itmax'] = c['itmax']

        self['mine']=c.get('mine',0.05)
        self['maxe']=c.get('maxe',0.8)
        self['nume']=c.get('nume',20)

        self['verbose'] = verbose


        # keywords for the ConvolvedImage
        convkeys = {'debug':debug, 'verbose':verbose}
        if 'sim_conv' in c:
            convkeys['conv'] = c['sim_conv']
        for k in ['minres','image_nsub','fft_nsub','eps','forcegauss']:
            if k in c:
                convkeys[k] = c[k]

        rgkeys = {'debug':debug, 'verbose':verbose}
        if 'regauss_conv' in c:
            convkeys['conv'] = c['regauss_conv']

        for k in ['image_nsub','admom_nsub']:
            if k in c:
                rgkeys[k] = c[k]

        self.convkeys=convkeys
        self.rgkeys=rgkeys
        self['debug']=debug


    def run_many_ellip(self):
        for ellip in self.ellipvals():
            self.run_ellip(ellip)
        wlog("Done many_ellip")


    def run_ellip(self, ellip):

        wlog("ellip:",ellip)
        outfile = self.outfile(ellip)
        wlog("outfile:",outfile)

        # get a n ew RandomConvolvedImage with this ellip
        wlog("getting convolved image")

        #if 's2n' in self:
        if True:
            trials_file = self.outfile(ellip,type='trials')
            #output, trials = self.do_regauss_trials(ci)
            output, trials = self.do_regauss_trials(ellip)
            eu.io.write(trials_file, trials, clobber=True)
        else:
            ci = self.new_convolved_image(ellip)
            output = self.do_regauss(ci)

        eu.io.write(outfile, output, clobber=True)

        if self['debug']:
            self.write_images(ellip, ci)


    def do_regauss(self, ci, verbose_local=True):
        """
        ci is a convolved image
        """
        rgkeys = self.rgkeys

        rgkeys['guess'] = (ci['cov_admom'][0] + ci['cov_admom'][2])/2
        rgkeys['guess_psf'] = (ci['cov_psf_admom'][0] + ci['cov_psf_admom'][2])/2
        if 's2n' in self:
            rgkeys['sigsky'] = self['sigsky']
        
        if verbose_local:
            wlog("running regauss")
        rgkeys['verbose'] = verbose_local

        rg = admom.ReGauss(ci.image, ci['cen'][0], ci['cen'][1], ci.psf, **rgkeys)
        rg.do_all()
        self.add_unweighted_truth(rg, ci.image0)
        if self['verbose']:
            wlog("uwcorrstats")
            wlog(rg['uwcorrstats'])


        if rg['rgstats'] == None or rg['rgcorrstats'] == None:
            raise RuntimeError("Failed to run regauss")
        if rg['rgstats']['whyflag'] != 0:
            raise RuntimeError("regauss failed: '%s'" % rg['rgstats']['whystr'])

        # copy out the data
        output = self.copy_output(ci, rg)
        return output

    def do_regauss_trials(self, ellip):
        """
        Run multiple trials of the noise
        """
        #ci.image_nonoise = ci.image
        ntrial = self['ntrial']
        dotper=10
        if 's2n' in self:
            print("doing",self['ntrial'],"trials at S/N:",self['s2n'])
        else:
            print("doing",self['ntrial'],"trials")
        print("one dot per",dotper,"trials")
        trials = numpy.zeros(ntrial, dtype=self.out_dtype())
        itrial = 0
        iter=0
        while (itrial < ntrial) and (iter < self['itmax']):
            if ((itrial+1) % dotper) == 0:
                stderr.write('.')
            theta = 45.0*numpy.random.random()
            if itrial == 0:
                verbose_local=True
            else:
                verbose_local=False

            ci = self.new_convolved_image(ellip, 
                                          verbose_local=verbose_local,
                                          theta=theta)
            if 's2n' in self:
                ci.image_nonoise = ci.image
                self.add_noise(ci)
            try:
                trial = self.do_regauss(ci, verbose_local=False)
                trials[itrial] = trial
                itrial+=1
            except RuntimeError:
                pass
            iter+=1
        stderr.write('\n')
        if iter >= self['itmax']:
            raise ValueError("reached itmax")

        print("total iterations:",iter)
        output = self.average_trials(trials)
        return output, trials

    def new_convolved_image(self, ellip, verbose_local=True, theta=45):
        pcov=fimage.conversions.ellip2mom(2*self['psf_sigma']**2,
                                          e=self['psf_ellip'],theta=0)
        if self['psfmodel'] == 'dgauss':
            pcov1=pcov
            pcov2=pcov*self['psf_sigrat']**2
            b=self['psf_cenrat']
            psfpars = dict(model = 'dgauss',
                           cov1 = pcov1,
                           cov2 = pcov2,
                           cenrat=b)
        else:
            psfpars = dict(model = 'gauss',
                           cov = pcov)

        sigma = self['psf_sigma']/sqrt(self['s2'])
        cov=fimage.conversions.ellip2mom(2*sigma**2,e=ellip,theta=theta)
        objpars = dict(model = self['objmodel'],
                       cov=cov)

        # default in convolved image is 16
        #ci = fimage.convolved.ConvolvedImage(objpars,psfpars, **self.convkeys)
        self.convkeys['verbose'] = verbose_local
        ci = fimage.convolved.ConvolvedImageFFT(objpars,psfpars, **self.convkeys)

        return ci
        
    def add_noise(self, ci):
        """
        Add noise based on requested S/N.  We only add background noise
        so the S/N is

              sum(pix)
        -------------------   = S/N
        sqrt(npix*skysig**2)

        thus
            
            sum(pix)
        ----------------   = sigsky
        sqrt(npix)*(S/N)

        
        We use an aperture of 3sigma but if no pixels
        are returned we increase the aperture until some
        are returned.

        Side effects:
            ci.image is set to the noisy image
            self['sigsky'] is set to the noise level
        """
        s2n = self['s2n']

        cen = ci['cen_uw']
        T = ci['cov_uw'][0] + ci['cov_uw'][2]
        sigma = fimage.mom2sigma(T)
        imagenn = ci.image_nonoise
        shape = imagenn.shape

        row,col=ogrid[0:shape[0], 0:shape[1]]
        rm = array(row - cen[0], dtype='f8')
        cm = array(col - cen[1], dtype='f8')

        radpix = sqrt(rm**2 + cm**2)
        # sigfac is 4 in admom
        sigfac = 4.0
        step=0.1
        w = where(radpix <= sigfac*sigma)
        npix = w[0].size + w[1].size
        while npix == 0:
            sigfac += step
            w = where(radpix <= sigfac*sigma)
            npix = w[0].size + w[1].size

        pix = imagenn[w]
        wt = sqrt(pix)
        wsignal = (wt*pix).sum()
        wsum = wt.sum()

        sigsky = wsignal/sqrt(wsum)/s2n

        noise_image = \
            sigsky*standard_normal(imagenn.size).reshape(shape)
        ci.image = imagenn + noise_image

        out = admom.admom(ci.image, cen[0], cen[1], guess=T/2., sigsky=sigsky)

        # fix up sigsky based on measurement
        sigsky = out['s2n']/s2n*sigsky
        noise_image = \
            sigsky*standard_normal(imagenn.size).reshape(shape)
        ci.image = imagenn + noise_image
        self['sigsky'] = sigsky


        if False:
            out = admom.admom(ci.image, cen[0], cen[1], guess=T/2., sigsky=sigsky)
            print("    target S/N:            ",s2n)
            print("    meas S/N after noise:  ",out['s2n'])

    def add_unweighted_truth(self, rg, image0):
        mom0=fimage.statistics.fmom(image0)

        uw = rg['uwcorrstats']

        uw['cen_image0'] = mom0['cen']
        uw['cov_image0'] = mom0['cov']

        T0 = mom0['cov'][0]+mom0['cov'][2]
        uw['e1_image0']=(mom0['cov'][2]-mom0['cov'][0])/T0
        uw['e2_image0']=2*mom0['cov'][1]/T0
        uw['e_image0'] = sqrt( uw['e1_image0']**2 + uw['e2_image0']**2 )

        rg['uwcorrstats'] = uw


    def copy_output(self, ci, rg):
        st = numpy.zeros(1, dtype=self.out_dtype())
        ims = rg['imstats']
        psfs = rg['psfstats']

        corrs=rg['corrstats']
        rgs = rg['rgstats']
        rgcorrs = rg['rgcorrstats']

        st['s2'] = self['s2']
        st['s2noweight'] = (ci['cov_psf_uw'][0]+ci['cov_psf_uw'][2])/(ci['cov_image0_uw'][0]+ci['cov_image0_uw'][2])
        st['s2admom'] = (ci['cov_psf_admom'][0]+ci['cov_psf_admom'][2])/(ci['cov_image0_admom'][0]+ci['cov_image0_admom'][2])

        st['e1true'] = ci['e1true']
        st['e2true'] = ci['e2true']
        st['etrue']  = ci['etrue']

        # "conv" here means "measured off the convolved image"
        st['e1conv'] = ims['e1']
        st['e2conv'] = ims['e2']
        st['econv'] = sqrt( ims['e1']**2 + ims['e2']**2 )
        st['uncerconv'] = ims['uncer']
        st['s2nconv'] = ims['s2n']

        st['R'] = corrs['R']
        st['R_rg'] = rgcorrs['R']
        st['e1corr'] = corrs['e1']
        st['e2corr'] = corrs['e2']
        st['ecorr'] = sqrt( corrs['e1']**2 + corrs['e2']**2 )

        st['e1corr_rg'] = rgcorrs['e1']
        st['e2corr_rg'] = rgcorrs['e2']
        st['ecorr_rg'] = sqrt( rgcorrs['e1']**2 + rgcorrs['e2']**2 )

        uwcorrs = rg['uwcorrstats']
        st['e1corr_uw'] = uwcorrs['e1']
        st['e2corr_uw'] = uwcorrs['e2']
        st['ecorr_uw'] = uwcorrs['e']

        if 's2n' in self:
            st['s2n'] = self['s2n']
            st['ntrial'] = self['ntrial']
            st['itmax'] = self['itmax']
        return st

    def out_dtype(self):
        dt = [('s2','f8'),
              ('s2noweight','f8'), # actual s2 of object
              ('s2admom','f8'),    # s2 from admom, generally different
              ('e1true','f8'),
              ('e2true','f8'),
              ('etrue','f8'),
              ('e1conv','f8'),
              ('e2conv','f8'),
              ('econv','f8'),
              ('uncerconv','f8'),
              ('s2nconv','f8'),
              ('R','f8'),
              ('R_rg','f8'),
              ('e1corr','f8'),
              ('e2corr','f8'),
              ('ecorr','f8'),
              ('e1corr_rg','f8'),
              ('e2corr_rg','f8'),
              ('ecorr_rg','f8'),
              ('e1corr_uw','f8'),
              ('e2corr_uw','f8'),
              ('ecorr_uw','f8')]
        if 'ntrial' in self:
            dt += [('s2n','f8'),('ntrial','i8'),('itmax','i8')]
        return dt

    def average_trials(self, trials):
        st = numpy.zeros(1, dtype=self.out_dtype())
        for n in st.dtype.names:
            st[n] = median(trials[n])
        return st
    def ellipvals(self):
        ellipvals = numpy.linspace(self['mine'], self['maxe'], self['nume'])
        return ellipvals

    def write_images(self, ellip, ci):
        f=self.outfile(ellip, 'image')
        wlog("writing image file:",f)
        eu.io.write(f, ci.image, clobber=True)

        f=self.outfile(ellip, 'image0')
        wlog("writing image0 file:",f)
        eu.io.write(f, ci.image0, clobber=True)

        f=self.outfile(ellip, 'psf')
        wlog("writing psf file:",f)
        eu.io.write(f, ci.psf, clobber=True)

    def outfile(self, ellip, type='res'):
        f = get_simfile(self['run'], 
                        self['objmodel'], self['psfmodel'], 
                        self['s2'], ellip, self['psf_ellip'],
                        type=type)

        dir=os.path.dirname(f)
        if not os.path.exists(dir):
            wlog("Making output dir:",dir)
            try:
                os.makedirs(dir)
            except:
                # probably a race condition
                pass
        return f


class RandomSDSSPSF(dict):
    """
    Read random PSFs from SDSS fields in filter and seeing range.

    If type is None, you must explicitly call next_kl() or next_dgauss()

    """
    def __init__(self, type, min_seeing, max_seeing, filter, verbose=False):
        import sdsspy
        if type not in ['kl','dgauss']:
            raise ValueError("type must be 'kl' or 'dgauss'")

        self['type'] = type
        self.verbose=verbose
        self['min_seeing'] = min_seeing
        self['max_seeing'] = max_seeing

        self['filter'] = filter
        filternum = sdsspy.FILTERNUM[filter]
        self['filternum'] = filternum

        self.cache_window()
        self.init_random()

        # immediately read some data
        self.next()

    def cache_window(self):
        import sdsspy
        w=sdsspy.window.Window()
        if self.verbose:
            wlog("Cacheing window flist")
        flist = w.read('flist')
        
        if self.verbose:
            wlog("Extracting matching fields")
        rfwhm = flist['psf_fwhm'][ :,self['filternum'] ] 
        w,=numpy.where(  (rfwhm > self['min_seeing']) 
                       & (rfwhm < self['max_seeing']) 
                       & (flist['score'] > 0.1) 
                       & (flist['rerun'] == '301') )
        if w.size == 0:
            raise ValueError("No runs found with seeing in [%0.3f,%0.3f]" % (min_seeing,max_seeing))

        if self.verbose:
            wlog("  Found:",w.size)
        self.flist = flist[w]

    def init_random(self):
        # get random index
        r = numpy.random.random(self.flist.size)
        self.randindex = r.argsort()
        self.i = 0
        self.current_i = 0

    def getmany(self, n):
        if self['type'] !=  'dgauss':
            raise ValueError("only support dgauss with getmany for now")
        res=numpy.zeros(n,dtype=[('sigma1','f8'),
                                 ('sigma2','f8'),
                                 ('b','f8')])
        for i in xrange(n):
            tmp = self.get()
            res['sigma1'][i] = tmp['sigma1']
            res['sigma2'][i] = tmp['sigma2']
            res['b'][i] = tmp['b']
            self.next()
        
        return res

    def get(self, meta=False):
        if self['type'] == 'dgauss':
            return self.get_dgauss(meta=meta)
        else:
            return self.get_kl(meta=meta)

    def get_dgauss(self, meta=False):
        m=self.meta
        fnum=self['filternum']
        out={'sigma1':m['psf_sigma1_2g'][fnum],
             'sigma2':m['psf_sigma2_2g'][fnum],
             'b':m['psf_b_2g'][fnum]}

        if meta:
            return out, self.meta
        else:
            return out

    def get_kl(self, meta=False):
        if meta:
            return self.kl, self.meta
        else:
            return self.kl

    def next(self):
        """
        Return the next random PSF model or image

        Make sure it is normalized.
        """

        if self['type'] == 'dgauss':
            self.next_dgauss()
        else:
            self.next_kl()


    def next_dgauss(self):

        i = self.randind()
        self.load_meta(i)

    def next_kl(self, meta=False):
        i = self.randind()

        self.load_meta(i)
        self.load_kl(i)

    def load_kl(self, i):
        import sdsspy
        self.kl = sdsspy.read('psField', 
                              self.flist['run'][i], 
                              self.flist['camcol'][i], 
                              self.flist['field'][i],
                              filter=self['filter'],
                              verbose=self.verbose)

    def load_meta(self, i):
        import sdsspy
        meta =  sdsspy.read('psField', 
                            self.flist['run'][i], 
                            self.flist['camcol'][i], 
                            self.flist['field'][i],
                            lower=True,
                            verbose=self.verbose)
        self.meta = meta[0]

    def randind(self):
        i = self.randindex[self.i % self.randindex.size] 
        self.current_i = i
        self.i += 1
        return i

class RegaussSDSSSimulator(dict):
    """
    Simulate images using a random SDSS PSF
    """
    def __init__(self, type, run, s2, objmodel, 
                 psf_fwhm_range=[1.3,1.5],
                 psf_filter='r',
                 mine=0.05, maxe=0.8, nume=20, 
                 nrand=100, trialfac=10,
                 conv='fast',
                 debug=False):

        if str(type) not in ['dgauss','kl']:
            raise ValueError("type must be 'kl' or 'dgauss'")
        self['type'] = type
        self['run']=run
        self['s2']=s2
        self['objmodel']=objmodel
        self['psfmodel']='sdss'

        self['mine']=mine
        self['maxe']=maxe
        self['nume']=nume
        self['nrand']=nrand
        self['trialfac'] = trialfac
        self['ntrial']=nrand*trialfac

        self['psf_fwhm_range'] = psf_fwhm_range
        self['psf_filter'] = psf_filter

        self['conv'] = conv

        self['debug'] = debug

        self.rsp = RandomSDSSPSF(type,psf_fwhm_range[0],psf_fwhm_range[1],psf_filter)

    def run_many_ellip(self):
        for ellip in self.ellipvals():
            self.run_ellip(ellip)
        wlog("Done many_ellip")
    def run_ellip(self, ellip):
        """

        Do nrand realizations for each ellipticity value

        If convergence fails, retry ntrial times 
        """

        wlog("ellip:",ellip)
        outfile = self.outfile(ellip)
        wlog("outfile:",outfile)

        dir=os.path.dirname(outfile)
        if not os.path.exists(dir):
            wlog("Making output dir:",dir)
            try:
                os.makedirs(dir)
            except:
                # probably race condition
                pass


        robj=None

        # get a n ew RandomConvolvedImage with this ellip
        rci = self.new_random_convolved_image(ellip)


        # do all the randoms, allowing for a certain number of failures to
        # retry

        randi=0
        trial = 0
        nrand = self['nrand']
        ntrial = self['ntrial']
        while randi < nrand and trial < ntrial:

            # moments of the object pre-convolution
            Tguess0 = rci.objpars['Irr_meas'] + rci.objpars['Icc_meas']

            # moments after convolution
            Tguess=rci['Irr'] + rci['Icc']

            # moments of the psf
            Tguess_psf = rci['Irr_psf'] + rci['Icc_psf']

            # get moments before convolution
            amtrue = admom.admom(rci.image0, rci['cen'][0], rci['cen'][1],
                                 Tguess=Tguess0)
            if amtrue['whyflag'] == 0:
                # do regauss here
                rg = admom.ReGauss(rci.image, rci['cen'][0], rci['cen'][1],
                                   rci.psf, Tguess=Tguess,Tguess_psf=Tguess_psf)
                rg.do_all()

                if rg['rgstats'] != None and rg['rgcorrstats'] != None:
                    if rg['rgstats']['whyflag'] == 0:
                        # copy out the data
                        output = self.copy_output(amtrue, rg, rci['theta'])
                        if robj is None:
                            robj = eu.sfile.Open(outfile, 'w')
                        robj.write(output)

                        # only now do we increment randi
                        randi += 1

            trial += 1
            if randi < nrand and trial < ntrial:
                rci = self.new_random_convolved_image(ellip)

        if robj is not None:
            robj.close()

        wlog("ntrial:",trial," nfail:",trial-nrand)
        wlog("randi/nrand: %s/%s" % (randi,nrand))
        if randi != nrand:
            wlog("Exceeded max trials, failed to get all",nrand," realizations")

        #sys.stdout.flush()
        
    def copy_output(self, amtrue, rg, theta):
        st = numpy.zeros(1, dtype=self.out_dtype())
        ims = rg['imstats']
        psfs = rg['psfstats']

        corrs=rg['corrstats']
        rgs = rg['rgstats']
        rgcorrs = rg['rgcorrstats']

        st['theta'] = theta
        st['s2'] = self['s2']
        st['s2prepsf'] = (psfs['Irr']+psfs['Icc'])/(amtrue['Irr'] + amtrue['Icc'])
        st['s2postpsf'] = (psfs['Irr']+psfs['Icc'])/(ims['Irr']+ims['Icc'])

        st['e1true'] = amtrue['e1']
        st['e2true'] = amtrue['e2']
        st['etrue'] = sqrt( amtrue['e1']**2 + amtrue['e2']**2 )

        st['e1conv'] = ims['e1']
        st['e2conv'] = ims['e2']
        st['econv'] = sqrt( ims['e1']**2 + ims['e2']**2 )

        st['R'] = corrs['R']
        st['e1corr'] = corrs['e1']
        st['e2corr'] = corrs['e2']
        st['ecorr'] = sqrt( corrs['e1']**2 + corrs['e2']**2 )

        st['e1corr_rg'] = rgcorrs['e1']
        st['e2corr_rg'] = rgcorrs['e2']
        st['ecorr_rg'] = sqrt( rgcorrs['e1']**2 + rgcorrs['e2']**2 )

        return st

    def out_dtype(self):
        dt = [('theta','f8'),
              ('s2','f8'),
              ('s2prepsf','f8'),
              ('s2postpsf','f8'),
              ('e1true','f8'),
              ('e2true','f8'),
              ('etrue','f8'),
              ('e1conv','f8'),
              ('e2conv','f8'),
              ('econv','f8'),
              ('R','f8'),
              ('e1corr','f8'),
              ('e2corr','f8'),
              ('ecorr','f8'),
              ('e1corr_rg','f8'),
              ('e2corr_rg','f8'),
              ('ecorr_rg','f8')]
        return dt

    def ellipvals(self):
        ellipvals = numpy.linspace(self['mine'], self['maxe'], self['nume'])
        return ellipvals

    def outfile(self, ellip):
        return get_simfile(self['run'], 
                           self['objmodel'], 
                           self['psfmodel'], 
                           self['s2'], ellip)

    def new_random_convolved_image(self, ellip):
        p = self.rsp.get()
        psfpars={'model':'dgauss',
                 'Irr1':p['sigma1']**2,
                 'Irc1':0.0,
                 'Icc1':p['sigma1']**2,
                 'sizerat': p['sigma2']/p['sigma1'],
                 'cenrat': p['b'],
                 'dims': [31,31],
                 'cen':[15,15]}

        # parameters for a model galaxy of the right size

        # just relative to the first, dominant one
        Tobj = (psfpars['Irr1']+psfpars['Icc1'])/self['s2']

        # random angle
        theta = numpy.random.random()*360
        Irr,Irc,Icc = fimage.ellip2mom(Tobj, e=ellip, theta=theta)

        objpars={'model':self['objmodel'],
                 'Irr':Irr,
                 'Irc':Irc,
                 'Icc':Icc}

        ci = fimage.ConvolvedImage(objpars, psfpars, conv=self['conv'])
        ci['theta'] = theta

        if self['debug']:
            self.rsp.load_kl(self.rsp.current_i)
            psf = self.rsp.kl.rec(1000,1000, trim=True)
            images.compare_images(ci.psf, psf)
            ci.show()
            key=raw_input('hit a key: ')

        self.rsp.next()

        return ci


def configdir():
    d=os.environ['ESPY_DIR']
    return path_join(d,'lensing','config')
def configfile(run):
    d=configdir()
    name='rgsim-%s.json' % run
    return path_join(d, name)
def read_config(run):
    f=configfile(run)
    return eu.io.read(f)

def get_wq_dir(run):
    dir=os.environ.get('REGAUSSIM_DIR')
    dir=path_join(dir, run, 'wq')
    return dir
def get_wq_url(run,objmodel,psfmodel,extra):
    dir=get_wq_dir(run)
    f='rgsim-%s-%s-%s-%s.yaml' % (run,objmodel,psfmodel,extra)
    f=os.path.join(dir,f)
    return f


def get_simfile(run, objmodel, psfmodel, s2, ellip, psf_ellip=None, type='res'):
    dir=get_simdir(run)

    if ellip != '*':
        ellip = '%0.3f' % ellip
    
    f = 'rgsim-%sobj-%s-%spsf' % (objmodel, ellip, psfmodel)
    if psf_ellip is not None:
        f += '-%0.3f' % psf_ellip
    f += '-s2-%0.3f' % s2

    if type == 'res':
        f += '.rec'
    elif type == 'trials':
        f += '-%s.rec' % type
    elif type in ['image','image0','psf']:
        f += '-%s.fits' % type
    else:
        raise ValueError("type must be 'rec','image','image0','psf'")

    #f = 'rgsim-%sobj-%s-%spsf-%0.3f-s2-%0.3f.rec' % (objmodel, ellip, psfmodel, psf_ellip, s2)
    f = path_join(dir, f)
    return f


def get_plotdir(run):
    dir=os.environ.get('REGAUSSIM_DIR')
    dir=path_join(dir, run, 'plots')
    return dir

def get_simdir(run):
    #dir=os.environ.get('REGAUSSIM_HDFS_DIR')
    dir=os.environ.get('REGAUSSIM_DIR')
    dir=path_join(dir, run, 'outputs')
    return dir

class RandomConvolvedImage(dict):
    """

    Generate a model convolved with a PSF using the imsim.ConvolvedImage class.

    The input 'psfmodel' parameter can be a string like 'gauss','dgauss' 
    or a generator with the method .next()
    
    random orientations for objects and PSF are chosen.

    The elliticity of the object pre-convolution is entered and if the psf
    is a model the psf_ellip can also be entered, with default 0.

    The optional psf_fwhm parameter is used for model PSFs as a way to scale
    all images.  The value is FWHM in *pixels*, so using a large number >10 would
    reduce pixelization sampling effects; seeing of 4'' in SDSS is that big.


    The galaxy will be scaled to the *measured* psf size by 1/s2.  s2 is the
    ratio (sizepsf/sizeobj)**2  This is only approximate as it depends somewhat
    on the model.

    For double gaussians, you need to send the Tratio and fluxfrac1 keywords,
    which determines the size of the second gaussian relative to the first

        T2 = T1*Tratio

    And how much of the flux is in the first gaussian:

        flux = fluxfrac1*gauss1 + (1-fluxfrac1)*guass2

    """
    def __init__(self, s2, objmodel, obj_ellip, psfmodel,
                 psf_ellip=0.0, Tratio_psf=None, fluxfrac1=None, 
                 psf_fwhm=10.0,  # in pixels.  4'' seeing in SDSS is this big
                 verbose=False):

        self.verbose=verbose
        self['s2'] = s2
        self['objmodel'] = objmodel
        self['obj_ellip'] = obj_ellip
        self['psfmodel'] = psfmodel
        
        if isinstance(psfmodel, (str,unicode)):
            self['psf_ellip'] = psf_ellip
            self['Tratio_psf'] = Tratio_psf
            self['fluxfrac1'] = fluxfrac1
            self['psf_fwhm'] = psf_fwhm
            self['Tpsf'] = admom.fwhm2mom(psf_fwhm)
        else:
            self['psf_ellip'] = None
            self['Tratio_psf'] = None
            self['fluxfrac1'] = None
            self['psf_fwhm'] = None
            self['Tpsf'] = None

        self['Irr_psf'] = None
        self['Irc_psf'] = None
        self['Icc_psf'] = None

        self.randomize()

    def randomize(self):
        """
        You can run generate multiple times to get more
        realizations
        """

        self.make_psf_stats()
        self.make_object_stats()
        self.make_images()


    def make_psf_stats(self):
        """
        Either generate moments from inputs or read an SDSS PSF
        """
        psfmodel = self['psfmodel']
        if isinstance(psfmodel, (str,unicode)):
            # generate a random orientation
            theta = 360.0*numpy.random.random()
            Irr, Irc, Icc = admom.ellip2mom(self['Tpsf'], e=self['psf_ellip'], theta=theta)
            self['Irr_psf'] = Irr
            self['Irc_psf'] = Irc
            self['Icc_psf'] = Icc
            self.psf = None
        else:
            # this is a psf generator.  We assume the psf has the center at
            # the image middle, is normalized
            psf = psfmodel.next()
            cen = [(psf.shape[0]-1)/2., (psf.shape[1]-1)/2.]
            out = admom.admom(psf, cen[0], cen[1])

            if out['whyflag'] != 0:
                raise RuntimeError("failure measuring psf admom")
            self.psf=psf
            self['Irr_psf'] = out['Irr']
            self['Irc_psf'] = out['Irc']
            self['Icc_psf'] = out['Icc']
            self['Tpsf'] = out['Irr'] + out['Icc']
            
    def make_object_stats(self):
        """
        Generate moments for object from the inputs before convolution
        """

        Tobj = self['Tpsf']/self['s2']

        theta = 360.0*numpy.random.random()
        Irr, Irc, Icc = admom.ellip2mom(Tobj, e=self['obj_ellip'], theta=theta)
        self['Irr'] = Irr
        self['Irc'] = Irc
        self['Icc'] = Icc
 
    def make_images(self):
        if self.psf is None:
            # just a string descriptor
            self.cm = imsim.ConvolvedImage(self['objmodel'], 
                                           self['Irr'], self['Irc'], self['Icc'],
                                           self['psfmodel'],
                                           Irr_psf=self['Irr_psf'],
                                           Irc_psf=self['Irc_psf'],
                                           Icc_psf=self['Icc_psf'],
                                           Tratio=self['Tratio_psf'],
                                           fluxfrac1=self['fluxfrac1'],
                                           verbose=self.verbose)
        else:
            # an image
            self.cm = imsim.ConvolvedImage(self['objmodel'], 
                                           self['Irr'], self['Irc'], self['Icc'],
                                           self.psf,
                                           Irr_psf=self['Irr_psf'],
                                           Irc_psf=self['Irc_psf'],
                                           Icc_psf=self['Icc_psf'],
                                           verbose=self.verbose)
        self.image0 = self.cm.image0
        self.image = self.cm.image
        self['cen'] = self.cm['cen']
        self.psf = self.cm.psf
        self['psfcen'] = self.cm['psfcen']

    def show(self):
        self.cm.show()

_wqtemplate="""
command: |
    source ~/.bashrc
    python $ESPY_DIR/lensing/bin/regauss-sim.py %(opts)s %(run)s

group: gen45
job_name: %(job_name)s
priority: med\n"""

def create_sim_wq_bys2(run):
    import pbs

    c = read_config(run)
    objmodel = c['objmodel']
    psfmodel = c['psfmodel']


    s2vals=sample_s2(c['mins2'],c['maxs2'],c['ns2'])
    for s2 in s2vals:
        extra = '%0.3f' % s2
        job_name='%s-%s' % (run,extra)

        wqurl = get_wq_url(run,objmodel,psfmodel,extra)
        opts = '-s2 %(s2)0.3f' % {'s2':s2}

        eu.ostools.makedirs_fromfile(wqurl,verbose=True)
        wlog("writing wq script:",wqurl)
        with open(wqurl,'w') as fobj:
            wqscript=_wqtemplate % {'run':run, 'job_name':job_name,'opts':opts}
            fobj.write(wqscript)

def create_sim_wq_byellip(run):
    import pbs

    c = read_config(run)
    objmodel = c['objmodel']
    psfmodel = c['psfmodel']

    s2vals=sample_s2(c['mins2'],c['maxs2'],c['ns2'])
    evals=numpy.linspace(c['mine'],c['maxe'],c['nume'])
    for s2 in s2vals:
        for ellip in evals:
            extra = '%0.3f-%0.3f' % (s2,ellip)
            job_name='%s-%s' % (run,extra)

            wqurl = get_wq_url(run,objmodel,psfmodel,extra)
            opts = '--s2 %(s2)0.3f -e %(ellip)0.3f' % {'s2':s2,'ellip':ellip}

            eu.ostools.makedirs_fromfile(wqurl,verbose=True)
            wlog("writing wq script:",wqurl)
            with open(wqurl,'w') as fobj:
                wqscript=_wqtemplate % {'run':run, 'job_name':job_name,'opts':opts}
                fobj.write(wqscript)



def create_sim_pbs(run):
    import pbs

    c = read_config(run)
    objmodel = c['objmodel']
    psfmodel = c['psfmodel']

    dir=pbsdir(run)

    name='%s-%s' % (objmodel,psfmodel)
    setups="""
setup esutil -r ~/exports/esutil-work
setup admom -r ~/exports/admom-work
setup fimage -r ~/exports/fimage-work
setup sdsspy -r ~/exports/sdsspy-work
setup espy -r ~/exports/espy-work
    """

    s2vals=sample_s2(c['mins2'],c['maxs2'],c['ns2'])
    plist=[]

    if psfmodel == 'sdss':
        commands="""
import lensing
rs = lensing.regauss_sim.RegaussSDSSSimulator('dgauss','{run}', {s2}, '{objmodel}',conv='{conv}',nrand={nrand},trialfac={trialfac})
rs.run_many_ellip()
        """
    else:
        commands="""
import lensing
rs = lensing.regauss_sim.RegaussSimulatorRescontrol('{run}', {s2})
rs.run_many_ellip()
        """

    for s2 in s2vals:
        p={'run':run, 's2':s2}
        plist.append(p)
    
    queue='fast'
    pbs.create_many(dir, name, commands, plist, python=True,
                    queue=queue, setups=setups)


def psfield_compare_model(rsp=None, generator='filter', next=False):
    """
    compare psf reconstructions with the best fit models in the 6th header.
    """

    import images
    import biggles
    from scipy.ndimage.filters import gaussian_filter
    import fimage

    filter='r'
    fnum=2
    if rsp is None:
        rsp=RandomSDSSPSF(1.3,1.5,filter,verbose=True)
        image,meta = rsp.next(meta=True)
    else:
        if rsp.psf is None:
            image,meta = rsp.current(meta=True)
        else:
            if next:
                image,meta = rsp.next(meta=True)
            else:
                image,meta = rsp.current(meta=True)

    cen = [(image.shape[0]-1)/2, (image.shape[1]-1)/2]

    extra='_2G'
    a = 1.0
    b = meta['psf_b'+extra][0,fnum]
    s1 = meta['psf_sigma1'+extra][0,fnum]
    s2 = meta['psf_sigma2'+extra][0,fnum]

    if generator == 'fimage':
        fake1 = fimage.makeimage('gauss',image.shape,cen,s1**2,0,s1**2,counts=1)
        fake2 = fimage.makeimage('gauss',image.shape,cen,s2**2,0,s2**2,counts=1)

        a /= (s1**2 + b*s2**2)
        fake = a*( s1**2*fake1 + b*s2**2*fake2 )
    elif generator == 'imsim':
        import imsim
        fake1 = imsim.mom2disk('gauss',s1**2,0,s1**2,image.shape,cen=cen,counts=1)
        fake2 = imsim.mom2disk('gauss',s2**2,0,s2**2,image.shape,cen=cen,counts=1)

        a /= (s1**2 + b*s2**2)
        fake = a*( s1**2*fake1 + b*s2**2*fake2 )

    elif generator == 'filter':
        a /= (s1**2 + b*s2**2)
        fake1 = numpy.zeros_like(image)
        # convolve delta function with gaussians
        fake1[cen[0],cen[1]] = 1
        fake = a * (s1**2 * gaussian_filter(fake1, (s1,s1)) + b*s2**2*gaussian_filter(fake1, (s2,s2)))
    else:
        raise ValueError("unknown generator type: '%s'" % generator)

    wlog("image counts:",image.sum(),"model counts:",fake.sum())

    resid = fake-image

    wlog("summed residuals:",resid.sum())

    maxval = max( image.max(), fake.max() )
    minval = 0.0

    levels=7
    tab=biggles.Table(2,3)
    #tab=biggles.Table(3,2)
    implt=images.view(image, levels=levels, show=False, min=minval, max=maxval)
    fakeplt=images.view(fake, levels=levels, show=False, min=minval, max=maxval)
    residplt=images.view(resid, show=False, min=minval, max=maxval)

    #sigma = numpy.sqrt((res['Irr']+res['Icc'])/2.0)
    #lab = biggles.PlotLabel(0.1,0.9,r'$\sigma$: %0.3f' % sigma, fontsize=4, halign='left')
    #fakeplt.add(lab)

    implt.title='original'
    fakeplt.title='gaussian '+generator
    residplt.title='residuals'


    # cross-sections
    imrows = image[:,cen[1]]
    imcols = image[cen[0],:]
    fakerows = fake[:,cen[1]]
    fakecols = fake[cen[0],:]
    resrows = resid[:,cen[1]]
    rescols = resid[cen[0],:]

    himrows = biggles.Histogram(imrows, color='blue')
    himcols = biggles.Histogram(imcols, color='blue')
    hfakerows = biggles.Histogram(fakerows, color='orange')
    hfakecols = biggles.Histogram(fakecols, color='orange')
    hresrows = biggles.Histogram(resrows, color='red')
    hrescols = biggles.Histogram(rescols, color='red')

    himrows.label = 'image'
    hfakerows.label = 'model'
    hresrows.label = 'resid'
    key = biggles.PlotKey(0.1,0.9,[himrows,hfakerows,hresrows]) 
    rplt=biggles.FramedPlot()
    rplt.add( himrows, hfakerows, hresrows,key )

    cplt=biggles.FramedPlot()
    cplt.add( himcols, hfakecols, hrescols )

    rplt.aspect_ratio=1
    cplt.aspect_ratio=1


    tab[0,0] = implt
    tab[0,1] = fakeplt
    tab[0,2] = residplt
    tab[1,0] = rplt
    tab[1,1] = cplt

    #tab[0,0] = implt
    #tab[0,1] = fakeplt
    #tab[1,0] = residplt
    #tab[1,1] = rplt
    #tab[2,0] = cplt


    tab.show()


    return rsp
