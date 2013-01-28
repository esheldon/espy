"""
    %prog [options]
"""

import sys
import os
import numpy
import lensing
import pprint
import numpy
from numpy import zeros, sqrt, where, array, sqrt

from cluster_step import sh1exp, sh2exp

import cluster_step
from cluster_step import files
from cluster_step import sh1exp, sh2exp

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')
parser.add_option('-s','--shnum',default=None,
                  help='The shear number, required')

parser.add_option('-p','--psfnum',default=None,
                  help='The psf number, required')
parser.add_option('-y','--yrange',default=None,
                  help='yrange for plot')


class Plotter(object):
    def __init__(self):
        options,args = parser.parse_args(sys.argv[1:])

        self.options=options

        if options.yrange is not None:
            self.yrange=[float(y) for y in options.yrange.split(',')]

        if (options.run is None 
                or options.psfnum is None
                or options.shnum is None):
            parser.print_help()
            sys.exit(1)

        self.run=options.run
        self.psfnum=int(options.psfnum)
        self.shnum=int(options.shnum)

        self.nsub=1

        self._read_data()

    def doplot(self):
        import biggles
        #self._measure_psf_uw()
        #self._measure_gals_uw()
        self._measure_psf_admom()
        self._measure_gals_admom()
        self._measure_gals_gmix()
        return
        
        g1=self._shear['e1']*0.5
        g1err=self._shear['e1err']*0.5
        g2=self._shear['e2']*0.5
        g2err=self._shear['e2err']*0.5
        s2n=self._shear['s2n']

        color1='red'
        color2='blue'
        e1_pts=biggles.Points(s2n,g1, type='filled circle',
                              color=color1)


        e1err_pts=biggles.SymmetricErrorBarsY(s2n,g1,g1err,
                                              color=color1)

        e2_pts=biggles.Points(s2n, g2, type='filled circle',
                              color=color2)

        e2err_pts=biggles.SymmetricErrorBarsY(s2n,g2,g2err,
                                              color=color2)

        plt=biggles.FramedPlot()

        plt.add(e1_pts,e1err_pts,e2_pts,e2err_pts)

        if self.shnum in sh1exp:
            g1exp=zeros(g1.size)+sh1exp[self.shnum]
            g2exp=zeros(g1.size)+sh2exp[self.shnum]
            g1exp_plt=biggles.Curve(s2n, g1exp)
            g2exp_plt=biggles.Curve(s2n, g2exp)
            plt.add(g1exp_plt)
            plt.add(g2exp_plt)


        plt.xlog=True
        plt.xrange=[0.8*s2n.min(), 1.2*s2n.max()]
        if self.yrange is not None:
            plt.yrange=self.yrange
        plt.title=self.run
        plt.show()



    def _measure_gals_admom(self):
        import admom
        nbin=self.im_stacks.size

        Tpsf=self.psf_ares['Irr'] + self.psf_ares['Icc']
        e1p=self.psf_ares['e1']
        e2p=self.psf_ares['e2']
        a4p=self.psf_ares['a4']

        st=zeros(nbin, dtype=[('s2n','f8'),
                              ('e1','f8'),
                              ('e1err','f8'),
                              ('e2','f8'),
                              ('e2err','f8'),
                              ('T','f8'),
                              ('R','f8')])
        self._ares_dicts=[]
        for i in xrange(nbin):
            print '-'*70
            s2n_min=self.im_stacks['s2n_min'][i]
            s2n_max=self.im_stacks['s2n_max'][i]

            print 's2n: [%.2f,%.2f]' % (s2n_min,s2n_max)


            im=self.im_stacks['images'][i,:,:]
            skysig=sqrt(self.im_stacks['skyvar'][i])

            cen=array(im.shape)/2.
            ares = admom.admom(im, cen[0], cen[1],
                               sigsky=skysig,
                               guess=4.0,
                               nsub=self.nsub)
            T=ares['Irr']+ares['Icc']

            corr=admom.correct(T, ares['e1'], ares['e2'],ares['a4'],
                               Tpsf, e1p, e2p, a4p)
            
            e1,e2,R,flags=corr
            uncer=ares['uncer']/R

            err2=( 0.32**2 + uncer**2)/self.im_stacks['nstack'][i]
            err=sqrt(err2)

            mess='e1: %s +/- %s e2: %s +/- %s R: %s  flags: %s'
            mess=mess % (e1,err,e2,err,R,flags)
            print mess

            st['s2n'][i] = (s2n_min+s2n_max)/2.
            st['e1'][i] = e1
            st['e1err'][i] = err
            st['e2'][i] = e2
            st['e2err'][i] = err
            st['R'][i] = R
            st['T'][i] = T

            self._ares_dicts.append(ares)

        self._admom_shear = st


    def _measure_psf_admom(self):
        import admom

        cen=array(self.psf_stack.shape)/2.
        ares = admom.admom(self.psf_stack, cen[0], cen[1],
                           sigsky=sqrt(self.psf_skyvar),
                           guess=4.0,
                           nsub=self.nsub)
        self.psf_ares=ares
        pprint.pprint(self.psf_ares)


    def _measure_gals_gmix(self):
        import gmix_image
        from gmix_image.util import print_pars
        from esutil.random import srandu
        import esutil as eu

        print '\n\n\n'
        nbin=self.im_stacks.size

        """
        psf_pars=numpy.array( [1.0, 
                               self.psf_ares['wrow'], 
                               self.psf_ares['wcol'],
                               self.psf_ares['Irr'], 
                               self.psf_ares['Irc'], 
                               self.psf_ares['Icc']] )

        gmix_psf=gmix_image.gmix.GMixCoellip(psf_pars)
        """


        psf_pixerr=sqrt(self.psf_skyvar)
        psfim =eu.numpy_util.to_native( self.psf_stack )
        psf_counts=psfim.sum()
        psf_Tguess=self.psf_ares['Irr'] + self.psf_ares['Icc']
        while True:
            e1guess=0.0
            e2guess=0.0
            if False:
                psf_prior=numpy.array( [psfim.shape[0]/2.*(1.+0.1*srandu()), 
                                        psfim.shape[1]/2.*(1.+0.1*srandu()), 
                                        e1guess+0.05*srandu(),
                                        e2guess+0.05*srandu(),
                                        psf_Tguess*0.6*(1.+0.1*srandu()),
                                        psf_Tguess*0.4*(1.+0.1*srandu()),
                                        psf_Tguess*0.2*(1.+0.1*srandu()),
                                        psf_counts*0.3*(1.+0.1*srandu()), 
                                        psf_counts*0.4*(1.+0.1*srandu()), 
                                        psf_counts*0.3*(1.+0.1*srandu())])
            elif True:
                psf_prior=numpy.array( [psfim.shape[0]/2.*(1.+0.1*srandu()), 
                                        psfim.shape[1]/2.*(1.+0.1*srandu()), 
                                        e1guess+0.05*srandu(),
                                        e2guess+0.05*srandu(),
                                        psf_Tguess*(1.+0.1*srandu()),
                                        psf_counts*(1.+0.1*srandu())] )

            psf_width=numpy.abs(psf_prior)*1.e6
            fitter=gmix_image.gmix_fit.GMixFitCoellip(psfim, psf_pixerr, psf_prior, psf_width, 
                                                      verbose=False)
            flags=fitter.get_flags()
            if flags==0:
                break
        
        psf_pars=fitter.get_pars()
        psf_perr=fitter.get_perr()
        print_pars(psf_pars, front='psf pars: ')
        print_pars(psf_perr, front='psf perr: ')
        gmix_psf = fitter.get_gmix()
        print gmix_psf

        for i in xrange(nbin):
            print '-'*70
            print self.im_stacks['s2n_min'][i], self.im_stacks['s2n_max'][i]


            im=eu.numpy_util.to_native( self.im_stacks['images'][i,:,:] )
            pixerr=numpy.sqrt(self.im_stacks['skyvar'][i])
            ivar=1./self.im_stacks['skyvar'][i]
            
            counts=im.sum() 
            #im /= counts
            #im += 0.0001*numpy.random.randn(im.size).reshape(im.shape)
            #pixerr=0.0001
            #pixerr /= counts
            #counts=1.
            #import images
            #images.multiview(im)
            #stop


            #print 'counts: %s pixerr: %s' % (counts,pixerr)

            #pixerr=sqrt(counts)

            #e1guess=self._admom_shear['e1'][i]
            #e2guess=self._admom_shear['e2'][i]
            e1guess=0.0
            e2guess=0.0
            #e1guess=0.2
            #e2guess=0.2
            Tguess = self._admom_shear['T'][i]

            row_guess=self._ares_dicts[i]['wrow']
            col_guess=self._ares_dicts[i]['wcol']

            while True:
                if False:
                    prior=numpy.array( [row_guess,
                                        col_guess,
                                        e1guess+0.05*srandu(),
                                        e2guess+0.05*srandu(),
                                        Tguess*0.6*(1.+0.1*srandu()),
                                        Tguess*0.4*(1.+0.1*srandu()),
                                        Tguess*0.2*(1.+0.1*srandu()),
                                        counts*0.3*(1.+0.1*srandu()), 
                                        counts*0.4*(1.+0.1*srandu()), 
                                        counts*0.3*(1.+0.1*srandu())])
                elif True:
                    prior=numpy.array( [row_guess,
                                        col_guess,
                                        e1guess+0.05*srandu(),
                                        e2guess+0.05*srandu(),
                                        Tguess*0.5*(1.+0.1*srandu()),
                                        Tguess*0.5*(1.+0.1*srandu()),
                                        counts*0.5*(1.+0.1*srandu()), 
                                        counts*0.5*(1.+0.1*srandu())])
                elif False:
                    prior=numpy.array( [row_guess,
                                        col_guess,
                                        e1guess+0.05*srandu(),
                                        e2guess+0.05*srandu(),
                                        Tguess*(1.+0.05*srandu()), 
                                        counts*(1.0+0.05*srandu())] )


                width=numpy.abs(prior)*1.e6
                #width[0] = 0.01
                #width[1] = 0.01

                #fitter=gmix_image.gmix_fit.GMixFitCoellip(im, pixerr, prior, width, 
                #                                          psf=gmix_psf, verbose=False)
                #flags=fitter.get_flags()

                fitter=gmix_image.gmix_fit.GMixFitSimple(im, 1./pixerr**2, gmix_psf, 'gexp',
                                                         self._ares_dicts[i])
                res=fitter.get_result()
                flags=res['flags']


                if flags==0:
                    break
                else:
                    print flags

            res=fitter.get_result()
            print_pars(prior,front='prior: ')
            print_pars(res['pars'], front='pars: ')
            print_pars(res['perr'], front='perr: ')
            #print 'chi2per:',res['chi2per']




    def _measure_psf_uw(self):
        import fimage
        
        res=fimage.fmom(self.psf_stack)
        pprint.pprint(res)

        self.psf_res=res

    def _measure_gals_uw(self):
        import fimage
        nbin=self.im_stacks.size

        for i in xrange(nbin):
            print '-'*70
            print self.im_stacks['s2n_min'][i], self.im_stacks['s2n_max'][i]


            im=self.im_stacks['images'][i,:,:]

            res=fimage.fmom(im)

            irr=res['cov'][0]-self.psf_res['cov'][0]
            irc=res['cov'][1]-self.psf_res['cov'][1]
            icc=res['cov'][2]-self.psf_res['cov'][2]

            T=irr+icc

            e1=(icc-irr)/T
            e2=2*irc/T
            
            print e1,e2


    def _read_data(self):
        import fitsio
        path=files.get_output_path(ftype='shear-stack',
                                   run=self.run,
                                   psfnum=self.psfnum,
                                   shnum=self.shnum)
        print path
        with fitsio.FITS(path) as fits:
            #res=fits[1].read()
            self.psf_stack=fits[2].read()

            self.psf_h=fits[2].read_header()
            self.psf_skyvar=self.psf_h['skyvar']

            self.im_stacks=fits[3].read()

def main():
    c=Plotter()
    c.doplot()

main()
