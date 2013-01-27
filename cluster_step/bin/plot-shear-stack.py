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

class Plotter(object):
    def __init__(self):
        options,args = parser.parse_args(sys.argv[1:])

        self.options=options

        if (options.run is None 
                or options.psfnum is None
                or options.shnum is None):
            parser.print_help()
            sys.exit(1)

        self.run=options.run
        self.psfnum=int(options.psfnum)
        self.shnum=int(options.shnum)

        self.nsub=16

        self._read_data()

    def doplot(self):
        import biggles
        #self._measure_psf_uw()
        #self._measure_gals_uw()
        self._measure_psf_admom()
        self._measure_gals_admom()
        
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
                              ('R','f8')])
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

        self._shear = st


    def _measure_psf_admom(self):
        import admom

        cen=array(self.psf_stack.shape)/2.
        ares = admom.admom(self.psf_stack, cen[0], cen[1],
                           sigsky=sqrt(self.psf_skyvar),
                           guess=4.0,
                           nsub=self.nsub)
        self.psf_ares=ares
        pprint.pprint(self.psf_ares)

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
