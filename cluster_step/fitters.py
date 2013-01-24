from numpy import zeros, sqrt, linspace, array, arange, where
import cluster_step
from cluster_step import files, stats


from esutil.stat import wmom

from fitting import LineFitter
from mcmc import PolyFitter

class BiasFitter(object):
    def __init__(self, data, order=1, s2n_field='s2n_w', fitter='mcmc'):
        #from esutil.numpy_util import aprint

        self.data=data
        self.order=order
        self.s2n_field=s2n_field
        self.fitter=fitter

        self.organize_data()
        self.set_averages()

        #aprint(self.avg, fancy=True)

        self.g1ind=arange(8)
        self.g2ind=arange(8)

        if self.fitter not in ['lm','mcmc']:
            raise ValueError("bad fitter type: '%s'" % self.fitter)
        self.do_fits()

    def organize_data(self):

        self.datadict={}
        for shnum in files.SHNUMS:
            w,=where(self.data['shnum'] == shnum)

            self.datadict[shnum]=self.data[w]
     
    def set_averages(self):
        nsh=len(files.SHNUMS)
        dt=[('s2n','f8'),
            ('g1true','f8'),
            ('g2true','f8'),
            ('g1','f8'),('g1_err','f8'),
            ('g2','f8'),('g2_err','f8'),
            ('g1sens','f8'),('g1sens_err','f8'),
            ('g2sens','f8'),('g2sens_err','f8')]
        out=zeros(nsh, dtype=dt)
        for i,shnum in enumerate(files.SHNUMS):
            data=self.datadict[shnum]
            shres=stats.get_mean_shear_shstruct(data)

            wts=stats.get_weights(data['gcov'])
            out['s2n'][shnum-1],err = wmom(data[self.s2n_field], wts)

            out['g1true'][shnum-1] = cluster_step.sh1exp[shnum]
            out['g2true'][shnum-1] = cluster_step.sh2exp[shnum]

            for n in shres:
                out[n][shnum-1] = shres[n]

        self.avg=out

    def get_g1poly(self):
        return self.g1fit.get_poly()
    def get_g2poly(self):
        return self.g2fit.get_poly()

    def do_fits(self):

        if self.order > 1 and self.fitter=='lm':
            raise ValueError("implement poly lm fitter")

        g1true=self.avg['g1true'][self.g1ind]
        g1vals=self.avg['g1'][self.g1ind]
        g1err=self.avg['g1_err'][self.g1ind]

        g2true=self.avg['g2true'][self.g2ind]
        g2vals=self.avg['g2'][self.g2ind]
        g2err=self.avg['g2_err'][self.g2ind]

        g1diff=g1vals-g1true
        g2diff=g2vals-g2true

        if self.fitter=='lm':
            self.g1fit=LineFitter(g1true, g1diff, g1err)
            self.g2fit=LineFitter(g2true, g2diff, g2err)
        else:
            nwalkers=200
            burnin=100
            nstep=100

            self.g1fit=PolyFitter(self.order, g1true,g1diff,nwalkers,burnin,nstep,
                                  yerr=g1err, guess=[0.01,0.001])
            self.g2fit=PolyFitter(self.order, g2true,g2diff,nwalkers,burnin,nstep,
                                  yerr=g2err, guess=[0.01,0.001])


        self.g1true=g1true
        self.g2true=g2true

        self.g1=g1vals
        self.g2=g2vals
        self.g1err=g1err
        self.g2err=g2err

        self.g1diff=g1diff
        self.g2diff=g2diff

