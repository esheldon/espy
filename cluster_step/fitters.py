from numpy import zeros, sqrt, linspace, array, arange, where
import cluster_step
from cluster_step import files, stats


from esutil.numpy_util import aprint
from esutil.stat import wmom

from fitting import LineFitter

class BiasFitter(object):
    def __init__(self, data, run, **keys):

        self.keys=keys
        self.data=data
        self.run     = run

        self.s2n_field=keys.get('s2n_field','s2n_w')

        # starting point for labels
        self.lab1_loc=[1.-0.075, 0.1]
        self.lab1_halign='right'
        self.lab1_yshift=0.075

        self.lab2_loc=[0.075,0.075]
        self.lab2_halign='left'
        self.lab2_yshift=+0.075

        self.lab3_loc=[1-0.075,1-0.075]
        self.lab3_halign='right'
        self.lab3_yshift=-0.075

        self.organize_data()
        self.set_averages()

        aprint(self.avg, fancy=True)

        self.g1ind=arange(8)
        self.g2ind=arange(8)
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

    def do_fits(self):

        g1true=self.avg['g1true'][self.g1ind]
        g1vals=self.avg['g1'][self.g1ind]
        g1err=self.avg['g1_err'][self.g1ind]

        g2true=self.avg['g2true'][self.g2ind]
        g2vals=self.avg['g2'][self.g2ind]
        g2err=self.avg['g2_err'][self.g2ind]

        g1diff=g1vals-g1true
        g2diff=g2vals-g2true

        self.g1fit=LineFitter(g1true, g1diff, g1err)
        self.g2fit=LineFitter(g2true, g2diff, g2err)


        print self.g1fit
        print self.g2fit



