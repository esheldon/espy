"""
    %prog [options]

shearnums 1-8
psfnums 1-6
"""

import sys
import os
from numpy import zeros, sqrt, where, median

import cluster_step
from cluster_step import files, stats
from cluster_step import sh1exp, sh2exp

import esutil as eu
from esutil.numpy_util import aprint
import biggles
from biggles import FramedPlot, FramedArray, Points, \
        Curve, SymmetricErrorBarsX, SymmetricErrorBarsY

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-m','--model',default=None,
                  help="which model to plot, required")
parser.add_option('-r','--run',default=None,
                  help='The run id, required, required')
parser.add_option('-s','--shnum',default=None,
                  help='The shear number, required, required')
parser.add_option('-p','--psfnum',default=None,
                  help='The psf num, required')

def get_grid(ntot):
    sq=int(sqrt(ntot))
    if ntot==sq*sq:
        return (sq,sq)
    elif ntot <= sq*(sq+1):
        return (sq,sq+1)
    else:
        return (sq+1,sq+1)

class PSFPlotter(object):
    def __init__(self):
        biggles.configure( 'default', 'fontsize_min', 1)
        options,args = parser.parse_args(sys.argv[1:])

        self.options=options

        if (options.run is None 
                or options.shnum is None
                or options.psfnum is None
                or options.model is None):
            parser.print_help()
            sys.exit(1)

        self.run=options.run
        self.shnum=int(options.shnum)
        self.psfnum=int(options.psfnum)
        self.model=options.model

        self._load_data()
        
        self.parsname=self.model+'_pars'
        if self.parsname not in self._data.dtype.names:
            raise ValueError("tag '%s' not found" % self.parsname)

    def doplot(self):
        data=self._data
        nobj=data.size

        if self.model=='admom':
            names=['e1','e2','T']
            npars=len(names)
            pars=zeros( (nobj, npars) )
            irr=data['admom_pars'][:,3]
            irc=data['admom_pars'][:,4]
            icc=data['admom_pars'][:,5]

            T=irr+icc
            e1=(icc-irr)/T
            e2=2*irc/T

            pars[:,0] = e1
            pars[:,1] = e2
            pars[:,2] = T

        else:
            names=['cen1','cen2','e1','e2','T','flux']
            pars=data[self.parsname]
            npars=pars.shape[1]


        nrow,ncol=get_grid(npars)
        plt=biggles.Table(nrow,ncol)
        plt.aspect_ratio=float(nrow)/ncol

        for i in xrange(npars):
            row=i/ncol
            col=i % ncol

            std=pars[:,i].std()
            binsize=0.2*std

            plti,hi=eu.plotting.bhist(pars[:,i], binsize=binsize, show=False,
                                      gethist=True)

            
            medval=median( pars[:,i] )
            err=std/sqrt(nobj)

            hmax=hi['hist'].max()
            
            medc=biggles.Curve([medval,medval],[0,hmax],color='blue',width=2)
            plti.add(medc)

            mess='%.3g +/- %.3g' % (medval,err)
            lab=biggles.PlotLabel(0.1,0.9,mess,halign='left',color='blue')
            plti.add(lab)

            plti.xlabel=names[i]

            plti.yrange=[-10,1.2*hmax]

            plt[row,col] = plti

        tit='{run} p{psfnum} s{shnum}'
        tit=tit.format(run=self.run,
                       psfnum=self.psfnum,
                       shnum=self.shnum)
        plt.title=tit
        plt.show()

    def _load_data(self):
        from esutil.numpy_util import combine_arrlist
        alldata=[]
        for ccd in xrange(1,62+1):
            datai=files.read_fits_output(run=self.run,
                                         ftype='psf',
                                         psfnum=self.psfnum,
                                         shnum=self.shnum,
                                         ccd=ccd)
            if datai is None:
                fname=files.get_output_path(run=self.run,
                                            ftype='psf',
                                            psfnum=self.psfnum,
                                            shnum=self.shnum,
                                            ccd=ccd)
                raise ValueError("error reading: %s")
            alldata.append(datai)

        data=combine_arrlist(alldata)
        flags_name=self.model+'_flags'
        w,=where(data[flags_name]==0)
        if w.size==0:
            raise ValueError("none with flags==0")
        data=data[w]
        self._data=data


def main():
    sp=PSFPlotter()
    sp.doplot()

main()
