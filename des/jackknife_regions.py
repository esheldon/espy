"""
code to work with jackknife regions
"""
from __future__ import print_function
import numpy
from numpy import where, log10
from .files import *

class JackknifeMaker(dict):
    """
    make jackknife regions from ra,dec in the specified file
    """
    def __init__(self, ra, dec, des_region, njack, maxiter=100):
        self._des_region=des_region
        self._tol=1.0e-5
        self._njack=njack
        self._maxiter=maxiter

        self.set_data(ra, dec)

    def set_data(self, ra, dec):
        """
        set the data, trimmed and randomized
        """

        radec=numpy.zeros( (ra.size, 2) )
        radec[:,0] = ra
        radec[:,1] = dec

        radec=self._trim_cat(radec)
        radec=self._randomize_cat(radec)

        print("data size:",radec.shape[0])

        self._radec=radec

    def run_kmeans(self):
        """
        run the kmeans algorithm, checking for convergence
        """
        import kmeans_radec
        self._km=kmeans_radec.kmeans_sample(self._radec,
                                            self._njack,
                                            maxiter=self._maxiter,
                                            tol=self._tol)
        if not self._km.converged:
            raise RuntimeError("did not converge")


    def _trim_cat(self, radec):
        """
        trim the data to the region of intereste
        """
        from esutil.numpy_util import between
        print("trimming to region:",self._des_region)

        conf=read_config('constants')
        reginfo=conf['regions'][self._des_region]

        ra_range=reginfo['ra_range']
        dec_range=reginfo['dec_range']

        w,=numpy.where(  between(radec[:,0], ra_range[0],ra_range[1])
                       & between(radec[:,1], dec_range[0], dec_range[1]) )

        radec=radec[w,:]
        return radec

    def _randomize_cat(self, radec):
        """
        randomize the catalog to help kmeans algorithm
        """

        print("randomizing catalog")
        r=numpy.random.random(radec.shape[0])
        s=r.argsort()
        radec=radec[s,:]

        return radec

    def write_centers(self, orig_file=None):
        """
        write out the centers
        """
        import fitsio

        d=get_jackknife_centers_dir(self._des_region)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        fname=get_jackknife_centers_file(self._des_region,
                                         self._njack)
        print("writing:",fname)
        h={'tol':self._tol,
           'njack':self._njack,
           'des_region':self._des_region,
           'maxiter':self._maxiter}
        if orig_file is not None:
            h['file']=orig_file

        fitsio.write(fname, self._km.centers, clobber=True, header=h)

    def make_plots(self):
        """
        plot regions and histogram
        """
        import biggles
        import esutil as eu
        import converter
        import pcolors
        # now make a plot
        ncolor=100
        all_colors=pcolors.rainbow(ncolor)


        points=['filled circle','diamond','square','circle','filled diamond',
               'cross','star','octagon','filled square']

        npoint=len(points)

        h,rev=eu.stat.histogram(self._km.labels, rev=True)

        ri=eu.random.random_indices(ncolor, h.size)
        colors=[all_colors[i] for i in ri]

        nbin=h.size

        plt=biggles.FramedPlot()
        size=0.7

        for i in xrange(nbin):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                color=colors[i]
                ptype=points[i % npoint]
                pts=biggles.Points(self._radec[w,0],
                                   self._radec[w,1],
                                   type=ptype,
                                   color=color,
                                   size=size)

                plt.add(pts)

        cpts=biggles.Points(self._km.centers[:,0],
                            self._km.centers[:,1],
                            type='filled circle',
                            color='black')
        plt.add( cpts )

        ra=self._radec[:,0]
        dec=self._radec[:,1]

        rarange=ra.max()-ra.min()
        decrange=dec.max()-dec.min()
        plt.aspect_ratio = decrange/rarange
        plt.xlabel='RA'
        plt.ylabel='DEC'

        epsfile=get_jackknife_centers_epsfile(self._des_region,
                                              self._njack,
                                              extra='regions')

        print("writing:",epsfile)
        plt.write_eps(epsfile)

        #converter.convert(epsfile,dpi=150,verbose=True)


        hpts=biggles.Histogram(h, x0=self._km.labels.min(), binsize=1)

        hplt=biggles.FramedPlot()
        hplt.add(hpts)
        hplt.xlabel='kmeans cluster'

        epsfile=get_jackknife_centers_epsfile(self._des_region,
                                              self._njack,
                                              extra='hist')

        print("writing:",epsfile)
        hplt.write_eps(epsfile)
        #converter.convert(epsfile,dpi=100,verbose=True)

