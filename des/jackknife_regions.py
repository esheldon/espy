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

        ra,dec=trim_radec(ra, dec, self._des_region)
        ra,dec=randomize_radec(ra, dec)


        radec=numpy.zeros( (ra.size, 2) )
        radec[:,0] = ra
        radec[:,1] = dec

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

    def make_plots(self, write_eps=True, show=False):
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



        hpts=biggles.Histogram(h, x0=self._km.labels.min(), binsize=1)

        hplt=biggles.FramedPlot()
        hplt.add(hpts)
        hplt.xlabel='kmeans cluster'

        hepsfile=get_jackknife_centers_epsfile(self._des_region,
                                              self._njack,
                                              extra='hist')



        if write_eps:

            print("writing:",epsfile)
            plt.write_eps(epsfile)
            converter.convert(epsfile,dpi=150,verbose=True)

            print("writing:",hepsfile)
            hplt.write_eps(hepsfile)
            converter.convert(hepsfile,dpi=100,verbose=True)

        if show:
            plt.show()
            hplt.show()

def get_jackknife_regions(ra, dec, des_region, njack):
    """
    get the jackknife regions associated with the input ra,dec
    """
    import kmeans_radec

    radec=numpy.zeros( (ra.size, 2) )
    radec[:,0] = ra
    radec[:,1] = dec

    centers=read_jackknife_centers(des_region, njack)

    km=kmeans_radec.KMeans(centers)
        
    labels=km.find_nearest(radec)

    return labels

def add_jackknife_regions(data, des_region, njacks):
    """
    add a number of sets of jackknifes to the input recarray

    the column names are 'jack_region%d'

    parameters
    ----------
    data: array
        array with 'ra' and 'dec' columns
    des_region: string
        des region, e.g. 'spte'
    njacks: scalar, list or array
        Number of jackknife regions e.g. 100 or [60,80,100]
        There should be a corresponding center file
    """
    import esutil as eu

    njacks=numpy.array(njacks, ndmin=1, dtype='i4')

    dt=[]
    for njack in njacks:
        regname='jack_region%d' % njack
        dt.append( (regname, 'i4') )

    ndata=eu.numpy_util.add_fields(data, dt)

    radec=numpy.zeros( (data.size, 2) )
    radec[:,0] = data['ra']
    radec[:,1] = data['dec']

    for njack in njacks:
        print("njack:",njack)

        regname='jack_region%d' % njack

        labels=get_jackknife_regions(data['ra'],
                                     data['dec'],
                                     des_region,
                                     njack)

        ndata[regname] = labels

    return ndata


def trim_radec(ra, dec, des_region, get_index=False):
    """
    trim the data to the region of intereste
    """
    from esutil.numpy_util import between
    print("trimming to region:",des_region)

    conf=read_config('constants')
    reginfo=conf['regions'][des_region]

    ra_range=reginfo['ra_range']
    dec_range=reginfo['dec_range']

    w,=numpy.where(  between(ra, ra_range[0],ra_range[1])
                   & between(dec, dec_range[0], dec_range[1]) )

    if get_index:
        return w
    else:
        return ra[w], dec[w]

def randomize_radec(ra, dec):
    """
    randomize the catalog to help kmeans algorithm
    """

    print("randomizing catalog")
    r=numpy.random.random(ra.size)
    s=r.argsort()

    return ra[s], dec[s]
