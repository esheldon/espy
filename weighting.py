"""
Codes to match the distributions of two populations

weight_match:
    derive weights so that the distributiins of two populations
    are proportional in N-dimensional space. Uses WeightNearest

WeightNearest:
    A class wrapping the nearest neighbor weighting code
    derive weights so that the distributiins of two populations
    are proportional in N-dimensional space

hist_match: 
    derive weights so that the boxcar histograms are proportional
    1-D only

hist_match_remove:
    remove objects until the boxcar histograms match.
    1-D only
"""

from __future__ import print_function

import os
import tempfile

import numpy
from numpy import where, zeros, ones, array, diag


try:
    import recfile
except:
    try:
        from esutil import recfile
    except:
        pass


def weight_match(data1, data2, n_near, npass=1, n_near2=None, cleanup=True):
    """
    Use the nearest neighbor method to calculate weights that force two
    distributions to be equal

    parameters
    ----------
    data1: numpy array
        The first dataset, to be matched to the second
    data2: numpy array
        The second dataset
    n_near: int
        Number of nearest neighbors to use for the calculation
    npass: int, optional
        Number of passes, 1 or 2.  Default 1.
    n_near2: int, optional
        Number of nearest neighbors to use in the second pass
    cleanup: bool, optional
        If True, clean up the temporary files.  default True

    returns
    -------
    A WeightNearest object
    """
    wn=WeightNearest(data1, data2, npass=npass)
    wn.go(n_near, n_near2=n_near2, cleanup=cleanup)

    return wn

class WeightNearest(dict):
    """
    Use the nearest neighbor method to calculate weights that force two
    distributions to be equal

    The first data set is matched to the second

    For 1 pass , you may want to ensure the domains of the two inputs overalap
    in a sensible way.  In principle the 2 pass method can do this for you, but
    it needs to be tuned.

    examples
    --------
    npass=1
    n_near=10
    wn=WeightNearest(data1, data2, npass)
    wn.go(n_near)
    wn.plot_results(show=True)

    # weights for the first data set that match the 2nd
    weights=wn.get_weights()

    # now try two passes
    n_near1=5
    n_near2=20
    wn=WeightNearest(data1, data2, npass)
    wn.go(n_near, n_near2)
    wn.plot_results(show=True)


    """
    def __init__(self, data1, data2, npass=1):

        if npass != 1 and npass != 2:
            raise ValueError("npass should be 1 or 2")

        self.data1=data1
        self.data2=data2
        self.npass=npass

        self.ndim=self._get_ndim(data1)
        ndim2=self._get_ndim(data2)
        mess="dims must be same, but %d %d" % (self.ndim,ndim2)
        assert self.ndim==ndim2,mess

    def get_weight_data(self):
        """
        get the full data structure for the weights
        """
        if not hasattr(self, 'weight_data'):
            raise RuntimeError("run the calcultion first")

        return self.weight_data

    def get_weights(self):
        """
        get the weights for dataset1 that match it to dataset 2
        """

        wdata=self.get_weight_data()
        return wdata['weight']

    def get_used_indices(self):
        """
        For 2 pass, the indices used for the second pass
        """
        if not hasattr(self, 'used_indices'):
            raise RuntimeError("run a 2 pass calculation first")

        return self.used_indices

    def get_num_data(self):
        """
        data on the number of neighbors
        """
        return self.num_data

    def go(self, n_near, n_near2=None, cleanup=True):
        """
        Generate a set of weights for data set 1 such that the distribution of
        observables are matched to dataset 2.

        If the number of passes was set to 2 during construction, an initial
        run will be made to remove zero-weight objects from the first set.
        In that case make sure to get the set of used indices with
        get_used_indices()

        Often it is useful to use a lower n_near for the first pass than the
        second
        
        parameters
        ----------
        n_near: int
            Number of nearest neighbors to use for the calculation
        n_near2: int, optional
            Number of nearest neighbors to use in the second pass
        cleanup: bool, optional
            If True, clean up the temporary files.  default True
        """

        if self.npass==2:
            if n_near2 is None:
                raise ValueError("send n_near2= for a 2 pass calculation")
            self._calc_2pass(n_near, n_near2, cleanup=cleanup)
        else:
            self._calc_1pass(n_near, cleanup=cleanup)

    def _calc_1pass(self, n_near, cleanup=True):
        
        data1=self.data1
        data2=self.data2

        # for input weight is not used
        data1_dtype = self._get_dtype1()
        data2_dtype = self._get_dtype2()

        input1 = numpy.zeros(data1.shape[0], dtype=data1_dtype)
        input2 = numpy.zeros(data2.shape[0], dtype=data2_dtype)

        #input1['id'] = numpy.arange(data1.shape[0])
        input1['data'] = data1

        input2['id'] = numpy.arange(data2.shape[0])
        input2['data'] = data2

        front='wts{ndim}'.format(ndim=self.ndim)
        f1=tempfile.mktemp(prefix=front+'-input-data1-', suffix='.dat')
        f2=tempfile.mktemp(prefix=front+'-input-data2-', suffix='.dat')

        wfile   = tempfile.mktemp(prefix=front+'-weights-', suffix='.dat')
        numfile = tempfile.mktemp(prefix=front+'-num-', suffix='.dat')

        try:
            print("writing temp file for data1:",f1)
            with recfile.Open(f1, mode='w', delim=' ') as fobj1:
                fobj1.write(input1)
            print("writing temp file for data2:",f2)
            with recfile.Open(f2, mode='w', delim=' ') as fobj2:
                fobj2.write(input2)


            command="calcweights{ndim} {train_file} {photo_file} {n_near} {weights_file} {num_file}"
            command=command.format(ndim         = self.ndim,
                                   train_file   = f1,
                                   photo_file   = f2,
                                   n_near       = n_near,
                                   weights_file = wfile,
                                   num_file     = numfile)

            if os.system(command) != 0:
                raise RuntimeError("command '%s' failed" % command)

            self.weight_data = self._read_type1(wfile)
            self.num_data = self._read_num(numfile)

        finally:
            if cleanup:
                for f in [f1,f2,wfile,numfile]:
                    if os.path.exists(f):
                        os.remove(f)


    def _calc_2pass(self, n_near1, n_near2, cleanup=True):
        
        data1=self.data1
        data2=self.data2

        # for input weight is not used
        data1_dtype = self._get_dtype1()
        data2_dtype = self._get_dtype2()

        input1 = numpy.zeros(data1.shape[0], dtype=data1_dtype)
        input2 = numpy.zeros(data2.shape[0], dtype=data2_dtype)

        input1['data'] = data1

        input2['id'] = numpy.arange(data2.shape[0])
        input2['data'] = data2

        front='wts{ndim}'.format(ndim=self.ndim)
        f1=tempfile.mktemp(prefix=front+'-input-data1-', suffix='.dat')
        f2=tempfile.mktemp(prefix=front+'-input-data2-', suffix='.dat')

        wfile_pass1        = tempfile.mktemp(prefix=front+'-weights-pass1-', suffix='.dat')
        wfile_nozero_pass1 = tempfile.mktemp(prefix=front+'-weights-nozero-pass1-', suffix='.dat')
        numfile_pass1      = tempfile.mktemp(prefix=front+'-num-pass1-', suffix='.dat')

        wfile_pass2        = tempfile.mktemp(prefix=front+'-weights-pass2-', suffix='.dat')
        numfile_pass2      = tempfile.mktemp(prefix=front+'-num-pass2-', suffix='.dat')


        script_file = tempfile.mktemp(prefix=front+'-script-', suffix='.sh')

        try:
            print("writing temp file for data1:",f1)
            with recfile.Open(f1, mode='w', delim=' ') as fobj1:
                fobj1.write(input1)
            print("writing temp file for data2:",f2)
            with recfile.Open(f2, mode='w', delim=' ') as fobj2:
                fobj2.write(input2)

            script = _script_2pass.format(ndim=self.ndim,
                                          train_file           = f1,
                                          photo_file           = f2,
                                          n_near1              = n_near1,
                                          n_near2              = n_near2,
                                          weights_file1        = wfile_pass1,
                                          num_file1            = numfile_pass1,
                                          weights_file_nozero1 = wfile_nozero_pass1,
                                          weights_file2        = wfile_pass2,
                                          num_file2            = numfile_pass2)

            print("writing script file:",script_file)
            with open(script_file,'w') as sobj:
                sobj.write(script)

            command="bash %s" % script_file
            if os.system(command) != 0:
                raise RuntimeError("command '%s' failed" % command)

            self.weight_data_pass1 = self._read_type1(wfile_pass1)

            self.used_indices, = where(self.weight_data_pass1['weight'] > 0)
            self.weight_data = self._read_type1(wfile_pass2)
            self.num_data = self._read_num(numfile_pass2)


        finally:
            if cleanup:
                for f in [
                          f1,f2,
                          wfile_pass1,
                          numfile_pass1,
                          wfile_nozero_pass1,
                          wfile_pass2,
                          numfile_pass2,
                          script_file,
                          ]:
                    if os.path.exists(f):
                        os.remove(f)

    def _get_dtype1(self):
        """
        data type for input 1
        """
        # first col is not used, we will make it an id
        #return [('id','i8'),('junk2','f8'),('weight','f8'), ('data','f8',self.ndim)]
        return [('junk1','f8'),('junk2','f8'),('weight','f8'), ('data','f8',self.ndim)]

    def _get_dtype2(self):
        """
        data type for input 2
        """
        return [('id','i8'), ('data','f8', self.ndim)]

    def _read_type1(self, fname):
        """
        read the type 1 file
        """
        dt = self._get_dtype1()
        with recfile.Open(fname,'r',dtype=dt, delim=' ') as robj:
            print("reading: ",fname)
            data=robj.read()
            return data

    def _read_type2(self, fname):
        """
        read the type 2 file
        """
        dt = self._get_dtype2()
        with recfile.Open(fname,'r',dtype=dt, delim=' ') as robj:
            print("reading: ",fname)
            data=robj.read()
            return data

    def _read_num(self, fname):
        """
        read the number file
        """
        dt=[('id','i8'),('num','i8')]
        with recfile.Open(fname,'r',dtype=dt, delim=' ') as robj:
            print("reading: ",fname)
            data=robj.read()
        return data



    def plot_results(self, binsizes=None, show=False):
        """
        plot a table with entries for each dimension, plus
        the weights histogram
        """
        import biggles
        import esutil as eu

        data1=self.data1
        data2=self.data2
        weights1=self.get_weights()
        if self.npass==1:
            ind=numpy.arange(weights1.size)
        else:
            ind=self.get_used_indices()

        nplots = self.ndim+1
        nrow,ncol=eu.plotting.get_grid(nplots)

        tab=biggles.Table(nrow,ncol)

        plist=[]
        for i in xrange(nplots):
            if binsizes is None:
                binsize=None
            else:
                binsize=binsizes[i]

            row = i/ncol
            col = i % ncol

            if i==(nplots-1):
                plt=plothist_weights(weights1,show=show)
            else:
                plt=plot_results1d(data1[ind,i],
                                   data2[:,i],
                                   weights1,
                                   binsize=binsize,
                                   show=show)
            plist.append(plt)

        return plist

    def _get_ndim(self, data):
        if len(data.shape) == 1:
            ndim=1
        elif len(data.shape)==2:
            ndim=data.shape[1]
        else:
            raise ValueError("data should have shape [npts] or [npts, ndim]")
        return ndim




def plothist_weights(weights, show=False):
    import biggles
    normweights = weights/weights.max()

    w,=numpy.where(normweights > 0)
    print("plothist_weights: greater than zero: %d/%d" % (w.size,normweights.size))

    lw = numpy.log10(normweights[w])
    #binsize = 0.1*lw.std()
    plt=biggles.plot_hist(lw,
                          nbin=50,
                          #binsize=binsize, 
                          xlabel=r'$log_{10}(relative weight)$',
                          visible=False)

    if show:
        plt.show()

    return plt


def plot_results1d(data1, data2, weights1, binsize=None, 
                   xmin=None, xmax=None, xlabel=None, title=None,
                   epsfile=None, pngfile=None, show=True,
                   label1='dataset 1',
                   label2='dataset 2'):
    """
    compare the histograms at the input binsize

    Unless the domains are exactlyl the same, you should restrict xmin,xmax so
    that the normalizations will match correctly.

    """
    import biggles
    from esutil.stat import histogram

    #if xmin is None:
    #    xmin = data2.min()
    #if xmax is None:
    #    xmax = data2.max()
    if xmin is None:
        xmin = min([data1.min(), data2.min()])
    if xmax is None:
        xmax = max([data1.max(), data2.max()])

    if binsize is None:
        w,=where( (data2 < xmax) & (data2 > xmin) )
        binsize=0.2*data2[w].std()

    nw=weights1/weights1.max()
    effnum = nw.sum()
    effperc = effnum/data1.size*100
    plabtext='effnum: %d/%d = %0.1f%%' % (effnum,data1.size,effperc)

    print("    plotting hist match results")
    print("    histogramming data set 1")
    h1dict = histogram(data1, binsize=binsize, more=True, 
                       min=xmin, max=xmax)
    print("    histogramming data set 1 with weights")
    h1wdict = histogram(data1, binsize=binsize, 
                        min=xmin, max=xmax,
                        weights=weights1,
                        more=True)

    print("    histogramming data set 2")
    h2dict = histogram(data2, binsize=binsize, more=True,
                       min=xmin, max=xmax)

    h1=h1dict['hist']/float(h1dict['hist'].sum())
    h1w=h1wdict['whist']/float(h1wdict['whist'].sum())
    h2=h2dict['hist']/float(h2dict['hist'].sum())

    hdiff = h2-h1w


    #arr=biggles.FramedArray(2,1)
    tab=biggles.Table(2,1)


    ph1 = biggles.Histogram(h1, binsize=binsize, x0=h1dict['low'][0],color='blue')
    ph1.label = label1

    ph1w = biggles.Histogram(h1w, binsize=binsize, x0=h1dict['low'][0], color='red', width=2)
    ph1w.label = label1+' weighted'

    ph2 = biggles.Histogram(h2, binsize=binsize, x0=h2dict['low'][0], width=2)
    ph2.label = label2

    #plt=arr[0,0]
    plt=biggles.FramedPlot()
    plt.title=title

    plt.add(ph1)
    plt.add(ph2)
    plt.add(ph1w)
    plt.xlabel=xlabel

    key=biggles.PlotKey(0.1,0.90,[ph1,ph2,ph1w],halign='left')
    plt.add(key)

    tab[0,0]=plt

    #pltdiff=arr[1,0]
    pltdiff=biggles.FramedPlot()

    phdiff = biggles.Points(h1dict['center'], hdiff)

    zero=biggles.Curve([xmin,xmax],[0,0])

    plab=biggles.PlotLabel(0.05,0.9,plabtext,halign='left')
    
    pltdiff.add(phdiff, zero, plab)
    pltdiff.xlabel = xlabel
    pltdiff.ylabel = '%s-%s weighted' % (label2, label1)

    pltdiff.title=title
    tab[1,0] = pltdiff

    #arr.xlabel=xlabel
    #arr.title=title

    if epsfile is not None:
        tab.write_eps(epsfile)
    if pngfile is not None:
        tab.write_img(800,800,pngfile)

    if show:
        tab.show()

    return tab

def hist_match(data1, data2, binsize, extra_weights1=None):
    """
    Generate a set of weights for data set 1 such that the distribution of
    observables are matched to dataset 2.  

    This is the simplest method for histogram matching and just works
    in rectangular bins

    parameters
    ----------
    data1:
        This data set is to be matched by weighting to data2
    data2:
        The data to be matched against
    binsize:
        The binsize to use for the histogram
    extra_weights1:
        An extra set of weights to apply to data1.  The returned weights
        will include this weight
    """

    weights1 = zeros(data1.size)
    min2=data2.min()
    max2=data2.max()

    if extra_weights1 is not None:
        bs1 = histogram(data1, binsize=binsize, min=min2, max=max2, rev=True,
                        weights=extra_weights1)
        h1=bs1['whist']
        rev1=bs1['rev']
    else:
        h1,rev1=histogram(data1, binsize=binsize, min=min2, max=max2, rev=True)

    h2 = histogram(data2, min=min2, max=max2, binsize=binsize)

    if h1.size != h2.size:
        raise ValueError("histogram sizes don't match: %d/%d" % (h1.size,h2.size))

    ratio = zeros(h1.size)
    w,=where(h1 > 0)
    ratio[w] = (h2[w]*1.0)/h1[w]

    # this is the weight for each object in the bin
    ratio /= ratio.max()

    for i in xrange(h1.size):
        if rev1[i] != rev1[i+1]:
            w1 = rev1[ rev1[i]:rev1[i+1] ]

            weights1[w1] = ratio[i]

    if extra_weights1 is not None:
        weights1 *= extra_weights1
    return weights1

def hist_match_remove(data1, data2, binsize, extra_weights1=None):
    """

    Similar to hist_match but instead of returning the weights, actually remove
    a random subset from data set 1

    """
    import esutil as eu
    min2=data2.min()
    max2=data2.max()

    if extra_weights1 is not None:
        bs1 = histogram(data1, binsize=binsize, min=min2, max=max2, rev=True,
                        weights=extra_weights1)
        h1=bs1['whist']
        rev1=bs1['rev']
    else:
        h1,rev1=histogram(data1, binsize=binsize, min=min2, max=max2, rev=True)

    h2 = histogram(data2, min=min2, max=max2, binsize=binsize)

    if h1.size != h2.size:
        raise ValueError("histogram sizes don't match: %d/%d" % (h1.size,h2.size))

    ratio = zeros(h1.size)
    w,=where(h1 > 0)
    ratio[w] = (h2[w]*1.0)/h1[w]

    # this is the weight for each object in the bin
    ratio /= ratio.max()

    keep=[]
    for i in xrange(h1.size):
        if rev1[i] != rev1[i+1]:
            w1 = rev1[ rev1[i]:rev1[i+1] ]

            # get a random subsample
            nkeep = int(w1.size*ratio[i])
            if nkeep > 0:
                # sort method is faster here.
                indices = eu.random.random_indices(w1.size, nkeep)
                keep.append(w1[indices])

    return eu.numpy_util.combine_arrlist(keep)




_script_2pass = """#!/bin/bash
photo_file="{photo_file}"
train_file="{train_file}"

n_near1={n_near1}
n_near2={n_near2}

weights_file1="{weights_file1}"
weights_file_nozero1="{weights_file_nozero1}"
weights_file2="{weights_file2}"

num_file1="{num_file1}"
num_file2="{num_file2}"

# first run

echo "
Running calcweights
"

calcweights{ndim}         \\
    "$train_file"    \\
    "$photo_file"    \\
    "$n_near1"       \\
    "$weights_file1" \\
    "$num_file1"

if [ "$?" != "0" ]; then
    echo Halting
    exit 45
fi


# remove training set objects with zero weight at first
# n_near
echo "
removing zero weight training objects to
    "$weights_file_nozero1"
"
awk '$3>0' < "$weights_file1" 1> "$weights_file_nozero1"
#wc -l $weights_file1 $weights_file_nozero1

if [ "$?" != "0" ]; then
    echo Error running awk.  Halting
    exit 45
fi

# second run

calcweights{ndim}                \\
    "$weights_file_nozero1" \\
    "$photo_file"           \\
    "$n_near2"              \\
    "$weights_file2"        \\
    "$num_file2"

if [ "$?" != "0" ]; then
    echo Halting
    exit 45
fi
    """

def test_gaussians(npass=1,
                   n_near1=5,
                   n_near2=100,
                   n1=100000,
                   n2=10000,
                   show=False, epsfile=None):
    """
    test matching a set of 5-d gaussians to another
    """
    from numpy.random import multivariate_normal
    from esutil.random import srandu
    ndim=5

    # we want to weight this set to look like the second
    ngauss1=3
    npergauss1=n1/ngauss1
    ntot1=npergauss1*ngauss1
    cen1=[ 1.0+0.1*srandu(ndim),
           0.9+0.1*srandu(ndim),
           1.1+0.1*srandu(ndim)]
    sig1=[ [0.75]*ndim,
           [1.5]*ndim,
           [2]*ndim ]

    data1=zeros( (ntot1, ndim) )

    ngauss2=3
    npergauss2=n1/ngauss1
    ntot2=npergauss2*ngauss2

    cen2=[ 1.0+0.1*srandu(ndim),
           1.0+0.1*srandu(ndim),
           1.0+0.1*srandu(ndim)]
    sig2=[ [0.5]*ndim,
           [0.75]*ndim,
           [1.2]*ndim ]

    data2=zeros( (ntot2, ndim) )

    for i in xrange(ngauss1):
        beg=i*npergauss1
        end=(i+1)*npergauss1

        c=array(cen1[i])
        s2=array(sig1[i])**2
        cov=diag(s2)
        pts = multivariate_normal( c, cov, npergauss1)
        data1[beg:end,:] = pts

    for i in xrange(ngauss2):
        beg=i*npergauss2
        end=(i+1)*npergauss2

        c=array(cen2[i])
        s2=array(sig2[i])**2
        cov=diag(s2)
        pts = multivariate_normal( c, cov, npergauss2)
        data2[beg:end,:] = pts

    if npass==1:
        wn=WeightNearest(data1, data2, npass=npass)
        wn.go(n_near1)
    else:
        wn=WeightNearest(data1, data2, npass=npass)
        wn.go(n_near1, n_near2)

    plist=wn.plot_results(show=show)
    return plist

def test_hist_match(binsize, binnum=9, l=None, r=None, remove=False):
    import lensing
    import biggles
    biggles.configure('screen','width', 1140)
    biggles.configure('screen','height', 1140)


    if l is None or r is None:
        l,r = load_test_data()

    binner=lensing.binning.N200Binner(12)
    print("selecting ",binner.bin_label(binnum))
    w=binner.select_bin(l, binnum)
    print("    kept %d/%d" % (w.size,l.size))

    if remove:
        keep = hist_match_remove(r['z'], l['z'][w], binsize)
        perc = keep.size/(1.*r.size)
        print("used number: %d/%d = %0.2f" % (keep.size,r.size, perc))

        weights = ones(keep.size)
        plot_results1d(r['z'][keep], l['z'][w], weights, binsize)
    else:
        weights = hist_match(r['z'], l['z'][w], binsize)

        effnum = weights.sum()
        effperc = effnum/r.size
        print("effective number: %d/%d = %0.2f" % (effnum,r.size, effperc))

        plot_results1d(r['z'], l['z'][w], weights, binsize)
    return weights
    

