"""

generate inputs for the weighting-cdim code, run it, and return the results.
Note this is only for the "calcweights" code

"""

from __future__ import print_function

import os

import numpy
from numpy import where, zeros

import esutil as eu
from esutil.numpy_util import where1
from esutil.stat import histogram

import tempfile

import recfile

def plothist_weights(weights):
    import biggles
    import esutil as eu
    normweights = weights/weights.max()

    w,=numpy.where(normweights > 0)
    print("plothist_weights: greater than zero: %d/%d" % (w.size,normweights.size))

    lw = numpy.log10(normweights[w])
    binsize = 0.1*lw.std()
    eu.plotting.bhist(lw, binsize=binsize, xlabel='log(relative weight)')


def plot_results1d(data1, data2, weights1, binsize, 
                   xmin=None, xmax=None, xlabel=None, title=None, epsfile=None, show=True):
    """
    compare the histograms at the input binsize

    Unless the domains are exactlyl the same, you should restrict xmin,xmax so
    that the normalizations will match correctly.

    """
    import biggles
    import esutil as eu

    if xmin is None:
        xmin = data2.min()
    if xmax is None:
        xmax = data2.max()

    nw=weights1/weights1.max()
    effnum = nw.sum()
    effperc = effnum/data1.size*100
    plabtext='effnum: %d/%d = %0.1f%%' % (effnum,data1.size,effperc)

    print("plotting hist match results")
    print("    histogramming data set 1")
    h1dict = eu.stat.histogram(data1, binsize=binsize, more=True, 
                               min=xmin, max=xmax)
    print("    histogramming data set 1 with weights")
    h1wdict = eu.stat.histogram(data1, binsize=binsize, 
                                min=xmin, max=xmax,
                                weights=weights1,
                                more=True)

    print("    histogramming data set 2")
    h2dict = eu.stat.histogram(data2, binsize=binsize, more=True,
                               min=xmin, max=xmax)

    h1=h1dict['hist']/float(h1dict['hist'].sum())
    h1w=h1wdict['whist']/float(h1wdict['whist'].sum())
    h2=h2dict['hist']/float(h2dict['hist'].sum())

    h2err = numpy.sqrt(h2dict['hist'])/float(h2dict['hist'].sum())
    hdiff = h2-h1w
    hdifferr = h2err


    arr=biggles.FramedArray(2,1)


    ph1 = biggles.Histogram(h1, binsize=binsize, x0=h1dict['low'][0],color='blue')
    ph1.label = 'dataset 1'

    ph1w = biggles.Histogram(h1w, binsize=binsize, x0=h1dict['low'][0], color='red', type='longdashed', width=2)
    ph1w.label = 'weighted 1'

    ph2 = biggles.Histogram(h2, binsize=binsize, x0=h2dict['low'][0], width=2)
    ph2.label = 'dataset 2'

    plt=arr[0,0]

    plt.add(ph1)
    plt.add(ph2)
    plt.add(ph1w)
    plt.xlabel=xlabel

    key=biggles.PlotKey(0.1,0.90,[ph1,ph2,ph1w],halign='left')
    plt.add(key)


    pltdiff=arr[1,0]

    phdiff = biggles.Points(h1dict['center'], hdiff)
    phdifferr = biggles.SymmetricErrorBarsY(h1dict['center'], hdiff, hdifferr)

    zero=biggles.Curve([xmin,xmax],[0,0])

    plab=biggles.PlotLabel(0.05,0.9,plabtext,halign='left')
    
    pltdiff.add(phdiff, phdifferr, zero, plab)
    pltdiff.xlabel = xlabel
    pltdiff.ylabel = 'hist2-weight hist1'

    arr.xlabel=xlabel
    arr.title=title

    if epsfile is not None:
        print("writing eps file:",epsfile)
        arr.write_eps(epsfile)

    if show:
        arr.show()

def hist_match(data1, data2, binsize):
    """
    The simplest method for histogram matching

    Generate a set of weights for data set 1 such that the distribution of
    observables are matched to dataset 2.  

    """

    weights1 = zeros(data1.size)
    min2=data2.min()
    max2=data2.max()

    h1,rev1 = histogram(data1, binsize=binsize, min=min2, max=max2, rev=True)
    h2 = histogram(data2, min=min2, max=max2, binsize=binsize)

    if h1.size != h2.size:
        raise ValueError("histogram sizes don't match: %d/%d" % (h1.size,h2.size))

    ratio = zeros(h1.size)
    #w=where1(h2 > 0)
    #ratio[w] = (h1[w]*1.0)/h2[w]
    w=where1(h1 > 0)
    ratio[w] = (h2[w]*1.0)/h1[w]

    # this is the weight for each object in the bin
    ratio /= ratio.max()

    for i in xrange(h1.size):
        if rev1[i] != rev1[i+1]:
            w1 = rev1[ rev1[i]:rev1[i+1] ]

            weights1[w1] = ratio[i]

    return weights1

class WeightCalculator(dict):
    def __init__(self, data1, data2):
        self.data1=data1
        self.data2=data2

        self.weight_data=None
        self.num_data=None

    def dtype1(self):
        return [('junk1','f8'),('junk2','f8'),('weight','f8'), ('data','f8')]
    def dtype2(self):
        return [('id','i8'), ('data','f8')]

    def read_type1(self, fname):
        dt = self.dtype1()
        with recfile.Open(fname,'r',dtype=dt, delim=' ') as robj:
            print("reading: ",fname)
            data=robj.read()
            return data

    def read_type2(self, fname):
        dt = self.dtype2()
        with recfile.Open(fname,'r',dtype=dt, delim=' ') as robj:
            print("reading: ",fname)
            data=robj.read()
            return data

    def read_num(self, fname):
        dt=[('id','i8'),('num','i8')]
        with recfile.Open(fname,'r',dtype=dt, delim=' ') as robj:
            print("reading: ",fname)
            data=robj.read()
        return data


    def calc_1pass(self, n_near, cleanup=True):
        """

        Generate a set of weights for data set 1 such that the distribution of
        observables are matched to dataset 2.  Use "nnear1" nearest neighbors
        on the first pass to remove zero weight objects from dataset 1, nnear2
        for the second.

        For 1 pass like this, you man want to ensure the domains of the two
        inputs overalap in a sensible way.  In principle the 2 pass method can
        do this for you, but it needs to be tuned.

        Dataset 1 is called the "training" set for photozs and dataset 2 is
        called the "photometric" sample. 

        The inputs for the calcweights code are wierd for the dataset 1, 
        consisting of 

            zspec extra weight n-dim-point

        For the weights calculation, the first three are not necessary at all.

        Inputs for dataset 2 is

            id n-dim-point

        Limitations:
            for now only support 1-d
        """

        
        data1=self.data1
        data2=self.data2
        
        ndim=len(data1.shape)
        if len(data1.shape) > 1 or len(data2.shape) > 1:
            raise ValueError("only support 1-d now")

        # for input weight is not used
        data1_dtype = self.dtype1()
        data2_dtype = self.dtype2()

        input1 = numpy.zeros(data1.size, dtype=data1_dtype)
        input2 = numpy.zeros(data2.size, dtype=data2_dtype)

        input1['data'] = data1

        input2['id'] = numpy.arange(data2.size)
        input2['data'] = data2

        front='wts{ndim}'.format(ndim=ndim)
        f1=tempfile.mktemp(prefix=front+'-input-data1-', suffix='.dat')
        f2=tempfile.mktemp(prefix=front+'-input-data2-', suffix='.dat')

        wfile   = tempfile.mktemp(prefix=front+'-weights-', suffix='.dat')
        numfile = tempfile.mktemp(prefix=front+'-num-', suffix='.dat')

        script_file = tempfile.mktemp(prefix=front+'-script-', suffix='.sh')

        try:
            print("writing temp file for data1:",f1)
            with recfile.Open(f1, mode='w', delim=' ') as fobj1:
                fobj1.write(input1)
            print("writing temp file for data2:",f2)
            with recfile.Open(f2, mode='w', delim=' ') as fobj2:
                fobj2.write(input2)


            command="calcweights{ndim} {train_file} {photo_file} {n_near} {weights_file} {num_file}"
            command=command.format(ndim         = ndim,
                                   train_file   = f1,
                                   photo_file   = f2,
                                   n_near       = n_near,
                                   weights_file = wfile,
                                   num_file     = numfile)

            if os.system(command) != 0:
                raise RuntimeError("command '%s' failed" % command)

            self.weight_data = self.read_type1(wfile)
            self.num_data = self.read_num(numfile)

        finally:
            if cleanup:
                for f in [f1,f2,wfile,numfile]:
                    if os.path.exists(f):
                        os.remove(f)


    def calc_2pass(self, n_near1, n_near2, cleanup=True):
        """

        Generate a set of weights for data set 1 such that the distribution of
        observables are matched to dataset 2.  Use "nnear1" nearest neighbors
        on the first pass to remove zero weight objects from dataset 1, nnear2
        for the second.

        Dataset 1 is called the "training" set for photozs and dataset 2 is
        called the "photometric" sample. 

        The inputs for the calcweights code are wierd for the dataset 1, 
        consisting of 

            zspec extra weight n-dim-point

        For the weights calculation, the first three are not necessary at all.

        Inputs for dataset 2 is

            id n-dim-point

        Limitations:
            for now only support 1-d
        """

        
        data1=self.data1
        data2=self.data2
        
        ndim=len(data1.shape)
        if len(data1.shape) > 1 or len(data2.shape) > 1:
            raise ValueError("only support 1-d now")

        # for input weight is not used
        data1_dtype = self.dtype1()
        data2_dtype = self.dtype2()

        input1 = numpy.zeros(data1.size, dtype=data1_dtype)
        input2 = numpy.zeros(data2.size, dtype=data2_dtype)

        input1['data'] = data1

        input2['id'] = numpy.arange(data2.size)
        input2['data'] = data2

        front='wts{ndim}'.format(ndim=ndim)
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

            script = _script_2pass.format(ndim=ndim,
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

        finally:
            if cleanup:
                if os.path.exists(f1):
                    os.remove(f1)
                if os.path.exists(f2):
                    os.remove(f2)
                if os.path.exists(script_file):
                    os.remove(script_file)




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

First run with n_near=$n_near1
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

if [ "$?" != "0" ]; then
    echo Error running awk.  Halting
    exit 45
fi

# second run
echo "
Second run with n_near=$n_near2
"

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

def test_nearest(n_near, limit=False, binnum=9, l=None, r=None):
    import lensing
    if l is None or r is None:
        l,r = load_test_data()

    binner=lensing.binning.N200Binner(12)
    print("selecting ",binner.bin_label(binnum))
    w=binner.select_bin(l, binnum)
    print("    kept %d/%d" % (w.size,l.size))

    # in this case, we know r covers a larger range
    if limit:
        print("limiting random range to match lenses")
        w,=where((r['z'] >= l['z'][w].min()) & (r['z'] <= l['z'][w].max()))
        wc=WeightCalculator(r['z'][w],l['z'][w])
    else:
        print("no limits applied")
        wc=WeightCalculator(r['z'],l['z'][w])

    wc.calc_1pass(n_near)
    return wc
    
def load_test_data():
    import lensing
    run='06'
    rrun='07'
    l=lensing.files.sample_read('collated',run)
    r=lensing.files.sample_read('collated',rrun)

    return l,r

def test_hist_match(binsize, binnum=9, l=None, r=None):
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

    weights = hist_match(r['z'], l['z'][w], binsize)

    effnum = weights.sum()
    effperc = effnum/r.size
    print("effective number: %d/%d = %0.2f" % (effnum,r.size, effperc))

    plot_results1d(r['z'], l['z'][w], weights, binsize)
    return weights
    

