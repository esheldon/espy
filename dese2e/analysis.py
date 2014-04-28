from __future__ import print_function
import os
from os import path

import numpy
from numpy import where
from . import files

class Namer(object):
    def __init__(self, front):
        self.front=front
    def __call__(self, suffix):
        return '%s_%s' % (self.front, suffix)

def select_bad_center(data, model):
    """
    Select centers with large offsets from (0,0)
    """

    n=Namer(model)
    flags=data[n("flags")]
    cen=data[n("pars")][:,0:0+2]

    w,=where(  (flags == 0)
             & (  (numpy.abs(cen[:,0]) > 2)
                | (numpy.abs(cen[:,1]) > 2) ) )

    return w

def plot_bad_center(**keys):
    """
    Plot a random subset of objects with bad centers.

    parameters
    ----------
    run: keyword, string
        The run identifier
    model: keyword, string
        Model name, e.g. 'exp'

    num: keyword, int
        Number to make images for.  
    """
    import images

    model=keys['model']
    num=keys['num']
    seed=keys.get('seed',35)
    numpy.random.seed(seed)

    data=files.read_output(**keys)
    m=files.open_meds()

    w=select_bad_center(data, model)

    r=numpy.random.random(w.size)
    s=r.argsort()

    wplot=w[s[0:num]]

    d=files.get_plot_dir(**keys)
    pf=path.join(d,'mosaic-%(run)s-%(index)06d.png')
    if not os.path.exists(d):
        os.makedirs(d)

    for i in xrange(num):
        index=wplot[i]


        im=m.get_mosaic(index)
        im=im.transpose()

        plt=images.view(im, title='%s' % index,show=False)

        keys['index'] = index
        fname=pf % keys

        print(fname)

        # view transposes image to match IDL
        xsize=1800
        ysize=xsize*float(im.shape[0])/im.shape[1]
        plt.write_img(xsize,ysize,fname)
