from __future__ import print_function
import os
from os import path

import numpy
from numpy import where
from . import files


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

def plot_images(**keys):
    """
    Plot a random subset of objects, potentially with 
    some selection

    parameters
    ----------
    run: keyword, string
        The run identifier
    model: keyword, string
        Model name, e.g. 'exp'
    imtype: keyword, string
        'mosaic' or 'composite'
    num: keyword, int
        Number to make images for.  

    type: keyword, string
        'bad-center' or 'any'
    """
    import images

    model=keys['model']
    num=keys['num']
    seed=keys.get('seed',35)
    keys['type']=keys.get('type','any')
    keys['imtype']=keys.get('imtype','mosaic')

    numpy.random.seed(seed)

    data=files.read_output(**keys)
    m=files.open_meds()

    if keys['type']=='bad-center':
        w=select_bad_center(data, model)
    else:
        w=numpy.arange(data.size)

    r=numpy.random.random(w.size)
    s=r.argsort()

    wplot=w[s[0:num]]

    d=files.get_plot_dir(**keys)
    pf=path.join(d,'%(imtype)s-%(type)s-%(run)s-%(index)06d.png')
    if not os.path.exists(d):
        os.makedirs(d)

    for i in xrange(num):
        index=wplot[i]


        im=m.get_mosaic(index)

        if keys['imtype']=='composite':
            cw=m.get_cseg_mosaic(index)
            im *= cw
        elif keys['imtype'] == 'interp':
            iseg=m.interpolate_coadd_seg_mosaic(index)
            w=where( (iseg != (index+1)) & (iseg != 0) )
            im[w] = 0.0


        plt=images.view(im, title='%s' % index,show=False)

        keys['index'] = index
        fname=pf % keys

        print(fname)

        xsize=1800
        ysize=xsize*float(im.shape[1])/im.shape[0]
        plt.write_img(xsize,ysize,fname)

class Namer(object):
    def __init__(self, front):
        self.front=front
    def __call__(self, suffix):
        return '%s_%s' % (self.front, suffix)

