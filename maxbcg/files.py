from __future__ import print_function
import os
import sdsspy
import es_sdsspy

import esutil as eu
from esutil.ostools import path_join

import maxbcg

def read_catalog(type):
    f = catalog_file(type)
    return eu.io.read(f, verbose=True, lower=True)

def catalog_file(type):
    if type == 'public-v1':
        name = 'maxbcg_public_catalog.fit'
    else:
        raise ValueError("don't support catalog type '%s' yet" % type)
    d = catalog_dir()
    f = path_join(d, name)
    return f

def catalog_dir():
    """
    location of the public catalog and others
    """
    return os.environ['MAXBCG_CATDIR']

def open_input_columns():
    import columns
    d=input_coldir()
    return columns.Columns(d)

def input_coldir():
    basedir = os.environ['MAXBCG_INPUT']
    prefix = maxbcg_prefix()

    coldir='maxbcg-input-%s.cols' % prefix
    coldir = path_join(basedir, coldir)
    return coldir
    
def maxbcg_prefix():
    photo_sweep = os.environ['PHOTO_SWEEP']
    if photo_sweep.find('dr8') != -1:
        return 'dr8'
    else:
        raise ValueError("currently only know prefix for dr8 sweeps")


def compare_public_catalog_to_input(data=None):
    import biggles
    retdata = False
    if data is None:
        retdata=True
        data={}
        cat = read_catalog('public-v1')

        data['cat'] = cat
        data['unmasked'] = maxbcg.select.do_radcut(cat)
        data['cat_ra_wrap'] = sdsspy.sdss_wrap(cat['ra'])

        c=open_input_columns()
        print('reading ra')
        ra = c['ra'][:]
        print('reading dec')
        dec = c['dec'][:]
        print("wrapping ra")
        ra_wrap = sdsspy.sdss_wrap(ra)

        data['ra'] = ra
        data['ra_wrap'] = ra_wrap
        data['dec'] = dec

        print("random subset")
        data['irand'] = eu.numpy_util.random_subset(data['ra'].size, data['ra'].size*0.005, verbose=True)

    plt = biggles.FramedPlot()

    ip = biggles.Points(data['ra_wrap'][data['irand']],
                        data['dec'][data['irand']],
                        type='dot')
    cp = biggles.Points(data['cat_ra_wrap'], data['cat']['dec'], type='dot', color='red')
    w=data['unmasked']
    cp_um = biggles.Points(data['cat_ra_wrap'][w], data['cat']['dec'][w], type='dot', color='magenta')

    # create labels
    ipfake = biggles.Points([-1000],[-1000], type='filled circle')
    ipfake.label = 'input'
    cpfake = biggles.Points([-1000],[-1000], type='filled circle', color='red')
    cpfake.label = 'MaxBCG'
    cpfake_um = biggles.Points([-1000],[-1000], type='filled circle', color='magenta')
    cpfake_um.label = 'MaxBCG Edgecut'

    key = biggles.PlotKey(0.1,0.9,[ipfake,cpfake, cpfake_um])

    plt.xrange = [-80,270] 
    plt.yrange = [-15,75]
    plt.xlabel = 'RA'
    plt.ylabel = 'DEC'

    plt.add(ip, cp, cp_um, key)

    plt.show()

    if retdata:
        return data
    else:
        return plt
