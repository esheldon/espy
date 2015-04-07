from __future__ import print_function, division
import os
import numpy
from numpy import sqrt

GMAX_HIGHS2N=0.985

def fit_g_prior(run, type='great-des', show=False, gmax_data=GMAX_HIGHS2N, gmax=1.0):
    import esutil as eu
    import ngmix

    g=cache_data(run, gmax_data, 'g',
                 min_s2n=100.0,
                 min_Ts2n=30.0,
                 max_g=gmax_data)
    hd=eu.stat.histogram(g, nbin=100, min=0, max=gmax,more=True)

    xdata=hd['center']
    ydata=hd['hist'].astype('f8')

    if type=='no-exp':
        p=ngmix.priors.GPriorGreatDESNoExp(gmax=gmax)
    elif type=='great-des':
        p=ngmix.priors.GPriorGreatDES(gmax=gmax)
    elif type=='great-des2':
        p=ngmix.priors.GPriorGreatDES2()
    elif type=='ba':
        p=ngmix.priors.GPriorBA(0.2)
    elif type=='merf':
        p=ngmix.priors.GPriorMErf()
    else:
        raise ValueError("bad type '%s'" % type)

    p.dofit(hd['center'], hd['hist'].astype('f8'), show=show)

def fit_log_TF_prior(run, show=False, gmax_data=GMAX_HIGHS2N):
    import esutil as eu
    import ngmix

    log_TF=cache_data(run, gmax_data, 'log_TF',
                 min_s2n=10.0,
                 min_Ts2n=10.0,
                 max_g=gmax_data)



def cache_data(run, gmax_data, type, **keys):
    """
    type should be g or log_TF
    """
    import esutil as eu
    import fitsio
    from . import analyze, files

    tmpdir=os.environ['TMPDIR']
    cachename=os.path.join(tmpdir,'%s-%s.fits' % (run,type))

    if os.path.exists(cachename):
        print("reading:",cachename)
        tmp=fitsio.read(cachename)
    else:

        conf=files.read_config(run=run)

        dlist=[]
        for i in xrange(conf['ng']):
            fname=files.get_collated_file(run=run, gnum=i)

            print("reading:",fname)
            data=fitsio.read(fname)

            print("    selecting")
            data=analyze.select_good(data, **keys)

            dlist.append(data)

        data=eu.numpy_util.combine_arrlist(dlist)

        if type=='g':
            tmp=numpy.zeros(data.size, dtype=[('g','f4')])
            g=sqrt(data['g'][:,0]**2 + data['g'][:,1]**2)
            tmp['g'] = g
        else:
            tmp=numpy.zeros(data.size, dtype=[('log_TF','f4',2)])
            tmp['log_TF'][:,0] = data['log_T']
            tmp['log_TF'][:,1] = data['log_F']

        print("writing:",cachename)
        fitsio.write(cachename, tmp, clobber=True)

    return tmp[type]
