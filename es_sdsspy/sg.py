"""

Find objects with > N observations.  Start by limiting to the stripe82 runs,
and then just read through them. Save a collated file of 

    combined flux in each band
    single epoch fluxes.
    seeing

Based on the combined flux, we can do excellent s/g separation.  Then bin the
objects by their single-epoch flux and seeing and look at the distribution of
concentration.

"""
from __future__ import print_function

import numpy
from numpy import log10
import esutil as eu
from esutil.numpy_util import where1
import sdsspy
from . import stomp_maps

def stripe82_run_list():
    """
    We need to get these fcombined runs
    """
    win = sdsspy.window.Window()
    fl = win.read('flist')

    map=stomp_maps.load('boss','basic')
    inmap = map.Contains(fl['ra'],fl['dec'],'eq')

    w=where1(  (fl['dec'] >= -1.25) 
             & (fl['dec'] <= 1.25) 
             & ((fl['ra'] > 300) | (fl['ra'] < 60)) 
             & (fl['score'] >= 0.1) 
             & (inmap==1))
    uruns = numpy.unique(fl['run'][w])
    print("number of unique runs:",uruns.size)
    return uruns


def avg_gri(flux_g, ivar_g, 
            flux_r, ivar_r, 
            flux_i, ivar_i):
    """

    NOTE: you should use the same weights (ivar) for the psf and modelfluxes
    for things like the concentration.

    if you don't want to use an object in a particular band,  you should
    set the ivar to zero *before* calling this function

    fluxes are clipped between [0.001,1.e4], which is magnitude [12.5,30]
    """

    # clip the fluxes on the high and low end
    # this is mag between 12.5 and 30
    flux_g = flux_g.clip(0.001, 1.e4)
    flux_r = flux_r.clip(0.001, 1.e4)
    flux_i = flux_i.clip(0.001, 1.e4)

    # clip ivar as well, although this should not really be a problem as ivar
    # seems to always be well behaved

    ivar_g = ivar_g.clip(0.0, 70)
    ivar_r = ivar_r.clip(0.0, 70)
    ivar_i = ivar_i.clip(0.0, 70)

    ivarsum = ivar_g + ivar_r + ivar_i

    fsum = flux_g*ivar_g + flux_r*ivar_r + flux_i*ivar_i

    flux = fsum/ivarsum

    return flux, ivarsum


def calc_c(modelflux, psfflux, log=False):
    if log:
        return -log10(psfflux/modelflux)
    else:
        return 1.0-psfflux/modelflux


def read_test_data():
    import esutil as eu
    from . import select
    gal=eu.io.read('~/data/boss/calibObj-000756-3-gal.fits',lower=True)
    star=eu.io.read('~/data/boss/calibObj-000756-3-star.fits',lower=True)

    # this is about 21 mag

    g_rflux_logic = gal['modelflux'][:,2] > 4.0
    s_rflux_logic = star['modelflux'][:,2] > 4.0

    gsel = select.Selector(gal)
    ssel = select.Selector(star)

    gfl = gsel.flag_logic()
    gol = gsel.object1_logic()
    gbl = gsel.binned_logic()

    sfl = ssel.flag_logic()
    sol = ssel.object1_logic()
    sbl = ssel.binned_logic()

    gw = where1(g_rflux_logic & gfl & gol & gbl)
    sw = where1(s_rflux_logic & sfl & sol & sbl)

    ntot = gw.size + sw.size
    dt=[('origtype','i4'),
        ('modelflux','f4',5),
        ('modelflux_ivar','f4',5),
        ('psfflux','f4',5),
        ('psfflux_ivar','f4',5)]

    data = numpy.zeros(ntot, dtype=dt)
    data['origtype'][0:gw.size] = 3
    data['origtype'][gw.size:] = 6
    for n in ['modelflux','modelflux_ivar','psfflux','psfflux_ivar']:
        data[n][0:gw.size] = gal[n][gw]
        data[n][gw.size:] = star[n][sw]

    return data

def test(data=None, logc=False):
    from biggles import Histogram, FramedPlot, PlotKey, Table
    if data is None:
        data = read_test_data()

    
    modelflux, modelflux_ivar = avg_gri(data['modelflux'][:,1],data['modelflux_ivar'][:,1],
                                        data['modelflux'][:,2],data['modelflux_ivar'][:,2],
                                        data['modelflux'][:,3],data['modelflux_ivar'][:,3])
    #psfflux, psfflux_ivar = avg_gri(data['psfflux'][:,1],data['psfflux_ivar'][:,1],
    #                                data['psfflux'][:,2],data['psfflux_ivar'][:,2],
    #                                data['psfflux'][:,3],data['psfflux_ivar'][:,3])
    psfflux, psfflux_ivar = avg_gri(data['psfflux'][:,1],data['modelflux_ivar'][:,1],
                                    data['psfflux'][:,2],data['modelflux_ivar'][:,2],
                                    data['psfflux'][:,3],data['modelflux_ivar'][:,3])

    fmin_log=1.e-3
    fmax_log=4.


    tab = Table(2,2)
    binsize=0.05
    col=0
    for type in ['modelflux','psfflux']:

        flux_plt = FramedPlot()

        h = eu.stat.histogram(log10(data[type][:,1]), min=fmin_log, max=fmax_log, binsize=binsize)
        gmod_h = Histogram(h, x0=fmin_log, binsize=binsize, color='green')
        gmod_h.label = 'g '+type

        h = eu.stat.histogram(log10(data[type][:,2]), min=fmin_log, max=fmax_log, binsize=binsize)
        rmod_h = Histogram(h, x0=fmin_log, binsize=binsize, color='red')
        rmod_h.label = 'r '+type

        h = eu.stat.histogram(log10(data[type][:,3]), min=fmin_log, max=fmax_log, binsize=binsize)
        imod_h = Histogram(h, x0=fmin_log, binsize=binsize, color='magenta')
        imod_h.label = 'i '+type


        if type == 'modelflux':
            h = eu.stat.histogram(log10(modelflux), min=fmin_log, max=fmax_log, binsize=binsize)
        else:
            h = eu.stat.histogram(log10(psfflux), min=fmin_log, max=fmax_log, binsize=binsize)

        mod_h = Histogram(h, x0=fmin_log, binsize=binsize, width=2)
        mod_h.label = 'combined '+type

        key = PlotKey(0.5,0.9,[gmod_h, rmod_h, imod_h, mod_h])
        
        flux_plt.add(gmod_h, rmod_h, imod_h, mod_h, key)
        flux_plt.xlabel = 'flux'

        tab[0,col] = flux_plt

        col += 1



    col=0
    for logc in [False,True]:
        if logc:
            xmin=-0.1
            #xmax=1
            xmax=0.6
            binsize = 0.01
        else:
            xmin=-0.1
            xmax=1.0
            binsize = 0.01

        gc = calc_c(data['modelflux'][:,1], data['psfflux'][:,1], log=logc)
        rc = calc_c(data['modelflux'][:,2], data['psfflux'][:,2], log=logc)
        ic = calc_c(data['modelflux'][:,3], data['psfflux'][:,3], log=logc)
        allc = calc_c(modelflux, psfflux, log=logc)
        
        c_plt = FramedPlot()

        h = eu.stat.histogram(gc, min=xmin, max=xmax, binsize=binsize)
        gch = Histogram(h, x0=xmin, binsize=binsize, color='green')
        gch.label = 'g'

        h = eu.stat.histogram(rc, min=xmin, max=xmax, binsize=binsize)
        rch = Histogram(h, x0=xmin, binsize=binsize, color='red')
        rch.label = 'r'

        h = eu.stat.histogram(ic, min=xmin, max=xmax, binsize=binsize)
        ich = Histogram(h, x0=xmin, binsize=binsize, color='magenta')
        ich.label = 'i'


        h = eu.stat.histogram(allc, min=xmin, max=xmax, binsize=binsize)
        ch = Histogram(h, x0=xmin, binsize=binsize, width=2)
        ch.label = 'combined '

        key = PlotKey(0.7,0.9,[gch, rch, ich, ch])
        
        c_plt.add(gch, rch, ich, ch, key)
        if logc:
            c_plt.xlabel = r'$-log_{10}(psf/model)$'
        else:
            c_plt.xlabel = '1-psf/model'

        tab[1,col] = c_plt
        col += 1

    tab.show()


