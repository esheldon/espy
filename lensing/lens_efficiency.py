from __future__ import print_function
import cosmology
import lensing
import esutil as eu

from numpy import sqrt, zeros

def plot_lens_s2n_bysample(lens_sample, pzrun, cumulative=True):
    import zphot
    conf=lensing.files.read_config('lcat', lens_sample)
    zconf=zphot.cascade_config(pzrun)

    #
    # generate the expected n(z)
    #

    # read the weights and create weighted histogram
    wstruct = zphot.weighting.read_weights(zconf['pofz']['wrun'], 2)

    binsize=0.0314
    zmax=1.1
    hdict = eu.stat.histogram(wstruct['z'], binsize=binsize, weights=wstruct['weight'],max=zmax)

    zs = hdict['center']
    pzs = hdict['whist']

    # now read the zl
    data = lensing.files.read_original_catalog(type='lens',sample=lens_sample)

    return plot_lens_s2n(data['z'], hdict['low'], hdict['high'], pzs)

def plot_lens_s2n(zl, zslow, zshigh, pzs, cumulative=True):
    """

    Calculate the cumulative expected usefulness of lenses as a function of
    redshift for the given source N(z).

    """

    import biggles
    from biggles import FramedPlot, Curve, PlotLabel, PlotKey, Table, Histogram

    biggles.configure('screen','width',1140)
    biggles.configure('screen','height',1140)
    dz = zshigh[0]-zslow[0]
    zlmin = 0.0
    zlmax = max([zl.max(), zshigh.max()])

    # first get the mean inverse critical density as a function
    # of lens redshift
    zsmid = (zshigh+zslow)/2.
    sc = lensing.sigmacrit.ScinvCalculator(zsmid, dz, zlmin, zlmax)

    mean_scinv = sc.calc_mean_scinv(pzs)


    # histogram the lens redshifts in a binning corresponding
    # to the mean_scinv

    hdict = eu.stat.histogram(zl, min=zlmin, max=zlmax, binsize=dz, more=True)

    # get distances to the bin centers
    cosmo=cosmology.Cosmo()
    Dlens = cosmo.Da(0.0, hdict['center'])

    # interpolate the mean_scinv onto the centers
    mean_scinv = eu.stat.interplin(mean_scinv, sc.zlvals, hdict['center'])

    # this is the S/N at a given lens redshift
    s2n = sqrt(hdict['hist']/Dlens**2)*mean_scinv

    hd = hdict['hist']/Dlens**2
    hdcum = hd.cumsum()
    s2ncum = (mean_scinv*hd).cumsum()/sqrt(hdcum.clip(1.,hdcum.max()))


    s2n /= s2n.max()
    s2ncum /= s2ncum.max()

    
    tab = Table(2,2)
    tab.aspect_ratio=1

    # redshift histograms
    tab[0,0] = FramedPlot()
    zl_histp = Histogram(hdict['hist']/float(hdict['hist'].sum()), x0=hdict['low'][0], binsize=dz, color='blue')
    zl_histp.label = 'lenses'
    zs_histp = Histogram(pzs/float(pzs.sum()), x0=zslow[0], binsize=dz, color='red')
    zs_histp.label = 'sources'

    hkey=PlotKey(0.9,0.9,[zl_histp,zs_histp],halign='right')
    tab[0,0].add(zl_histp, zs_histp, hkey)
    tab[0,0].xlabel = 'z'

    # mean inverse critical density
    tab[0,1] = FramedPlot()
    tab[0,1].xlabel = 'lens redshift'
    tab[0,1].ylabel = r'$\langle \Sigma_{crit}^{-1} \rangle$'

    mc = Curve(hdict['center'], mean_scinv)
    tab[0,1].add(mc)


    # cumulative volume
    tab[1,0] = FramedPlot()
    tab[1,0].xlabel = 'z'
    tab[1,0].ylabel = 'cumulative volume'
    v=zeros(hdict['center'].size)
    for i in xrange(v.size):
        v[i] = cosmo.V(0.0, hdict['center'][i])
    v /= v.max()
    cv=Curve(hdict['center'], v)
    tab[1,0].add(cv)


    # S/N
    tab[1,1] = FramedPlot()
    tab[1,1].xlabel = 'lens redshift'
    tab[1,1].ylabel = 'S/N'

    cs2n = Curve(hdict['center'], s2n, color='blue')
    cs2ncum = Curve(hdict['center'], s2ncum, color='red')

    cs2n.label = 'S/N'
    cs2ncum.label = 'Cumulative S/N'
    ckey=PlotKey(0.9,0.9,[cs2n, cs2ncum],halign='right')

    tab[1,1].add(cs2n, cs2ncum, ckey)


    tab.show()

    return tab
