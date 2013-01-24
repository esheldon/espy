"""
Don't get fooled!  kappa is not a linear function of n!

The lookup table is faster and more accurate
"""
import os
import numpy

def get_kappa(n):
    """
    I=A*exp( -kappa*( (r/re)^{1/n} - 1 ) )

    kappa is chosen such that re is the half life radius.
    """
    ktable=KappaTable()
    return ktable(n)

class KappaTable(object):
    _table=None
    def __init__(self):
        if KappaTable._table is None:
            npoints=100
            nmin=0.25
            nmaxa=8.0
            nvals=numpy.logspace(numpy.log10(0.25), 
                                 numpy.log10(8.0), 
                                 npoints)
            table=numpy.zeros(npoints, dtype=[('n','f8'),('kappa','f8')])

            kfinder=KappaFinder()
            for i in xrange(npoints):
                n=nvals[i]
                table['n'][i] = n

                table['kappa'][i] = kfinder(n)
            KappaTable._table=table

    def __call__(self, n):
        kappa=numpy.interp(n, 
                           KappaTable._table['n'], 
                           KappaTable._table['kappa'])
        return kappa

class KappaFinder(object):
    """
    This uses a solver to find the minimum.

    I=A*exp( -kappa*( (r/re)^{1/n} - 1 ) )

    kappa is chosen such that re is the half life radius.

    It will be slow, so use this to make a lookup table
    """
    def __init__(self):
        pass

    def __call__(self, n):
        import scipy.optimize

        self._n=n
        kappa_start=numpy.array([3.67])
        kappa = scipy.optimize.fsolve(self._kappa_func, kappa_start)

        return kappa[0]

    def _kappa_func(self, kappa):
        from scipy.special import gamma,gammainc
        n=self._n
        val=gamma(2.*n)*(1-2*gammainc(2.*n, kappa))
        return val

    def __repr__(self):
        return '%.16g' % self._kappa


class KappaApprox(object):
    _ply=numpy.poly1d([9.06733525e-08,-4.62271203e-06,1.04365771e-04,
                       -1.37472817e-03, 1.17184236e-02,-6.77862447e-02, 
                       2.71475361e-01, -7.54487848e-01, 1.43751841e+00, 
                       -1.82717720e+00,  1.47845721e+00,  1.29710871e+00,
                       -1.66880136e-01])
    def __init__(self):
        pass
    def __call__(self, n):
        return KappaApprox._ply(n)


def _find_approx_kappa_vs_n(nvals, kappas, order):
    pcoeff=numpy.polyfit(nvals, kappas, order) 
    return pcoeff

def plot_kappa_vs_n(show=False, order=5):
    import biggles
    np=1000

    nmin=0.25
    nmax=8.0
    nvals=numpy.logspace(numpy.log10(nmin), numpy.log10(nmax), np)
    kappas=numpy.zeros(np)

    kfinder=KappaFinder()
    for i in xrange(np):
        kappas[i] = kfinder(nvals[i])
    
    tab=biggles.Table(2,1)

    plt=biggles.FramedPlot()
    plt.xlabel='sersic n'
    plt.ylabel=r'$\kappa$'
    plt.xlog=True
    plt.ylog=True

    pts=biggles.Points(nvals, kappas, type='filled circle',
                       size=0.5, color='blue')

    pcoeff=_find_approx_kappa_vs_n(nvals, kappas, order)
    print pcoeff
    ply=numpy.poly1d(pcoeff)

    yfit=ply(nvals)
    c=biggles.Curve(nvals, yfit)

    plt.add(pts,c)

    tab[0,0]=plt

    fplt=biggles.FramedPlot()
    fplt.xlabel='sersic n'
    fplt.ylabel=r'$\Delta \kappa/\kappa$'
    fplt.xlog=True
    fdiff=(yfit-kappas)/kappas
    fdiff_pts=biggles.Points(nvals, fdiff, type='filled circle',
                            size=0.5)
    lab=biggles.PlotLabel(0.9,0.9,'order: %s' % order,
                          halign='right')
    fplt.add(fdiff_pts,lab)

    tab[1,0] = fplt

    if show:
        tab.show()

    f=os.path.expanduser('~/tmp/kappa-vs-n.eps')
    print f
    plt.write_eps(f)

def test_speed():
    import time
    np=10000

    nmin=0.25
    nmax=8.0
    nvals= nmin + (nmax-nmin)*numpy.random.random(np)

    kapprox=KappaApprox()
    kfinder=KappaFinder()
    ktable=KappaTable()

    t0=time.time()
    for i in xrange(np):
        k=kfinder(nvals[i])
    tm=time.time()-t0
    print 'finder:',tm,'per:',tm/float(np)

    t0=time.time()
    for i in xrange(np):
        k=kapprox(nvals[i])
    tm=time.time()-t0
    print 'approx:',tm,'per:',tm/float(np)

    t0=time.time()
    for i in xrange(np):
        k=ktable(nvals[i])
    tm=time.time()-t0
    print 'table:',tm,'per:',tm/float(np)

def test_accuracy():
    import biggles
    np=100

    nmin=0.25
    nmax=8.0
    nvals= nmin + (nmax-nmin)*numpy.random.random(np)

    kapprox=KappaApprox()
    kfinder=KappaFinder()
    ktable=KappaTable()

    kfvals=numpy.zeros(np)
    kavals=numpy.zeros(np)
    ktvals=numpy.zeros(np)

    for i in xrange(np):
        kfvals[i]=kfinder(nvals[i])

    for i in xrange(np):
        kavals[i]=kapprox(nvals[i])

    for i in xrange(np):
        ktvals[i]=ktable(nvals[i])

    plt=biggles.FramedPlot()

    #fpts=biggles.Points(nvals, kfvals, color='red', type='filled circle', size=0.5)
    apts=biggles.Points(nvals, kavals/kfvals-1, color='blue', type='filled diamond', size=0.5)
    tpts=biggles.Points(nvals, ktvals/kfvals-1, color='magenta', type='circle', size=0.5)

    plt.add(apts, tpts)

    plt.show()
