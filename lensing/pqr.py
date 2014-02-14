import esutil as eu
from math import sqrt, atan2, tanh, atanh
import copy
from . import util
from .util import ShapeRangeError
import numpy

def get_shear_pqr_sums(P,Q,R):
    """
    Create the sums used to calculate shear from BA13 PQR
    """
    from numpy import zeros,where

    w,=where(P > 0)
    if w.size != P.size:
        print 'P > 0: %s/%s' (w.size,P.size)

    QQ = zeros( (w.size,2,2), dtype=Q.dtype)
    Cinv_all = zeros( (w.size,2,2), dtype=Q.dtype)
    QbyP = zeros( (w.size,2), dtype=Q.dtype)

    # outer product
    QQ[:,0,0] = Q[w,0]*Q[w,0]
    QQ[:,0,1] = Q[w,0]*Q[w,1]
    QQ[:,1,0] = Q[w,1]*Q[w,0]
    QQ[:,1,1] = Q[w,1]*Q[w,1]

    Pinv = 1/P[w]
    P2inv = Pinv*Pinv

    #Cinv_all = QQ/P**2 - R/P
    Cinv_all[:,0,0] = QQ[:,0,0]*P2inv - R[w,0,0]*Pinv
    Cinv_all[:,0,1] = QQ[:,0,1]*P2inv - R[w,0,1]*Pinv
    Cinv_all[:,1,0] = QQ[:,1,0]*P2inv - R[w,1,0]*Pinv
    Cinv_all[:,1,1] = QQ[:,1,1]*P2inv - R[w,1,1]*Pinv

    P_sum = P[w].sum()

    Cinv_sum = Cinv_all.sum(axis=0)

    QbyP[:,0] = Q[w,0]*Pinv
    QbyP[:,1] = Q[w,1]*Pinv
    Q_sum = QbyP.sum(axis=0)

    return P_sum, Q_sum, Cinv_sum


def get_shear_pqr(P,Q,R, get_sums=False):
    """
    Extract a shear estimate from the p,q,r values from
    Bernstein & Armstrong

    parameters
    ----------
    P: array[nobj]
        Prior times jacobian
    Q: array[nobj,2]
        gradient of P with respect to shear
    R: array[nobj,2,2]
        gradient of gradient

    output
    ------
    [g1,g2]: array

    notes
    -----
    If done on a single object, the operations would look simpler

    QQ = numpy.outer(Q,Q)
    Cinv = QQ/P**2 - R/P
    C = numpy.linalg.inv(Cinv)
    g1g2 = numpy.dot(C,Q/P)

    """

    P_sum, Q_sum, Cinv_sum = get_shear_pqr_sums(P,Q,R)

    C = numpy.linalg.inv(Cinv_sum)
    g1g2 = numpy.dot(C,Q_sum)

    if get_sums:
        return g1g2, C, Q_sum, Cinv_sum
    else:
        return g1g2, C

def pqr_jackknife(P, Q, R,
                  get_sums=False,
                  verbose=False,
                  show=False,
                  eps=None,
                  png=None):
    """
    Get the shear covariance matrix using jackknife resampling.

    The trick is that this must be done in pairs
    """

    if verbose:
        import progressbar
        pg=progressbar.ProgressBar(width=70)

    ntot = P.size
    if ( (ntot % 2) != 0 ):
        raise  ValueError("expected factor of two, got %d" % ntot)

    npair = ntot/2

    P_sum, Q_sum, Cinv_sum = get_shear_pqr_sums(P,Q,R)
    C = numpy.linalg.inv(Cinv_sum)
    shear = numpy.dot(C,Q_sum)

    shears = numpy.zeros( (npair, 2) )
    for i in xrange(npair):

        if verbose:
            frac=float(i+1)/npair
            pg.update(frac=frac)

        ii = i*2

        Ptmp = P[ii:ii+2]
        Qtmp = Q[ii:ii+2,:]
        Rtmp = R[ii:ii+2,:,:]

        P_sum, Q_sum_tmp, Cinv_sum_tmp = get_shear_pqr_sums(Ptmp,Qtmp,Rtmp)
        
        Q_sum_tmp    = Q_sum - Q_sum_tmp
        Cinv_sum_tmp = Cinv_sum - Cinv_sum_tmp

        Ctmp = numpy.linalg.inv(Cinv_sum_tmp)
        shear_tmp = numpy.dot(C,Q_sum_tmp)

        shears[i, :] = shear_tmp

    shear_cov = numpy.zeros( (2,2) )
    fac = (npair-1)/float(npair)

    shear = shears.mean(axis=0)

    shear_cov[0,0] = fac*( ((shear[0]-shears[:,0])**2).sum() )
    shear_cov[0,1] = fac*( ((shear[0]-shears[:,0]) * (shear[1]-shears[:,1])).sum() )
    shear_cov[1,0] = shear_cov[0,1]
    shear_cov[1,1] = fac*( ((shear[1]-shears[:,1])**2).sum() )

    if show or eps or png:
        _plot_shears(shears, show=show, eps=eps, png=png)

    if get_sums:
        return shear, shear_cov, Q_sum, Cinv_sum
    else:
        return shear, shear_cov

def pqr_bootstrap(P, Q, R, nsamples, verbose=False, show=False, eps=None, png=None):
    """
    Get the shear covariance matrix using boot resampling.

    The trick is that this must be done in pairs
    """

    if verbose:
        import progressbar
        pg=progressbar.ProgressBar(width=70)

    ntot = P.size
    if ( (ntot % 2) != 0 ):
        raise  ValueError("expected factor of two, got %d" % ntot)

    npair = ntot/2

    Pboot = P.copy()
    Qboot = Q.copy()
    Rboot = R.copy()

    rind1 = numpy.zeros(npair, dtype='i8')
    rind2 = numpy.zeros(npair, dtype='i8')
    rind = numpy.zeros(ntot, dtype='i8')

    shears = numpy.zeros( (nsamples, 2) )

    for i in xrange(nsamples):
        if verbose:
            frac=float(i+1)/nsamples
            pg.update(frac=frac)

        # first of the pair
        rind1[:] = 2*numpy.random.randint(low=0,high=npair,size=npair)
        # second of the pair
        rind2[:] = rind1[:]+1

        rind[0:npair] = rind1
        rind[npair:]  = rind2

        Pboot[:] = P[rind]
        Qboot[:,:] = Q[rind,:]
        Rboot[:,:,:] = R[rind,:,:]

        sh, C_not_used =  get_shear_pqr(Pboot, Qboot, Rboot)

        shears[i, :] = sh

    shear = shears.mean(axis=0)

    shear_cov = numpy.zeros( (2,2) )

    shear_cov[0,0] = ( (shears[:,0]-shear[0])**2 ).sum()/(nsamples-1)
    shear_cov[0,1] = ( (shears[:,0]-shear[0])*(shears[:,1]-shear[1]) ).sum()/(nsamples-1)
    shear_cov[1,0] = shear_cov[0,1]
    shear_cov[1,1] = ( (shears[:,1]-shear[1])**2 ).sum()/(nsamples-1)

    if show or eps or png:
        _plot_shears(shears, show=show, eps=eps, png=png)

    return shear, shear_cov



def _plot_shears(shears, show=True, eps=None, png=None):
    import biggles
    tab=biggles.Table(2,1)
    std=shears.std(axis=0)

    plt1=eu.plotting.bhist(shears[:,0], binsize=0.2*std[0],
                           color='blue',show=False,
                           xlabel=r'$\gamma_1$')
    plt2=eu.plotting.bhist(shears[:,1], binsize=0.2*std[1],
                           color='red',show=False,
                           xlabel=r'$\gamma_2$')
    tab[0,0] = plt1
    tab[1,0] = plt2

    if png is not None:
        print png
        tab.write_img(800,800,png)
    if eps is not None:
        print eps
        tab.write_eps(eps)


    if show:
        tab.show()


