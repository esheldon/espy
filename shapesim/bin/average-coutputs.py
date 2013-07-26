"""
    %prog run
"""

import sys
import os
import numpy

import esutil as eu
from shapesim import shapesim
import lensing

import time

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--skip',default=None,
                  help="is2n elements to skip")

def pqr_jackknife(P, Q, R, verbose=False):
    """
    Get the shear covariance matrix using bootstrap resampling.
    We use this for errors when doing ring tests

    The trick is that this must be done in pairs

    This is "unweighted", although there are built-in weights
    """

    if verbose:
        import progressbar
        pg=progressbar.ProgressBar(width=70)

    ntot = P.size
    if ( (ntot % 2) != 0 ):
        raise  ValueError("expected factor of two, got %d" % ntot)

    npair = ntot/2

    Q_sum, Cinv_sum = lensing.shear.get_shear_pqr_sums(P,Q,R)
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

        Q_sum_tmp, Cinv_sum_tmp = lensing.shear.get_shear_pqr_sums(Ptmp,Qtmp,Rtmp)
        
        Q_sum_tmp    = Q_sum - Q_sum_tmp
        Cinv_sum_tmp = Cinv_sum - Cinv_sum_tmp

        Ctmp = numpy.linalg.inv(Cinv_sum_tmp)
        shear_tmp = numpy.dot(C,Q_sum_tmp)

        shears[i, :] = shear_tmp

    shear_cov = numpy.zeros( (2,2) )
    fac = (npair-1)/float(npair)

    shear_cov[0,0] = fac*( ((shear[0]-shears[:,0])**2).sum() )
    shear_cov[0,1] = fac*( ((shear[0]-shears[:,0]) * (shear[1]-shears[:,1])).sum() )
    shear_cov[1,0] = shear_cov[0,1]
    shear_cov[1,1] = fac*( ((shear[1]-shears[:,1])**2).sum() )
    return shear, shear_cov, Q_sum, Cinv_sum

def pqr_bootstrap(P, Q, R, shear, nsamples=100, verbose=False):
    """
    Get the shear covariance matrix using bootstrap resampling

    The trick is that this must be done in pairs
    """

    ntot = P.size
    if ( (ntot % 2) != 0 ):
        raise  ValueError("expected factor of two, got %d" % ntot)

    npair = ntot/2

    if verbose:
        import progressbar
        pg=progressbar.ProgressBar(width=70)
    shears = numpy.zeros( (nsamples, 2) )

    Pboot = P.copy()
    Qboot = Q.copy()
    Rboot = R.copy()

    rind1 = numpy.zeros(npair, dtype='i8')
    rind2 = numpy.zeros(npair, dtype='i8')
    rind = numpy.zeros(ntot, dtype='i8')

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

        sh, C =  lensing.shear.get_shear_pqr(Pboot, Qboot, Rboot)

        shears[i, :] = sh

    shear_cov = numpy.zeros( (2,2) )

    shear_cov[0,0] = ( (shears[:,0]-shear[0])**2 ).sum()/(nsamples-1)
    shear_cov[0,1] = ( (shears[:,0]-shear[0])*(shears[:,1]-shear[1]) ).sum()/(nsamples-1)
    shear_cov[1,0] = shear_cov[0,1]
    shear_cov[1,1] = ( (shears[:,1]-shear[1])**2 ).sum()/(nsamples-1)

    return shear_cov

def get_averaged(data, s2n_matched):
    dt =data.dtype.descr
    dt += [('s2n_matched','f8'),

           ('nsum','i8'),

           ('Q_sum','f8',2),
           ('Cinv_sum','f8',(2,2)),
           ('shear','f8',2),
           ('shear_cov','f8',(2,2)),
           ('shear_cov_inv_sum','f8',(2,2)),
    
    ]


    print 'jackknifing'
    t0=time.time()
    shear, shear_cov, Q_sum, Cinv_sum = pqr_jackknife(data['P'],data['Q'],data['R'],
                                                      verbose=True)
    shear_cov_inv = numpy.linalg.inv(shear_cov)

    print 'time:',time.time()-t0
    d=numpy.zeros(1, dtype=dt)

    d['s2n_matched'] = s2n_matched
    d['nsum'] = data.size

    d['Q_sum'][0] = Q_sum
    d['Cinv_sum'][0] = Cinv_sum

    d['shear'][0] = shear
    d['shear_cov'][0] = shear_cov
    d['shear_cov_inv_sum'][0] = shear_cov_inv

    print 'shear1: %.16g +/- %.16g' % (shear[0],numpy.sqrt(shear_cov[0,0]))

    return d

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]

    if options.skip is None:
        skip=[]
    else:
        skip = [int(v) for v in options.skip.split(',')]

    c = shapesim.read_config(run)
    s2n_vals    = c['s2n_vals']

    n_s2n = len(s2n_vals)

    dlist=[]
    for is2n in reversed(xrange(n_s2n)):

    #for is2n in xrange(n_s2n):
        s2n_matched = s2n_vals[is2n]
        fname=shapesim.get_output_url(run, 0, is2n)
        print fname
        data=eu.io.read(fname)

        d = get_averaged(data, s2n_matched)
        dlist.append(d)

    output = eu.numpy_util.combine_arrlist(dlist)
    out_fname=shapesim.get_averaged_url(run, 0)
    print 'writing:',out_fname
    eu.io.write(out_fname, output, clobber=True)

main()
