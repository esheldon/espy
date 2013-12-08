"""
    %prog run
"""

import sys
import os
import numpy
from numpy import sqrt

import esutil as eu
from shapesim import shapesim
import lensing

import time

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--skip',default=None,
                  help="is2n elements to skip")
parser.add_option('--verbose',action='store_true',
                  help="show progress")
parser.add_option('--gerror',action='store_true',
                  help="Use scatter in g to get errors")
parser.add_option('--show',action='store_true',
                  help="plot histogram of jackknife shears")


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

def get_averaged_gerror(data, s2n_matched, verbose=False):
    dt =data.dtype.descr
    dt += [('s2n_matched','f8'),

           ('nsum','i8'),

           ('Q_sum','f8',2),
           ('Cinv_sum','f8',(2,2)),
           ('shear','f8',2),
           ('shear_cov','f8',(2,2)),
           ('shear_cov_inv_sum','f8',(2,2)),
    
    ]

    if 'g' in data.dtype.names:
        do_lensfit=True
        dt += [('g_sum','f8',2),
               ('gsens_sum','f8',2),
               ('shear_lensfit','f8',2),
               ('shear_lensfit_cov','f8',(2,2)),
               ('shear_lensfit_cov_inv_sum','f8',(2,2)),
              ]

    else:
        do_lensfit=False

    do_simple=False
    if data['pars'].shape[1] == 6:
        do_simple=True
        dt += [('flux_sum','f8'),
               ('flux_err2invsum','f8'),
               ('flux_s2n_sum','f8'),
               ('flux','f8'),
               ('flux_err','f8'),
               ('flux_s2n','f8'),

               ('T_sum','f8'),
               ('T_err2invsum','f8'),
               ('T_s2n_sum','f8'),
               ('T','f8'),
               ('T_err','f8'),
               ('T_s2n','f8')]


    Q_sum, Cinv_sum = lensing.pqr.get_shear_pqr_sums(data['P'],
                                                     data['Q'],
                                                     data['R'])
    C = numpy.linalg.inv(Cinv_sum)
    shear = numpy.dot(C,Q_sum)

    cov = data['pcov'][:,2:2+2, 2:2+2]
    cov_inv = cov.copy()

    det = cov[:,0,0]*cov[:,1,1] - cov[:,0,1]*cov[:,1,0]
    cov_inv[:,0,0] = cov[:,1,1]
    cov_inv[:,1,1] = cov[:,0,0]
    cov_inv[:,0,1] = - cov[:,0,1]
    cov_inv[:,1,0] = - cov[:,1,0]

    idet = 1.0/det
    cov_inv[:,0,0] *= idet
    cov_inv[:,0,1] *= idet
    cov_inv[:,1,0] *= idet
    cov_inv[:,1,1] *= idet

    shear_cov_inv = cov_inv.sum(axis=0)
    shear_cov = numpy.linalg.inv(shear_cov_inv)


    d=numpy.zeros(1, dtype=dt)


    d['s2n_matched'] = s2n_matched
    d['nsum'] = data.size

    d['Q_sum'][0] = Q_sum
    d['Cinv_sum'][0] = Cinv_sum

    d['shear'][0] = shear

    d['shear_cov'][0,:,:] = shear_cov
    d['shear_cov_inv_sum'][0,:,:] = shear_cov_inv

    sherr=numpy.sqrt(shear_cov[0,0])
    print 'shear1:      %.16g +/- %.16g' % (shear[0],sherr)

    if do_simple:
        flux=data['pars'][:,5]
        flux_var=data['pcov'][:,5,5]
        flux_err=sqrt(flux_var)

        d['flux_sum'] = flux.sum()
        d['flux_err2invsum'] = (1.0/flux_var).sum()

        d['flux'] = d['flux_sum']/d['nsum']
        d['flux_err'] = sqrt(1.0/d['flux_err2invsum'])

        flux_s2n_vals = flux/flux_err
        d['flux_s2n_sum'] = flux_s2n_vals.sum()
        d['flux_s2n'] = d['flux_s2n_sum']/d['nsum']


        T=data['pars'][:,4]
        T_var=data['pcov'][:,4,4]
        T_err=sqrt(T_var)

        d['T_sum'] = T.sum()
        d['T_err2invsum'] = (1.0/T_var).sum()

        d['T'] = d['T_sum']/d['nsum']
        d['T_err'] = sqrt(1.0/d['T_err2invsum'])

        T_s2n_vals = T/T_err
        d['T_s2n_sum'] = T_s2n_vals.sum()
        d['T_s2n'] = d['T_s2n_sum']/d['nsum']


    if do_lensfit:
        gmean = data['g'].mean(axis=0)
        gsens = data['gsens'].mean(axis=0)
        
        shear_cov[0,0] /= (gsens[0]*gsens[0])
        shear_cov[0,1] /= (gsens[0]*gsens[1])
        shear_cov[1,0] /= (gsens[0]*gsens[1])
        shear_cov[1,1] /= (gsens[1]*gsens[1])
        shear_cov_inv = numpy.linalg.inv(shear_cov)

        shear=gmean/gsens
        sherr=numpy.sqrt(shear_cov[0,0])

        d['g_sum'] = data['g'].sum(axis=0)
        d['gsens_sum'] = data['gsens'].sum(axis=0)

        d['shear_lensfit'] = shear
        d['shear_lensfit_cov'][0,:,:] = shear_cov
        d['shear_lensfit_cov_inv_sum'][0,:,:] = shear_cov_inv

        print 'shear1 (lensfit): %.16g +/- %.16g' % (shear[0],sherr)
    return d


def get_averaged_jackknife(data, s2n_matched, verbose=False, show=False, fname=None):
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
    shear, shear_jack, shear_cov, Q_sum, Cinv_sum = \
            lensing.pqr.pqr_jackknife(data['P'],data['Q'],data['R'],
                                      verbose=verbose,
                                      show=show,
                                      fname=fname)
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

    print 'shear1:      %.16g +/- %.16g' % (shear[0],numpy.sqrt(shear_cov[0,0]))

    return d

def get_image_file(plot_dir, fname):
    bname=os.path.basename(fname)
    bname += '-jackknife.png'

    return os.path.join(plot_dir, bname)

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
    #for is2n in reversed(xrange(n_s2n)):

    for is2n in xrange(n_s2n):
        if is2n in skip:
            continue

        s2n_matched = s2n_vals[is2n]
        fname=shapesim.get_output_url(run, 0, is2n)
        print fname
        data=eu.io.read(fname)

        if options.gerror:
            d = get_averaged_gerror(data, s2n_matched, verbose=options.verbose)
        else:
            plot_dir=shapesim.get_plot_dir(run)
            plot_name=get_image_file(plot_dir, fname)
            d = get_averaged_jackknife(data, s2n_matched,
                                       verbose=options.verbose,
                                       show=options.show,
                                       fname=plot_name)
        dlist.append(d)

    output = eu.numpy_util.combine_arrlist(dlist)
    out_fname=shapesim.get_averaged_url(run, 0)
    print 'writing:',out_fname
    eu.io.write(out_fname, output, clobber=True)

main()
