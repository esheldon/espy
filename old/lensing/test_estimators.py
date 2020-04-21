"""
    usage: %prog [options] n ntrial 

 Two simple estimators where we think we know the p(zsource) for each
 source galaxy.  We will try two methods:

   1) we assign each source a z drawn from the expected distribution and
       treat it as truth in the calculations.
   2) we use the p(z) to calculate the sigma_crit^{-1}(zlens) for each source
       and store that, interpolating as required.

 we will want to test a few different source distributions and lens 
 redshifts too.
"""

import os
import sys
import numpy
import lensing
import esutil as eu
from esutil import sfile

from time import time


from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-e","--extra",dest="extra", default="",
                  help="extra string for file name")


def plot_results(files, type, plotfile=None):
    hdr=eu.sfile.read_header(files[0])
    dtrue = hdr['deltasig_true']

    t=eu.io.read(files, combine=True)

    nbin=20

    plt=eu.plotting.setuplot()
    plt.clf()

    if type == 'sample':
        deltasig = t['deltasig_sample']
        mean_err = t['deltasig_sample_err'].mean()
        label="sampled P(z)"
    else:
        deltasig = t['deltasig_meanscinv']
        mean_err = t['deltasig_meanscinv_err'].mean()
        label="mean P(z)"

    hist, xhist, patches = plt.hist(deltasig, nbin, 
                                    normed=True, align='center', 
                                    label='realizations')

    # overplot the expected gaussian
    xvals = numpy.arange(0.0,20.0,0.1)
    norm=1.0/numpy.sqrt(2.0*numpy.pi*mean_err**2)
    gauss=norm*numpy.exp(-0.5*(xvals-dtrue)**2/mean_err**2 )

    l=plt.plot(xvals,gauss,label='expected result')
    leg=plt.legend()
    leg.draw_frame(False)

    plt.xlim([5,15])
    plt.xlabel('$\Delta\Sigma$')
    #plt.ylabel('P\left( \Delta\Sigma \right)')
    plt.ylabel('$P(\Delta\Sigma)$')

    plt.annotate( label, (6,0.35) )
    ax = plt.axes()


    eu.plotting.set_minor_ticks(ax)
    if plotfile is not None:
        plt.savefig(plotfile)
    else:
        plt.show()


def generate_pbs(n, ntrial, npbs, start=0):
    pbsdir = '~/pbs/estimators'
    pbsdir=os.path.expanduser(pbsdir)
    if not os.path.exists(pbsdir):
        os.mkdir(pbsdir) 

    for i in range(start,npbs+start):
        extra='%02i' % i
        
        name="sample-meanscinv-%s-%s-%s" % (n,ntrial,extra)
        fname=os.path.join(pbsdir,name+'.pbs')
        logf=os.path.join(pbsdir,name+'.pbs.log')
        shortname="est-%s" % (extra,)

        pbs="""

#!/bin/bash
#PBS -S /bin/bash
#PBS -N %s
#PBS -j oe
#PBS -o %s.pbslog
#PBS -m a
#PBS -V
#PBS -r n
#PBS -W umask=0022
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1

source /global/data/products/eups/bin/setups.sh

setup python
setup esutil -r ~/exports/esutil-work

logf="%s"
python ~/python/lensing/test_estimators.py -e %s %s %s &> "$logf"
        """ % (shortname,name,logf,extra,n,ntrial)

        print fname
        fobj=open(fname,'w')
        fobj.write(pbs)
        fobj.close()

def generate_lognormal_zsource(m, sigma, n):
    """
    relationshnip between mean,sigma and the actuall mean and stdev of the
    resulting distribution:

    mu = log(m)
    mean = exp(mu + 0.5*sigma**2)
    sdev = sqrt( (exp(sigma**2) - 1.0)*exp(2*mu + sigma**2) )

    example:
        m=0.5 (mu=-0.69)
        sigma=0.5
    
        mean = 0.57
        sdev = 0.30


    """
    mu = numpy.log( m )
    s = numpy.random.lognormal(mu, sigma, n)
    return s

def shear_noise():
    # keep error down a bit
    return numpy.sqrt( 0.16**2 )

def generate_shears(shear_true):
    print 'generating shears'

    noise = shear_noise()
    shears = shear_true + noise*numpy.random.normal(size=shear_true.size)

    return shears

def estimate_deltasigma_sample(shear, shear_err, scinv):
    # note I'm cancelling the top sigma_crit
    eweight = 1.0/shear_err
    eweight = eweight**2

    numer_vals = shear*eweight*scinv
    denom_vals = eweight*scinv**2
    denom = denom_vals.sum()
    numer = numer_vals.sum()
    deltasig = numer/denom
    deltasig_err = numpy.sqrt( 1.0/denom )

    return deltasig, deltasig_err

def main(n, ntrial, extra=''):
    
    # parameters of source redshift distribution and lens
    m=0.5
    sigma = 0.5
    zlens = 0.2
    deltasig = 10.0

    if extra != '':
        if extra[0] == '-':
            extra = extra[1:]

    outfile = '~/oh/estimators/samp-meanscinv-%s-%s-%s.rec' % (n, ntrial, extra)
    print 'writing to file:',outfile
    sf = sfile.Open(outfile,'w+', delim=',')

    hdr = {'deltasig_true':deltasig, 
           'zlens':zlens,
           'nsource': n, 
           'm':m,
           'sigma':sigma,
           'ntrial':ntrial}

    dtype=[('deltasig_sample','f4'),('deltasig_sample_err','f4'),
           ('deltasig_meanscinv','f4'),('deltasig_meanscinv_err','f4')]

    # append as we go in case of a crash or I  kill it
    tm0 = time()
    for i in range(ntrial):

        tm0_i = time()
        st = numpy.zeros(1, dtype=dtype)


        print 'generating zs and zphot'
        ztrue = generate_lognormal_zsource(m, sigma, n)
        zphot = generate_lognormal_zsource(m, sigma, n)

        print 'generating scinv and scinv_zphot'
        scinv_true = lensing.tools.sigmacritinv(zlens, ztrue)
        scinv_zphot = lensing.tools.sigmacritinv(zlens, zphot)

        # for this will will just plug in the mean expectation value
        # for the inverse critical density
        scinv_zphot_mean = numpy.zeros(n)
        scinv_zphot_mean[:] = scinv_zphot.mean()


        shear_true = deltasig*scinv_true
        shears_noisy = generate_shears(shear_true)

        shear_err = numpy.zeros(n)
        shear_err[:] = shear_noise()

        # first sampling method

        print 'sample: estimating delta sigma'
        deltasig_est, deltasig_err = \
                estimate_deltasigma_sample(shears_noisy, 
                                           shear_err, 
                                           scinv_zphot)
        print 'sample: ',deltasig,deltasig_est, '+/-', deltasig_err
        st['deltasig_sample'][0] = deltasig_est
        st['deltasig_sample_err'][0] = deltasig_err


        # now using mean scinv
        print 'meanscinv: estimating delta sigma'
        deltasig_est, deltasig_err = \
                estimate_deltasigma_sample(shears_noisy, 
                                           shear_err, 
                                           scinv_zphot_mean)
        print 'meanscinv: ',deltasig,deltasig_est, '+/-', deltasig_err
        st['deltasig_meanscinv'][0] = deltasig_est
        st['deltasig_meanscinv_err'][0] = deltasig_err


        sf.write(st, header=hdr)
        sf.flush()
        sys.stdout.flush()

        tm1_i = time()
        eu.misc.ptime(tm1_i-tm0_i)

    tm1 = time()
    eu.misc.ptime(tm1-tm0)

if __name__=='__main__':
    options, args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    n=long(args[0])
    ntrial=long(args[1])

    extra=options.extra

    main(n,ntrial,extra=extra)
