"""
    %prog gmix_run scinv_sample

Description:

    Create mean inverse critical density as a function of lens redshift
    for regauss objects with criteria for the input sample

"""
from __future__ import print_function

import sys
import numpy
import gmix_sdss
import lensing
import zphot
import esutil as eu

from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])
parser.add_option('--clobber',action='store_true',
                  help="clobber existing column")
parser.add_option('--chunksize',default=100000,
                  help="chunksize to process")


def add_scinv(gmix_run, scinv_sample, 
              clobber=False,
              chunksize=100000):
    """

    Add a column to the db that is the mean scinv as a function of lens
    redshift.

    scinv are created using the *corrected* p(z).

    The dzl,zlmin,zlmax,zlvals will be in the meta data

    Procedure
    ---------
    - get the correction factor N(z)/sumpofz(z)
        calls zphot.weighting.pofz_correction(pzrun)

    - for each p(z) multiply by that factor
    - generate mean inverse critical density.


    """

    scinv_conf = lensing.files.read_config('scinv', scinv_sample)
    cosmo = lensing.files.read_config('cosmo',scinv_conf['cosmo_sample'])
    pzrun=scinv_conf['pzrun']


    cols = gmix_sdss.collate.open_columns(gmix_run)

    # the column which will hold the inverse critical density.
    # depending on keywords, we might want to raise an error
    colname = lensing.scat.scinv_colname(scinv_sample)
    print("Writing to column:\n",colname)
    if colname in cols:
        if not clobber:
            raise ValueError("Column already exists")
        else:
            print("  removing column")
            cols[colname].delete()

    # get the matches
    zphot_matchname = lensing.scat.zphot_matchname(pzrun)

    if zphot_matchname not in cols:
        raise ValueError("zphot match column not "
                         "found: '%s'" % zphot_matchname)

    #
    # correction factor to apply to all p(z)
    #
    print("getting p(z) correction function\n")
    corrfile=zphot.weighting.pofz_correction_file(pzrun)
    corrstruct = eu.io.read(corrfile)
    zs = (corrstruct['zmax']+corrstruct['zmin'])/2.
    corr = corrstruct['corr']


    print("")
    scalc = lensing.sigmacrit.ScinvCalculator(scinv_conf['zlmin'],
                                              scinv_conf['zlmax'],
                                              scinv_conf['nzl'], 
                                              zs[0],
                                              zs[-1],
                                              H0=cosmo['H0'],
                                              omega_m=cosmo['omega_m'])


    zlvals = scalc.zlvals

    meta={'scinv_sample':scinv_sample,
          'zlmin':scinv_conf['zlmin'],
          'zlmax':scinv_conf['zlmax'],
          'nzl':scinv_conf['nzl'],
          'zlvals':zlvals}

    print("opening corresponding p(z) columns: '%s'\n" % pzrun)
    pzcols = zphot.weighting.open_pofz_columns(pzrun)

    # work on chunks
    print("using chunksize:",chunksize)
    ntot = cols['run'].size
    nchunk = ntot/chunksize
    if (ntot % chunksize) > 0:
        nchunk += 1

    print("nchunk:",nchunk)
    for i in xrange(nchunk):
        imin = i*chunksize
        imax = (i+1)*chunksize
        print("  %d %d (%d)" % (imin,imax,ntot))
        match = cols[zphot_matchname][imin:imax]
        
        nrows = match.size
        scinv = numpy.zeros((nrows, zlvals.size), dtype='f8') -9999.

        w,=numpy.where(match >= 0)
        print("    ",w.size,"with good matches")
        if w.size > 0:
            # read the p(z) for the matches
            pofz = pzcols['pofz'][match[w]]

            pofz = pofz[:]*corr[:]

            for j in xrange(w.size):
                scinv[w[j], :] = scalc.calc_mean_scinv(pofz[j])
        # metadata only gets written once
        cols.write_column(colname, scinv, meta=meta)

def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]
    scinv_sample=args[1]

    chunksize=int(options.chunksize)

    add_scinv(gmix_run, scinv_sample,
              clobber=options.clobber,
              chunksize=chunksize)

main()
