from __future__ import print_function
import os
import sys
import numpy
from numpy import log10,sqrt,linspace,where
import lensing
import esutil as eu
from esutil.numpy_util import where1

#
# collated outputs
#

def make_collated_lensred(run):
    """
    
    run
        collate_lensout_and_catalog

    end write out the file

    """
    data = collate_lensred_and_catalog(run)
    f = lensing.files.sample_file('lensred-collate', run)
    print(f)
    eu.io.write(f, data, verbose=True)

def collate_lensred_and_catalog(run):
    """

    Collate the reduced lensout (lens-by-lens) file created with

        make_reduced_lensout

    with the original input catalog.

    """

    conf = lensing.files.cascade_config(run)
    lsample = conf['lens_config']['sample']
    cat = lensing.files.read_original_catalog('lens',lsample)
    lout = lensing.files.lensred_read(run)

    # trim down to the ones we used
    print("Extracting by zindex")
    cat = cat[ lout['zindex'] ]

    # create a new struct with the lensing outputs at the end
    print('collating')
    data = eu.numpy_util.add_fields(cat, lout.dtype)
    eu.numpy_util.copy_fields(lout, data)

    return data

#
# Code to reduce the split lensums (we split the source
# catalog into chunks and run all lenses on each chunk)
#

def make_reduced_lensout(run):
    """

    reduce the lensout splits into a single, summed lensum (lensout) file.
    This is still lens-by-lens, but the sums from different source splits are
    added.

    """
    file = lensing.files.sample_file('lensred', run)
    print("Will combine into file:",file)
    print("Reading data\n")
    data = lensing.files.lensout_read(run=run)

    eu.io.write(file, data, verbose=True)


def add_lensums(l1, l2):
    """

    Add the sums from l2 to l1.  Note this is different than reduce_lensums,
    which produces averages

    The rows of l1 must correspond to those of l2; zindex must match

    """

    w=where1(l1['zindex'] != l2['zindex'])
    if w.size > 0:
        raise ValueError("zindex do not line up")
    for n in ['weight','npair','rsum','wsum','dsum','osum']:
        l1[n] += l2[n]




#
# Codes for combining the lensout "lensum" data and
# getting averages
#

def reduce_lensums(lout):
    """

    Reduce the lens-by-lens lensums by summing over
    all the individual sums and producing averages

    """

    nlens = lout.size
    nbin = lout['rsum'][0].size
    dt = lensred_dtype()

    comb=numpy.zeros(nbin,dtype=dt)

    comb['weight'] = lout['weight'].sum()

    for i in xrange(nbin):
        npair = lout['npair'][:,i].sum()
        rsum  = lout['rsum'][:,i].sum()


        wsum = lout['wsum'][:,i].sum()
        dsum = lout['dsum'][:,i].sum()
        osum = lout['osum'][:,i].sum()

        comb['npair'][i] = npair
        comb['rsum'][i] = rsum
        comb['wsum'][i] = wsum
        comb['dsum'][i] = dsum
        comb['osum'][i] = osum

        # averages
        comb['r'][i] = rsum/npair
        comb['dsig'][i] = dsum/wsum
        comb['osig'][i] = osum/wsum

        comb['dsigerr'][i] = numpy.sqrt(1.0/wsum)

    return comb

def lensred_dtype():
    dt=[('r','f8'),
        ('dsig','f8'),
        ('dsigerr','f8'),
        ('osig','f8'),
        ('npair','i8'),
        ('weight','f8'),
        ('rsum','f8'),
        ('wsum','f8'),
        ('dsum','f8'),
        ('osum','f8')]
    return numpy.dtype(dt)


