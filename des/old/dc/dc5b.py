import numpy
from sys import stdout
import os
import esutil as eu
import deswl

def check_astrometry(band):
    flist = deswl.files.collated_redfiles_read('dc5b', band)['flist']
    #flist = flist[0:500]

    ntot = len(flist)

    nbad=0
    i=0
    badstd = 0.05
    for finfo in flist:
        if (i+1) % 100 == 0 or i == 0:
            print '%s/%s (%0.2f%%)' % (i+1, ntot, 100*float(i+1)/ntot)
        imf = os.path.expandvars(finfo['imfile'])
        catf = os.path.expandvars(finfo['catfile'])
        #print imf

        cat = eu.io.read(catf, ext=2, lower=True)
        h = eu.pyfits.getheader(imf, ext=1)
        wcs = eu.wcsutil.WCS(h)

        ra,dec = wcs.image2sky(cat['x_image'], cat['y_image'])

        dist = eu.coords.sphdist(cat['alphamodel_j2000'], cat['deltamodel_j2000'], ra, dec)

        std_arcsec = 3600.0*dist.std()

        if std_arcsec > badstd:
            #print 'For image:',imf
            #print '    large std: %s arcsec' % std_arcsec
            nbad += 1

        i+=1


    print '\nBad rate: %s/%s (%0.2f%%)' % (nbad,ntot,nbad/float(ntot))

def compare_one():
    wcs = eu.wcsutil.WCS(eu.pyfits.getheader('/global/data/DES/red/20100423103614_20100324/red/decam--23--42-i-5/decam--23--42-i-5_01.fits.fz',ext=1))
    cat = eu.io.read('/global/data/DES/red/20100423103614_20100324/red/decam--23--42-i-5/decam--23--42-i-5_01_cat.fits',lower=True,ext=2)

    ra,dec = wcs.image2sky(cat['x_image'], cat['y_image'])
    dist = eu.coords.sphdist(cat['alphamodel_j2000'], cat['deltamodel_j2000'], ra, dec, 
                             units=['deg','deg'])
    dist *= 3600.0
    
    fname='~/tmp/compare-one.dat'
    print 'writing to file:',fname
    eu.misc.colprint(cat['x_image'], cat['y_image'], cat['alphamodel_j2000'], cat['deltamodel_j2000'], ra, dec, 
                     (cat['alphamodel_j2000']-ra)*3600, (cat['deltamodel_j2000']-dec)*3600, dist, 
                    names=['x_image','y_image','alphamodel_j2000', 'deltamodel_j2000', 'ra_wcs','dec_wcs', 
                           'radiff(arcsec)','decdiff(arcsec)','sphdist(arcsec)'], file=fname, 
                    format='%0.9f')

