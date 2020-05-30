from __future__ import print_function
import cosmology

def relative_npix(ngal, da, pixel_scale, nobs):
    return ngal/da**2/pixel_scale**2 * nobs

c = cosmology.Cosmo(h=0.72)

des_area = 5000
lsst_area = 20000

des_med_z = 0.7
lsst_med_z = 1.5

des_Da = c.Da(0.0, des_med_z)
lsst_Da = c.Da(0.0, lsst_med_z)

des_raw_density = 19.5
des_raw_number = des_raw_density*3600.*des_area

lsst_raw_number = 10e9
lsst_raw_density = lsst_raw_number/lsst_area

des_pixel_scale = 0.27
lsst_pixel_scale = 0.2

des_nobs = 2*10
lsst_nobs = 2*230

des_rel_npix = relative_npix(des_raw_number, des_Da, des_pixel_scale, des_nobs)
lsst_rel_npix = relative_npix(lsst_raw_number, lsst_Da, lsst_pixel_scale, lsst_nobs)


npixratio = lsst_rel_npix/des_rel_npix

mess="""Relavent quantity is the number of pixels processed

    npixels = Ngalaxy * npixels/galaxy * Nobservations
    npixels = Ngalaxy * area(arcsec^2)/galaxy * (pixels/arcsec)^2 * Nobservations

Assume galaxies are roughly the same physical size for LSST and DES, then the
second term is proportional to the inverse angular diameter distance squared.
      
Morgan desired that we use numbers from the white papers

The DES numbers are from 
    http://www.darkenergysurvey.org/reports/proposal-standalone.pdf
The LSST numbers are from 
    Ivezic et al. astro-ph/0805.2366

                DES       LSST

<z>         {des_z}    {lsst_z}
Da (Mpc)    {des_da}    {lsst_da}
ngal/10^9   {des_ngal}    {lsst_ngal}
pixsize''   {des_pixsize}    {lsst_pixsize}
nobs (r+i)  {des_nobs}    {lsst_nobs}

thus the ratio of number of pixels is

LSST/DES = {lsst_ngal} * {des_da}^2 * {des_pixsize}^2 * {lsst_nobs}
           -------------------------------------------
           {des_ngal} * {lsst_da}^2 * {lsst_pixsize}^2 * {des_nobs}

    
         ~ {npixratio}

""".format(des_z='%7.1f' % des_med_z,
           lsst_z='%7.1f' % lsst_med_z,
           des_da='%7.2f' % des_Da,
           lsst_da='%7.2f' % lsst_Da,
           des_ngal='%7.2f' % (des_raw_number/1.e9,),
           lsst_ngal='%7.2f' % (lsst_raw_number/1.e9,),
           des_pixsize='%7.2f' % des_pixel_scale,
           lsst_pixsize='%7.2f' % lsst_pixel_scale,
           des_nobs='%7d' % des_nobs,
           lsst_nobs='%7d' % lsst_nobs,
           npixratio='%d' % npixratio)


print(mess)
