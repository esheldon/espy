"""
Tools for selecting a good set of galaxies
"""

# this package
import des

import os
from sys import stdout
import deswl
import columns
import numpy
from numpy import where

import esutil as eu
from esutil.ostools import path_join

class SizeMagSelector(dict):
    def __init__(self, run):
        self._run=run
        coldir = deswl.files.coldir(run)
        if not os.path.exists(coldir):
            raise ValueError("No such coldir: %s" % coldir)
        self['cols'] = columns.Columns(coldir)

    @property
    def cols(self):
        return self._cols

    def cache(self, force=False):

        do_wgood = False
        names = \
            ['shear_flags','size_flags','sigma0','interp_psf_sigma','imag']
        for n in names:
            if (n not in self) or force:
                if n in self:
                    del self[n]

                stdout.write("Cacheing %s\n" % n)
                self[n] = self['cols'][n][:]

                if n in ['shear_flags','size_flags'] or force:
                    do_wgood = True

        if 'wgood' not in self or do_wgood:
            if 'wgood' in self:
                del self['wgood']

            stdout.write("Cacheing good indices\n")
            wgood, = \
                where((self['shear_flags'] == 0) & (self['size_flags'] == 0))
            self['wgood'] = wgood

       
    def dc4_region_select(self, ra, dec, region):
        return des.dc4.region_select(ra,dec,region)

    def plot(self, flip=True):

        collate_dir= deswl.files.wlse_collated_dir(self._run)
        plotdir=path_join(collate_dir,'plots')

        band='i'
        eps_file= path_join(plotdir, '%s-%s-sizemag.eps' % (self._run,band))

        import biggles

        self.cache()
        
        size_mag_all = self.get_size_mag_hist2d()
        dens = image_norm(asinh_scale(size_mag_all['hist']), 
                          reverse=True )
        size_mag_allp = biggles.Density(dens, size_mag_all['ranges_reverse'])

        # for a contour plot
        size_mag_allcp = biggles.Contours(size_mag_all['hist'], 
                                          size_mag_all['ycenter'],
                                          size_mag_all['xcenter'],
                                          color='white')
        size_mag_allcp.levels=15


        size_mag = self.get_size_mag_hist2d(all=False)
        dens = image_norm(asinh_scale(size_mag['hist']), 
                          reverse=True )
        size_magp = biggles.Density(dens, size_mag['ranges_reverse'])

        # for a contour plot
        size_magcp = biggles.Contours(size_mag['hist'], 
                                      size_mag['ycenter'],
                                      size_mag['xcenter'],
                                      color='white')

        size_magcp.levels=15



        sizeratio_mag_all = self.get_sizeratio_mag_hist2d()
        dens = image_norm(asinh_scale(sizeratio_mag_all['hist']), 
                          reverse=True )
        sizeratio_mag_allp = biggles.Density(dens, sizeratio_mag_all['ranges_reverse'])

        # for a contour plot
        sizeratio_mag_allcp= biggles.Contours(sizeratio_mag_all['hist'], 
                                              sizeratio_mag_all['ycenter'],
                                              sizeratio_mag_all['xcenter'],
                                              color='white')

        sizeratio_mag_allcp.levels=50


        sizeratio_mag = self.get_sizeratio_mag_hist2d(all=False)
        dens = image_norm(asinh_scale(sizeratio_mag['hist']), 
                          reverse=True )
        sizeratio_magp = biggles.Density(dens, sizeratio_mag['ranges_reverse'])

        # for a contour plot
        sizeratio_magcp= biggles.Contours(sizeratio_mag['hist'], 
                                          sizeratio_mag['ycenter'],
                                          sizeratio_mag['xcenter'],
                                          color='white')
        sizeratio_magcp.levels=15


        stdout.write("making size ratio mag plot\n")
        t = biggles.Table(2,2)

        p00 = biggles.FramedPlot()
        p00.add( size_mag_allp )
        p00.add( size_mag_allcp )

        p00.xlabel = 'i-mag (uncalibrated)'
        p00.ylabel = 'FWHM'
        p00.title = 'size_flags==0'
        t[0,0] = p00

        p01 = biggles.FramedPlot()
        p01.add( size_magp )
        p01.add( size_magcp )
        p01.xlabel = 'i-mag (uncalibrated)'
        p01.ylabel = 'FWHM'
        p01.title = 'size_flags==0 && shear_flags==0'
        t[0,1] = p01

        p10 = biggles.FramedPlot()
        p10.add( sizeratio_mag_allp )
        p10.add( sizeratio_mag_allcp )
        p10.xlabel = 'i-mag (uncalibrated)'
        p10.ylabel = r'$\sigma_{PSF}/\sigma_{Object}'
        p10.title = 'size_flags==0'
        t[1,0] = p10

        p11 = biggles.FramedPlot()
        p11.add( sizeratio_magp )
        p11.add( sizeratio_magcp )
        p11.xlabel = 'i-mag (uncalibrated)'
        p11.ylabel = r'$\sigma_{PSF}/\sigma_{Object}'
        p11.title = 'size_flags==0 && shear_flags==0'
        t[1,1] = p11


        #t.show()

        stdout.write("Writing eps file: %s\n" % eps_file)
        t.write_eps(eps_file)
        stdout.write('Converting to png\n')
        eu.ostools.exec_process('converter -d 90 %s' % eps_file,verbose=True)

    def get_size_mag_hist2d(self, all=True):
        if all:
            stdout.write("Creating 'all' size mag 2d histogram\n")
            wgood = where( self['size_flags'] == 0 )
        else:
            stdout.write("Creating 'restricted' size mag 2d histogram\n")
            wgood = self['wgood']

        imag = self['imag'][wgood]
        fwhm = self['sigma0'][wgood]*2.35
        imagmin = 12.0
        imagmax = 17.0
        sizemin = 0.0
        sizemax = 3.5

        # note flip
        hist_dict = eu.stat.histogram2d(fwhm, imag, 
                                        nx=100, ny=100,
                                        xmin=sizemin,xmax=sizemax,
                                        ymin=imagmin,ymax=imagmax,
                                        more=True)

        return hist_dict

    def get_sizeratio_mag_hist2d(self, all=True):
        if all:
            stdout.write("Creating 'all' sizeratio mag 2d histogram\n")
            wgood = where( self['size_flags'] == 0 )
        else:
            stdout.write("Creating 'restricted' sizeratio "
                         "mag 2d histogram\n")
            wgood = self['wgood']
        imag = self['imag'][wgood]
        size_ratio = \
            self['interp_psf_sigma'][wgood]/self['sigma0'][wgood]

        imagmin = 12.0
        imagmax = 17.0
        sizemin = 0.0
        sizemax = 1.2

        # note flip
        hist_dict = eu.stat.histogram2d(size_ratio, imag, 
                                        nx=100, ny=100,
                                        xmin=sizemin,xmax=sizemax,
                                        ymin=imagmin,ymax=imagmax,
                                        more=True)

        return hist_dict


def asinh_scale(image, alpha=0.02, nonlinearity=8.0):
    image_out=numpy.array(image, dtype='f8', copy=True)

    image_out[:] = \
        numpy.arcsinh( alpha*nonlinearity*image )/nonlinearity

    return image_out

def image_norm(image, reverse=False):
    image_out=numpy.array(image, dtype='f8', copy=True)
    image_out /= image_out.max()

    if reverse:
        image_out = 1.0 - image_out

    return image_out
