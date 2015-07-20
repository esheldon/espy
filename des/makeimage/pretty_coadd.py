from __future__ import print_function

import sys
import os
import glob
import esutil as eu
from esutil import wcsutil
import numpy
from numpy import array, zeros, flipud, where, sqrt

import des
import fitsio

NOMINAL_EXPTIME=900.0

NONLINEAR=.12
DEFAULT_RELEASE='y1a1'

def make_coadd_jpg(run, rebin=None, release=DEFAULT_RELEASE):
    """
    make a color jpeg for the specified run

    parameters
    ----------
    run: string
        DES coadd run
    rebin: int, optional
        Amount to rebin image
    release: string
        Release; scaling differs for different releases
        currently y1a1 and sva1
    """
    cj=CoaddJPGFiles(run, rebin=rebin)

    make_jpg(cj['gfile'],
             cj['rfile'],
             cj['ifile'],
             cj['jpg_file'],
             rebin=rebin,
             release=release)

def make_jpg(gfile, rfile, ifile, outfile,
             rebin=None,
             release=DEFAULT_RELEASE):
    """
    make a color jpeg for the specified files

    parameters
    ----------
    gfile, rfile, ifile: string
        DES g,r,i coadd files
    outfile: string
        Output jpeg file
    rebin: int, optional
        Amount to rebin image
    release: string
        Release; scaling differs for different releases
        currently y1a1 and sva1
    """

    image_maker=RGBImageMaker(gfile, rfile, ifile,
                              rebin=rebin,
                              release=release)

    image_maker.make_image()
    image_maker.write_image(outfile)

class RGBImageMaker(object):
    def __init__(self, gfile, rfile, ifile,
                 rebin=None,
                 release=DEFAULT_RELEASE):

        self.gfile=gfile
        self.rfile=rfile
        self.ifile=ifile

        self.rebin=rebin
        self.release=release

        self.satval=1.0e9

    def _make_imlist(self):
        imlist=[]
        for fname in [self.gfile, self.rfile, self.ifile]:
            print(fname)
            im = ImageTrans(fname)
            imlist.append(im)

        for im in imlist:
            #im.scale_image()
            im.flip_ud()
            #im.transpose()
            if self.rebin is not None:
                im.rebin(rebin)
            
            print()

        self.imlist=imlist

    def make_image(self):
        import images

        import images

        self._make_imlist()

        scales=self._get_scales()

        print("using satval:",self.satval)
        print('getting color image')
        imlist=self.imlist
        colorim=images.get_color_image(imlist[2].image,
                                       imlist[1].image,
                                       imlist[0].image,
                                       scales=scales,
                                       nonlinear=NONLINEAR,
                                       satval=self.satval)

        print('bytescaling')
        colorim = images.bytescale(colorim)

        self.colorim=colorim

    def write_image(self, outfile):
        from PIL import Image
        print('writing:',outfile)

        outfile=os.path.expanduser(outfile)
        outfile=os.path.expandvars(outfile)

        pim=Image.fromarray(self.colorim)
        make_dir(outfile)
        pim.save(outfile, quality=90)

    def _get_scales(self):
        # this will be i,r,g -> r,g,b
        if self.release=='y1a1':
            print('getting scaled color for y1')
            # smaller scale means darker, so noise is more suppressed
            # compared to the peak. remember it is all scaled below
            # one, so we are also cutting off some point in the image
            #SCALE=.010
            #SCALE=.010*sqrt(2.0)
            SCALE=.010*sqrt(2.0)
            relative_scales= array([1.00, 1.2, 2.0])
        elif self.release=='sva1':
            # SVA seems to want a different scaling
            print('getting scaled color for sv')
            #SCALE=.050*0.666
            SCALE=.050*0.88
            #relative_scales= array([1.1, 1.1, 2.0])
            relative_scales= array([1.00, 1.2, 2.5])
        else:
            raise ValueError("bad release: '%s'" % self.release)

        scales= SCALE*relative_scales

        for i in xrange(3):
            im=self.imlist[i]
            print("    scaling",im.band,im.exptime)
            scales[i] *= sqrt(NOMINAL_EXPTIME/im.exptime)
        return scales



class ImageTrans(object):
    def __init__(self, filename):
        with fitsio.FITS(filename) as fits:
            image=fits[1].read()
            header=fits[1].read_header()
            wt=fits[2].read()

        self.image=image
        self.header=header
        self.weight=wt

        self.band=header['FILTER'].split()[0]
        self.exptime=header['exptime']
        self.satval=header['saturate']

    def zero_bad_weightmap(self, minval=0.001):
        print("    zeroing bad weight map")

        wt=self.weight
        w=where(wt < minval)

        if w[0].size > 0:
            print("        max val from image:",self.image[w].max())
            self.image[w] = 0.0

    def flip_ud(self):
        print("    flipping",self.band)
        self.image = flipud(self.image)

    def transpose(self):
        print("    transposing",self.band)
        self.image = self.image.transpose()

    def rebin(self, rebin):

        print("    rebinning",self.band)
        image=self.image

        nrows,ncols=image.shape

        # pad nrows,ncols for rebin
        row_remain=(nrows % rebin)
        if row_remain != 0:
            nrows += (rebin-row_remain)
        col_remain=(ncols % rebin)
        if col_remain != 0:
            ncols += (rebin-col_remain)

        imrebin=zeros( (nrows, ncols), dtype='f4' )

        imrebin[0:image.shape[0], 0:image.shape[1]] = image[:,:]

        imrebin = images.rebin(imrebin, rebin)

        del self.image
        self.image=imrebin

    def scale_image(self, exptime=NOMINAL_EXPTIME):
        print('    scaling',self.band, self.exptime,"to",exptime)

        self.image *= (exptime/self.exptime)
        #self.image *= (exptime/self.exptime)



def make_dir(fname):
    dname=os.path.dirname(fname)
    if dname=='':
        return

    if not os.path.exists(dname):
        print('making dirs:',dname)
        os.makedirs(dname)


class CoaddJPGFiles(dict):
    def __init__(self, run, rebin=None):
        self['run']=run
        self['tile']=run.split('_')[-1]

        self.rebin=rebin

        self._set_files()

    def _set_files(self):
        d=self.get_coadd_dir()
        odir=self.get_output_dir()

        jpg_name='%(tile)s_gri' % self

        if self.rebin is not None:
            jpg_name='%s_rebin%02d' % (jpg_name,int(self.rebin))
        jpg_name='%s.jpg' % jpg_name
        jpg_name=os.path.join(odir, jpg_name)

        self['jpg_file'] = jpg_name
        for band in ['g','r','i']:
            fname='%s_%s.fits.fz' % (self['tile'], band)

            fname=os.path.join(d, fname)
            self[band+'file'] = fname

    def get_coadd_dir(self):
        d='$DESDATA/OPS/coadd/%(run)s/coadd' % self
        d=os.path.expandvars(d)
        return d

    def get_output_dir(self):
        d='$DESDATA/jpg/OPS/coadd/%(run)s/coadd' % self
        d=os.path.expandvars(d)
        return d


