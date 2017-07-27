from __future__ import print_function

import sys
import os
import subprocess
import numpy
from numpy import array, zeros, flipud, where, sqrt
import fitsio

NOMINAL_EXPTIME=900.0

NONLINEAR=.12
DEFAULT_CAMPAIGN='y3a1_coadd'

def make_coadd_jpg(campaign, tilename, rebin=None, cleanup=True):
    """
    make a color jpeg for the specified run

    parameters
    ----------
    campaign: string
        campaign, e.g. y3a1_coadd
    tilename: string
        DES coadd run
    rebin: int, optional
        Amount to rebin image
    """
    files=Files(campaign, tilename, rebin=rebin)
    files.sync()

    image_maker=RGBImageMaker(
        files,
        rebin=rebin,
    )

    image_maker.make_image()
    image_maker.write_image()

    if cleanup:
        files.clean()

class RGBImageMaker(object):
    """
    class to actually make the color image and write it
    """
    def __init__(self,
                 files,
                 rebin=None):

        self.files=files
        self.rebin=rebin

        self.satval=1.0e9

    def _make_imlist(self):
        imlist=[]
        files=self.files
        for fname in [files['gfile'], files['rfile'], files['ifile']]:
            print(fname)
            im = ImageTrans(fname)
            imlist.append(im)

        for i,im in enumerate(imlist):
            #im.scale_image()
            im.flip_ud()
            #im.zero_bad_weightmap()
            #im.transpose()
            if self.rebin is not None:
                im.rebin(rebin)

            """
            minval=0.001
            if i==0:
                bad_logic = im.weight < minval
            else:
                bad_logic = bad_logic | (im.weight < minval)
            """
            print()

        """
        wbad=numpy.where(bad_logic)
        if wbad[0].size != 0:
            print("zeroing",wbad[0].size,"pixels")
            for im in imlist:

                image=im.image
                print("max in bad:",image[wbad].max())
                image[wbad]=0.0
                print("max in bad after:",image[wbad].max())
        """

        self.imlist=imlist

    def make_image(self):
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

    def write_image(self):
        from PIL import Image

        outfile=self.files['jpg_file']
        print('writing:',outfile)

        pim=Image.fromarray(self.colorim)
        make_dir(outfile)
        pim.save(outfile, quality=90)

    def _get_scales(self):
        # this will be i,r,g -> r,g,b
        campaign=self.files['campaign'].upper()
        print('getting scaled color for',campaign)
        if campaign=='Y3A1_COADD':
            # smaller scale means darker, so noise is more suppressed
            # compared to the peak. remember it is all scaled below
            # one, so we are also cutting off some point in the image
            #SCALE=.010
            #SCALE=.010*sqrt(2.0)
            SCALE=.010*sqrt(2.0)
            relative_scales= array([1.00, 1.2, 2.0])
        elif campaign=='Y1A1':
            print('getting scaled color for y1')
            # smaller scale means darker, so noise is more suppressed
            # compared to the peak. remember it is all scaled below
            # one, so we are also cutting off some point in the image
            #SCALE=.010
            #SCALE=.010*sqrt(2.0)
            SCALE=.010*sqrt(2.0)
            relative_scales= array([1.00, 1.2, 2.0])
        elif campaign=='SVA1':
            # SVA seems to want a different scaling
            print('getting scaled color for sv')
            #SCALE=.050*0.666
            SCALE=.050*0.88
            #relative_scales= array([1.1, 1.1, 2.0])
            relative_scales= array([1.00, 1.2, 2.5])
        else:
            raise ValueError("bad campaign: '%s'" % campaign)

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
            print("        zeroing",w[0].size,"bad pixels")
            print("        max val from image:",self.image[w].max())
            self.image[w] = 0.0

    def flip_ud(self):
        print("    flipping",self.band)
        self.image = flipud(self.image)

    def transpose(self):
        print("    transposing",self.band)
        self.image = self.image.transpose()

    def rebin(self, rebin):
        import images
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


class Files(dict):
    """
    deal with files, including syncing
    """
    def __init__(self, campaign, tilename, rebin=None, cleanup=True):
        self['campaign'] = campaign.upper()
        self['tilename'] = tilename
        self._bands=['g','r','i']

        self._rebin=rebin
        self.cleanup=cleanup

        self._set_files()


    def get_remote_coadd_file(self, band):
        """
        local location of coadd fits file
        """
        return os.path.join(
            os.path.expandvars('$DESREMOTE_RSYNC'),
            self._get_coadd_path(band),
        )

    def get_coadd_file(self, band):
        """
        local location of coadd fits file
        """
        return os.path.join(
            self.get_temp_dir(),
            os.path.basename(self._get_coadd_path(band)),
        )

    def _get_coadd_path(self, band):
        """
        get local location for coadd file
        """
        flist=self.get_flist()

        key = '%s-%s' % (self['tilename'], band)
        return flist[key]


    def get_flist(self):
        """
        get the coadd file list
        """
        if not hasattr(self, '_flist'):
            flist_file=self.get_flist_file()
            print("reading:",flist_file)
            data=fitsio.read(flist_file)

            flist={}
            for i in xrange(data.size):
                key=data['key'][i].strip()
                path=data['path'][i].strip()
                flist[key] = path

            self._flist = flist

        return self._flist

    def get_list_dir(self):
        """
        we keep lists here
        """
        return os.path.expandvars('$DESDATA/jpg/lists')

    def get_flist_file(self):
        """
        holds paths to coadds
        """
        dir=self.get_list_dir()
        fname='coadd-flist-%(campaign)s.fits' % self
        return os.path.join(dir, fname)


    def get_temp_dir(self):
        """
        temporary location for the input fits files
        """
        return os.path.join(
            self.get_output_dir(),
            'sources',
        )

    def get_output_dir(self):
        """
        location for the image and temp files
        """
        d='$DESDATA/jpg/%(campaign)s/%(tilename)s' % self
        d=os.path.expandvars(d)
        return d

    def sync(self):
        """
        sync the coadd images
        """
        odir=self.get_temp_dir()
        if not os.path.exists(odir):
            os.makedirs(odir)

        remote_url=self.get_remote_coadd_file('g')
        remote_url = remote_url.replace('_g.fits','_[g,r,i].fits')
        cmd = r"""
    rsync                                   \
        -aP                                 \
        --password-file $DES_RSYNC_PASSFILE \
        %(remote_url)s \
        %(local_dir)s/
        """ % dict(
            remote_url=remote_url,
            local_dir=odir,
        )

        print(cmd)
        subprocess.check_call(cmd,shell=True)

    def clean(self):
        """
        clean up the source files
        """
        odir=self.get_temp_dir()
        if os.path.exists(odir):
            print("removing sources:",odir)
            shutil.rmtree(odir)

    def _set_files(self):
        odir=self.get_output_dir()

        jpg_name='%(tilename)s_gri' % self

        if self._rebin is not None:
            jpg_name='%s_rebin%02d' % (jpg_name,int(self._rebin))
        jpg_name='%s.jpg' % jpg_name
        self['jpg_file'] =os.path.join(odir, jpg_name)

        for band in self._bands:
            self['%sfile' % band] = self.get_coadd_file(band)
            self['%sfile_remote' % band] = self.get_remote_coadd_file(band)

