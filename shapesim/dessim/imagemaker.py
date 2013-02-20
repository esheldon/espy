import numpy
import lensing
import esutil as eu
import gmix_image

from . import files
from . import noise
from .util import FILTERNUM, FILTERCHAR

PADDING=5.0

# sigma ~ fwhm/TURB_SIGMA_FAC
TURB_SIGMA_FAC=1.68

class ImageMaker(dict):
    def __init__(self, simname, pointing):
        """

        Create an image for the input sim and pointing number. The catalog
        should already exist.

        If poisson noise is added, the data type will be int64, otherwise
        float64

        """

        conf=files.read_config(simname)
        self.update(conf)

        self['pointing_id']=pointing
        self['fnum']=FILTERNUM[self['filter']]

        # use 2*seed for images, seed for catalogs
        numpy.random.seed(2*self['seed'])

        self._load_pointing()
        self._load_catalog()

    def go(self):

        self._make_image()

        nobj=self._data.size
        for i in xrange(nobj):
            if ((i % 1000)==0) or (i==0):
                print '%d/%d' % (i+1,nobj)
            #print '%d/%d' % (i+1,nobj)

            gmix=self._get_gmix(i)
            self._put_object(gmix)

        if self['noise_type'] is not None:
            self._add_noise()
        self._write()

    def _add_noise(self):
        if self['noise_type']=='poisson':
            self._add_poisson_noise()
        else:
            raise ValueError("no other noises implemented")

    def _add_poisson_noise(self):
        print 'adding poisson noise'
        sky=noise.get_sky(self['filter'], self['exptime'],units='e')
        im=self._image
        im += sky

        pim = numpy.random.poisson(im)

        self._image = pim
        del im

    def _put_object(self, gmix):
        from gmix_image.render import _render
        cen=gmix.get_cen()
        T=gmix.get_T()

        sigma=numpy.sqrt(T/2)
        rad=PADDING*sigma

        row_low  = int(cen[0]-rad)
        row_high = int(cen[0]+rad)
        col_low  = int(cen[1]-rad)
        col_high = int(cen[1]+rad)

        _render.fill_model_bbox(self._image, 
                                gmix, 
                                self['nsub'],
                                row_low,
                                row_high,
                                col_low,
                                col_high)
                                    
    def _make_image(self):
        image=numpy.zeros( (self['nrow'], self['ncol']) )
        self._image=image

    def _write(self):

        url=files.get_image_url(self['name'], self['pointing_id'])
        self._makedir(url)
        print url
        eu.io.write(url, self._image, clobber=True, 
                    header=self._header)

    def _makedir(self, url):
        try:
            eu.ostools.makedirs_fromfile(url)
        except:
            pass



    def _get_gmix(self, i):
        obj=self._data[i]
        model=obj['model']
        if model=='gauss':
            model='coellip'

        if model=='full':
            raise ValueError("implement full gmix")
        else:
            shape=lensing.shear.Shear(g1=obj['g1'],
                                      g2=obj['g2'])
            T=(obj['sigma']*2)**2

            pars=numpy.array([obj['row'],
                              obj['col'],
                              shape.e1,
                              shape.e2,
                              T,
                              obj['tflux'] ])

            gmix0=gmix_image.gmix.GMix(pars, type=model)

        if self['psf_type'] is not None:
            psf_gmix=self._get_psf_gmix(i)
            gmix=gmix0.convolve(psf_gmix)
            return gmix
        else:
            return gmix0

    def _get_psf_gmix(self, i):
        if self['psf_type'] != 'constant':
            raise ValueError("implement non-constant PSF")


        if self['psf_model'] == 'gauss':
            sigma = self['psf_sigma']
            if not hasattr(self, '_psf_gmix_cache'):
                fwhm_fac=get_fwhm_fac()
                fwhm_arcsec = sigma*fwhm_fac*self['pixscale']

                Tpsf = 2*sigma**2

                mess='gauss psf, fwhm: %s arcsec sigma: %s pixels'
                print mess % (fwhm_arcsec, sigma)

                pars=numpy.array([1.,1.,
                                  self['psf_e1'],
                                  self['psf_e2'],
                                  Tpsf,
                                  1.0])
                gmix=gmix_image.gmix.GMixCoellip(pars)
                self._psf_gmix_cache=gmix
            return self._psf_gmix_cache
        else:
            raise ValueError("implement other PSFs")

    def _load_catalog(self):
        url=files.get_catalog_url(self['name'], self['pointing_id'])
        print url
        self._data=eu.io.read(url)

        for i in xrange(self._data.size):
            self._data['model'][i] = self._data['model'][i].strip()

    def _load_pointing(self):
        pointings=files.read_pointings(self['name'])

        w,=numpy.where(pointings['index']==self['pointing_id'])
        if w.size==0:
            raise ValueError("bad pointing: %d" % self['pointing_id'])

        self._pointing=pointings[w[0]]
        h={}
        for n in self._pointing.dtype.names:
            val=self._pointing[n]
            if n == 'index':
                h['pointing'] = val
            else:
                h[n] = val

        h['simname'] = self['name']

        self._header=h

def get_fwhm_fac():
    from math import sqrt,log
    fac = 2*sqrt(2*log(2)) # ~2.35
    return fac
