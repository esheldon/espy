import numpy
import lensing
import esutil as eu
import gmix_image

from . import files
from .util import FILTERNUM, FILTERCHAR

PADDING=5.0

class ImageMaker(dict):
    def __init__(self, simname, pointing):
        """
        Create an for the input sim and pointing number. The catalog
        should already exist.

        For now just set the ellipticities.  Will want to scale the
        fluxes later when we add noise.
        """

        conf=files.read_config(simname)
        self.update(conf)

        self['pointing_id']=pointing
        self['fnum']=FILTERNUM[self['filter']]

        self._load_catalog()

    def go(self):

        self._make_image()

        nobj=self._data.size
        for i in xrange(nobj):
            if ((i % 100)==0) or (i==0):
                print '%d/%d' % (i+1,nobj)
            #print '%d/%d' % (i+1,nobj)

            gmix=self._get_gmix(i)
            self._put_object(gmix)

        self._write()

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
        eu.io.write(url, self._image, clobber=True)

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

        if self['psf_model_type'] is not None:
            psf_gmix=self._get_psf_gmix(i)
            gmix=gmix0.convolve(psf_gmix)
            return gmix
        else:
            return gmix0

    def _get_psf_gmix(self, i):
        raise ValueError("implement psf")

    def _load_catalog(self):
        url=files.get_catalog_url(self['name'], self['pointing_id'])
        print url
        self._data=eu.io.read(url)

        for i in xrange(self._data.size):
            self._data['model'][i] = self._data['model'][i].strip()

