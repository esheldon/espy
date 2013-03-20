import pprint
import numpy
import esutil as eu
from . import files
from . import gprior
from . import noise

from .util import FILTERNUM, FILTERCHAR

def get_seed(overall_seed, pointing):
    return overall_seed + pointing

class SimpleCatalogMaker(dict):
    def __init__(self, simname, pointing):
        """
        Create a catalog for the input sim and pointing number
        """

        conf=files.read_config(simname)
        self.update(conf)

        self['pointing_id']=pointing
        self['fnum']=FILTERNUM[self['filter']]

        seed = get_seed(self['seed'],pointing)
        numpy.random.seed(seed)

        pprint.pprint(self)

    def go(self):
        self._load_pointing()
        self._load_cols()
        self._load_orig_data()
        self._set_output()
        self._set_tmag_tflux()
        self._set_row_col()
        self._generate_models()

        self._load_ellip_generator()
        self._generate_ellip()

        self._load_size_generator()
        self._generate_sizes()

        self._add_shear()

    def write_fits(self):
        url=files.get_catalog_url(self['name'], self['pointing_id'])
        self._makedir(url)
        print url
        eu.io.write(url, self._data, clobber=True)

    def write_ascii(self):
        self.write_gal_ascii()
        self.write_psf_ascii()

    def write_gal_ascii(self):
        from esutil import recfile
        import lensing
        url=files.get_catalog_url(self['name'], self['pointing_id'],
                                  ftype='ascii')
        self._makedir(url)
        data=self._get_ascii_struct(self._data.size)

        data['model'] = self._data['model']
        w,=numpy.where( data['model'] == 'gdev' )
        if w.size != 0:
            data['model'][w] = 'dev'
        w,=numpy.where( data['model'] == 'gexp' )
        if w.size != 0:
            data['model'][w] = 'exp'

        data['row'] = self._data['row']
        data['col'] = self._data['col']

        g1=self._data['g1']
        g2=self._data['g2']
        e1,e2 = lensing.util.g1g2_to_e1e2(g1,g2)

        data['e1']    = e1
        data['e2']    = e2
        data['sigma'] = self._data['sigma']
        data['tflux'] = self._data['tflux']

        data['psf_model'] = self['psf_model']
        data['psf_e1']    = self['psf_e1']
        data['psf_e2']    = self['psf_e2']
        data['psf_sigma'] = self['psf_sigma']

        print url
        with recfile.Recfile(url,mode='w',delim=' ',ignorenull=True) as fobj:
            fobj.write(data)

    def write_psf_ascii(self):
        """
        Write a grid of psfs at different central sub-pixel
        locations

        Give them all the flux of the brightest object in the field.
        """
        from esutil import recfile
        import lensing

        if self['psf_type'] is None:
            print 'no psf, not writing psf catalog'
            return

        flux = self._data['tflux'].max()

        url=files.get_catalog_url(self['name'],
                                  self['pointing_id'],
                                  type='psf',
                                  ftype='ascii')
        self._makedir(url)

        nrow=self['psf_nrow']
        ncol=self['psf_ncol']
        nstar_row = self['psf_nstar_row']
        nstar_col = self['psf_nstar_col']
        nstars=nstar_row*nstar_col

        data=self._get_ascii_struct(nstars)

        rows_per_star=nrow/float(nstar_row)
        cols_per_star=nrow/float(nstar_col)

        irow,icol = numpy.mgrid[0:nstar_row, 0:nstar_col]
        row = rows_per_star*(0.5 + irow)
        col = cols_per_star*(0.5 + icol)

        row += numpy.random.random(nstars).reshape(nstar_row,nstar_col)
        col += numpy.random.random(nstars).reshape(nstar_row,nstar_col)

        data['model'] = 'star'

        data['row'] = row.ravel()
        data['col'] = col.ravel()

        data['e1']    = -0.0999999 # not used
        data['e2']    = -0.0999999 # not used
        data['sigma'] =  0.0999999 # not used
        data['tflux'] = flux

        data['psf_model'] = self['psf_model']
        data['psf_e1']    = self['psf_e1']
        data['psf_e2']    = self['psf_e2']
        data['psf_sigma'] = self['psf_sigma']

        print url
        with recfile.Recfile(url,mode='w',delim=' ',ignorenull=True) as fobj:
            fobj.write(data)

    def _makedir(self, url):
        try:
            eu.ostools.makedirs_fromfile(url)
        except:
            pass

    
    def _load_pointing(self):
        purl = files.get_pointings_url(self['name'])

        pointings = eu.io.read(purl,verbose=True)
        w,=numpy.where(pointings['index']==self['pointing_id'])
        if w.size==0:
            raise ValueError("no such pointing found: %d" % self['pointing_id'])

        self._pointing=pointings[w]
        self._wcs = eu.wcsutil.WCS(self._pointing)

    def _load_cols(self):
        self._cols=files.open_columns(self['orig_vers'])

    def _generate_models(self):
        if len( self['models'] ) > 1:
            models = self._draw_random_models()
            self._data['model'] = models
        else:
            self._data['model'] = self['models'][0]

    def _draw_random_models(self):
        """
        For now only allow discreet models
        """
        if len( self['models'] ) != 2:
            raise ValueError("only two models for now")
        fracs = numpy.array(self['model_fracs'])
        fracs /= fracs.sum()

        rand=numpy.random.random(self._data.size)
        w0,=numpy.where(rand < fracs[0])
        w1,=numpy.where(rand > fracs[0])

        models=numpy.zeros(self._data.size, dtype='S5')
        if w0.size > 0:
            models[w0] = self['models'][0]
        if w1.size > 0:
            models[w1] = self['models'][1]
        return models

    def _load_ellip_generator(self):
        print 'getting ellipticity generator'
        if self['model_ellip_type'] == 'cluster-step':
            self._ellip_generator = gprior.GPriorVsMag()
        else:
            raise ValueError("Bad model ellip type: '%s'" % self['model_ellip_type'])

    def _generate_ellip(self):
        print 'generate intrinsic ellipticities'
        fnum=self['fnum']
        mags=self._orig_data['tmag'][:,fnum]
        g1, g2 = self._ellip_generator.sample2d(mags, self._data['model'])

        self._data['g1'] = g1
        self._data['g2'] = g2

    def _set_tmag_tflux(self):
        tmag=self._orig_data['tmag'][:,self['fnum']]
        tflux=noise.get_flux(self['filter'], 
                             tmag, 
                             self['exptime'])

        self._data['tmag'] = tmag
        self._data['tflux'] = tflux

    def _set_row_col(self):
        self._data['row'] = self._row
        self._data['col'] = self._col

    def _load_size_generator(self):
        print 'getting size generator'
        if self['model_size_type'] == 'cluster-step':
            raise ValueError("implement cluster step size")
        elif self['model_size_type'] == 'catalog':
            self._size_generator=None
        elif self['model_size_type'] == 'fixed':
            self._size_generator=None
        else:
            raise ValueError("Bad model ellip type: '%s'" % self['model_size_type'])

    def _generate_sizes(self):
        if self['model_size_type']=='fixed':
            self._data['sigma'] = self['model_sigma'] 
        elif self['model_size_type'] == 'catalog':
            flux_radius=self._orig_data['tsize']/self['pixscale']
            self._data['sigma'] = self._flux_radius_to_sigma(flux_radius)
        else:
            raise ValueError("implement other size generators")

    def _flux_radius_to_sigma(self, flux_radius):
        """
        2*flux_radius=FWHM for gaussians
        """
        return flux_radius*2/2.3548

    def _add_shear(self):
        from lensing.shear import Shear

        data=self._data
        if self['shear_type'] == "constant":
            # over-writing existing cosmolical shear
            print 'adding constant shear:',self['shear']
            data['gamma1'] = self['shear'][0]
            data['gamma2'] = self['shear'][1]
        elif self['shear_type'] == "catalog":
            # these are already in the gamma1,gamma2 fields
            print 'using cosmological shear'
        else:
            raise ValueError("bad shear type: '%s'" % self['shear_type'])

        for i in xrange(data.size):
            shear=Shear(g1=data['gamma1'][i],
                        g2=data['gamma2'][i])
            shape=Shear(g1=data['g1'][i],
                        g2=data['g2'][i])

            sheared_shape = shape + shear

            data['g1'][i] = sheared_shape.g1
            data['g2'][i] = sheared_shape.g2



    def _load_orig_data(self):
        wp=self._get_objects()

        colnames=self._get_colnames()
        self._orig_data=self._cols.read_columns(colnames, 
                                                rows=wp,
                                                verbose=True)

    def _get_colnames(self):
        return ['id','ra','dec','tmag','tsize','gamma1','gamma2']

    def _get_objects(self):

        cols=self._cols

        field_size_pix=max(self['nrow'], self['ncol'])
        field_size_deg=self['pixscale']/3600*field_size_pix

        # buffer by a field size around the center, meaning
        # half a field size on either side
        margin = field_size_deg

        pnt=self._pointing
        print 'selecting objects around pointing'
        wbox=( cols['ra'].between(pnt['crval1'][0] - margin, pnt['crval1'][0]+margin) 
              &
               cols['dec'].between(pnt['crval2'][0] - margin, pnt['crval2'][0]+margin)  )

        if wbox.size==0:
            raise ValueError("no objects found in box")

        print '    found:',wbox.size
        print 'getting image locations'
        ra=cols['ra'][wbox]
        dec=cols['dec'][wbox]
        col,row=self._wcs.sky2image(ra, dec, find=False)

           
        print '    trimming'
        wp,=numpy.where(  (col > 0) 
                        & (col <= self['ncol'])
                        & (row > 0)
                        & (row <= self['nrow']) )
        if wp.size == 0:
            raise ValueError("no objects found in pointing")

 
        self._row=row[wp]
        self._col=col[wp]

        wp=wbox[wp]

        print '    kept:',wp.size,'objects'
        return wp

    def _set_output(self):
        data = self._get_struct(self._orig_data.size)

        for n in data.dtype.names:
            if n in self._orig_data.dtype.names:
                if n not in ['tmag','tflux','row','col']:
                    data[n] = self._orig_data[n]
        self._data=data

    def _get_ascii_struct(self, n):
        dt=[('model','S10'),
            ('row','f8'),
            ('col','f8'),
            ('e1','f8'),
            ('e2','f8'),
            ('sigma','f8'),
            ('tflux','f8'),
            ('psf_model','S10'),
            ('psf_e1','f8'),
            ('psf_e2','f8'),
            ('psf_sigma','f8')]

        return numpy.zeros(n, dtype=dt)

    def _get_struct(self, n):
        dt=[('id','i4'),
            ('ra','f8'),
            ('dec','f8'),
            ('row','f8'),
            ('col','f8'),
            ('tmag','f8'),
            ('tflux','f8'),
            ('sigma','f8'),
            ('model','S10'),
            ('g1','f8'),
            ('g2','f8'),
            ('gamma1','f8'),
            ('gamma2','f8')]

        return numpy.zeros(n, dtype=dt)
