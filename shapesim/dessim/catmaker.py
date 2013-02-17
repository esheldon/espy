import pprint
import numpy
import esutil as eu
from . import files
from . import gprior
from . import noise

from .util import FILTERNUM, FILTERCHAR

class SimpleCatalogMaker(dict):
    def __init__(self, simname, pointing):
        """
        Create a catalog for the input sim and pointing number

        For now just set the ellipticities.  Will want to scale the
        fluxes later when we add noise.
        """

        conf=files.read_config(simname)
        self.update(conf)

        self['pointing_id']=pointing
        self['fnum']=FILTERNUM[self['filter']]

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

        self._write()

    def _write(self):
        url=files.get_catalog_url(self['name'], self['pointing_id'])
        self._makedir(url)
        print url
        eu.io.write(url, self._data, clobber=True)

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
            raise ValueError("implement multiple models")
        self._data['model'] = self['models'][0]

    def _load_ellip_generator(self):
        print 'getting ellipticity generator'
        if self['model_ellip_type'] == 'cluster-step-nosplit':
            self._ellip_generator = gprior.GPriorVsMagNoSplit()
        else:
            raise ValueError("Bad model ellip type: '%s'" % self['model_ellip_type'])

    def _generate_ellip(self):
        print 'generate intrinsic ellipticities'
        fnum=self['fnum']
        mags=self._orig_data['tmag'][:,fnum]
        g1, g2 = self._ellip_generator.sample2d(mags)

        self._data['g1'] = g1
        self._data['g2'] = g2

    def _set_tmag_tflux(self):
        tmag=self._orig_data['tmag'][:,self['fnum']]
        if self['nexp'] > 1:
            raise ValueError("make sure gain stuff is right")
        tflux=noise.get_flux(tmag, self['exptime']*self['nexp'])

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
            flux_radius=self._orig_data['tsize']
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

    def _get_struct(self, n):
        dt=[('id','i4'),
            ('ra','f8'),
            ('dec','f8'),
            ('row','f8'),
            ('col','f8'),
            ('tmag','f8'),
            ('tflux','f8'),
            ('model','S10'),
            ('g1','f8'),
            ('g2','f8'),
            ('gamma1','f8'),
            ('gamma2','f8')]

        if 'full' in self['models']:
            raise ValueError("implement full gaussian mixture models")
        else:
            dt += [('sigma','f8')]

        return numpy.zeros(n, dtype=dt)
