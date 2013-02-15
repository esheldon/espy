import pprint
import numpy
import esutil as eu
from . import files
from . import gprior

from .util import FILTERNUM, FILTERCHAR, MAGLIMS, TEFF

class SimpleCatalogMaker(dict):
    def __init__(self, simname, pointing):
        """
        Create a catalog for the input sim and pointing number

        This is "simple" becuase a constant shear is used, and shapes are drawn
        from a simple shape distribution with random orientations.  All galaxies
        are, for now, the same size.
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
        self._load_ellip_generator()

        self._generate_ellip()
        self._add_shear()


    
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

    def _generate_ellip(self):
        print 'generate intrinsic ellipticities'
        fnum=self['fnum']
        mags=self._orig_data['tmag'][:,fnum]
        g1, g2 = self._ellip_generator.sample2d(mags)

        self._data['g1'] = g1
        self._data['g2'] = g2


    def _add_shear(self):
        from lensing.shear import Shear

        data=self._data
        if self['shear'] is not None:
            print 'adding constant shear:',self['shear']
            data['gamma1'] = self['shear'][0]
            data['gamma2'] = self['shear'][1]
        else:
            print 'using cosmological shear'

        for i in xrange(data.size):
            shear=Shear(g1=data['gamma1'][i],
                        g2=data['gamma2'][i])
            shape=Shear(g1=data['g1'][i],
                        g2=data['g2'][i])

            sheared_shape = shape + shear

            data['g1'][i] = sheared_shape.g1
            data['g2'][i] = sheared_shape.g2


    def _load_ellip_generator(self):
        print 'getting ellipticity generator'
        if self['model_ellip'] == 'cluster-step-nosplit':
            self._ellip_generator = gprior.GPriorVsMagNoSplit()
        else:
            raise ValueError("Bad model ellip type: '%s'" % self['model_ellip'])

    def _load_orig_data(self):
        wp=self._get_objects()

        colnames=self._get_colnames()
        self._orig_data=self._cols.read_columns(colnames, rows=wp,
                                           verbose=True)

    def _get_colnames(self):
        return ['ra','dec','flux','tmag','tsize','gamma1','gamma2']

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
                data[n] = self._orig_data[n]
        self._data=data

    def _get_struct(self, n):
        dt=[('ra','f8'),
            ('dec','f8'),
            ('row','f8'),
            ('col','f8'),
            ('tmag','f8',5),
            ('flux','f8',5),
            ('g1','f8'),
            ('g2','f8'),
            ('gamma1','f8'),
            ('gamma2','f8')]

        return numpy.zeros(n, dtype=dt)
