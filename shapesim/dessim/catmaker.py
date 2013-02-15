import pprint
import numpy
import esutil as eu
from . import files
from . import gprior

class SimpleCatalogMaker(dict):
    def __init__(self, simname, pointing):
        """
        Create a catalog for the input sim and pointing number

        This is "simple" becuase a constant shear is used, and shapes are drawn
        from a simple shape distribution with random orientations.  All galaxies
        are, for now, the same size.
        """

        self['pointing_id']=pointing

        conf=files.read_config(simname)
        self.update(conf)
        pprint.pprint(self)

    def go(self):
        self._load_pointing()
        self._load_cols()
        self._load_orig_data()
        self._set_output()
        self._load_ellip_generator()

        self._generate_ellips()


    
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

    def _generate_ellips(self):
        print 'generate intrinsic ellipticities'
        mags=self._orig_data['tmag'][:,self['model_mag']]
        g1, g2 = self._ellip_generator.sample2d(mags)

        self._data['g1'] = g1
        self._data['g2'] = g2

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
        self._data = self._get_struct(self._orig_data.size)

    def _get_struct(self, n):
        dt=[('ra','f8'),
            ('dec','f8'),
            ('row','f8'),
            ('col','f8'),
            ('tmag','f8'),
            ('flux','f8'),
            ('g1','f8'),
            ('g2','f8'),
            ('gamma1','f8'),
            ('gamma2','f8')]

        return numpy.zeros(n, dtype=dt)
