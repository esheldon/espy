"""
code to create xshear input files

    first run the reformat-and-trim.py to add jackknife regioins, randomize
    ra,dec.  currently this is just a script in one of the directories

    then run

        /bin/make-lcat-wq lcat-vers

    to make wq scripts for doing the lcat creation, which involves the quad
    checks. These scripts run /bin/make-xshear-lcat

"""
from __future__ import print_function
import numpy
from esutil.numpy_util import between
try:
    import cosmology
except ImportError:
    pass

from .files import *

QUADEQ_ALL_OK = 31

def make_xshear_input(lcat_vers, chunk=None):
    """
    write xshear input files

    parameters
    ----------
    chunk: optional
        Optional chunk, default do all
    """

    xi=XShearInput(lcat_vers)
    if chunk is None:
        xi.write_all()
    else:
        xi.write_chunk(chunk)

def make_wq(lcat_vers):
    """
    make the wq scripts
    """
    from .wqscripts import MakeLcatWQJob

    conf=read_config(lcat_vers)
    for chunk in xrange(conf['nchunk']):
        job=MakeLcatWQJob(lcat_vers, chunk) 
        job.write()


class XShearInput(dict):
    """
    write out ascii files with maskflags
    """
    def __init__(self, lcat_vers):
        self['lcat_vers']=lcat_vers
        conf=read_config(lcat_vers)
        self.update(conf)

        
        cconf=read_config(self['cosmo_vers'])
        self.cosmo = cosmology.Cosmo(omega_m=cconf['omega_m'], H0=cconf['H0'])

        self.mask_info=read_config(self['mask_vers'])

    def write_all(self):
        """
        write all chunks
        """
        for chunk in xrange(self['nchunk']):
            print("-"*70)
            print("processing chunk: %d/%d" % (chunk+1,self['nchunk']))
            self.write_chunk(chunk)

    def write_chunk(self, chunk):
        """
        select, cal
        """
        data=self.read_original()
        beg,end=self.get_chunk_range(data.size, chunk)
        print("    working on chunk: %d:%d" % (beg,end))
        data=data[beg:end]

        w=self.select(data)
        data=data[w]

        output=self.get_output(data)
        output['maskflags'] = self.get_maskflags(output['ra'],
                                                 output['dec'],
                                                 output['z'])
        w,=numpy.where(output['maskflags'] > 0)
        print("        keeping %d/%d maskflags" % (w.size,output.size))
        output=output[w]

        fname=get_lcat_file(self['lcat_vers'], chunk)
        write_lcat(fname, output)

    def get_output(self, data):
        """
        copy into output struct
        """
        ndata=self.get_struct(data.size)

        ndata['index'] = data[self['index_col']]

        if self['ra_col'] == 'ra_cent':
            cen_index = self['cen_index']
            print("    using ra_cent[%d]" % cen_index)

            ra = data[self['ra_col']][:,cen_index]
            dec = data[self['dec_col']][:,cen_index]
        else:
            ra = data[self['ra_col']]
            dec = data[self['dec_col']]

        ndata['ra'] = ra
        ndata['dec'] = dec
        ndata['z'] = data[self['z_col']]
        
        return ndata

    def select(self, data):
        """
        apply any cuts
        """
        z_logic=self.get_z_logic(data[self['z_col']])

        logic = z_logic

        w,=numpy.where(logic)
       
        return w

    def get_radec_range_logic(self, ra, dec):
        """
        not currently used

        apply cuts on ra,dec range if set
        """

        ra_range=self.get('ra_range',None)
        dec_range=self.get('dec_range',None)
        if ra_range is not None:
            print("    Cutting ra to [%g, %g]" % tuple(ra_range))
            print("    Cutting dec to [%g, %g]" % tuple(dec_range))

            logic = (  between(ra, ra_range[0], ra_range[1])
                     & between(dec, dec_range[0], dec_range[1])  )

            w,=numpy.where(logic)
            print("        remaining %d/%d" % (w.size,ra.size))
            if w.size == 0:
                raise ValueError("No objects passed ra,dec range cut")

        else:
            logic = numpy.ones(ra.size, dtype='bool')


        return logic

    def get_z_logic(self, z):
        """
        apply cuts in redshift
        """
        print("    Cutting z to [%g, %g]" % (self['zmin'],self['zmax']))
        logic = (z > self['zmin']) & (z < self['zmax']) 

        w,=numpy.where(logic)
        print("        remaining %d/%d" % (w.size,z.size))
        if w.size == 0:
            raise ValueError("No objects passed z cut")

        return logic

    def get_maskflags(self, ra, dec, z):
        """
        get the maskflags
        """

        mask_info=self.mask_info
        mask_type=mask_info['mask_type']

        if mask_type is None:
            print("        setting all maskflags to QUADEQ_ALL_OK")
            return numpy.zeros(ra.size, dtype='i8') + QUADEQ_ALL_OK
        else:
            import healpix_util as hu


            if mask_type != 'healpix':
                raise ValueError("only healpix supported for now")

            # Da and rmax in Mpc
            print("    max radius for maskflags: %0.1f" % self['rmax'])

            Da = self.cosmo.Da(0.0, z)

            rmax = self['rmax']
            radius_degrees = rmax/Da*180./numpy.pi

            print("    reading healpix map:",mask_info['mask_file'])
            hmap=hu.readDensityMap(mask_info['mask_file'])
            
            maskflags=numpy.zeros(ra.size, dtype='i8')
            print("    getting maskflags")
            for i in xrange(ra.size):
                if (i % 1000) == 0:
                    print("        %d/%d" % (i,ra.size))

                maskflags[i] = hmap.check_quad(ra[i],
                                               dec[i],
                                               radius_degrees[i],
                                               mask_info['mask_ellip_max'])

            w,=numpy.where(maskflags > 1)
            print("    %d/%d had good quadrant pairs" % (w.size, ra.size))

        return maskflags


    def get_chunk_range(self, nobj, chunk):
        """
        get the index range for a chunk
        """
        if chunk < 0 or chunk >= self['nchunk']:
            raise ValueError("chunk %s outside of "
                             "range [0,%d]" % (chunk,self['nchunk']-1))

        nper=nobj//self['nchunk']

        beg=chunk*nper

        if chunk==(self['nchunk']-1):
            end=nobj
        else:
            end=(chunk+1)*nper

        return beg,end

    def read_original(self):
        """
        read the original catalog
        """
        return read_lcat_original(self['lcat_name'])
    
    def read(self, chunk=0):
        """
        read the xshear input file
        """
        from esutil.recfile import Recfile
        dt=self.get_dtype()

        with Recfile(fname,'r',dtype=dt, delim=' ') as robj:
            data=robj.read()

        return data

    def get_struct(self, n):
        """
        output structure
        """
        dt=self.get_dtype()
        return numpy.zeros(n, dtype=dt)

    def get_dtype(self):
        dt=[('index','i8'),
            ('ra','f8'),
            ('dec','f8'),
            ('z','f8'),
            ('maskflags','i4')]
        return dt




