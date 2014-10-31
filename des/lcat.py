"""
code to create xshear input files
"""
from __future__ import print_function
import numpy
import cosmology

from .files_common import *

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
    conf=read_config(lcat_vers)
    d=get_wq_dir(lcat_vers)

    if not os.path.exists(d):
        print("making dir:",d)
        os.makedirs(d)

    for chunk in xrange(conf['nchunk']):

        fname=get_wq_file(lcat_vers, chunk)
        job_name="%s-%06d" % (lcat_vers, chunk)

        text=_wq_template.format(lcat_vers=lcat_vers,
                                 chunk=chunk,
                                 job_name=job_name)

        print("writing:",fname)
        with open(fname,'w') as fobj:
            fobj.write(text)


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
        data=data[beg:end]

        w=self.select(data)
        data=data[w]

        output=self.get_output(data)
        output['maskflags'] = self.get_maskflags(output['ra'],
                                                 output['dec'],
                                                 output['z'])

        fname=get_lcat_file(self['lcat_vers'], chunk)
        write_lcat(fname, output)

    def get_output(self, data):
        """
        copy into output struct
        """
        ndata=self.get_struct(data.size)

        ndata['index'] = data[self['index_col']]
        ndata['ra'] = data[self['ra_col']]
        ndata['dec'] = data[self['dec_col']]
        ndata['z'] = data[self['z_col']]
        
        return ndata

    def select(self, data):
        """
        apply any cuts
        """
        z_logic=self.get_z_logic(data[self['z_col']])
        
        w,=numpy.where(z_logic)
        return w

    def get_z_logic(self, z):
        """
        apply cuts in redshift
        """
        print("    Cutting z to [%g, %g]" % (self['zmin'],self['zmax']))
        logic = (z > self['zmin']) & (z < self['zmax']) 

        w,=numpy.where(logic)
        print("    Keeping %d/%d" % (w.size,z.size))
        if w.size == 0:
            raise ValueError("No objects passed z cut")

        return logic

    def get_maskflags(self, ra, dec, z):
        """
        get the maskflags
        """
        mask_type=self.get('mask_type',None)
        if mask_type is None:
            return numpy.zeros(ra.size)
        else:
            import healpix_util as hu
            if mask_type != 'healpix':
                raise ValueError("only healpix supported for now")

            # Da and rmax in Mpc
            print("    max radius for maskflags: %0.1f" % self['rmax'])

            Da = self.cosmo.Da(0.0, z)

            rmax = self['rmax']
            radius_degrees = rmax/Da*180./numpy.pi

            print("    reading healpix map:",self['mask_file'])
            hmap=hu.readDensityMap(self['mask_file'])
            
            maskflags=numpy.zeros(ra.size, dtype='i8')
            print("    getting maskflags")
            for i in xrange(ra.size):
                if (i % 1000) == 0:
                    print("        %d/%d" % (i,ra.size))

                maskflags[i] = hmap.check_quad(ra[i],
                                               dec[i],
                                               radius_degrees[i],
                                               self['mask_ellip_max'])

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
        return read_lcat_original_file(self['lcat_name'])
    
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



def get_wq_dir(lcat_vers):
    dir=os.environ['TMPDIR']
    dir=os.path.join(dir,'des-lcat-make',lcat_vers)
    return dir

def get_wq_file(lcat_vers, chunk):
    dir=get_wq_dir(lcat_vers)
    name='%s-%06d.yaml' % (lcat_vers, chunk)
    return os.path.join(dir, name)

_wq_template="""
command: |
    lcat_vers={lcat_vers}
    chunk={chunk}
    $ESPY_DIR/des/bin/make-xshear-lcat --chunk $chunk $lcat_vers

job_name: {job_name}
"""


