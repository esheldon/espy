from __future__ import print_function
import os
import numpy

def get_base_dir():
    return '/astro/u/astrodat/data/DES/EXTRA/photoz/combined'

def get_version_dir(version):
    d=get_base_dir()
    return os.path.join(d, version)

def get_version_info_file():
    dir=get_base_dir()
    name='version-info.yaml'
    return os.path.join(dir, name)

def read_version_info():
    import yaml
    fname=get_version_info_file()
    with open(fname) as fobj:
        data=yaml.load(fobj)
    return data

def get_h5_file(version):
    dir=get_version_dir(version)
    name='DES_photoz_PDFS_%s.h5' % version
    return os.path.join(dir, name)

def get_scinv_dir(version, type):
    dir=get_version_dir(version)
    return os.path.join(dir, '%s-scinv' % type)

def get_scinv_file(version, type, chunk=None):
    dir=get_scinv_dir(version,type)
    name='DES_scinv_%s_%s' % (version, type)

    if chunk is not None:
        name='%s_%06d' % (name,chunk)
    name='%s.fits' % name
    return os.path.join(dir, name)

class DESPofz(object):
    """
    Wrapper for hdf5 file

    dz=DESPofz('0.1', 'skynet')

    print dz.zvals
    print dz[35]
    print dz[35:40]
    """

    def __init__(self, version, type):
        self.version=version
        self.type=type
        self._load()

    def __len__(self):
        return self.table.size

    def __getitem__(self, arg):
        data=self.table[arg]
        data.dtype.names=('index','pofz')
        return data

    def _load(self):
        import h5py

        self.key, self.zvals = get_info(self.version, self.type)
        self.fname=get_h5_file(self.version)

        print("loading:",self.fname)
        self.h5=h5py.File(self.fname)

        self.table=self.h5[self.key]['table']
        self.size=self.table.size

class SCinv(object):
    def __init__(self, version, type, chunksize):
        self.version=version
        self.type=type
        self.pz=DESPofz(version,type)
        self.chunksize=chunksize

        tsize=self.pz.size
        nchunk = tsize/chunksize
        if (tsize % chunksize) != 0:
            nchunk += 1
        
        self.nchunk=nchunk
   
    def make_scinv_chunk(self, chunk, zlmin, zlmax, nzl):
        from lensing.sigmacrit import ScinvCalculator
        outfile=get_scinv_file(self.version, self.type, chunk=chunk)
        print("will write to:",outfile)

        beg=chunk*self.chunksize
        end=(chunk+1)*self.chunksize

        data=self.pz[beg:end]

        out=numpy.zeros(data.size, dtype=[('index','i8'),
                                          ('scinv','f8',nzl)])
        out['index'] = data['index'].astype('i8')

        zs=self.pz.zvals
        scalc=ScinvCalculator(zlmin, zlmax, nzl, zs[0], zs[-1])

        for i in xrange(data.size):
            out['scinv'][i,:]=scalc.calc_mean_scinv(zs, data['pofz'][i,:])

        self._write_data(outfile, scalc.zlvals, out)

    def _write_data(self, outfile, zlvals, out):
        import fitsio
        d=os.path.dirname(outfile)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        print("writing:",outfile)
        with fitsio.FITS(outfile,'rw',clobber=True) as fits:
            fits.write(zlvals, extname='zlvals')
            fits.write(out, extname='scinv')

    def write_wq(self):
        """
        write a wq script for each chunk
        """
        pass

    def get_chunks(self):
        beglist=[]
        endlist=[]

        for chunk in xrange(self.nchunk):
            beg=chunk*self.chunksize
            end=(chunk+1)*self.chunksize

            beglist.append(beg)
            endlist.append(end)

        return beglist, endlist

def get_info(version, type):

    allinfo=read_version_info()
    if version not in allinfo:
        raise KeyError("bad version: '%s'" % version)

    vinfo = allinfo[version]
    if type not in vinfo:
        raise KeyError("bad type : '%s'" % type)
    
    info = vinfo[type]

    z_max = info['z_max']
    nbins = info['nbins']
    key   = info['key']

    z_values = numpy.linspace(0.0, z_max, nbins + 1)
    z_values = (z_values[1:] + z_values[:-1]) / 2.0
    return key, z_values


