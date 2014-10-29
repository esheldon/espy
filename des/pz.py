from __future__ import print_function
import os
import numpy

class DESPofz(object):
    """
    Wrapper for hdf5 file

    dz=DESPofz('0.1', 'skynet')

    print dz.zvals
    print dz[35]
    print dz[35:40]
    """

    def __init__(self, version, type, store='pytables'):
        self.version=version
        self.type=type

        if store=='pytables':
            self._load_pytables()
        else:
            self._load_h5py()

    def __len__(self):
        return self.table.size

    def __enter__(self):
        return self
    def __exit__(self, exception_type, exception_value, traceback):
        self.h5.close()

    def __getitem__(self, arg):
        data=self.table[arg]
        data.dtype.names=('index','pofz')
        return data

    def _load_h5py(self):
        import h5py

        self.key, self.zvals = get_info(self.version, self.type)
        self.fname=get_h5_file(self.version)

        print("loading:",self.fname)
        self.h5=h5py.File(self.fname)

        self.table=self.h5[self.key]['table']
        self.size=self.table.size

    def _load_pytables(self):
        import tables

        self.key, self.zvals = get_info(self.version, self.type)
        self.fname=get_h5_file(self.version)

        print("loading:",self.fname)
        self.h5=tables.open_file(self.fname)

        node=self.h5.getNode('/'+self.key)
        self.table=node.table
        self.size=self.table.shape[0]


def make_scinv_wq(version, type, chunksize):
    """
    make the wq scripts
    """
    chunk=0
    sc=SCinv(version, type, chunksize, chunk)
    sc.write_wq()

def combine_scinv(version, type, chunksize):
    """
    combine all the chunks into one big file
    """
    import fitsio
    outfile=get_scinv_file(version, type)
    print("will write to:",outfile)

    chunk=0
    sc=SCinv(version, type, chunksize, chunk)

    nchunk=sc.nchunk

    with fitsio.FITS(outfile,'rw',clobber=True) as fits:
        for chunk in xrange(nchunk):
            infile=get_scinv_file(version, type, chunk=chunk)
            print(infile)

            with fitsio.FITS(infile) as fin:
                data=fin['scinv'][:]
                if chunk==0:
                    zvals=fin['zlvals'][:]
                    fits.write(zvals,extname='zlvals')
                    fits.write(data,extname='scinv')
                else:
                    fits['scinv'].append(data)


class SCinv(object):
    """
    create scinv in chunks
    """
    def __init__(self, version, type, chunksize, chunk,
                 zlmin=0.095, zlmax=0.95, nzl=57):
        self.version=version
        self.type=type
        self.zlmin=zlmin
        self.zlmax=zlmax
        self.nzl=nzl

        self.chunksize=chunksize
        self.chunk=chunk

        self._load()

    def _load(self):
        with DESPofz(self.version,self.type) as pofz:
            tsize=pofz.size

            nchunk = get_nchunk(self.chunksize, tsize)

            if self.chunk > (nchunk-1):
                raise RuntimeError("chunk %d out of bounds: "
                                   "[0,%d)" % (chunk,nchunk))

            beg=self.chunk*self.chunksize
            end=(self.chunk+1)*self.chunksize

            self.data = pofz[beg:end]
            self.zsvals = pofz.zvals
        
        self.nchunk=nchunk

   
    def make_chunk(self):
        """
        make the chunk
        """
        from lensing.sigmacrit import ScinvCalculator
        outfile=get_scinv_file(self.version, self.type, chunk=self.chunk)
        print("will write to:",outfile)

        data=self.data
        out=self._get_output()

        zs=self.zsvals
        scalc=ScinvCalculator(self.zlmin, self.zlmax, self.nzl, zs[0], zs[-1])

        nobj=data.size
        printstep=nobj//10
        for i in xrange(nobj):
            if (i==0) or ((i+1) % printstep) == 0:
                print("    %d/%d" % (i+1,nobj))
            out['scinv'][i,:]=scalc.calc_mean_scinv(zs, data['pofz'][i,:])

        self._write_data(outfile, scalc.zlvals, out)

    def _get_output(self):
        out=numpy.zeros(self.data.size, dtype=[('index','i8'),
                                               ('scinv','f8',self.nzl)])
        out['index'] = self.data['index'].astype('i8')
        return out


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

        dir=get_scinv_wq_dir(self.version, self.type)
        if not os.path.exists(dir):
            print("making dir:",dir)
            os.makedirs(dir)

        beglist,endlist=get_chunks(self.chunksize, self.nchunk)
        for chunk in xrange(len(beglist)):
            beg=beglist[chunk]
            end=endlist[chunk]
            fname=get_scinv_wq_file(self.version,self.type,chunk)

            job_name='scinv-%s-%06d' % (self.type, chunk)
            text=_wq_scinv.format(version=self.version,
                                  type=self.type,
                                  chunksize=self.chunksize,
                                  chunk=chunk,
                                  job_name=job_name)

            print(fname)
            with open(fname,'w') as fobj:
                fobj.write(text)

def get_nchunk(chunksize, nobj):
    nchunk = nobj/chunksize
    if (nobj % chunksize) != 0:
        nchunk += 1
    return nchunk

def get_chunks(chunksize, nchunk):
    beglist=[]
    endlist=[]

    for chunk in xrange(nchunk):
        beg=chunk*chunksize
        end=(chunk+1)*chunksize

        beglist.append(beg)
        endlist.append(end)

    return beglist, endlist

_wq_scinv="""
command: |
    version={version}
    type={type}
    chunk={chunk}
    chunksize={chunksize}
    $ESPY_DIR/des/bin/make-scinv $version $type $chunksize $chunk

job_name: {job_name}
"""

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
        dir=os.path.join(dir,'chunks')
        name='%s_%06d' % (name,chunk)
    name='%s.fits' % name
    return os.path.join(dir, name)

def get_scinv_wq_dir(version, type):
    dir=os.environ['TMPDIR']
    dir=os.path.join(dir,'des-scinv',version,type)
    return dir

def get_scinv_wq_file(version, type, chunk):
    dir=get_scinv_wq_dir(version, type)
    name='DES_scinv_%s_%s_%06d.yaml' % (version, type, chunk)
    return os.path.join(dir, name)



