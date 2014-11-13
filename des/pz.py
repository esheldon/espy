"""
code for p(zs) and scinv(zl)
"""

from __future__ import print_function
import os
import numpy

from .files import *

class DESPofz(object):
    """
    Wrapper for hdf5 file

    dz=DESPofz('0.1', 'skynet')

    print dz.zvals
    print dz[35]
    print dz[35:40]
    """

    def __init__(self, pz_vers, pz_type, store='pytables'):
        self.pz_vers=pz_vers
        self.pz_type=pz_type

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
        """
        ANNZ is currently broken
        """
        data=self.table[arg]

        if self.key=='TPZ':
            data=self._extract_tpz(data)
        else:
            data=self._extract(data)

        return data

    def _extract_tpz(self, data):
        """
        first two elements of "pofz" are mean and mode
        or something
        """

        if len(data.shape)==0:
            data=data.reshape([1])
            is_scalar=True
        else:
            is_scalar=False

        data.dtype.names = ['index','values']
        nobj=data.size

        nbins=data['values'].shape[1]-2

        print("nobj:",nobj)

        dt=[('index','i8'),('pofz','f2',nbins)]
        ndata=numpy.zeros(nobj, dtype=dt)

        ndata['index'] = data['index']
        ndata['pofz']  = data['values'][:,2:]

        if is_scalar:
            ndata=ndata[0]
        return ndata

    def _extract(self, data):
        """
        pull out what we want
        """

        data.dtype.names = ['index','pofz','z_mean','z_peak']
        return data



    def _load_h5py(self):
        import h5py

        res = get_info(self.pz_vers, self.pz_type)
        self.key, self.zvals = res

        self.fname=get_pz_h5_file(self.pz_vers)

        print("loading:",self.fname)
        self.h5=h5py.File(self.fname)

        self.table=self.h5[self.key]['table']
        self.size=self.table.size

    def _load_pytables(self):
        import tables

        res = get_info(self.pz_vers, self.pz_type)
        self.key, self.zvals = res

        self.fname=get_pz_h5_file(self.pz_vers)

        print("loading:",self.fname)
        self.h5=tables.open_file(self.fname)

        node=self.h5.getNode('/'+self.key)
        self.table=node.table
        self.size=self.table.shape[0]


def make_scinv_wq(cosmo_vers, pz_vers, pz_type, chunksize):
    """
    make the wq scripts
    """
    chunk=0
    sc=SCinv(cosmo_vers, pz_vers, pz_type, chunksize, chunk)
    sc.write_wq()

def combine_scinv(cosmo_vers, pz_vers, pz_type, chunksize):
    """
    combine all the chunks into one big file
    """
    import fitsio
    outfile=get_scinv_file(pz_vers, pz_type)
    print("will write to:",outfile)

    chunk=0
    sc=SCinv(cosmo_vers, pz_vers, pz_type, chunksize, chunk)

    nchunk=sc.nchunk

    with fitsio.FITS(outfile,'rw',clobber=True) as fits:
        for chunk in xrange(nchunk):
            infile=get_scinv_file(pz_vers, pz_type, chunk=chunk)
            print(infile)

            with fitsio.FITS(infile) as fits_in:
                data=fits_in['scinv'][:]
                if chunk==0:
                    h=fits_in['scinv'].read_header()
                    h.delete('chunk')
                    fits.write(data,extname='scinv',header=h)
                else:
                    fits['scinv'].append(data)

                if chunk==(nchunk-1):
                    print("creating zlvals extension")
                    zvals=fits_in['zlvals'][:]
                    fits.write(zvals,extname='zlvals')

class SCinv(object):
    """
    create scinv in chunks
    """
    def __init__(self, cosmo_vers, pz_vers, pz_type, chunksize, chunk,
                 zlmin=0.0, zlmax=0.95, nzl=63):

        self.cosmo_vers=cosmo_vers
        self.pz_vers=pz_vers
        self.pz_type=pz_type
        self.zlmin=zlmin
        self.zlmax=zlmax
        self.nzl=nzl

        self.chunksize=chunksize
        self.chunk=chunk

        self._load()

    def _load(self):

        self.cosmo_conf = read_config(self.cosmo_vers)
        with DESPofz(self.pz_vers,self.pz_type) as pofz:
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
        outfile=get_scinv_file(self.pz_vers, self.pz_type, chunk=self.chunk)
        print("will write to:",outfile)

        data=self.data
        out=self._get_output()

        zs=self.zsvals

        conf=self.cosmo_conf
        print("using H0:",conf['H0'],"omega_m:",conf['omega_m'])
        scalc=ScinvCalculator(self.zlmin, self.zlmax, self.nzl, zs[0], zs[-1],
                              H0=conf['H0'],
                              omega_m=conf['omega_m'])

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

        header={'cosmo_vers':self.cosmo_vers,
                'pz_vers':self.pz_vers,
                'pz_type':self.pz_type,
                'chunksize':self.chunksize,
                'chunk':self.chunk}


        print("writing:",outfile)
        with fitsio.FITS(outfile,'rw',clobber=True) as fits:
            fits.write(out, extname='scinv', header=header)
            fits.write(zlvals, extname='zlvals')

    def write_wq(self):
        """
        write a wq script for each chunk
        """

        dir=get_scinv_wq_dir(self.pz_vers, self.pz_type)
        if not os.path.exists(dir):
            print("making dir:",dir)
            os.makedirs(dir)

        beglist,endlist=get_chunks(self.chunksize, self.nchunk)
        for chunk in xrange(len(beglist)):
            beg=beglist[chunk]
            end=endlist[chunk]
            fname=get_scinv_wq_file(self.pz_vers,self.pz_type,chunk)

            job_name='scinv-%s-%06d' % (self.pz_type, chunk)
            text=_wq_scinv.format(cosmo_vers=self.cosmo_vers,
                                  pz_vers=self.pz_vers,
                                  pz_type=self.pz_type,
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
    cosmo_vers={cosmo_vers}
    pz_vers={pz_vers}
    pz_type={pz_type}
    chunk={chunk}
    chunksize={chunksize}
    $ESPY_DIR/des/bin/make-scinv $cosmo_vers $pz_vers $pz_type $chunksize $chunk

job_name: {job_name}
"""

def get_info(pz_vers, pz_type):

    allinfo=read_pz_vers_info()
    if pz_vers not in allinfo:
        raise KeyError("bad pz_vers: '%s'" % pz_vers)

    vinfo = allinfo[pz_vers]
    if pz_type not in vinfo:
        raise KeyError("bad pz_type : '%s'" % pz_type)
    
    info = vinfo[pz_type]

    z_max = info['z_max']
    nbins = info['nbins']
    key   = info['key']

    if key=='TPZ':
        dz = 0.007475

        edge_low  = 0.005-dz/2
        edge_high = 1.5-dz/2
    else:
        edge_low  = 0.0
        edge_high = z_max

    z_edges = numpy.linspace(edge_low, edge_high, nbins+1)

    z_values = (z_edges[1:] + z_edges[:-1]) / 2.0
    return key, z_values



def get_scinv_wq_dir(pz_vers, pz_type):
    dir=os.environ['TMPDIR']
    dir=os.path.join(dir,'des-scinv',pz_vers,pz_type)
    return dir

def get_scinv_wq_file(pz_vers, pz_type, chunk):
    dir=get_scinv_wq_dir(pz_vers, pz_type)
    name='DES_scinv_%s_%s_%06d.yaml' % (pz_vers, pz_type, chunk)
    return os.path.join(dir, name)


