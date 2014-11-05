"""
tools to match the scinv catalogs to shape catalogs

tools to select and write the xshear input (these we call scat)
"""
from __future__ import print_function
import numpy
from .files_common import *
from . import pz
from . import sg

SHAPENOISE=0.16
SHAPENOISE2=SHAPENOISE**2


def make_xshear_input(scat_vers, tilename=None):
    """
    write xshear input files

    parameters
    ----------
    tilename: optional
        Optional tilename, default do all
    """

    if tilename is None:
        conf=read_config(scat_vers)
        tilenames=get_tilenames(conf['scat_name'])
    else:
        tilenames=[tilename]

    ntile=len(tilenames)

    xi=XShearInput(scat_vers)
    for i,tilename in enumerate(tilenames):
        print("-"*70)
        print("processing %d/%d: %s" % (i+1,ntile,tilename))

        xi.write(tilename)

class XShearInput(dict):
    """
    Convert big fits files into stripped down ascii files for
    xshear.  make cuts and selections.
    """
    def __init__(self, scat_vers):
        self['scat_vers']=scat_vers
        conf=read_config(scat_vers)
        self.update(conf)

    def write(self, tilename):
        """
        select objects and write the ascii inputs
        """

        data0=self.read_input_data(tilename)
        if data0 is None:
            return

        w=self.select(data0)
        if w.size == 0:
            print("skipping write: no objects passed cuts")
        else:
            print("%d/%d passed all cuts" % (w.size,data0.size))
            data=data0[w]
            output=self.extract_cols(data)
            self.fix_data(output)
            self.write_data(output, tilename)

    def read_input_data(self, tilename):
        """
        read the input catalog
        """
        import fitsio

        if self['scinv_style']=='interp':
            ext='model_fits'
            inf=get_scinv_matched_file(self['scat_name'],
                                       self['pz_vers'],
                                       self['pz_type'],
                                       tilename)
        else:
            ext=1
            inf=get_dg_scat_file(self['scat_name'],
                                 tilename)
            
        if not os.path.exists(inf):
            print("skipping missing file:",inf)
            return None

        print("reading:",inf)
        data=fitsio.read(inf,ext=ext,lower=True)

        return data

    def extract_cols(self, data):
        """
        get the actual columns we want to write
        """
        from esutil.numpy_util import copy_fields

        dt = [(self['ra_col'], 'f8'),
              (self['dec_col'],'f8'),
              (self['e1_col'], 'f8'),
              (self['e2_col'], 'f8')]

        if self['shear_style']=='lensfit':
            dt += [(self['e1sens_col'],'f8'),
                   (self['e2sens_col'],'f8')]

        dt += [('weight','f8')]

        if self['scinv_style']=='interp':
            nzl=data[self['scinv_col']].shape[1]
            dt += [(self['scinv_col'],'f8',nzl)]
        else:
            dt += [(self['z_col'],'f8')]

        newdata=numpy.zeros(data.size, dtype=dt)

        copy_fields(data, newdata)

        var=(2*SHAPENOISE2
             +     data[self['e_cov_11_col']]
             + 2.0*data[self['e_cov_12_col']]
             +     data[self['e_cov_22_col']])
        weight=1.0/var
        newdata['weight'] = weight

        return newdata

    def fix_data(self, data):
        """
        fix conventions and ranges
        """
        if self['minus_e1']:
            print("        multiplying e1 by -1")
            data[self['e1_col']] *= -1

        if self['minus_e2']:
            print("        multiplying e2 by -1")
            data[self['e2_col']] *= -1

        if self['shear_style']=='lensfit':
            print("        clipping sensitivities")
            c1=self['e1sens_col']
            c2=self['e2sens_col']
            data[c1].clip(min=0.0,max=None,out=data[c1])
            data[c2].clip(min=0.0,max=None,out=data[c2])
            print("        min1:",data[c1].min(),"min2:",data[c2].min())

    def write_data(self, data, tilename):
        """
        utility function to write the data
        """
        d=get_scat_dir(self['scat_vers'])
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        fname=get_scat_file(self['scat_vers'], tilename)
        print("writing:",fname)
        write_scat(fname, data)

    def select(self, data):
        """
        select based on flags and other cuts
        """
        flag_logic=self.get_flag_logic(data)
        sg_logic=self.get_sg_logic(data)
        cut_logic=self.get_cut_logic(data)

        logic=flag_logic & sg_logic & cut_logic
        w,=numpy.where(logic)
        return w

    def get_flag_logic(self, data):
        """
        get logic for sx flags and other flags
        """
        logic=(  (data[self['sxflags_col']]==0)
               & (data[self['flags_col']]==0) )
        return logic

    def get_sg_logic(self, data):
        """
        get logic for some s/g separation
        """
        if self['sg_type']=="modest":
            logic=sg.get_modest(data[self['mag_auto_col']],
                                data[self['mag_psf_col']],
                                data[self['class_star_col']],
                                data[self['spread_model_col']],
                                data[self['spreaderr_model_col']])
        else:
            raise ValueError("bad sg_type: '%s'" % self['sg_type'])

        return logic

    def get_cut_logic(self, data):
        """
        get logic for various cuts
        """
        logic=(  (data[self['s2n_col']] > self['s2n_min'])
               & (data[self['Ts2n_col']] > self['Ts2n_min'])
               & (data[self['arate_col']] > self['arate_range'][0])
               & (data[self['arate_col']] < self['arate_range'][1]) )
        return logic

def match_scinv(scat_vers, tilenames=None):
    """
    match the dg scat to scinv using a Matcher

    parameters
    ----------
    scat_vers: string
        e.g. scat-001, implying a config/scat-001.yaml. This
        holds the scat_name, pz_vers, pz_type etc.
    """
    import fitsio

    conf=read_config(scat_vers)
        
    if tilenames is None:
        tilenames=get_tilenames(conf['scat_name'])

    matcher=Matcher(conf['scat_name'],
                    conf['pz_vers'],
                    conf['pz_type'])

    cvers=matcher.header['cosmo_vers'].strip()
    if cvers != conf['cosmo_vers']:
        raise RuntimeError("mismatch in cosmology: '%s' vs "
                           "'%s'" % (cvers,conf['cosmo_vers']))

    ntile=len(tilenames)
    nuse=0
    for i,tilename in enumerate(tilenames):
        print("-"*70)
        print("processing tile %d/%d: %s" % (i+1,ntile,tilename))
        # not all tiles are in Daniel's match file
        fname=get_dg_scat_file(conf['scat_name'], tilename)
        if not os.path.exists(fname):
            print("skipping missing file:",fname)
            continue

        nuse += 1
        matcher.match(tilename)

    print("finally used: %d/%d" % (nuse, ntile))

class Matcher(object):
    """
    match the scat to the sigma crit

    e.g. Matcher('ngmix009'
    """

    def __init__(self, scat_name, pz_vers, pz_type):
        self.scat_name=scat_name
        self.pz_vers=pz_vers
        self.pz_type=pz_type

        res = pz.read_scinv(self.pz_vers,
                            self.pz_type,
                            get_header=True)
        self.data, self.zlvals, self.header=res
        self.nz=self.zlvals.size

        self.ids = self.data['index']

    def match(self, tilename):
        """
        match and write the output file
        """
        from esutil import numpy_util as nu
        dg_data=read_dg_scat(self.scat_name, tilename)

        print("matching")
        mdg, msc = nu.match(dg_data['coadd_objects_id'], self.ids)
        print("    matched: %d/%d" % (mdg.size, dg_data.size))

        if mdg.size == 0:
            print("    skipping")
            return

        dg_data=dg_data[mdg]
        newdata=nu.add_fields(dg_data, [('scinv','f8',self.nz)])
        newdata['scinv'] = self.data['scinv'][msc,:]

        self._write_data(newdata, tilename)

    def _write_data(self, data, tilename):
        """
        write the data, making dir if necessary
        """
        import fitsio
        d=get_scinv_matched_dir(self.scat_name,
                                self.pz_vers,
                                self.pz_type)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        outf=get_scinv_matched_file(self.scat_name,
                                    self.pz_vers,
                                    self.pz_type,
                                    tilename)
        print("writing to:",outf)
        with fitsio.FITS(outf,'rw',clobber=True) as fits:
            fits.write(data, extname='model_fits', header=self.header)
            fits.write(self.zlvals, extname='zlvals')


def get_tilenames(scat_name):
    """
    this assumes the scat_name corresponds to a database table
    """
    print("getting tile list for:",scat_name)

    tilenames=read_tilenames_cache(scat_name)
    print("found",len(tilenames),"tiles")
    return tilenames

def read_tilenames_cache(scat_name):
    """
    write the tilenames to a file
    """
    import yaml
    import desdb
    fname=get_tilename_cache_file(scat_name)

    if not os.path.exists(fname):
        print("    caching from database")
        dir=get_tilename_cache_dir()
        if not os.path.exists(dir):
            print("making dir:",dir)
            os.makedirs(dir)

        with desdb.Connection() as conn:
            q="select distinct(tilename) from %s" % scat_name
            res=conn.quick(q)
            tilenames=[r['tilename'] for r in res]

        
        print("    writing to:",fname)
        with open(fname,'w') as fobj:
            yaml.dump(tilenames, fobj)
    else:
        print("    reading:",fname)
        with open(fname) as fobj:
            tilenames=yaml.load(fobj)

    return tilenames

def get_tilename_cache_dir():
    """
    directory to hold source catalog files
    """
    d=get_des_lensdir()
    return os.path.join(d, 'tilename_cache')

def get_tilename_cache_file(scat_name):
    """
    file to hold tilenames
    """
    dir=get_tilename_cache_dir()
    fname='%s-tilenames.yaml' % scat_name
    return os.path.join(dir, fname)


