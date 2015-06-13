"""
tools to match the scinv catalogs to shape catalogs

tools to select and write the xshear input (these we call scat)

example steps:

    scat_path=$LENSDIR/catalogs/ngmix011-v15/des_sv_wl_ngmix.fits
    info_path=$LENSDIR/catalogs/ngmix011-v15/des_sv_wl_info.fits
    scat_name='ngmix011-v15b'
    tablename='ngmix011'
    create_matt_cat(scat_path, info_path, tablename, scat_name)

    # split everything by coadd tile
    split_orig_scat_by_tile(scat_name, tablename)

    # you can run this from bin/match-scinv as well
    cosmo_vers='cosmo-01'
    pz_vers='v0.1.7'
    pz_type: 'skynet_mag_auto'
    scat_table='ngmix011'
    match_scinv(scat_name, cosmo_vers, pz_vers, pz_type, scat_table=scat_table)

    # then run script from /bin
    make-xshear-scat ngmix011-v15b

precursor will be creating the scinv for the p(z) catalog.  See the .pz module
"""
from __future__ import print_function
import numpy
from .files import *
from . import pz, sg

import esutil as eu


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
        # scat_table is original table the scat was drawn from
        tilenames=get_tilenames(conf['scat_table'])
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

        if 'mask_vers' in self:
            self.mask_info=read_config(self['mask_vers'])

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
                                       self['cosmo_vers'],
                                       tilename)
        else:
            ext=1
            inf=get_orig_scat_file(self['scat_name'],
                                   tilename)
            
        if not os.path.exists(inf):
            print("skipping missing file:",inf)
            return None

        print("reading:",inf)
        data=fitsio.read(inf,ext=ext)

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

        var=(2*self['shapenoise']**2
             +     data[self['e_cov_11_col']]
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

        exists_logic=self.get_exists_logic(data)
        flag_logic=self.get_flag_logic(data)
        sg_logic=self.get_sg_logic(data)
        scinv_logic=self.get_scinv_logic(data)
        cut_logic=self.get_cut_logic(data)
        mask_logic = self.get_mask_logic(data)

        logic=(  exists_logic & flag_logic & sg_logic
               & scinv_logic & cut_logic
               & mask_logic)

        w,=numpy.where(logic)
        return w

    def get_mask_logic(self, data):
        """
        get the maskflags
        """

        logic = numpy.ones(data.size, dtype=bool)
        if 'mask_vers' in self:

            hmap = self._get_density_map()
            weight=hmap.get_weight(data[self['ra_col']], data[self['dec_col']])

            logic = logic & (weight > 0)

            w,=numpy.where(logic)
            print("    %d/%d in mask" % (w.size, data.size))

        return logic

    def _get_density_map(self):
        if not hasattr(self, '_hmap'):
            import healpix_util as hu

            mask_info=self.mask_info
            mask_type=mask_info['mask_type']

            if mask_type != 'healpix':
                raise ValueError("only healpix supported for now")

            print("    reading healpix map:",mask_info['mask_file'])
            self._hmap=hu.readDensityMap(mask_info['mask_file'])

        return self._hmap

    def get_exists_logic(self, data):
        """
        some catalogs have multiple catalogs in them and an exists flag
        """

        logic = numpy.ones(data.size, dtype=bool)
        if 'exists_col' in self:
            exists_logic = data[self['exists_col']] == 1
            w,=numpy.where(exists_logic)
            print("    kept %d/%d that 'exists'" % (w.size, data.size))

            logic = logic & exists_logic

        return logic

    def get_flag_logic(self, data):
        """
        get logic for sx flags and other flags
        """

        logic = numpy.ones(data.size, dtype=bool)

        if 'sxflags_col' in self:
            sxlogic = data[self['sxflags_col']]==0
            w,=numpy.where(sxlogic)
            print("    kept %d/%d sxflags" % (w.size, data.size))

            logic = logic & sxlogic

        if 'flags_col' in self:
            flogic  = data[self['flags_col']]==0
            w,=numpy.where(flogic)
            print("    kept %d/%d source flags" % (w.size, data.size))

            logic = logic & flogic

        return logic

    def get_sg_logic(self, data):
        """
        get logic for some s/g separation
        """
        logic = numpy.ones(data.size, dtype=bool)
        
        if 'sg_type' in self:
            if self['sg_type']=="modest":
                sg_logic=sg.get_modest(data[self['mag_auto_col']],
                                       data[self['mag_psf_col']],
                                       data[self['class_star_col']],
                                       data[self['spread_model_col']],
                                       data[self['spreaderr_model_col']])
                w,=numpy.where(sg_logic)
                print("    kept %d/%d galaxies" % (w.size, data.size))

                logic = logic & sg_logic
            else:
                raise ValueError("bad sg_type: '%s'" % self['sg_type'])

        return logic

    def get_scinv_logic(self, data):
        """
        cut bad p(z) or scinv
        """
        logic=data[self['scinv_flags_col']]==0
        w,=numpy.where(logic)
        print("    kept %d/%d good scinv" % (w.size, data.size))
        return logic

    def get_cut_logic(self, data):
        """
        get logic for various cuts
        """

        logic = numpy.ones(data.size, dtype=bool)

        if 's2n_min' in self:
            elogic = (data[self['s2n_col']] > self['s2n_min'])

            w,=numpy.where(elogic)
            print("        kept %d/%d s2n cuts" % (w.size, data.size))

            logic = logic & elogic

        if 'Trat_min' in self:
            Tg = data[self['T_col']]
            Tp = data[self['psf_T_col']]

            Trat = Tg/Tp
            elogic = (Trat > self['Trat_min'])
            w,=numpy.where(elogic)
            print("        kept %d/%d Trat cuts" % (w.size, data.size))

            logic = logic & elogic


        if 'arate_range' in self:
            elogic = (  (data[self['arate_col']] > self['arate_range'][0])
                      & (data[self['arate_col']] < self['arate_range'][1]) )

            w,=numpy.where(elogic)
            print("        kept %d/%d arate cuts" % (w.size, data.size))

            logic = logic & elogic

        if 'keep_photoz_bins' in self:
            pbins=self['keep_photoz_bins']
            elogic = (data['photoz_bin'] == pbins[0])
            for i in xrange(1,len(pbins)):
                elogic = elogic | (data['photoz_bin'] == pbins[i])

            w,=numpy.where(elogic)
            print("        kept %d/%d p(z) bins cuts" % (w.size, data.size))

            logic = logic & elogic
                
        w,=numpy.where(logic)
        print("    kept %d/%d all cuts" % (w.size, data.size))
        return logic


 
def match_scinv(scat_name, cosmo_vers, pz_vers, pz_type, tilenames=None, scat_table=None):
    """
    match the dg scat to scinv using a Matcher

    parameters
    ----------
    scat_name: string
        name of a source catalog, e.g. ngmix011-v14-mb which implies a
        catalog dir

    cosmo_vers: string
        the cosmology version
    pz_vers: string
        p(z) version in the hdf5 file from Bonnett
    pz_type: string
        p(z) type in the hdf5 file from Bonnett
    tilenames: list, optional
        subset to process

    scat_table: string
        table from which the tiles should be drawn.  Default is
        to use the scat_name as the table name
    """
    import fitsio


    if tilenames is None:
        if scat_table is None:
            scat_table = scat_name
        tilenames=get_tilenames(scat_table)

    matcher=Matcher(scat_name,
                    pz_vers,
                    pz_type,
                    cosmo_vers)

    ntile=len(tilenames)
    nuse=0
    for i,tilename in enumerate(tilenames):
        print("-"*70)
        print("processing tile %d/%d: %s" % (i+1,ntile,tilename))
        # not all tiles are in Daniel's match file
        fname=get_orig_scat_file(scat_name, tilename)
        if not os.path.exists(fname):
            print("skipping missing file:",fname)
            continue

        nuse += 1
        matcher.match(tilename)

    print("finally used: %d/%d" % (nuse, ntile))


class Matcher(object):
    """
    match the scat to the sigma crit

    e.g. Matcher('ngmix009-v14-mb','v0.1.7','tpz')
    """

    def __init__(self, scat_name, pz_vers, pz_type, cosmo_vers):
        self.scat_name=scat_name
        self.pz_vers=pz_vers
        self.pz_type=pz_type
        self.cosmo_vers=cosmo_vers

        if '-mb' in scat_name:
            self.scat_type='mb'
        else:
            self.scat_type='dg'

        res = read_scinv(self.pz_vers,
                         self.pz_type,
                         self.cosmo_vers,
                         get_header=True)
        self.data, self.zlvals, self.header=res
        self.nz=self.zlvals.size

        self.ids = self.data['index']

    def match(self, tilename):
        """
        match and write the output file
        """
        from esutil import numpy_util as nu

        scat_data = read_orig_scat(self.scat_name, tilename)

        print("matching")
        mdg, msc = nu.match(scat_data['coadd_objects_id'], self.ids)
        print("    matched: %d/%d" % (mdg.size, scat_data.size))

        if mdg.size == 0:
            print("    skipping")
            return

        scat_data=scat_data[mdg]
        add_dt=[('scinv_flags','i4'),('scinv','f8',self.nz)]
        newdata=nu.add_fields(scat_data, add_dt)
        newdata['scinv'] = self.data['scinv'][msc,:]
        newdata['scinv_flags'] = self.data['flags'][msc]

        self._write_data(newdata, tilename)

    def _write_data(self, data, tilename):
        """
        write the data, making dir if necessary
        """
        import fitsio
        d=get_scinv_matched_dir(self.scat_name,
                                self.pz_vers,
                                self.pz_type,
                                self.cosmo_vers)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

        outf=get_scinv_matched_file(self.scat_name,
                                    self.pz_vers,
                                    self.pz_type,
                                    self.cosmo_vers,
                                    tilename)
        print("writing to:",outf)
        with fitsio.FITS(outf,'rw',clobber=True) as fits:
            fits.write(data, extname='model_fits', header=self.header)
            fits.write(self.zlvals, extname='zlvals')



def create_matt_cat(scat_path, info_path, tablename, scat_name, remove_cat=None):
    """
    convert the info and data file into a combined file, 
    and write to the standard file path according to
    the scat_name.  get tilenames from the indicated table

    add photoz_bin to the output catalog

    parameters
    ----------
    scat_path: string
        original path
    info_path: string
        original info path
    tablename: string
        name of db table to get tilenames from
    scat_name: string
        name for output cat
    remove_cat: string
        name of catalog to match by id and remove members.

    e.g.

    scat_path='$LENSDIR/catalogs/ngmix011-v15/des_sv_wl_ngmix.fits'
    info_path='$LENSDIR/catalogs/ngmix011-v15/des_sv_wl_info.fits'
    scat_name='ngmix011-v15'
    tablename='ngmix011'
    remove_cat=$LENSDIR/catalogs/redmapper-6.3.3-lgt05-ubermem/redmapper-6.3.3-lgt05-ubermem.fits
    create_matt_cat(scat_path, info_path, tablename, scat_name, remove_cat=remove_cat)
    """

    tile_data = get_tilenames(tablename, full=True)

    print("reading:",info_path)
    info  = fitsio.read(info_path)
    print("reading:",scat_path)
    scat = fitsio.read(scat_path)


    print("selecting")
    w, = numpy.where(  (info['sva1_flags'] == 0)
                     & (info['ngmix_flags'] == 0))

    print("keeping %d/%d" % (w.size, scat.size))
    info = info[w]
    scat = scat[w]

    if remove_cat is not None:
        print("reading remove cat:",remove_cat)
        rmcat=fitsio.read(remove_cat,lower=True)

        # remove with probability equal to membership probability
        r=numpy.random.random(rmcat.size)
        w,=numpy.where(rmcat['p'] > r)
        print("    %d/%d of rmcat passed random "
              "membership probability test" % (w.size,rmcat.size))

        if w.size > 0:
            rmcat=rmcat[w]

            ms,mr = eu.numpy_util.match(scat['coadd_objects_id'],rmcat['id'])
            print("    removing %d/%d" % (ms.size, scat.size))

            if ms.size > 0:
                ind0 = numpy.arange(scat.size)
                ind = numpy.delete(ind0, ms)
                print("    keeping %d/%d" % (ind.size, scat.size))
                scat=scat[ind]
                info=info[ind]

    add_dt = [('ra','f8'), ('dec','f8'),
              ('tilename','S12'),('photoz_bin','i4')]
    print("adding fields:",add_dt)
    output = eu.numpy_util.add_fields(scat, add_dt)
    del scat

    output['ra'] = info['ra']
    output['dec'] = info['dec']
    output['photoz_bin'] = info['photoz_bin'].astype('i4')

    del info

    print("matching to get tiles")
    md, mt = eu.numpy_util.match(output['coadd_objects_id'],
                                 tile_data['coadd_objects_id'])

    print("matched %d/%d" % (md.size, output.size))
    if md.size != output.size:
        raise RuntimeError("not all matched")

    output['tilename'][md] = tile_data['tilename'][mt]

    del tile_data

    fname = get_orig_scat_file_full(scat_name)
    eu.ostools.makedirs_fromfile(fname,verbose=True)

    print("writing:",fname)
    fitsio.write(fname, output, clobber=True)


def split_orig_scat_by_tile(scat_name, tablename):
    """
    match the input scat to the associated table by coadd_objects_id, and
    write out a file for each tile.  cache the tilenames.
    """
    data = read_orig_scat_full(scat_name)

    utiles = numpy.unique(data['tilename'])
    ntile=utiles.size

    for i,tilename in enumerate(utiles):
        w,=numpy.where(data['tilename'] == tilename)
        print("%d/%d %s: %d" % (i+1,ntile,tilename, w.size))

        tdata = data[w]

        fname=get_orig_scat_file(scat_name, tilename)
        print("writing:",fname)
        fitsio.write(fname, tdata, clobber=True)
