from __future__ import print_function
from .files_common import *
from . import pz

def match_scat(scat_vers, tilenames=None):
    """
    match the dg scat to scinv using a Matcher

    parameters
    ----------
    scat_vers: string
        e.g. scat-001, implying a config/scat-001.yaml. This
        holds the scat_name, pz_vers, pz_type etc.
    """

    conf=read_config(scat_vers)
        
    if tilenames is None:
        tilenames=get_tilenames(conf['scat_name'])

    matcher=Matcher(conf['scat_name'],
                    conf['pz_vers'],
                    conf['pz_type'])

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

def get_tilenames(scat_name):
    """
    this assumes the scat_name corresponds to a database table
    """
    import desdb
    print("getting tile list")
    with desdb.Connection() as conn:
        q="select distinct(tilename) from %s" % scat_name
        res=conn.quick(q)
        tilenames=[r['tilename'] for r in res]

    print("found",len(tilenames),"tiles")
    return tilenames

class Matcher(object):
    """
    match the scat to the sigma crit

    e.g. Matcher('ngmix009'
    """

    def __init__(self, scat_name, pz_vers, pz_type):
        self.scat_name=scat_name
        self.pz_vers=pz_vers
        self.pz_type=pz_type

        self.data, self.zlvals = pz.read_scinv_file(self.pz_vers,
                                                    self.pz_type,
                                                    get_zlvals=True)
        self.nz=self.zlvals.size

        self.ids = self.data['index']

    def match(self, tilename):
        """
        match and write the output file
        """
        from esutil import numpy_util as nu
        dg_data=read_dg_scat_file(self.scat_name, tilename)

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
            fits.write(data, extname='model_fits')
            fits.write(self.zlvals, extname='zlvals')



# e.g "ngmix009" for scat_name
dg_name={'ngmix009':'{tilename}_{scat_name}m.fits.gz'}

# daniel's files
def get_dg_scat_file(scat_name, tilename):
    d=get_cat_dir(scat_name+'-dg')
    pattern=dg_name[scat_name]

    fname=pattern.format(scat_name=scat_name, tilename=tilename)
    return os.path.join(d, fname)

def read_dg_scat_file(scat_name, tilename):
    import fitsio
    fname=get_dg_scat_file(scat_name, tilename)
    print("reading:",fname)
    return fitsio.read(fname, lower=True)


def get_scinv_matched_dir(scat_name, pz_vers, pz_type):
    d=get_cat_basedir()
    pattern='{scat_name}-{pz_vers}-{pz_type}-match'
    sub_dir=pattern.format(scat_name=scat_name,
                           pz_vers=pz_vers,
                           pz_type=pz_type)
    return os.path.join(d, sub_dir)

def get_scinv_matched_file(scat_name, pz_vers, pz_type, tilename):
    d=get_scinv_matched_dir(scat_name, pz_vers, pz_type)

    pattern='{tilename}-{scat_name}-{pz_vers}-{pz_type}-match.fits'
    fname=pattern.format(tilename=tilename,
                         scat_name=scat_name,
                         pz_vers=pz_vers,
                         pz_type=pz_type)

    return os.path.join(d, fname)

