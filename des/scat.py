from .files_common import *
from . import pz

def match_scat(scat_vers, tilenames=None):
    import desdb

    conf=read_config(scat_vers)
        
    if tilenames is None:
        tilenames=get_tilenames(conf['scat_name'])

    matcher=Matcher(conf['scat_name'],
                    conf['pz_vers'],
                    conf['pz_type'])

    for tile in tilenames:
        # not all tiles are in Daniel's match file
        fname=get_dg_scat_file(conf['scat_name'], tilename)
        if not os.path.exists(fname):
            print("skipping missing file:",fname)
            continue

        matcher.match(tilename)

def get_tilenames(scat_name):
    """
    this assumes the scat_name corresponds to a database table
    """
    print("getting tile list")
    with desdb.Connection() as conn:
        q="select distinct(tilename) from %(scat_name)s" % scat_name
        res=conn.quick(q)
        tilenames=res['tilename']
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

        self.zlvals, self.data = read_scinv_file(self.pz_vers,
                                                 self.pz_type,
                                                 get_zlvals=True)
        self.nz=self.zlvals.size

        self.ids = self.data['index']

    def match(self, tilename):
        """
        match and write the output file
        """
        import fitsio
        from esutil import numpy_util as nu
        dg_data=read_dg_scat_file(self.scat_name, tilename)

        mdg, msc = nu.match(dg_data['coadd_objects_id'], self.ids)
        print("matched: %d/%d" % (mdg.size, dg_data.size))

        newdata=nu.add_fields(dg_data, [('scinv','f8',self.nz)])
        newdata['scinv'][mdg,:] = self.data['scinv'][msc,:]

        outf=get_scinv_matched_file(self.scat_name,
                                    self.pz_vers,
                                    self.pz_type,
                                    tilename)
        print("writing to:",outf)
        return
        with fitsio.FITS(outf,'rw',clobber=True) as fits:
            fits.write(newdata, extname='model_fits')
            fits.write(self.zlvals, extname='zlvals')



dg_name={'ngmix009':'{tilename}_{catname}m.fits.gz'}

# daniel's files
def get_dg_scat_file(scat_name, tilename):
    d=get_cat_dir(scat_name+'-dg')
    pattern=dg_name[scat_name]

    fname=pattern.format(tilename=tilename)
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

