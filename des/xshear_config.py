from __future__ import print_function
from .files import *

class XShearConfig(dict):
    """
    For writing the xshear config file
    """
    def __init__(self, run):
        self['run']=run
        conf=cascade_config(run)
        self.update(conf)

        self.comb={}
        self.comb.update(conf)
        self.comb.update(conf['lens_conf'])
        self.comb.update(conf['source_conf'])
        self.comb.update(conf['cosmo_conf'])

    def write(self):
        """
        write the config
        """
        fname=get_xshear_config_file(self['run'])

        d=get_xshear_config_dir(self['run'])
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)
        
        text=self.get_text()
        print("writing:",fname)
        with open(fname,'w') as fobj:
            fobj.write(text)

    def get_text(self):
        """
        get the text of the config file
        """

        comb=self.comb
        text = [_config_top % comb]

        scinv_text=self.get_scinv_text()
        text.append(scinv_text)

        return '\n'.join(text)

    def get_scinv_text(self):
        """
        the config file has different members for interpolation vs point z
        """

        comb=self.comb
        text=[]
        if comb['scinv_style']=="point":
            text.append(_config_point)
            if 'zdiff_min' in comb:
                tmp=_config_zdiff_min % comb
                text.append(tmp)
        elif comb['scinv_style']=="interp":
            zlvals_string=self.get_zlvals_string()
            tmp=_config_interp % {'zlvals':zlvals_string}
            text.append(tmp)
        else:
            raise ValueError("bad scinv_style: '%s'" % comb['scinv_style'])

        return '\n'.join(text)

    def get_zlvals(self):
        """
        get the zlvals for the scinv interpolation
        """
        import fitsio
        from . import pz

        conf=self['source_conf']
        chunk=0
        fname=pz.get_scinv_file(conf['pz_vers'],
                                conf['pz_type'],
                                conf['cosmo_vers'],
                                chunk=chunk)
        zlvals=fitsio.read(fname, ext='zlvals')
        return zlvals


    def get_zlvals_string(self):
        """
        get string version of the zlvals for the scinv interpolation
        """

        zlvals=self.get_zlvals()
        return str(zlvals)


_config_top="""
# cosmology parameters
H0                = %(H0)g
omega_m           = %(omega_m)g

# nside for healpix
healpix_nside     = %(healpix_nside)d

# masking style, for quadrant cuts. "none" "sdss", "equatorial"
mask_style        = "%(mask_style)s"

# shear style
#  "reduced": ordinary reduced shear, source catalog rows are like
#      ra dec g1 g2 weight ...
#  "lensfit": for lensfit with sensitivities
#      ra dec g1 g2 g1sens g2sens weight ...

shear_style       = "%(shear_style)s"

# number of logarithmically spaced radial bins to use
nbin              = %(nbin)d

# min and max radius in Mpc
rmin              = %(rmin)g
rmax              = %(rmax)g

# sigma crit style
#  "point": using source z as truth. Implies the last column in source cat is z
#  "interp": Interpolate 1/sigmacrit (created from full P(z). 
#     Implies last N columns in source cat are \Sigma_{crit}(zlens)_i
"""


_config_point="""
sigmacrit_style   = "point"
"""

_config_zdiff_min="""
# demand zs > zl + zdiff_min
# optional, only makes sense for point z
zdiff_min         = %(zdiff_min)g
"""

_config_interp="""
sigmacrit_style   = "interp"

# zlens values for the \Sigam_{crit}(zlens) values tabulated for each source
# note the elements of arrays can be separated by either spaces or commas

zlvals = %(zlvals)s
"""



