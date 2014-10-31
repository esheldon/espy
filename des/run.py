"""

runs are a combination of source and lens sample with configuration details

"""

from __future__ import print_function
from .files_common import *
from .import scat
from . import lcat
from . import output

from .wqscripts import XShearWQJob, RedshearWQJob

class Run(dict):
    """
    for writing run configs and wq submit files
    """
    def __init__(self, run):
        conf=cascade_config(run)
        self.update(conf)

    def write_all(self):
        """
        write config and yaml files
        """
        self.write_config()
        self.write_wq()

    def write_config(self):
        """
        write the xshear config
        """
        xshear_conf=XShearConfig(self['run'])
        xshear_conf.write()

    def write_wq(self):
        """
        write the cfg file
        """
        self.write_xshear_wq()
        self.write_redshear_wq()
        return
        self.write_collate_wq()
    
    def write_xshear_wq(self):
        """
        write all the chunks and tiles
        """
        lens_nchunk=self['lens_conf']['nchunk']
        tilenames=scat.get_tilenames(self['source_conf']['scat_name'])

        ntile=len(tilenames)
        for lens_chunk in xrange(lens_nchunk):
            for i,tilename in enumerate(tilenames):
                print("    %d/%d: %s" % (i,ntile,tilename))
                # first check if this source catalog exists
                if self._scat_exists(tilename):
                    job=XShearWQJob(self['run'],
                                    lens_chunk,
                                    tilename)
                    job.write()
                else:
                    print("    skipping due to missing scat")

    def write_redshear_wq(self):
        """
        write all the chunks
        """
        lens_nchunk=self['lens_conf']['nchunk']

        for lens_chunk in xrange(lens_nchunk):
            job=RedShearWQJob(self['run'], lens_chunk)
            job.write()

    def _scat_exists(self, tilename):
        """
        check if the source catalog exists
        """
        fname=scat.get_scat_file(self['scat_vers'],
                                 tilename)
        return os.path.exists(fname)

                                    
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

        d=get_run_dir(self['run'])
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
        fname=pz.get_scinv_file(conf['pz_vers'], conf['pz_type'], chunk)
        zlvals=fitsio.read(fname, ext='zlvals')
        return zlvals


    def get_zlvals_string(self):
        """
        get string version of the zlvals for the scinv interpolation
        """

        zlvals=self.get_zlvals()
        return str(zlvals)

def get_run_basedir():
    """
    lensdir/run
    """
    d=get_des_lensdir()
    return os.path.join(d, 'run')

def get_run_dir(run):
    """
    lensdir/run/run_name
    """
    d=get_run_basedir()
    return os.path.join(d, run)

def get_xshear_config_file(run):
    """
    lensdir/run/{run_name}/{run_name}.cfg

    parameters
    ----------
    run:
        e.g. run-001
    """
    d=get_run_dir(run)
    fname='%s.cfg' % run
    return os.path.join(d, fname)

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


_xshear_wq_template="""
command: |
    source ~esheldon/.bashrc

    hostname
    module load xshear/work

    config=%(config_file)s
    scat=%(scat_file)s
    lcat=%(lcat_file)s

    outf=%(output_file)s

    tmp_dir=$TMPDIR/sobjshear-$RANDOM-$RANDOM
    mkdir -vp $tmp_dir

    dname=$(dirname $outf)
    mkdir -vp $dname

    tmp_outf=$tmp_dir/$(basename $outf)
    tmp_lcat=$tmp_dir/$(basename $lcat)
    tmp_scat=$tmp_dir/$(basename $scat)
    tmp_config=$tmp_dir/$(basename $config)

    echo -e "staging lcat file to local disk:\\n  $lcat\\n  $tmp_lcat"
    cp $lcat $tmp_lcat
    echo -e "staging scat file to local disk:\\n  $scat\\n  $tmp_scat"
    cp $scat $tmp_scat

    echo -e "staging config file to local disk:\\n  $config\\n  $tmp_config"
    cp $config $tmp_config

    xshear $tmp_config $tmp_lcat < $tmp_scat 1> $tmp_outf

    err=$?
    if [[ $err != "0" ]]; then
        echo "Error running xshear: $err"
    fi
    if [[ -e $tmp_outf ]]; then
        echo -e "pushing temp file to\\n  $outf"
        mv -fv $tmp_outf $outf
    fi
    rm -rvf $tmp_dir 2>&1

    date

job_name: "%(job_name)s"
"""

_redshear_template="""
command: |
    source ~/.bashrc
    module load xshear/work

    dir=%(output_dir)s
    conf=%(config_file)s
    outf=%(reduced_file)s

    cat $dir/%(pattern)s | redshear $conf > $outf

job_name: "%(job_name)s"
"""
