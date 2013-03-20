"""
    %prog simname

Make the configuration files and the wq scripts

Both go in the same directory
"""

import sys
import os
import esutil as eu
import shapesim
from shapesim.dessim import files, noise, catmaker

from optparse import OptionParser
parser=OptionParser(__doc__)

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-g','--groups',default=None,
                  help='groups for wq, csv. Default is unconstrained')
parser.add_option('-n','--notgroups',default=None,
                  help='not groups for wq, csv. Default is unconstrained')
parser.add_option('-p','--priority',default='med',
                  help='priority for queue')
parser.add_option('--vers',default='work',
                  help="version gsim")

_wqtemplate="""
command: |
    source ~/.bashrc
    module unload gsim && module load gsim/%(vers)s
    config=%(config)s
    psf_config=%(psf_config)s

    %(obj_command)s
    %(psf_command)s

%(groups)s
%(notgroups)s
job_name: %(job_name)s
priority: %(priority)s
"""


_cfg_template="""
nrow = %(nrow)d
ncol = %(ncol)d

noise_type = "%(noise_type)s"
nsub = %(nsub)d

sky = %(sky).16g

seed = %(seed)d

# optional WCS
wcs = {
    cd1_1 = %(cd1_1).16g
    cd1_2 = %(cd1_2).16g
    cd2_1 = %(cd2_1).16g
    cd2_2 = %(cd2_2).16g

    crval1 = %(crval1).16g
    crval2 = %(crval2).16g

    crpix1 = %(crpix1).16g
    crpix2 = %(crpix2).16g

    ctype1 = "%(ctype1)s"
    ctype2 = "%(ctype2)s"
}
"""

def write_cfg(simname, pstruct, conf, type='objects'):
    pointing=pstruct['index']
    cfg_url=files.get_gsim_cfg_url(simname,pointing,type=type)

    seed = catmaker.get_seed(conf['seed'],pointing)

    sky = noise.get_sky(conf['filter'],
                        conf['exptime'])

    if type=='objects':
        nrow=conf['nrow']
        ncol=conf['ncol']
    else:
        nrow=conf['psf_nrow']
        ncol=conf['psf_ncol']

    d={'nrow':nrow,
       'ncol':ncol,
       'noise_type':conf['noise_type'],
       'nsub':conf['nsub'],
       'sky':sky,
       'seed':seed}

    for n in ['cd1_1','cd1_2','cd2_1','cd2_2',
              'crval1','crval2','crpix1','crpix2',
              'ctype1','ctype2']:
        d[n] = pstruct[n]

    text = _cfg_template % d

    print cfg_url
    with open(cfg_url,'w') as fobj:
        fobj.write(text)

def write_wq(simname, pointing, conf, options):

    groups=''
    if options.groups is not None:
        groups = 'group: [%s]' % options.groups
    notgroups=''
    if options.notgroups is not None:
        notgroups = 'notgroup: [%s]' % options.notgroups

    cat=files.get_catalog_url(simname,pointing,ftype='ascii')
    psf_cat=files.get_catalog_url(simname,pointing,type='psf',
                                  ftype='ascii')
    im=files.get_image_url(simname,pointing)
    psf_im=files.get_image_url(simname,pointing,type='psf')

    cfg_url=files.get_gsim_cfg_url(simname,pointing)
    psf_cfg_url=files.get_gsim_cfg_url(simname,pointing,type='psf')

    job_name_t=simname+'-'+files.get_pid_format()
    job_name=job_name_t % pointing

    obj_command="""
    catalog=%(catalog)s
    image=%(image)s
    gsim "$config" "$catalog" "$image"
    """ % {'catalog':cat, 'image':im}

    if conf['psf_type'] is not None:
        psf_command="""
    psf_catalog=%(catalog)s
    psf_image=%(image)s
    gsim "$psf_config" "$psf_catalog" "$psf_image"
        """ % {'catalog':psf_cat, 'image':psf_im}
    else:
        psf_command=""

    text=_wqtemplate % {'config':cfg_url,
                        'psf_config':psf_cfg_url,
                        'obj_command':obj_command,
                        'psf_command':psf_command,
                        'job_name':job_name,
                        'priority':options.priority,
                        'groups':groups,
                        'notgroups':notgroups,
                        'vers':options.vers}

    wq_url=files.get_gsim_wq_url(simname,pointing)
    print wq_url
    with open(wq_url,'w') as fobj:
        fobj.write(text)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    simname=args[0]
    conf=files.read_config(simname)
    pointings=files.read_pointings(simname)

    d=files.get_gsim_dir(simname)
    if not os.path.exists(d):
        os.makedirs(d)

    for pstruct in pointings:    
        write_cfg(simname, pstruct, conf)
        write_cfg(simname, pstruct, conf, type='psf')
        write_wq(simname,pstruct['index'], conf, options)


main()
