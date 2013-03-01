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
    module unload gsim && module load gsim/%(vers)s
    config=%(config)s
    catalog=%(catalog)s
    image=%(image)s
    gsim "$config" "$catalog" "$image"

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

sky = %(sky)f

seed = %(seed)d
"""

def write_cfg(simname, pointing, conf):
    cfg_url=files.get_gsim_cfg_url(simname,pointing)

    seed = catmaker.get_seed(conf['seed'],pointing)

    sky = noise.get_sky(conf['filter'],
                        conf['exptime'])

    text = _cfg_template % {'nrow':conf['nrow'],
                            'ncol':conf['ncol'],
                            'noise_type':conf['noise_type'],
                            'nsub':conf['nsub'],
                            'sky':sky,
                            'seed':seed}

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

    cat=files.get_catalog_url(simname,pointing,type='ascii')
    im=files.get_image_url(simname,pointing)

    cfg_url=files.get_gsim_cfg_url(simname,pointing)

    job_name_t=simname+'-'+files.get_pid_format()
    job_name=job_name_t % pointing
    text=_wqtemplate % {'config':cfg_url,
                        'catalog':cat,
                        'image':im,
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
        write_cfg(simname, pstruct['index'], conf)
        write_wq(simname,pstruct['index'], conf, options)

main()
