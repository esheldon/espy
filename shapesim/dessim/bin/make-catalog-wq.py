"""
    %prog simname
"""

import sys
import os
import esutil as eu
import shapesim
from shapesim.dessim import files

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
                  help="version for espy and gmix image")

_wqtemplate="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/%(vers)s
    module unload gmix_image && module load gmix_image/%(vers)s
    module unload esutil && module load esutil/%(vers)s

    python $ESPY_DIR/shapesim/dessim/bin/make-catalog.py %(simname)s %(pointing)d

%(groups)s
%(notgroups)s
job_name: %(job_name)s
priority: %(priority)s\n"""


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    simname=args[0]
    
    pointings=files.read_pointings(simname)

    groups=''
    if options.groups is not None:
        groups = 'group: [%s]' % options.groups
    notgroups=''
    if options.notgroups is not None:
        notgroups = 'notgroup: [%s]' % options.notgroups

    d=files.get_wq_dir(simname)
    if not os.path.exists(d):
        os.makedirs(d)

    job_name_t=simname+'-'+files.get_pid_format()
    for pointing in pointings:    
        url=files.get_wq_url(simname,pointing['index'])

        job_name=job_name_t % pointing['index']
        text=_wqtemplate % {'simname':simname,
                            'pointing':pointing['index'],
                            'job_name':job_name,
                            'priority':options.priority,
                            'groups':groups,
                            'notgroups':notgroups,
                            'vers':options.vers}

        print url
        with open(url,'w') as fobj:
            fobj.write(text)
main()
