"""
    %prog [options] run njob

Description
    Generate wq scripts for parallel collation
"""
import des
import sys
from sys import stdout


from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('--hosts',default=None,
                  help="A list of hosts")

options, args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

run=args[0]
njob=int(args[1])

hosts=options.hosts
if isinstance(hosts,basestring):
    hosts = hosts.split(',')

c = des.collate.CollateWQJob(run,njob,hosts=hosts)
c.write()
