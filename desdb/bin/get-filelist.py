"""
    %prog [options] release band types

Look up the listed file types and write out there local path information.  The
types should be a comma separated list.

By default the paths are relative to the $DESDATA environment variable, which
is expanded for your system. Send --noexpand to prevent expansion of this
variable, so you will see $DESDATA/{type}/{run}/....

You can also request the remote url be printed along with the local path.
"""
import os
import sys
from sys import stderr,stdout
import desdb

try:
    import cjson
    have_cjson=True
except:
    import json
    have_cjson=False

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-u","--user",default=None, help="Username.")
parser.add_option("-p","--password",default=None, help="Password.")
parser.add_option("-s","--show",action='store_true', help="Show query on stderr.")
parser.add_option("--url",action='store_true', 
                  help="Print remote url and local path in two columns.")
parser.add_option("--noexpand",action='store_true', 
                  help="Don't expand the $DESDATA environment variable.")
parser.add_option("-o","--orderby",default=None, 
                  help="A csv list of fields by which to order the results.")

def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(45)


    release=args[0]
    band=args[1]
    types=args[2]

    net_rootdir=desdb.files.des_net_rootdir()

    # need to put single quotes around the types
    types = types.split(',')
    types = ["'%s'" % t for t in types]
    types = ','.join(types)

    orderby=""
    if options.orderby is not None:
        orderby = 'order by %s' % options.orderby
    query="""
    select
        '$DESDATA/' || path as path,
        '%(netroot)s/' || path as url
    from
        %(release)s_files
    where
        filetype in (%(types)s)
        and band='%(band)s'
    %(orderby)s\n""" % {'types':types,
                        'netroot':net_rootdir,
                        'release':release,
                        'band':band,
                         'orderby':orderby}

    if not options.noexpand:
        query=os.path.expandvars(query)

    conn=desdb.Connection(user=options.user,password=options.password)
    res=conn.quick(query,show=options.show)

    for r in res:
        if options.url:
            print r['url'],r['path']
        else:
            print r['path']

if __name__=="__main__":
    main()
