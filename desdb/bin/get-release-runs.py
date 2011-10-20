"""
    %prog [options] release filetype

Print all runs for input release and file type. You might want to see red runs
or coadd runs for instance.  By default the output is one run per line. If you
request the urls, either a csv or json file are written depending on the
requested format.

"""

import os
import sys
from sys import stdout,stderr
import desdb

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-u","--user",default=None, help="Username.")
parser.add_option("-p","--password",default=None, help="Password.")
parser.add_option("-s","--show", action='store_true', help="Show the query on stderr.")
parser.add_option("--url", action='store_true', help="Show the URL for this run.")
parser.add_option("-f","--format",default='json',help=("File format when outputting urls.  csv, json. "
                                                              "Default %default."))


def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    release=args[0].strip()
    filetype=args[1].strip()

    if filetype == 'red':
        extra='order by nite'
    elif filetype == 'coadd':
        extra='order by tilename'
    else:
        extra='order by run'

    base=desdb.files.des_net_rootdir()

    query="""
    select
        distinct(run),
        '%(base)s/' || '%(filetype)s/' || run || '/%(filetype)s' as url,
        '$DESDATA/' || '%(filetype)s/' || run || '/%(filetype)s' as path,
        nite,
        tilename
    from
        %(release)s_files
    where
        filetype='%(filetype)s' %(extra)s\n""" % {'base':base,
                                                  'release':release,
                                                  'filetype':filetype,
                                                  'extra':extra}

    conn=desdb.Connection(user=options.user,password=options.password)

    if options.url:
        conn.quickWrite(query,type=options.format,show=options.show)
    else:
        res=conn.quick(query,show=options.show)
        for r in res:
            print r['run']


if __name__=="__main__":
    main()
