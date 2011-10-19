"""
    %prog [options] release filetype

Print all runs for input release and file type. You might want
to see red runs or coadd runs for instance.

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
parser.add_option("-f","--format",default='csv',help=("File format for output.  csv, json, "
                                                      "json-pretty. Default %default."))


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
        '{base}/' || '{filetype}/' || run || '/{filetype}' as url,
        '$DESDATA/' || '{filetype}/' || run || '/{filetype}' as path,
        nite,
        tilename
    from
        {release}_files
    where
        filetype='{filetype}' {extra}\n""".format(base=base,
                                                  release=release,
                                                  filetype=filetype,
                                                  extra=extra)

    conn=desdb.Connection(user=options.user,password=options.password)

    if options.url:
        conn.quickWrite(query,type=options.format,show=options.show)
    else:
        res=conn.quick(query,show=options.show)
        if options.format == 'csv':
            print 'run'
            for r in res:
                print r['run']
        else:
            res=[{'run':r['run']} for r in res]
            desdb.desdb.write_json(res,options.format) 


if __name__=="__main__":
    main()
