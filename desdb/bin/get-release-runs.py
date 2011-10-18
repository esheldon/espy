"""
    %prog [options] release filetype

Print all runs for input release and file type. You might want
to see red runs or coadd runs for instance.

"""

import os
import sys
from desdb import desdb
from desdb import files
import csv

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-u","--user",default=None, help="Username.")
parser.add_option("-p","--password",default=None, help="Password.")
parser.add_option("-s","--show", action='store_true', help="Show the query on stderr.")
parser.add_option("--url", action='store_true', help="Show the URL for this run.")


def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    release=args[0].strip()
    filetype=args[1].strip()
    show=options.show
    url=options.url

    if filetype == 'red':
        extra='order by nite'
    elif filetype == 'coadd':
        extra='order by tilename'
    else:
        dtype=None
        extra='order by run'


    # ugh, jython is still on 2.5, no nice string formatting
    query="""
    select
        distinct(run),nite,tilename
    from
        %s_files
    where
        filetype='%s' %s\n""" % (release,filetype,extra)

    conn=desdb.Connection(user=options.user,password=options.password)

    res = conn.execute(query,show=show)

    if url:
        dtype='%s_run' % filetype
        df=files.DESFiles(root='web')
        urls = [df.url(dtype,run=r['run']) for r in res]
        print 'run,url'
        for i in xrange(len(res)):
            print "%s,%s" % (res[i]['run'],urls[i])
    else:
        print 'run'
        for r in res:
            print r['run']


if __name__=="__main__":
    main()
