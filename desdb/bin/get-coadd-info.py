"""
    %prog [options] release band

Look up all coadd images in the input release and write out their file ids,
along with some other info. A release id is something like 'dr012' (dc6b)

"""
import os
import sys
from sys import stdout
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
parser.add_option("-f","--format",default='csv',help=("File format for output.  csv, json, "
                                                      "json-pretty. Default %default."))

def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)


    release=args[0].strip()
    band=args[1].strip()

    net_rootdir=desdb.files.des_net_rootdir()
    query="""
    select
        id,
        filetype,
        run,
        tilename,
        band,
        filename,
        '$DESDATA/' || path as path,
        '%s/' || path as url
    from
        %s_files
    where
        filetype='coadd'
        and band = '%s'
        order by tilename\n""" % (net_rootdir,release,band)

    conn=desdb.Connection(user=options.user,password=options.password)

    conn.quickWrite(query,type=options.format,show=options.show)

if __name__=="__main__":
    main()
