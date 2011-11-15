"""
    %prog [options] release band

Look up all red catalogs and images in the input release and write out their
file ids, path info, and external url.  A release id is something like 'dr012'
(dc6b)

"""
import os
import sys
from sys import stderr,stdout
import desdb

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-u","--user",default=None, help="Username.")
parser.add_option("-p","--password",default=None, help="Password.")
parser.add_option("-s","--show",action='store_true', help="Show query on stderr.")
parser.add_option("-f","--format",default='json',help=("File format for output.  csv, json, "
                                                       "cjson, pyobj. Default %default."))

def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)


    release=args[0].strip()
    band=args[1].strip()

    net_rootdir=desdb.files.des_net_rootdir()

    # Note the kludge on filename to remove dups associated with standard star
    # fields
    query="""
    select
        im.file_exposure_name as exposurename,
        im.band,
        im.ccd,
        im.id as image_id,
        '$DESDATA/' || im.path as image_url,
        '%(netroot)s/' || im.path as image_url_remote,
        cat.id as cat_id,
        '$DESDATA/' || cat.path as cat_url,
        '%(netroot)s/' || cat.path as cat_url_remote
    from
        %(release)s_files cat,
        %(release)s_files im
    where
        cat.filetype='red_cat'
        and cat.band='%(band)s'
        and cat.catalog_parentid = im.id
    order by 
        cat_id\n""" % {'netroot':net_rootdir,'release':release,'band':band}

    conn=desdb.Connection(user=options.user,password=options.password)
    conn.quickWrite(query,fmt=options.format,show=options.show)

if __name__=="__main__":
    main()
