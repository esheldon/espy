"""
    %prog [options] release band

Look up all red catalogs and images in the input release and write out their
path information on the local machine.   By default the paths are relative to
the $DESDATA environment variable, send --noexpand to prevent expansion
of this variable.

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
        im.nite as nite,
        '$DESDATA/' || im.path as image_path,
        '%(netroot)s/' || im.path as image_url,
        '$DESDATA/' || cat.path as cat_path,
        '%(netroot)s/' || cat.path as cat_url
    from
        %(release)s_files cat,
        %(release)s_files im
    where
        cat.filetype='red_cat'
        and cat.band='%(band)s'
        and im.filename not like 'decam%%-0-%%.fits%%'
        and cat.catalog_parentid = im.id
    order by 
        nite\n""" % {'netroot':net_rootdir,'release':release,'band':band}

    if not options.noexpand:
        query=os.path.expandvars(query)

    conn=desdb.Connection(user=options.user,password=options.password)
    res=conn.quick(query,show=options.show)

    for r in res:
        if options.url:
            print r['image_url'],r['image_path']
            print r['cat_url'],r['cat_path']
        else:
            print r['image_path']
            print r['cat_path']

if __name__=="__main__":
    main()
