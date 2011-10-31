"""
    %prog [options] release band

Look up all coadd images in the input release and write out their file ids,
along with some other info. A release id is something like 'dr012' (dc6b)

"""
import os
import sys
from sys import stdout
import desdb


from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-u","--user",default=None, help="Username.")
parser.add_option("-p","--password",default=None, help="Password.")
parser.add_option("-s","--show",action='store_true', help="Show query on stderr.")
parser.add_option("-f","--format",default='pyobj',help=("File format for output.  pyobj, json-pretty."
                                                        "Default %default."))

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
        im.id as image_id,
        im.filetype as image_filetype,
        cat.id as cat_id,
        cat.filetype as cat_filetype,
        im.run,
        im.tilename,
        im.band,
        '$DESDATA/' || im.path as image_url,
        '%(netroot)s/' || im.path as image_url_remote,
        '$DESDATA/' || cat.path as cat_url,
        '%(netroot)s/' || cat.path as cat_url_remote
    from
        %(release)s_files cat,
        %(release)s_files im
    where
        cat.filetype='coadd_cat'
        and cat.band = '%(band)s'
        and cat.catalog_parentid = im.id
    order by tilename\n""" % {'netroot':net_rootdir,
                              'release':release,
                              'band':band}

    conn=desdb.Connection(user=options.user,password=options.password)

    res=conn.quick(query,show=options.show)

    out={}
    for c in res:
        image_id=c['image_id']
        out[image_id] = c

    if options.format[0:4] == 'json':
        import json
        json.dump(out, stdout, indent=1, separators=(',', ':'))
    else:
        import pprint
        pprint.pprint(out)

if __name__=="__main__":
    main()
