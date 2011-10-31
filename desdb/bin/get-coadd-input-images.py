"""
    %prog [options] release band

Look up all coadd images for the requested release and bandpass, find the
single epoch 'red' images that were used as input, and write out a json file
with the coadd and red image info.  The json file is keyed by coadd_id.
Release is something like 'dr012'

Because the json file can't be written until the end, the results are also
written to stderr along the way to show progress.

"""
import os
import sys
from sys import stdout,stderr
import desdb
import csv

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-u","--user",default=None, help="Username.")
parser.add_option("-p","--password",default=None, help="Password.")
parser.add_option("-v","--verbose",action="store_true",default=False, 
                  help="Print out queries as they are executed.")
parser.add_option("-f","--format",default='pyobj',help=("File format for output.  pyobj, json-pretty."
                                                        "Default %default."))
def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)


    release=args[0].strip()
    band=args[1].strip()
    verbose=options.verbose

    # ugh, jython is still on 2.5, no nice string formatting
    query="""
    select
        im.id as image_id,
        im.run,
        im.tilename,
        im.band,
        '$DESDATA/' || im.path as image_url,
        '$DESDATA/' || cat.path as cat_url
    from
        %(release)s_files cat,
        %(release)s_files im
    where
        cat.filetype='coadd_cat'
        and cat.band = '%(band)s'
        and cat.catalog_parentid = im.id\n""" % {'release':release,
                                                 'band':band}


    conn=desdb.Connection(user=options.user,password=options.password)

    res = conn.quick(query, show=verbose)
    output={}

    for cdict in res:
        coadd_id = cdict['image_id']
        query="""
        SELECT
            image.parentid
        FROM
            image,coadd_src
        WHERE
            coadd_src.coadd_imageid = %d
            AND coadd_src.src_imageid = image.id\n""" % coadd_id

        res = conn.quick(query, show=verbose)

        idlist = [str(d['parentid']) for d in res]

        ftype=None
        itmax=5

        i=0 
        while ftype != 'red' and i < itmax:
            idcsv = ', '.join(idlist)

            query="""
            SELECT
                id,
                imagetype,
                parentid
            FROM
                image
            WHERE
                id in (%s)\n""" % idcsv

            res = conn.quick(query, show=verbose)
            idlist = [str(d['parentid']) for d in res]
            ftype = res[0]['imagetype']
            
            if verbose: stderr.write('ftype: %s' % ftype)
            i+=1

        if ftype != 'red':
            raise ValueError("Reach itmax=%s before finding 'red' images. last is %s" % (itmax, ftype))

        if verbose: stderr.write("Found %d red images after %d iterations" % (len(idlist),i))

        #net_rootdir=desdb.files.des_net_rootdir()
        query="""
        select
            id,
            run,
            file_exposure_name as exposurename,
            ccd,
            '$DESDATA/' || path as path
        from
            %(release)s_files
        where
            id in (%(idcsv)s)
        order by
            id\n""" % {'release':release, 
                       'idcsv':idcsv}
            

        res = conn.quick(query, show=verbose)

        # write to stderr so we can see progress
        stderr.write('coadd_id id path\n')
        for r in res:
            stderr.write('%d %s %s\n' % (coadd_id,r['id'],r['path']))

        thisone={'coadd_id':coadd_id,
                 'release':release,
                 'tilename':cdict['tilename'],
                 'run':cdict['run'],
                 'band':band,
                 'image_url':cdict['image_url'],
                 'cat_url':cdict['cat_url'],
                 'srclist':res}
        output[coadd_id] = thisone
    
    if options.format[0:4] == 'json':
        import json
        json.dump(output, stdout, indent=1, separators=(',', ':'))
    else:
        import pprint
        pprint.pprint(output)


if __name__=="__main__":
    main()
