"""
    %prog [options] release band

Look up all coadd images for the requested release and bandpass, find the
single epoch 'red' images that were used as input, and write out a json file
with the coadd and red image info.  The json file is keyed by coadd_id.
Release is something like 'dr012'

The results are also written to stderr as you go so you can monitor progress.

columns will be

coadd_id,red_id,filetype,run,exposurename,band,ccd,filename

"""
import os
import sys
from sys import stdout,stderr
import desdb
import csv
import json

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-u","--user",default=None, help="Username.")
parser.add_option("-p","--password",default=None, help="Password.")
parser.add_option("-v","--verbose",action="store_true",default=False, 
                  help="Print out queries as they are executed.")

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
    SELECT
        id
    FROM
        %s_files
    WHERE
        filetype='coadd'
        and band = '%s'\n""" % (release,band)

    conn=desdb.Connection(user=options.user,password=options.password)

    res = conn.quick(query, show=verbose)
    first=True

    json_output={}

    for iddict in res:
        coadd_id = iddict['id']
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

        # now the idlist comes from id instead of parentid
        query="""
        select
            %s as coadd_id,
            id as red_id,
            filetype,
            run,
            exposurename,
            band,
            ccd,
            filename
        from
            location
        where
            id in (%s)
        """ % (coadd_id,idcsv)

        #net_rootdir=desdb.files.des_net_rootdir()
        query="""
        select
            %s as coadd_id,
            id as red_id,
            '$DESDATA/' || path as path
        from
            %s_files
        where
            id in (%s)
        order by
            id\n""" % (coadd_id, release, idcsv)
            

        res = conn.quick(query, show=verbose)

        # write to stderr so we can see progress
        #keywords=['coadd_id','id','path','url']
        keywords=['coadd_id','red_id','path']
        w=csv.DictWriter(stderr,keywords)
        stderr.write(','.join(keywords))
        stderr.write('\n')
        for r in res:
            w.writerow(r)

        json_output[coadd_id] = res
        first=False
    
    json.dump(json_output, stdout, indent=1, separators=(',', ':'))

if __name__=="__main__":
    main()
