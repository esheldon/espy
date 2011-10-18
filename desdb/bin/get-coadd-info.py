"""
    %prog [options] release band

Look up all coadd images in the input release and write out their file ids,
along with some other info. A release id is something like 'dr012' (dc6b)

"""
import os
import sys
from desdb import desdb
import csv

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-u","--user",default=None, help="Username.")
parser.add_option("-p","--password",default=None, help="Password.")

def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)


    release=args[0].strip()
    band=args[1].strip()

    # ugh, jython is still on 2.5, no nice string formatting
    query="""
    select
        id,filetype,run,tilename,band,filename
    from
        %s_files
    where
        filetype='coadd' %s
        and band = '%s'\n""" % (release,band)

    conn=desdb.Connection(user=options.user,password=options.password)

    r = conn.executeWrite(query,show=True)

if __name__=="__main__":
    main()
