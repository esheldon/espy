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
parser.add_option("-f","--format", default='json',
                  help=("File format for output.  csv, json, "
                        "cjson, pyobj. Default %default."))

def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)


    release=args[0].strip()
    band=args[1].strip()

    desdb.files.get_red_info(release,band,
                             user=options.user,password=options.password,
                             show=options.show, 
                             doprint=True, fmt=options.format)
    return

if __name__=="__main__":
    main()
