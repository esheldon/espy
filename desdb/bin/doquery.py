"""
    %prog [options]

Take the query on standard input, write the results on standard output.

"""
import os
import sys
from sys import stdin,stdout
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

    query = stdin.read()

    conn=desdb.Connection(user=options.user,password=options.password)
    res=conn.quickWrite(query,show=options.show,type=options.format)

if __name__=="__main__":
    main()
