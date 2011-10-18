"""
    %prog [options] table [optional] rownum

Print the column names,typecode,typename,precision,scale,value.

value is for the input rownum, default 1
"""
import os
import sys
import csv
from desdb import desdb

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():

    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    if len(args) > 1:
        rownum = int(args[1])
    else:
        rownum = 1

    table=args[0].strip()
    c=desdb.Connection()
    s=c.conn.createStatement()
    rset=s.executeQuery('select * from %s where rownum = %s' % (table,rownum))

    meta = rset.getMetaData()
    ncols = meta.getColumnCount()

    rset.next()
    for col in xrange(1,ncols+1): 
        print   meta.getColumnName(col),\
                meta.getColumnType(col),\
                meta.getColumnTypeName(col),\
                meta.getPrecision(col),\
                meta.getScale(col),\
                rset.getString(col)

if __name__=='__main__':
    main()
