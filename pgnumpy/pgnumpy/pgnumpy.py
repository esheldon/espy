"""
Module:
    pgnumpy
Purpose:
    A numerical python interface to the Postgres database.
    

See docs for the Classes and convenience function for more info.


Classes:
    PgNumpy:
        The class used in all database interactions.  This class represents a
        database connection and facilitates executing queries and extracting
        results.  See docs for pgnumpy.PgNumpy for more details.
    PgInput:  
        A class for writing input files for use in a COPY into the database.
    ArrayWriter:  
        Write arrays to a file for input to postgres.  This slower version can
        be used if recfile is not available.
    ArrayStringifier:  
        Make a string from an array, possibly with brackets indicating
        dimensions.

Convenience Functions:
    connect:
        Create a database connection, returning a PgNumpy object. If conninfo
        is None or "" then the "default" connection based on the PGUSER and
        PGDATABASE environment variables is used.

    array2table:
        Write array with fields (a structure) to a postgres table.  If the
        table does not yet exist it is created with column definitions based on
        the input array. If it does exist the data are appended as new rows in
        the table. 



"""
import os
import sys
from sys import stdout,stderr
import numpy
from numpy import where
import re
import tempfile
import cpgnumpy
import pydoc

# try get get recfile for ascii input
# files.  
have_recfile=False
try:
    import recfile
    have_recfile=True
except:
    try:
        from esutil import recfile
        have_recfile=True
    except:
        pass

stringmatch = re.compile('^string.*')

adminconn_def = 'user=postgres'

# dictionary for type name translation
# the dict key is the numpy type name.  
# Numpy integers are of the form 
#       iX where X is bytes
# or
#       intY where Y is bits.  
# Numpy floating point are 
#       fX where X is in bytes
# or
#       floatY where Y is in bits
#
# Postgres integers or of the form
#       intX with X in bytes
# Note int2 == smallint, int4 == integer or int, int8 == bigint.  Postgres
# does not support unsigned integer types.
# Postgres floating point binary types are:
#       floatX with X in bytes
# Note, float4 == real, float8 == double precision


_npy2pg={}

# convert one-byte int to two byte int2 in pg
# also, no unsigned
_npy2pg['i1']    = 'int2'
_npy2pg['int8']  = 'int2'
_npy2pg['u1']    = 'int2'
_npy2pg['uint8'] = 'int2'

_npy2pg['i2']     = 'int2'
_npy2pg['int16']  = 'int2'
_npy2pg['u2']     = 'int2'
_npy2pg['uint16'] = 'int2'

_npy2pg['i4']     = 'int4'
_npy2pg['int32']  = 'int4'
_npy2pg['u4']     = 'int4'
_npy2pg['uint32'] = 'int4'

_npy2pg['i8']      = 'int8'
_npy2pg['int64']   = 'int8'
_npy2pg['u8']      = 'int8'
_npy2pg['uint64']  = 'int8'

_npy2pg['f4']      = 'float4'
_npy2pg['float32'] = 'float4'
_npy2pg['f8']      = 'float8'
_npy2pg['float64'] = 'float8'

#string_store = 'varchar'
string_store = 'character'
#string_store = 'text'
_npy2pg['string'] = 'character'


def connect(conninfo=None):
    """
    Name:
        connect
    Purpose:
        Get a PgNumpy connection to a postgres database.
    Calling Sequence:
        import pgnumpy
        pg = pgnumpy.connect(conninfo=None)
    Keywords:
        conninfo: 
            Connection info for the database. If blank, the defaults from
            PGUSER and PGDATABASE environment variables is used.

            Example: conninfo='dbname=mydb user=myname'

    """
    if conninfo is not None:
        return PgNumpy(conninfo)
    else:
        return PgNumpy()

def array2table(array, tablename, 
                conninfo=None, 
                primary_key=None, 
                unique=None, 
                serial_extra=None,
                tmpdir=None, 
                keep=False, 
                verbose=False):
    """
    Name:
        array2able
    Purpose:
        Stuff a numerical python array into a table, creating the table if it
        does not exist.  This is a wrapper for the PgInput class.

    Calling Sequence:
        array2table(array, tablename, 
                    conninfo=None, 
                    keep=False, 
                    primary_key=None, 
                    unique=None, 
                    serial_extra=None,
                    tmpdir='.', 
                    verbose=False):
    Inputs:
        array: A numpy array with fields, aka structure or records.
        tablename:  The table name where the data will be stored.
    Keywords:
        conninfo: 
            Connection info for the database. If blank, the defaults from
            PGUSER and PGDATABASE environment variables is used.

            Example: conninfo='dbname=mydb user=myname'

        primary_key:
            Name of the primary key in the table.
        unique:
            Names of columns that should be unique.

        tmpdir: 
            Where to put the temporary file.  Default is the current
            working directory, but remember that this must be visible
            to the postgresql server.
        keep: 
            If True, don't remove the temporary file.

        serial_extra:
            If serial_extra not None, then we specify all the columns that will
            be copied not including the serial column

    Limitations:
        Unsigned integers and 1-byte integers are not supported by postgres.
        Unsigned will be converted to signed, so be careful about ranges.
        1-byte are converted to 2-byte.

        Bool types and object types are not supported.

    """

    pgi = PgInput(array.dtype, tablename, 
                  conninfo=conninfo,
                  primary_key=primary_key,
                  unique=unique,
                  serial_extra=serial_extra,
                  tmpdir=tmpdir,
                  verbose=verbose)
    pgi.write(array, keep=keep)
    del pgi


class PgNumpy(cpgnumpy.cPgNumpy):
    """
    Class:
        PgNumpy
    Purpose:

        A class for working with a PostgreSql database.  The focus is on
        retrieving results as a numerical python array with fields (aka
        structure or records), although generic queries can be executed.

        The base class is written as a C++ extension module to maximumize
        efficiency.

        This class has limited functionality compared to the full Python
        database api specification.  It can execute arbitrary queries and
        extract results into numpy arrays.  However, cursors are not yet
        supported.  For getting results, only the fetchall() command is
        available, as the goal is always to extract all rows into a single
        numpy structure rather than work row by row.

    Construction:
        # there are a few ways
        import pgnumpy
        pg = pgnumpy.connect(conninfo=None)
        pg = pgnumpy.PgNumpy()
        pg = pgnumpy.PgNumpy(conninfo)

    Methods:
        # see the docs for each method for more info.  E.g.
        # pgnumpy.PgNumpy.execute? in ipython

        open(conninfo=None):
            Open a new connection.
            conninfo: 
                Connection info for the database. If blank, the defaults from
                PGUSER and PGDATABASE environment variables is used.

                Example: conninfo='dbname=mydb user=myname'


        execute(query): 
            Execute a query.
            Limitations:
                Some types don't map onto numpy types well.  The most common
                cases are boolean and dates.  It is recommended you convert
                boolean to int or character types, and dates to varchar.

                e.g. SELECT date::varchar FROM mytable 

        fetch(query=None):
        fetchall(query=None):
            Fetch all results from the query into a numerical python array with
            fields.  AKA a structure.  If query is sent it will be executed
            and the results returned.  The fetchall name is an alias for fetch
            and is available for compatibility with the pthon DB API
        
        ntuples():
            Return the number of tuples available in the result.
        nfields():
            Return the number of fields available in the result.
        status():
            Return the query status.


        show(query, page=False, ...other keywords):
            Execute the query and make a nicely formatted printout of the 
            results, optionally sending to a pager.

        describe(table=None):
            Print a description of the indicated table, or if no table name is
            sent, list all public tables.

        tables():
            Get a list of public table names.
        table_exists():
            True if it exists, False if not.


        create_index(tablename, column, unique=False):
            Create an index on the specified table and column.  This shortcut
            method automatically picks a name for the index.  The column spec
            can be a comma separated list for a multi-column index. Columns can
            also be array elements.


        set_field_lengths(flength_dict):

            Send a dictionary specifying the length of one or more
            fields.  e.g.:
                {'colname1':10, 'colname2':25}    

            Since numpy arrays are fixed length, the length must be determined
            before creating the array.  For variable length columns this cannot
            be determined without searching all the results, which could be
            slow in some cases. This gives a way of controlling the lengths
            without searching the result data.

        clear_field_lengths():
            Delete any previously set field length dictionary.
 

    Examples:

        >>> import pgnumpy

        # select results from a table
        >>> pg.execute('select * from mytable')
        >>> res = pg.fetchall()

        # you can also shortcut the execute/fetch process by sending a query to
        # fetchall
        >>> res = pg.fetchall('select * from mytable')


        # print a nice list of all public tables in the default database
        >>> pg=pgnumpy.connect()
        >>> pg.describe()
             Public Tables
         Schema |  Name   | Owner
        --------+---------+------
         public | mytable | des
         public | test    | des

        # describe a table
        >>> pg.describe('test')
                      Table "public.test"
             Column      |       Type       | Modifiers
        -----------------+------------------+----------
         run             | integer          | not null
         rerun           | character(3)     |
         camcol          | integer          | not null
         field           | integer          | not null
         id              | integer          | not null
         objc_type       | integer          | not null
         objc_flags      | integer          | not null
         objc_flags2     | integer          | not null
         colc            | real[]           | not null
         petror50        | real[]           | not null
        ....snip...


        # print out some results of a query in a nice format.  show()
        # has lots of options, including paging the output
        >>> pg.show('select run,rerun,camcol,field,id,ra,dec from test')
         run  | rerun | camcol | field | id  |      ra       |      dec
        ------+-------+--------+-------+-----+---------------+---------------
         1010 | 137   |      1 |    69 | 500 | 22.7893646604 | -3.24415193812
         1010 | 137   |      1 |    70 |  93 | 22.9322016922 | -3.42616284608
         1010 | 137   |      1 |    70 |  95 | 22.9372923859 | -3.42525311698
         1010 | 137   |      1 |    70 | 106 | 22.8226045865 | -3.35162787925
         1010 | 137   |      1 |    70 | 172 | 22.8960652826 | -3.38928962863
         1010 | 137   |      1 |    70 | 314 | 22.8921524021 | -3.29063488036
         1010 | 137   |      1 |    70 | 341 | 22.9106492235 | -3.33561534067
         1010 | 137   |      1 |    70 | 395 | 22.9536925023 | -3.33270891541
         1010 | 137   |      1 |    70 | 484 |   22.88095547 | -3.22553133121
         1010 | 137   |      1 |    71 |  95 | 23.0386121293 | -3.31202356311


        # get a list of public tables
        >>> print pg.tables()
        ['mytable' 'test']

        # see if a table exists
        if pg.table_exists(tablename):
            # do something


        # If you have a structure in memory and you want to write it to a 
        # table, you can use the array2table convenience function.  See
        # docs for pgnumpy.array2table for more info.
        pgnumpy.array2table(array, tablename)
 
    """
    def execute(self, query):
        """
        Class:
            PgNumpy
        Name:
            execute

        Purpose:
            Execute a query.

        Limitations:
            Some types don't map onto numpy types well.  The most common cases
            are boolean and dates.  It is recommended you convert boolean to
            integer or a character type, and dates to varchar.  You can cast
            using the :: operator

            For example: 
                SELECT date_column::varchar FROM mytable 
                SELECT bool_column::smallint FROM mytable 

        Calling Sequence:
            pg=pgnumpy.connect()

            # execute and retrieve results
            pg.execute('select * from mytable')
            if pg.ntuples() > 0:
                res = pg.fetchall()
                # process results

            # These is equivalent
            pg.execute('select * from mytable')
            res = pg.fetchall()
            if res.size:
                # process results

            # These is also equivalent using the fetchall
            # shortcut
            res = pg.fetchall('select * from mytable')
            if res.size:
                # process results
        
        """
        cpgnumpy.cPgNumpy.execute(self,query)


    def fetchall(self, query=None):
        """
        This is an alias for fetch
        """
        if query is not None:
            self.execute(query)
        return cpgnumpy.cPgNumpy.fetchall(self)

    def fetch(self, query=None):
        """
        Class:
            PgNumpy
        Name:
            fetchall

        Purpose:
            Fetch all available results in a numerical python array with
            records.

            If the optional query is sent, execute the query and then
            retrieve the results.

        Calling Sequence:
            pg=pgnumpy.connect()

            # execute the query and fetch the results
            pg.execute('select * from mytable')
            res = pg.fetchall()

            # equivalent of above
            res = pg.fetchall('select * from mytable')

        Keywords:
            query:  An optional query to execute.  This is equivalent of
                using execute(query) and then fetchall().
        
        """
        
        if query is not None:
            self.execute(query)
        return cpgnumpy.cPgNumpy.fetchall(self)

    def show(self, query, fancy=True, **keys):
        """
        Class:
            PgNumpy
        Name:
            show

        Purpose:
            Execute a query and show the results, possibly paging. By
            default a fancy print format is used but this can be turned
            off.

            If fancy printint is turned off, more control can be had over
            the format.

        Calling Sequence:
            If fancy=True
                show(query, fancy=True, page=False, file=None, fields=ALL, nlines=ALL,
                     title=None)

            Additional keywords if fancy=False
                show(query, ..., sep=' ', format=, nformat=)

        Inputs:
            query: A query to execute.
        Keywords:
            delim: 
                Delimiter between fields.  Default ','
            array_delim: 
                Delimiter between array elements.  Default ','
            fancy: 
                If True, print with a visually appealing, but less machine readable
                format
            page: 
                If True, send results to a pager.
            file: 
                Send results to this file rather than standard output. 

            nlines: 
                Print only the top N lines.  Default is to print all.
            fields: 
                Only print a subset of the fields.

            title:
                A title to place above the printout during fancy printing.

            altnames: 
                A list of names for each argument.  There must be an entry for
                each argument. The names are printed above each column when
                doing fancy printing.


        Examples:
            pg=pgnumpy.connect()

            # execute the query and fetch the results
            pg.show('select * from mytable')

            pg.show('select * from mytable', page=True, fancy=False)
            
        
        """
 
        self.execute(query)
        res = self.fetchall()
        if res.size > 0:
            aw = ArrayWriter(fancy=fancy, **keys)
            aw.write(res, **keys)
            aw.close()

    def database(self):
        """
        Class:
            PgNumpy
        Name:
            database
        Purpose:
            Return the name of the current database.
        """
        res = self.fetchall('select current_database()')
        if res.size == 0:
            raise ValueError("Could not get database name!")
        return res['current_database'][0]


    def describe(self, table=None, page=False):
        """
        Class:
            PgNumpy
        Name:
            describe

        Purpose:
            Print a description of the indicated table, or if no table name is
            sent, list all public tables.

        Calling Sequence:
            pg=pgnumpy.connect()
            # List all public tables
            pg.describe()

            # Describe a table
            pg.describe(tablename)

            # send results to a pager
            pg.describe(tablename,page=True)

        Keywords:
            table:  A table name in the database.  This is optional.
            page: If True, send the output to a pager.
        
        """
        
        # give a list of public tables
        if table is None:
            self._describe_public_tables(page=page)
        else:
            self._describe_table(table, page=page)

    def _describe_public_tables(self, page=False):
            q = """
            select
                schemaname,
                tablename,
                tableowner
            from
                pg_tables
            where
                schemaname = 'public'
            """
            #self.set_field_lengths({'tablename':100,
            #                        'schemaname':25,
            #                        'tableowner':25})
            self.execute(q)
            res = self.fetchall()

            self.clear_field_lengths()

            if res.size == 0:
                stdout.write("No tables found\n")
                return
            res.dtype.names = ['Schema','Name','Owner']
            #aprint(res, fancy=True, title='Public Tables')
            aw=ArrayWriter(fancy=True,page=page)
            aw.write(res, title='Public Tables')

    def _describe_table(self, table, page=False):
        # get the oid
        oid_query="""
        SELECT 
            c.oid,
            n.nspname,
            c.relname
        FROM 
            pg_catalog.pg_class c
        LEFT JOIN 
            pg_catalog.pg_namespace n 
        ON 
            n.oid = c.relnamespace
        WHERE 
            c.relname ~ '^({table})$'
            AND pg_catalog.pg_table_is_visible(c.oid)
        ORDER BY 
            2, 3
        """.format(table=table)

        #self.set_field_lengths({'nspname':25,'relname':25})
        self.execute(oid_query)
        oid_res = self.fetchall()
        if oid_res.size == 0:
            stderr.write("No such table found: '%s'" % table)
        oid = oid_res['oid'][0]


        # This gets basic column info for columns in the table
        # with this oid
        col_query="""
        SELECT 
            a.attname,
            pg_catalog.format_type(a.atttypid, a.atttypmod),
                (SELECT 
                    substring(pg_catalog.pg_get_expr(d.adbin, d.adrelid) for 128)
                FROM 
                    pg_catalog.pg_attrdef d
                WHERE 
                    d.adrelid = a.attrelid AND d.adnum = a.attnum AND a.atthasdef),
            a.attnotnull::varchar, a.attnum
        FROM 
            pg_catalog.pg_attribute a
        WHERE 
            a.attrelid = {oid}
            AND a.attnum > 0 
            AND NOT a.attisdropped
        ORDER BY 
            a.attnum
        """.format(oid=oid)
        #self.set_field_lengths({'attname':25,'format_type':25})
        self.execute(col_query)
        t=self.fetchall()

        # also try to get a row, this will help with array 
        # descriptors
        example=self.fetchall('select * from %s limit 1' % table)
        have_example=False
        if example.size > 0:
            have_example=True


        # this is the output array we will print
        output = numpy.zeros(t.size, dtype=[('Column','S100'),
                                            ('Type','S100'),
                                            ('Modifiers','S25')])
        output['Column'] = t['attname']
        w,=where(t['attnotnull'] == 'true')
        if w.size > 0:
            output['Modifiers'][w] = 'not null'

        # Get type representation.  If we have an example, include
        # that to get actual dimensions
        for i in xrange(t.size):
            cname = t['attname'][i]
            ctype = t['format_type'][i]
            if have_example:
                # replace [] with the actual dimensions
                if ctype.find('[]') != -1:
                    # get array dims
                    shape=example[cname][0].shape
                    if len(shape) > 0:
                        sstr = [str(s) for s in shape]
                        #sstr = '['+ ']['.join(sstr)+']'
                        sstr = '['+ ','.join(sstr)+']'
                        ctype = ctype.replace('[]',sstr)
            ctype = ctype.replace('double precision','double')
            output['Type'][i] = ctype

        title = 'Table "{nspname}.{table}"'.format(nspname=oid_res['nspname'][0],
                                                   table=table)

        # get info about indexes
        index_query = """
        SELECT 
            c2.relname, 
            i.indisprimary, 
            i.indisunique, 
            i.indisclustered, 
            i.indisvalid, 
            pg_catalog.pg_get_indexdef(i.indexrelid, 0, true), 
            c2.reltablespace
        FROM 
            pg_catalog.pg_class c, 
            pg_catalog.pg_class c2, 
            pg_catalog.pg_index i
        WHERE 
            c.oid = {oid} 
            AND c.oid = i.indrelid 
            AND i.indexrelid = c2.oid
        ORDER BY 
            i.indisprimary DESC, i.indisunique DESC, c2.relname
        """.format(oid=oid)
        #self.set_field_lengths({'pg_get_indexdef':100})
        self.execute(index_query)
        index_res = self.fetchall()
        if index_res.size > 0:
            trailer = ['Indexes:']
            for i in xrange(index_res.size):
                tmp=index_res['pg_get_indexdef'][i]
                trailer.append("    %s" % tmp)
            trailer = '\n'.join(trailer)
        else:
            trailer=None

        aw = ArrayWriter(fancy=True,page=page)
        aw.write(output, title=title, trailer=trailer)




    def tables(self):
        """
        Class:
            PgNumpy
        Name:
            tables
        Purpose:
            Return a list of public table names.  For more detail, use the
            describe() function.  The result is a numerical python array.
        Example:
            >>> pg=pgnumpy.connect()
            >>> print pg.tables()
            >>> ['mytable' 'data' 'fields']
        """
        query="""
        select 
            tablename 
        from 
            pg_tables 
        where 
            schemaname = 'public'
        """
        #self.set_field_lengths({'tablename':100})
        self.execute(query)
        tb = self.fetchall()
        if tb.size == 0:
            return tb
        return tb['tablename']

    def table_exists(self, tablename):
        """
        Class:
            PgNumpy
        Name:
            table_exists

        Calling Sequence:
            pg=pgnumpy.connect()
            if pg.table_exists(tablename):
                # do something
        """
        tb=self.tables()
        if tb.size == 0:
            return False

        w,=where(tb == tablename)
        if w.size > 0:
            return True
        else:
            return False

    def current_queries(self):
        """
        Class:
            PgNumpy
        Name:
            current_queries
        Purpose:
            Show all the current running queries.  Note some 
            queries the user may not have permission to view.
        """
        # get all the current activity, not including this query
        # and IDLE queries
        query="""
        SELECT 
            datname,usename,current_query,waiting::varchar,query_start::varchar 
        FROM 
            pg_stat_activity 
        WHERE 
            current_query not like '%datname,usename,current_query,waiting%'
            AND current_query not like '%IDLE%'
        """
        res=self.fetch(query)
        aw = ArrayWriter(fancy=True)
        names = list(res.dtype.names)
        names[0] = 'dbname'
        names[2] = 'query'
        aw.write(res,altnames=names)

    def create_index(self, table, column_spec, unique=False):
        """
        Class:
            PgNumpy
        Name:
            create_index
        Calling Sequence:
            pg=pgnumpy.connect()
            pg.create_index(table,column_spec)

        Purpose:

            Create an index on the specified table and column.  This shortcut
            method automatically picks a name for the index.  The column spec
            can be a comma separated list for a multi-column index. Columns can
            also be array elements.

        Inputs:
            table: 
                The table name

            column_spec:  

                The column specifier.  This can be a name like 'colname' or a
                comma-separated list for multi-column indices.  Also you can
                specify array elements, e.g. 'colname[2][3]'.  Remember
                postgres databases by default use 1-offset array indices.


                For scalar columns, e.g 'colname', the index will be named:
                    tablename_colname_idx

                If multiple columns are being indexed, for example 'col1,col2':
                    tablename_col1_col2_idx

                If the column specification is for an array element, then this
                is incorporated into the index name.  For example if the
                request is for colname[3][4], then the index name is:

                    tablename_colname_i3_i4_idx

                Note this could cause a clash with an unlikely column named
                'colname_i3_i4'.  In those cases you should not use this
                function for creating the index.

                LIMITATION:  You cannot make a multi-column index that involves
                array elements.

        """

        if not self.table_exists(table):
            dbname=self.database()
            raise ValueError("table '%s' does not exist in "
                             "database '%s'" % (table,dbname))


        idxname = '{table}_{column_spec}_idx'.format(table=table,
                                                     column_spec=column_spec)

        # For the query, we put () around all names.  This protects array
        # subscripts.  

        csplit = column_spec.split(',')
        csplit = ['(%s)' % c for c in csplit]
        cstr = ','.join(csplit)

        if column_spec.find('[') != -1:
            idxname = idxname.replace('[','_i')
            idxname = idxname.replace(']','')


        # deal with multiple columns by replacing ',' with '_'
        idxname = idxname.split(',')
        idxname = '_'.join(idxname)

        if unique:
            unique_string = ' UNIQUE'
        else:
            unique_string = ''

        # if arrays specified, simply replace '[' with 'i' and remove ']'
        query="CREATE{unique} INDEX {idxname} on {table} ({column_spec})"
        query=query.format(idxname=idxname,
                           table=table,
                           column_spec=cstr,
                           unique=unique_string)

        stdout.write(query+'\n')
        self.execute(query)
        
        query = "ANALYZE %s" % table
        stdout.write(query+'\n')
        self.execute(query)


def npy2pg(npy_typename):

    if stringmatch.match(npy_typename):
        return _npy2pg['string']

    if npy_typename not in _npy2pg:
        raise ValueError('Unsupported type: %s' % npy_typename)
    return _npy2pg[npy_typename]



def make_bigendian(array, each=False, verbose=False):
    if each:
        # do each field separately
        for n in array.dtype.names:
            if numpy.little_endian and array[n].dtype.isnative:
                if verbose: 
                    sys.stdout.write('\t\t** Byte swapping field: %s\n' % n)
                # True means in place
                array[n].byteswap(True)
    else:
        if numpy.little_endian and array.dtype.isnative:
            # True means in place
            array.byteswap(True)







class PgInput:
    """
    Class:
        PgInput
    Purpose:
        Write an array out to a file that can be read through
        a COPY command into a postgrs table.

        It also creates the table if necessary, and checks the
        datatype on subsequent writes using the same object.

    Calling Sequence:
        pginput = PgInput(dtype, tablename, 
                          conninfo=None, 
                          primary_key=None, 
                          unique=None, 
                          serial_extra=None,
                          tmpdir=None, 
                          verbose=False)
        pginput.write(array)

    Limitations:
        Unsigned integers and 1-byte integers are not supported by postgres.
        Unsigned will be converted to signed, so be careful about ranges.
        1-byte are converted to 2-byte.

        Bool types and object types are not supported.

    """
    def __init__(self, dtype, tablename, 
                 conninfo=None, 
                 primary_key=None, 
                 unique=None, 
                 serial_extra=None,
                 tmpdir=None, 
                 verbose=False):

        self.set_verbosity(verbose)

        self.tablename = tablename
        self.set_dtype(dtype)
        self.conninfo=conninfo
        self.primary_key=primary_key
        self.unique=unique
        self.serial_extra=serial_extra

        self.set_coldefs()
        self.set_create_table_statement()

        if tmpdir is None:
            self.tmpdir = '.'
        else:
            tmpdir = os.path.expanduser(tmpdir)
            tmpdir = os.path.expandvars(tmpdir)
            self.tmpdir=tmpdir

    def write(self, array, keep=False):
        """
        Write an array with fields to the table
        """

        if array.dtype != self.dtype:
            raise ValueError("Array does not have the right dtype")

        # get a random filename
        self.fname = self.get_filename()

        # This will only create if the table doesn't already exist
        self.create_table()

        # this writes out the ascii file for input to the table
        self.write_input(array)

        # this actually runs the postgres copy command
        self.copy_to_table()

        if not keep:
            self.cleanup()
        else:
            stdout.write("===> Not removing temporary file: %s\n" % self.fname)

    def write_input(self, array):
        # Now write the file for input use recfile if we have it, much faster
        if have_recfile:
            if self.verbose:
                stdout.write("Writing ascii file using recfile: %s\n" % self.fname)
            r=recfile.Open(self.fname,'w',delim='\t',bracket_arrays=True)
            r.write(array)
            r.close()
        else:
            if self.verbose:
                stdout.write("Writing ascii file using ArrayWriter: %s\n" % self.fname)
            aw = ArrayWriter(self.fname,delim='\t',bracket_arrays=True)
            aw.write(array)
            aw.close()

    def create_table(self):
        # does the table exist?  If not, create it
        pg = connect(self.conninfo)
        if not pg.table_exists(self.tablename):
            self.print_create_table_statement()
            pg.execute(self.create_table_statement)
        pg.close()


    def cleanup(self):
        if self.verbose:
            stdout.write("Removing temporary file: %s\n" % self.fname)
        os.remove(self.fname)
        self.fname=None

    def get_filename(self):
        pref='pgcopy-'+self.tablename+'-'
        fname=tempfile.mktemp(prefix=pref,dir=self.tmpdir, suffix='.pgsql')
        fname=os.path.expanduser(fname)
        return fname


    def copy_query(self):
        """
        This requires admin privledges.

        If serial_extra not None, then we specify all the columns that will be
        copied not including the serial column
        """

        columns = ''
        if self.serial_extra is not None:
            # the serial extra columns were put at the beginning
            columns = '(%s)' % ','.join( self.names )

        if have_recfile:
            q="""
            COPY
                %s %s
            FROM
                '%s'
            """ % (self.tablename, columns, self.fname)
        else:
            q="""
            COPY BINARY 
                %s %s
            FROM
                '%s'
            """ % (self.tablename, columns, self.fname)

        return q


    def copy_to_table(self):
        """
        This requires admin privledges.

        If serial_extra not None, then we specify all the columns that will be
        copied not including the serial column
        """

        if self.verbose:
            stdout.write("Running COPY command\n")
        q = self.copy_query()
        if self.verbose >= 2:
            sys.stdout.write('%s\n' % q)
        conn=self.conninfo
        if conn is None:
            conn = adminconn_def 
        else:
            # libpq is smart and will take the last keyword-value pair when
            # there are duplicates
            conn = conn + ' ' + adminconn_def 
        
        conn = conn.strip()
        if self.verbose >= 2:
            sys.stdout.write("connection info = '%s'\n" % conn)
        
        pg = connect(conn)
        pg.execute(q)


    def set_verbosity(self, verbose):
        self.verbose=verbose
        # allow granularity in verbosity
        if self.verbose is True:
            self.verbose = 1
        elif self.verbose is False:
            self.verbose = 0


    def set_dtype(self, dtype):
        # This can be list, dict, or already a dtype
        self.dtype = numpy.dtype(dtype)
        self.set_field_info()


    def calc_array_totlen(self, ndim, nel, nbytes):
        """
        this is the length of the output since each row contains info
        about sizes, lengths, etc
        """
        # this is for python2 python3 compatability
        ta = numpy.zeros(5, dtype='i4')
        ta[0] = 2*4
        ta[1] = 4
        ta[2] = 2*ndim*4
        ta[3] = nel*4
        ta[4] = nel*nbytes
        return ta.sum()


    def set_field_info(self):
        ncolumns = len(self.dtype)
        maxdim = numpy.MAXDIMS
        self.maxdim = maxdim

        # attributes for output structure
        dt = self.dtype
        self.names = list(dt.names)
        self.nbytes = numpy.zeros(ncolumns, dtype='i4')
        self.ndim = numpy.zeros(ncolumns, dtype='i4')
        self.dims = numpy.zeros((ncolumns, maxdim), dtype='i4')

        self.nel = numpy.zeros(ncolumns, dtype='i4')
        self.tname = []

        for i in range(ncolumns):

            nbytes = dt[i].base.itemsize
            self.nbytes[i] = nbytes

            ndim = len(dt[i].shape)
            self.ndim[i] = ndim
            if ndim > 0:
                self.dims[i][0:ndim] = dt[i].shape
                w,=where(self.dims[i] > 0)
                if w.size > 0:
                    nel = self.dims[i][w].prod()
            else:
                nel = 1

            self.nel[i] = nel
            tname = dt[i].base.name
            self.tname.append(tname)


            if self.verbose >= 2:
                sys.stdout.write("field='"+dt.names[i]+"'\n")
                sys.stdout.write('\tnbytes= %s\n' % nbytes)
                sys.stdout.write('\tndim= %s\n' % ndim)
                sys.stdout.write('\tdims= %s\n' % self.dims[i])
                sys.stdout.write('\tnel= %s\n' % self.nel[i])
                sys.stdout.write('\ttotlen= %s\n' % self.totlen[i])



    def get_dimstring(self,i):
        if self.ndim[i] == 0:
            return ''

        dims = self.dims[i,:].copy()
        w,=where(dims > 0)
        
        dims = list(dims[w])
        dimstring=''
        for dim in dims:
            dimstring += '[%d]' % dim
        return dimstring


    def set_coldefs(self):
        self.coldefs = []
        for i in range(len(self.dtype)):
            # get the postgres type name
            tname = self.tname[i]
            pgtype = npy2pg(tname)
            if pgtype == 'varchar' or pgtype == 'character':
                pgtype += '(%d)' % self.nbytes[i]

            # dimensions
            pgtype += self.get_dimstring(i)
            extra=''
            if not stringmatch.match(tname):
                extra = ' NOT NULL'
            coldef = '%s %s%s' % (self.names[i], pgtype, extra) 
            self.coldefs.append(coldef)
            if self.verbose >= 2:
                sys.stdout.write('\n\t=> %s\n' % coldef)


        if self.serial_extra is not None:
            sadd = ['%s SERIAL8' % s for s in self.serial_extra]
            self.serial_extra_coldefs = sadd
        else:
            self.serial_extra_coldefs = None

        if self.primary_key is not None:
            self.coldefs.append('PRIMARY KEY (%s)' % self.primary_key)

        if self.unique is not None:
            uadd = ['UNIQUE (%s)' % u for u in unique]
            self.coldefs = self.coldefs + uadd

    def print_create_table_statement(self):
        coldefs = self.coldefs
        if self.serial_extra_coldefs is not None:
            coldefs = self.serial_extra_coldefs + coldefs

        pcoldefs = ['          '+cd.lower() for cd in coldefs]
        pcoldefs = ',\n'.join(pcoldefs)
        printable = """
        CREATE TABLE 
            %s 
        (
%s
        )\n""" % (self.tablename, pcoldefs)
        sys.stdout.write(printable)

    def set_create_table_statement(self):
        coldefs = self.coldefs
        if self.serial_extra_coldefs is not None:
            coldefs = self.serial_extra_coldefs + coldefs
        coldefs = ',\n '.join(coldefs)
        tabledef = \
            'create table %s (%s)' % (self.tablename, coldefs)
        self.create_table_statement = tabledef
            


def compare_arrays(arr1, arr2, verbose=False):
    nfail = 0
    for n2 in arr2.dtype.names:
        n1 = n2
        if n1 not in arr1.dtype.names:
            n1 = n1.lower()
            if n1 not in arr1.dtype.names:
                n1 = n1.upper()
                if n1 not in arr1.dtype.names:
                    raise ValueError('field name %s not found in array 1' % n2)
            
        if verbose:
            sys.stdout.write("    testing field: '%s'\n" % n2)
            sys.stdout.write('        shape...........')
        if arr2[n2].shape != arr1[n1].shape:
            nfail += 1
            if verbose:
                sys.stdout.write('shapes differ\n')
        else:
            if verbose:
                sys.stdout.write('OK\n')
                sys.stdout.write('        elements........')
            w,=where(arr1[n1].ravel() != arr2[n2].ravel())
            if w.size > 0:
                nfail += 1
                if verbose:
                    sys.stdout.write('\n        Error: %s elements in field %s differ\n' % \
                                     (w.size,n2))
            else:
                if verbose:
                    sys.stdout.write('OK\n')

    if nfail == 0:
        if verbose:
            sys.stdout.write('All tests passed\n')
        return True
    else:
        if verbose:
            sys.stdout.write('%d errors found\n' % nfail)
        return False


def test_simple(verbose=False):
    """
    This is designed to test all the basic datatypes in scalar and array
    forms
    """
    dt = [('i4','i4'),('f4','f4')]

    narr = 2
    arr=numpy.zeros(narr, dtype=dt)

    # initialize fields
    # f4 fields
    arr['i4'] = [3,25]
    arr['f4'] = [1.2, 66.3]

    table = 'test_npy_input_output'

    if table_exists(table):
        sys.stdout.write('removing existing test table %s\n' % table)
        pgqueryln('drop table %s' % table)
    sys.stdout.write("Testing database input to table '%s'\n" % table)
    array2table(arr, table, verbose=verbose, keep=True)

    sys.stdout.write("\nTesting database extraction from table '%s'\n" % table)
    res = pgquery('select * from %s' % table)

    sys.stdout.write('Testing that input and output match\n')
    compare_arrays(arr, res, verbose=True)
    
    sys.stdout.write("removing temporary table '%s'\n" % table)
    pgqueryln('truncate %s' % table)
    pgqueryln('vacuum full %s' % table)
    pgqueryln('drop table %s' % table)



def test(verbose=False, keep=False, pad_strings=True):
    """
    This is designed to test all the basic datatypes in scalar and array
    forms
    """
    
    # use postgres user to avoid any permissions issues
    conninfo='user=postgres'

    f4array_shape = (3,2,2)
    f8array_shape = (2,)
    strarray_shape = (2,2)
    i2array_shape = (3,)
    i4array_shape = (2,3)
    i8array_shape = (1,2)
    dt = [('f4','f4'), ('f4array','f4',f4array_shape),
          ('f8','f8'), ('f8array','f8',f8array_shape),
          ('str','S5'),('strarray','S5',strarray_shape),
          ('i2','i2'), ('i2array','i2',i2array_shape),
          ('i4','i4'), ('i4array','i4',i4array_shape),
          ('i8','i8'), ('i8array','i8',i8array_shape)]

    narr = 5
    arr=numpy.zeros(narr, dtype=dt)

    # initialize fields
    # f4 fields
    arr['f4'] = numpy.arange(len(arr))
    for i in range(narr):
        arr[i]['f4array'][:] = numpy.random.random(f4array_shape)[:]

    # f8 fields
    arr['f8'] = numpy.arange(narr)
    for i in range(narr):
        arr[i]['f8array'][:] = numpy.random.random(f8array_shape)[:]

   
    # string fields
    # no nul characters are allowed.  array2table will pad the strings
    # with spaces automatically, so we add them here by hand in order
    # for the tests to succeed
    if pad_strings:
        arr['str'][0] = 'hello'
        arr['str'][1] = 'world'
        arr['str'][2] = 'test '
        arr['str'][3] = 'str  '
        arr['str'][4] = 'field'
    else:
        arr['str'][0] = 'hello'
        arr['str'][1] = 'world'
        arr['str'][2] = 'test'
        arr['str'][3] = 'str'
        arr['str'][4] = 'field'


    for i in range(narr):
        arr[i]['strarray'][:] = numpy.random.random(strarray_shape)[:]
    

    # i2 fields
    arr['i2'] = numpy.arange(narr)
    for i in range(narr):
        arr[i]['i2array'][:] = \
                numpy.random.randint(-32000,32000,arr[i]['i2array'].size)

    # i4 fields
    arr['i4'] = numpy.arange(narr)
    for i in range(narr):
        arr[i]['i4array'][:] = \
                numpy.random.random(i4array_shape)*100000

    # i8 fields
    arr['i8'] = numpy.arange(narr)
    for i in range(narr):
        arr[i]['i8array'][:] = \
                numpy.random.random(i8array_shape)*100000



    table = 'test_npy_input_output'

    if table_exists(table):
        sys.stdout.write('removing existing test table %s\n' % table)
        pgqueryln('truncate %s' % table, conninfo=conninfo)
        pgqueryln('vacuum full %s' % table, conninfo=conninfo)
        pgqueryln('drop table %s' % table, conninfo=conninfo)
    sys.stdout.write("Testing database input to table '%s'\n" % table)
    array2table(arr, table, verbose=verbose, keep=keep, conninfo=conninfo,
                pad_strings=pad_strings)

    sys.stdout.write("Testing database extraction from table '%s'\n" % table)
    res = pgquery('select * from %s' % table, conninfo=conninfo)

    sys.stdout.write('Testing that input and output match\n')
    compare_arrays(arr, res, verbose=True)
    
    if not keep:
        sys.stdout.write("\nremoving temporary table '%s'\n" % table)
        pgqueryln('truncate %s' % table, conninfo=conninfo)
        pgqueryln('vacuum full %s' % table, conninfo=conninfo)
        pgqueryln('drop table %s' % table, conninfo=conninfo)


def aprint(array, **keys):
    """
    Name:
      aprint

        This is a hacked version of esutil.numpy_util.aprint to avoid
        further dependencies.
    Purpose:

        Print out the rows and columns of the array with fields, aka structure
        or records.  The focus is on visualizing the data rather than speed or
        efficient file output.
        
        Optionally results can be sent to a pager or file.  
        
        By default, the columns are printed in a simple, machine readable
        format, but if fancy=True a more visually appealing format is used.

        Subsets of the fields can be chosen.  Also, printing can be 
        restricted to the top N lines.

        If not fancy printing, the user has more control over the format.

    Calling Sequence:
        aprint(array, fancy=False, page=False, nlines=ALL, fields=ALL, 
               file=None)

        For fancy=True, you can also specify a title.

        If fancy is not True, more keywords are available
            names=use_field_names, sep=' ', format=, names=field_names, nformat=None
    
    Inputs:
        array: A numpy array with fields.

    Keywords:
        fancy: 
            If True, print with a visually appealing, but less machine readable
            format
        page: 
            If True, send results to a pager.
        file: 
            Send results to this file rather than standard output. 
        nlines: 
            Print only the top N lines.  Default is to print all.
        fields: 
            Only print a subset of the fields.

        # keywords available only when fancy=True
        title:
            A title to place above the printout.

        # keywords available when fancy=False
        format: 
            A format string to apply to every argument.  E.g. format='%15s'
            Since every arg gets the same format, only %s type formats should
            be used unless the types are homogeneous.
        sep: 
            Separator, default is ' '
        names: 
            A list of names for each argument.  There must be an entry for
            each argument. The names are printed above each column.
        nformat:
            A Format to apply to the names.  By default, the same format used
            for the arguments is tried.  If formatting fails, a simple '%s' is
            used for the names.


    Example:



    """

    # if fancy, use the fancy array printer
    fancy = keys.get('fancy', False)

    if not fancy:
        # use colprint
        fields = keys.get('fields', array.dtype.names)
        if 'names' not in keys:
            keys['names'] = fields
        names = keys.get('names',fields)
        if len(names) != len(fields):
            raise ValueError("names= keyword must be same length as number of fields")

        ftup = split_fields(array, fields=fields)

        #colprint = esutil.misc.colprint

        arglist=[]
        for i in range(len(fields)):
            arglist.append('ftup[%s]' % i)

        arglist = ', '.join(arglist)
        command = 'colprint('+arglist+', **keys)'
        eval(command)

    else:
        # fancy printing, more memory usage and not good for machine reading
        #center_text = esutil.misc.center_text
        title = keys.get('title', None)
        page = keys.get('page',False)

        if not page:
            # should we print to a file?
            f = keys.get('file', stdout)
            if isinstance(f, file):
                fobj = f
            else:
                fobj = open(f,'w')

        # if we are paging, we will store the lines, otherwise this won't be used
        lines = []

        fields = keys.get('fields', array.dtype.names)
        nlines = keys.get('nlines', array.size)

        max_lens = {}
        for name in fields:
            max_lens[name] = len(name)

        # first pass through data to get lengths
        for i in xrange(nlines):
            for name in fields:
                max_lens[name] = max(max_lens[name], len(str(array[name][i])))

        forms = {}
        separator = ''
        i=0
        ntot = len(fields)
        for name in fields:
            if isinstance(array[name][0], numpy.string_):
                forms[name]    = ' %-'+str(max_lens[name])+'s '
            else:
                forms[name]    = ' %'+str(max_lens[name])+'s '

            pad = 2
            if i == (ntot-1):
                pad=1
            this_sep = '%s' % '-'*(max_lens[name]+pad)

            if i > 0:
                forms[name] = '|' + forms[name]
                this_sep = '+' + this_sep
            separator += this_sep
            i+=1

        header = ''
        for n in fields:
            header += forms[n] % center_text(n,max_lens[n])

        if title is not None:
            title = center_text(title, len(header))

        if page:
            if title is not None:
                lines.append(title)
            lines.append(header)
            lines.append(separator)

        else:
            if title is not None:
                fobj.write(title)
                fobj.write('\n')

            fobj.write(header)
            fobj.write('\n')

            fobj.write(separator)
            fobj.write('\n')

        for i in xrange(nlines):
            line = ''
            for name in fields:
                if page:
                    line += forms[name] % array[name][i]
                else:
                    fobj.write(forms[name] % array[name][i])
            if page:
                lines.append(line)
            else:
                fobj.write('\n')

        if page:
            lines = '\n'.join(lines)
            pydoc.pager(lines)
        else:
            # close if this is not stdout
            if fobj != stdout:
                fobj.close()


def colprint(*args, **keys):
    """
    Name:
        colprint
    Purpose:
        print the input sequences or arrays in columns.  All must be the
        same length.
    Calling Sequence:
        colprint(var1, var2, ..., nlines=all, sep=' ', format=None, 
                 names=None, nformat=None, file=None, page=False)

    Inputs:
        A set of python objects.  Each must be a sequence or array and all must
        be the same length.

    Optional Inputs:
        nlines:  
            Number of lines to print.  Default is all.
        sep: 
            Separator, default is ' '
        file:  
            A file path or file object.  Default is to print to standard
            output. Ignored if paging.

        format: 
            A format string to apply to every argument.  E.g. format='%15s'
            Since every arg gets the same format, only %s type formats should
            be used unless the types are homogeneous.

        names: 
            A list of names for each argument.  There must be an entry for
            each argument. The names are printed above each column.
        nformat:
            A Format to apply to the names.  By default, the same format used
            for the arguments is tried.  If formatting fails, a simple '%s' is
            used for the names.

        page: If True, run the output through a pager.

    Revision History:
        Create: 2010-04-05, Erin Sheldon, BNL
    """
    nargs = len(args)
    if nargs == 0:
        return

    try:
        n1 = len(args[0])
    except:
        raise ValueError("Could not get len() of argument 1")

    # Should we print only a subset?
    nlines = keys.get('nlines',n1)
    if nlines is None:
        nlines = n1
    elif nlines > n1:
        nlines = n1

    # what separator should be used?
    sep = keys.get('sep',' ')

    # should we page the results?
    page = keys.get('page', False)

    if not page:
        # should we print to a file?
        f = keys.get('file', stdout)
        if isinstance(f, file):
            fobj = f
        else:
            fobj = open(f,'w')

    # make sure all the arguments are the same length.
    for i in range(nargs):
        try:
            l=len(args[i])
        except:
            raise ValueError("Could not get len() of argument %s" % (i+1))
        if l != n1:
            e="argument %s has non-matching length.  %s instead of %s" \
                    % (i+1, l, n1)
            raise ValueError(e)

    # if we are paging, we will store the lines, otherwise this won't be used
    lines = []

    # print a header
    names = keys.get('names',None)
    if names is not None:
        if isinstance(names, basestring):
            names = [names]
        nnames = len(names)
        if len(names) != nargs:
            raise ValueError("Expected %s names, got %s" % (nargs,nnames))
        
        # see if explicit format has been requested.
        nformat =keys.get('nformat',None)

        if nformat is not None:
            nformat = [nformat]*nnames
        else:
            # try to use the other format
            fmt=keys.get('format','%s')
            if fmt is None:
                fmt='%s'
            nformat=[fmt]*nnames

        nformat = sep.join(nformat)
        try:
            line = nformat % tuple(names)
        except:
            nformat = ['%s']*nnames
            nformat = sep.join(nformat)
            line = nformat % tuple(names)

        if page:
            lines.append(line)
        else:
            fobj.write(line)
            fobj.write('\n')


    # format for columns.  Same is used for all.
    format = keys.get('format','%s')
    if format is not None:
        format = [format]*nargs
    else:
        format = ['%s']*nargs

    format = sep.join(format)

    # loop over and print columns
    for i in range(nlines):
        data = []
        for iarg in range(nargs):
            data.append(args[iarg][i])
        
        data = tuple(data)

        line = format % data
        line = line.replace('\n','')

        if page:
            lines.append(line)
        else:
            fobj.write(line)
            fobj.write('\n')

    if page:
        lines = '\n'.join(lines)
        pydoc.pager(lines)
    else:
        # close if this is not stdout
        if fobj != stdout:
            fobj.close()


def center_text(text, width):
    text = text.strip()
    space = width - len(text)
    return ' '*(space/2) + text + ' '*(space/2 + space%2)

def split_fields(data, fields=None, getnames=False):
    """
    Name:
        split_fields

    Calling Sequence:
        The standard calling sequence is:
            field_tuple = split_fields(data, fields=)
            f1,f2,f3,.. = split_fields(data, fields=)

        You can also return a list of the extracted names
            field_tuple, names = split_fields(data, fields=, getnames=True)

    Purpose:
        Get a tuple of references to the individual fields in a structured
        array (aka recarray).  If fields= is sent, just return those
        fields.  If getnames=True, return a tuple of the names extracted
        also.

        If you want to extract a set of fields into a new structured array
        by copying the data, see esutil.numpy_util.extract_fields

    Inputs:
        data: An array with fields.  Can be a normal numpy array with fields
            or the recarray or another subclass.
    Optional Inputs:
        fields: A list of fields to extract. Default is to extract all.
        getnames:  If True, return a tuple of (field_tuple, names)

    """

    outlist = []
    allfields = data.dtype.fields

    if allfields is None:
        if fields is not None:
            raise ValueError("Could not extract fields: data has "
                             "no fields")
        return (data,)
    
    if fields is None:
        fields = allfields
    else:
        if isinstance(fields, (str,unicode)):
            fields=[fields]

    for field in fields:
        if field not in allfields:
            raise ValueError("Field not found: '%s'" % field)
        outlist.append( data[field] )

    output = tuple(outlist)
    if getnames:
        return output, fields
    else:
        return output


class ArrayWriter:
    """
    Class:
        ArrayWriter
    Purpose:

        A python only class to write numpy arrays.  Can write arrays with
        fields in columns.  Can also do a "fancy" print of the array which
        is easy on the eyes but not good for machine reading.

        This is much slower than using the recfile package, but as it
        is python only it is more flexible.

    Constructor:
        aw = ArrayWriter(file=None, 
                         delim=',', 
                         array_delim=',',
                         bracket_arrays=False,
                         fancy=False,
                         page=False)

    Inputs:
        file:
            File name or file object to use for printing.
        delim: 
            The delimiter between fields.
        array_delim: 
            The delimiter between array elements.
        bracket_arrays:  
            
            Put brackets in place to delintate dimensional boundies.  e.g.
            {{a,b,c},{d,e,f}}

            Note if fancy=True, brackets are always used.

        fancy:  
            If True, use a fancy printing scheme instead of the simple machine
            readable format that is default.
        page:
            If True, send the output to a pager.

    Examples:
        aw = ArrayWriter()
        aw.write(arr)

        # if fancy=True, can sent a title to place above the printout
        aw = ArrayWriter(delim='\t', fancy=True)
        aw.write(array, title='My Table')


    """

    def __init__(self, **keys):
        self.open(**keys)

    def open(self,**keys):

        self._delim = keys.get('delim',',')
        self._array_delim = keys.get('array_delim',',')
        self._fancy=keys.get('fancy',False)

        # we only will page if fancy is set
        if self._fancy:
            self._page=keys.get('page',False)
            self._bracket_arrays=True
        else:
            self._page=False
            self._bracket_arrays = keys.get('bracket_arrays',False)

        self._fobj=None
        self._close_the_fobj = False

        # Only load a file object if page is False
        # which in turn can only be true if fancy
        # is also True
        if not self._page:
            fobj = keys.get('file',stdout)

            if isinstance(fobj,file):
                self._fobj = fobj
            else:
                self._close_the_fobj=True
                fname = os.path.expanduser(fobj)
                fname = os.path.expandvars(fname)
                self._fobj = open(fname,'w')


    def write(self, arrin, **keys):
        """
        Class:
            ArrayWriter
        Name:
            write
        Calling Sequence:
            aw=ArrayWriter(**keywords)
            aw.write(array, **keywords)
        Purpose:
            Write an array.
        
        Inputs:
            array: 
                The array to write
            nlines: 
                The number of lines to write
            fields: 
                Only print a subset of the fields.

            title:
                A title to place above the printout during fancy printing.
            trailer:
                Text to print after the array data.

            altnames: 
                A list of names for each argument.  There must be an entry for
                each argument. The names are printed above each column when
                doing fancy printing.


        """
        if self._fancy:
            self.fancy_write(arrin, **keys)
            return

        arr=arrin.view(numpy.ndarray)

        nlines = keys.get('nlines', arr.size)
        names = keys.get('fields', arr.dtype.names)
        if names is None:
            # simple arrays are easy
            arr.tofile(self._fobj, sep='\n')
            return
        
        nnames = len(names)

        # we have fields
        for i in xrange(arr.size):
            iname=0
            for n in names:
                data = arr[n][i]
                if data.ndim > 0:
                    self.write_array(data)
                else:
                    self._fobj.write(str(data))

                if iname < (nnames-1):
                    self._fobj.write(self._delim)
                else:
                    self._fobj.write('\n')
                iname += 1

        
    def fancy_write(self, arrin, **keys):
        array=arrin.view(numpy.ndarray)

        title = keys.get('title', None)

        # if we are paging, we will store the lines, otherwise this won't be used
        lines = []

        fields = keys.get('fields', array.dtype.names)
        printnames = keys.get('altnames', fields)

        if len(fields) != len(printnames):
            raise ValueError("altnames must correspond directly to the fields "
                             "being printed")

        nlines = keys.get('nlines', array.size)
 
        max_lens = {}
        for name in fields:
            max_lens[name] = len(name)


        # first pass through data to get lengths

        # for array fields
        astr = ArrayStringifier(delim=self._array_delim,
                                brackets=self._bracket_arrays)
        for i in xrange(nlines):
            for name in fields:
                if array[name][i].ndim > 0:
                    strval=astr.stringify(array[name][i])
                    max_lens[name] = max(max_lens[name], len(strval))
                else:
                    max_lens[name] = max(max_lens[name], len(str(array[name][i])))

        # now create the forms for writing each field
        forms = {}
        separator = ''
        i=0
        ntot = len(fields)
        for name in fields:
            if isinstance(array[name][0], numpy.string_)  \
                    or (array[name][0].ndim > 0):
                forms[name]    = ' %-'+str(max_lens[name])+'s '
            else:
                forms[name]    = ' %'+str(max_lens[name])+'s '

            pad = 2
            if i == (ntot-1):
                pad=1
            this_sep = '%s' % '-'*(max_lens[name]+pad)

            if i > 0:
                forms[name] = '|' + forms[name]
                this_sep = '+' + this_sep
            separator += this_sep
            i+=1

        # possible header and title
        header = ''
        for i in xrange(len(fields)): 
            n=fields[i]
            pname=printnames[i]
            header += forms[n] % center_text(pname,max_lens[n])

        if title is not None:
            title = center_text(title, len(header))

        if self._page:
            if title is not None:
                lines.append(title)
            lines.append(header)
            lines.append(separator)

        else:
            if title is not None:
                self._fobj.write(title)
                self._fobj.write('\n')

            self._fobj.write(header)
            self._fobj.write('\n')

            self._fobj.write(separator)
            self._fobj.write('\n')

        for i in xrange(nlines):
            line = ''
            for name in fields:
                val = array[name][i]

                if val.ndim > 0:
                    val = astr.stringify(val)

                if self._page:
                    line += forms[name] % val
                else:
                    self._fobj.write(forms[name] % val)

            if self._page:
                lines.append(line)
            else:
                self._fobj.write('\n')

        trailer=keys.get('trailer',None)
        if trailer is not None:
            if self._page:
                lines.append(trailer)
            else:
                self._fobj.write(trailer)
                self._fobj.write('\n')

        if self._page:
            lines = '\n'.join(lines)
            pydoc.pager(lines)



    def write_array(self, arr):
        """
        Write a simple array, possibly with brackets indicating the dimensions.
        """
        if self._bracket_arrays:
            self._fobj.write('{')

        i=0

        dimsize=arr.shape[0]

        for a in arr:
            if a.ndim > 0:
                self.write_array(a)
            else:
                self._fobj.write(str(a))
        
            if i < (dimsize-1):
                self._fobj.write(',')
            i+=1

        if self._bracket_arrays:
            self._fobj.write('}')


    def close(self):
        if self._close_the_fobj:
            self._fobj.close()

    def __del__(self):
        if self._close_the_fobj:
            self._fobj.close()

class ArrayStringifier:
    """
    Stringify a simple array using a delimiter and
    possibly brackets
    """
    def __init__(self, delim=',', brackets=False):
        self._delim=delim
        self._brackets=brackets
        self._values=[]

    def stringify(self, arr):
        self._values=[]
        if arr.dtype.names is not None:
            raise ValueError("array must be simple, not structured")
        self._process(arr)
        return ''.join(self._values)
    
    def _process(self, arr):

        if self._brackets:
            self._values.append('{')

        i=0

        dimsize=arr.shape[0]

        for a in arr:
            if a.ndim > 0:
                self._process(a)
            else:
                self._values.append(str(a))
        
            if i < (dimsize-1):
                self._values.append(self._delim)
            i+=1

        if self._brackets:
            self._values.append('}')




# this is the old one with all the binary stuff in it
class PgInputWithBinary:
    """
    input = PgInput(tablename=None, dtype=None,
                    conninfo='user=postgres', verbose=False)

    pgi = PgInput()
    pgi.array2table(array, tablename, tmpdir=None)
    # see also the array2table convenience function

    """
    def __init__(self, tablename=None, dtype=None, 
                 conninfo=None, 
                 pad_strings=True,
                 verbose=False):
        if dtype is not None:
            self.set_dtype(dtype)

        self.verbose=verbose
        if self.verbose is True:
            self.verbose = 1
        elif self.verbose is False:
            self.verbose = 0

        self.tablename = tablename
        self.pad_strings = pad_strings
        self.conninfo=conninfo

    def array2table(self, array_in, tablename, tmpdir=None, keep=False,
                    primary_key=None, unique=None, serial_extra=None):
        array = array_in.view(numpy.ndarray)

        pref='pgcopy-'+tablename+'-'
        self.fname=tempfile.mktemp(prefix=pref,dir=tmpdir, suffix='.pgsql')
        self.fname=os.path.expanduser(self.fname)

        self.tablename = tablename
        self.serial_extra=serial_extra
        self.set_dtype(array.dtype)

        self.set_coldefs(primary_key=primary_key, unique=unique)
        self.set_create_table_statement()

        # does the table exist?  If not, create it
        pg = connect(self.conninfo)
        if not pg.table_exists(tablename):
            self.print_create_table_statement()
            pg.execute(self.create_table_statement)
        pg.close()

        # Now write the file for input use recfile if we have it, much faster
        # and simpler

        if have_recfile:
            if self.verbose:
                stdout.write("Writing ascii file: %s\n" % self.fname)
            r=recfile.Open(self.fname,'w',delim='\t',bracket_arrays=True)
            r.write(array)
            r.close()
        else:
            if self.verbose:
                stdout.write("Writing ascii file using ArrayWriter: %s\n" % self.fname)
            aw = ArrayWriter(self.fname,delim='\t',bracket_arrays=True)
            aw.write(array)
            aw.close()

            # this is the old binary
            crap="""
            self.fname = self.fname
            self.open_for_writing()
            self.write_header()
            self.write_data(array)
            self.write_trailer()
            self.close()
            """

        self.copy_to_table()

        if self.verbose:
            stdout.write("Removing temporary file: %s\n" % self.fname)
        os.remove(self.fname)
        self.fname=None

    def copy_query(self):
        """
        This requires admin privledges.

        If serial_extra not None, then we specify all the columns that will be
        copied not including the serial column
        """

        columns = ''
        if self.serial_extra is not None:
            # the serial extra columns were put at the beginning
            columns = '(%s)' % ','.join( self.names )

        if have_recfile:
            q="""
            COPY
                %s %s
            FROM
                '%s'
            """ % (self.tablename, columns, self.fname)
        else:
            q="""
            COPY BINARY 
                %s %s
            FROM
                '%s'
            """ % (self.tablename, columns, self.fname)

        return q


    def copy_to_table(self):
        """
        This requires admin privledges.

        If serial_extra not None, then we specify all the columns that will be
        copied not including the serial column
        """

        if self.verbose:
            stdout.write("Running COPY command\n")
        q = self.copy_query()
        if self.verbose >= 2:
            sys.stdout.write('%s\n' % q)
        conn=self.conninfo
        if conn is None:
            conn = adminconn_def 
        else:
            # libpq is smart and will take the last keyword-value pair when
            # there are duplicates
            conn = conn + ' ' + adminconn_def 
        
        conn = conn.strip()
        if self.verbose >= 2:
            sys.stdout.write("connection info = '%s'\n" % conn)
        
        pg = connect(conn)
        pg.execute(q)


    def set_dtype(self, dtype):
        # This can be list, dict, or already a dtype
        self.dtype = numpy.dtype(dtype)
        self.set_field_info()
        self.set_output_dtype()

        # we could implement these in the future
        self.flags = numpy.array([0], dtype='>i4')
        self.hdext_len = numpy.array([0], dtype='>i4')
        self.trailer = numpy.array([-1], dtype='>i2')
        

    def open_for_writing(self):
        if self.fname is None:
            raise ValueError("Set self.fname")
        self.fobj = open(self.fname, 'w')

    def write_header(self):
        if self.verbose:
            stdout.write('writing header\n')
        self.write_sig()
        self.write_flags()
        self.write_hdext_len()

    def close(self):
        self.fobj.close()

    def write_sig(self):
        """
        Signature is always the same, might change in future?
        """
        self.fobj.write('PGCOPY\n\377\r\n\0')

    def write_flags(self):
        """
        just zero big endian if no oids are in the file
        """
        self.flags.tofile(self.fobj)

    def write_hdext_len(self):
        """
        Zero big endian if not an extended header
        """
        self.hdext_len.tofile(self.fobj)

    def write_trailer(self):
        """
        The number -1 as 16 bit big endian
        """
        if self.verbose:
            stdout.write('writing trailer\n')
        self.trailer.tofile(self.fobj)


    def write_data(self, array, mbperwrite=100):
        """
        This creates a temporary structure that contains field lengths 
        embedded, which matches the postgres output file structure
        """
    
        n=len(array)
        if self.verbose:
            stdout.write("creating output\n")
        output = self.init_output(n)
        if self.verbose:
            stdout.write("copying to output\n")
        self.copy_to_output(array, output)
        if self.verbose >= 2:
            sys.stdout.write('Input dtype:\n')
            for el in array.dtype.descr:
                sys.stdout.write('%s\n' % str(el))
            sys.stdout.write('Output dtype:\n')
            for el in output.dtype.descr:
                sys.stdout.write('%s\n' % str(el))
        make_bigendian(output, each=True, verbose=(self.verbose >= 2))
        if self.verbose:
            stdout.write("writing data: %s\n" % self.fobj.name)
        output.tofile(self.fobj)


    def calc_array_totlen(self, ndim, nel, nbytes):
        """
        this is the length of the output since each row contains info
        about sizes, lengths, etc
        """
        # this is for python2 python3 compatability
        ta = numpy.zeros(5, dtype='i4')
        ta[0] = 2*4
        ta[1] = 4
        ta[2] = 2*ndim*4
        ta[3] = nel*4
        ta[4] = nel*nbytes
        return ta.sum()

    def get_array_sig(self, typename):
        sig = numpy.zeros(4, dtype='u1')
        # don't support numpy int8 or uint8 types.
        if typename == 'int8' or typename == 'uint8':
            raise ValueError("1-byte integers not supported for binary input")
            # convert to 2 byte integer
            sig[:] = [0,0,0,21]
        if typename == 'int16' or typename == 'uint16':
            sig[:] = [0,0,0,21]
        elif typename == 'int32' or typename == 'uint32':
            sig[:] = [0,0,0,23]
        elif typename == 'int64' or typename == 'uint64':
            sig[:] = [0,0,0,20]
        elif typename == 'float32':
            sig[:] = [0,0,2,188]
        elif typename == 'float64':
            sig[:] = [0,0,2,189]
        elif stringmatch.match(typename):
            # This is for character(n) fields
            if _npy2pg['string'] == 'character':
                sig[:] = [0,0,4,18]
            elif _npy2pg['string'] == 'varchar':
                # This is for varchar(n) fields
                sig[:] = [0,0,4,19]
            elif _npy2pg['string'] == 'text':
                # This is for text fields with no length.  Not recommended
                # since numpy has fixed length string fields
                sig[:] = [0,0,0,25]
        else:
            raise ValueError('Unsupported type: '+typename)

        return sig

    def set_field_info(self):
        ncolumns = len(self.dtype)
        maxdim = numpy.MAXDIMS
        self.maxdim = maxdim

        # attributes for output structure
        dt = self.dtype
        self.names = list(dt.names)
        self.ncolumns = ncolumns
        self.nbytes = numpy.zeros(ncolumns, dtype='i4')
        self.ndim = numpy.zeros(ncolumns, dtype='i4')
        self.dims = numpy.zeros((ncolumns, maxdim), dtype='i4')
        self.shapes = []
        self.nel = numpy.zeros(ncolumns, dtype='i4')
        self.totlen = numpy.zeros(ncolumns, dtype='i4')
        #self.type = numpy.zeros(ncolumns, dtype='i2')
        self.tname = []
        self.tdescr = []
        self.arrsig = numpy.zeros((ncolumns, 4), dtype='u1')
        self.bytesperrow = 0 

        for i in range(ncolumns):

            nbytes = dt[i].base.itemsize
            self.nbytes[i] = nbytes

            self.shapes.append( dt[i].shape )

            ndim = len(dt[i].shape)
            self.ndim[i] = ndim
            if ndim > 0:
                self.dims[i][0:ndim] = dt[i].shape
                w,=where(self.dims[i] > 0)
                if w.size > 0:
                    nel = self.dims[i][w].prod()
            else:
                nel = 1

            self.nel[i] = nel
            self.totlen[i] = \
                self.calc_array_totlen(ndim, nel, nbytes)
            tname = dt[i].base.name
            self.tname.append(tname)
            self.arrsig[i][:] = self.get_array_sig(tname)

            tdescr = dt[i].base.descr[0][1]
            self.tdescr.append(tdescr)

            if ndim == 0:
                self.bytesperrow += 4 + self.nbytes[i]*self.nel[i]
            else:
                self.bytesperrow += self.totlen[i]

            if self.verbose >= 2:
                sys.stdout.write("field='"+dt.names[i]+"'\n")
                sys.stdout.write('\tnbytes= %s\n' % nbytes)
                sys.stdout.write('\tndim= %s\n' % ndim)
                sys.stdout.write('\tdims= %s\n' % self.dims[i])
                sys.stdout.write('\tnel= %s\n' % self.nel[i])
                sys.stdout.write('\ttotlen= %s\n' % self.totlen[i])
                sys.stdout.write('\tarrsig= %s\n' % self.arrsig[i])

        if self.verbose >= 2:
            sys.stdout.write('bytes per row: %s\n' % self.bytesperrow)

    def set_output_dtype(self):
        descr = []

        descr.append( ('nfields', 'i2') )

        dt=self.dtype
        for i in range(len(self.dtype)):

            # this creates the field definition for the output.  We will have
            # to include some descriptions fields also, see below
            npname = self.tdescr[i]
            # force to be big endian
            #if npname[1] != 'S':
                #    npname = '>' + npname[1:]
            # force signed integers since postgres has no unsigned
            if npname[1] == 'u':
                npname = npname[0] + 'i' + npname[2:]

            # now description fields
            if self.ndim[i] == 0:
                descr.append( (self.names[i]+'_len', 'i4') )
                descr.append( (self.names[i], npname) )
            else:
                descr.append( (self.names[i]+'_totlen', 'i4') )
                descr.append( (self.names[i]+'_ndim', 'i4') )
                descr.append( (self.names[i]+'_zero', 'i4') )
                descr.append( (self.names[i]+'_asig', 'u1',4) )
                for j in range(self.ndim[i]):
                    descr.append((self.names[i]+'_dim'+str(j)+'_size', 'i4'))
                    descr.append( (self.names[i]+'_dim'+str(j)+'_one', 'i4') )
                for j in range(self.nel[i]):
                    descr.append( (self.names[i]+'_'+str(j)+'_len', 'i4') )
                    descr.append( (self.names[i]+'_'+str(j)+'_data', npname) )

        self.output_dtype = numpy.dtype(descr)

    def init_output(self, num):
        """
        Set all the length fields
        """
        arr = numpy.zeros(num, dtype=self.output_dtype)

        arr['nfields'] = len(self.dtype)
        for i in range(len(self.dtype)):
            name=self.names[i]
            if self.ndim[i] == 0:
                arr[name+'_len'] = self.nbytes[i]
            else:
                arr[name+'_totlen'] = self.totlen[i]
                arr[name+'_ndim'] = self.ndim[i]
                arr[name+'_zero'] = 0
                for j in range(4):
                    ta = arr[name+'_asig']
                    #ta[:][j] = self.arrsig[i][j]
                    arr[name+'_asig'][:,j] = self.arrsig[i][j]
                for j in range(self.ndim[i]):
                    arr[name+'_dim'+str(j)+'_size'] = self.dims[i][j]
                    arr[name+'_dim'+str(j)+'_one'] = 1
                for j in range(self.nel[i]):
                    arr[name+'_'+str(j)+'_len'] = self.nbytes[i]
        return arr




    # methods for copying from the input data to the output structure
    # these deal properly with strings as well as byte ordering issues
    # which can occur during the copying process
    def space_pad_string(self, instring, length):
        return instring.ljust(length,' ')

    def copy_string(self, instring, length=None):
        """
        Since postgresql version 8.1 we must pad the strings with spaces
        for binary copy.
        """
        if self.pad_strings:
            if length is not None:
                return self.space_pad_string(instring, length)
            else:
                return instring
        else:
            return instring

    def copy_scalar_field_old(self, inarr, outarr, index):
        name = self.names[index]
        if stringmatch.match(self.tname[index]):
            for i in range(inarr.size):
                outarr[i][name] = self.copy_string(inarr[i][name],
                                                   self.nbytes[index])
        else:
            outarr[name] = inarr[name]

    def copy_array_field_old(self, inarr, outarr, index):
        name = self.names[index]
        # loop over rows
        for j in range(len(inarr)):
            # get a row, unraveled
            din = inarr[j][name].ravel()
            # loop over elements in this row/field
            for k in range( din.size ):
                oname = name+'_'+str(k)+'_data'

                # sometimes this copy gets native byte order, so check
                val = din[k]
                #if (val.dtype.byteorder != inarr[j][name].dtype.byteorder):
                #    val = val.byteswap()

                if stringmatch.match(self.tname[index]):
                    outarr[j][oname]=self.copy_string(val,self.nbytes[index])
                else:
                    outarr[j][oname] = val


    def copy_scalar_field(self, inarr, outarr, index):
        name = self.names[index]
        if stringmatch.match(self.tname[index]):
            for i in xrange(inarr.size):
                outarr[name][i] = self.copy_string(inarr[name][i],
                                                   self.nbytes[index])
        else:
            outarr[name] = inarr[name]

    def copy_array_field(self, inarr, outarr, index):
        name = self.names[index]
        # loop over rows
        for j in xrange(len(inarr)):
            # get a row, unraveled
            din = inarr[name][j].ravel()
            # loop over elements in this row/field
            for k in xrange( din.size ):
                oname = name+'_'+str(k)+'_data'

                # sometimes this copy gets native byte order, so check
                val = din[k]
                #if (val.dtype.byteorder != inarr[name][j].dtype.byteorder):
                #    val = val.byteswap()

                if stringmatch.match(self.tname[index]):
                    outarr[oname][j]=\
                        self.copy_string(val,self.nbytes[index])
                else:
                    outarr[oname][j] = val



    def copy_field(self, inarr, outarr, index):
        name = self.names[index]
        if self.ndim[index] == 0:
            if self.verbose >= 2:
                sys.stdout.write('%s: scalar\n' % name)
            self.copy_scalar_field(inarr, outarr, index)
        else:
            if self.verbose >= 2:
                sys.stdout.write('%s: array\n' % name)
            self.copy_array_field(inarr, outarr, index)

    def copy_to_output(self, inarr, outarr):
        # loop over fields
        for index in range(len(self.dtype)):
            self.copy_field(inarr, outarr, index)





    def get_dimstring(self,i):
        if self.ndim[i] == 0:
            return ''

        dims = self.dims[i,:].copy()
        w,=where(dims > 0)
        
        dims = list(dims[w])
        dimstring=''
        for dim in dims:
            dimstring += '[%d]' % dim
        return dimstring


    def set_coldefs(self, primary_key=None, unique=None):
        self.coldefs = []
        for i in range(len(self.dtype)):
            # get the postgres type name
            tname = self.tname[i]
            pgtype = npy2pg(tname)
            if pgtype == 'varchar' or pgtype == 'character':
                pgtype += '(%d)' % self.nbytes[i]

            # dimensions
            pgtype += self.get_dimstring(i)
            extra=''
            if not stringmatch.match(tname):
                extra = ' NOT NULL'
            coldef = '%s %s%s' % (self.names[i], pgtype, extra) 
            self.coldefs.append(coldef)
            if self.verbose >= 2:
                sys.stdout.write('\n\t=> %s\n' % coldef)

        if self.serial_extra is not None:
            sadd = ['%s SERIAL8' % s for s in self.serial_extra]
            self.serial_extra_coldefs = sadd
        else:
            self.serial_extra_coldefs = None

        if primary_key is not None:
            self.coldefs.append('PRIMARY KEY (%s)' % primary_key)

        if unique is not None:
            uadd = ['UNIQUE (%s)' % u for u in unique]
            self.coldefs = self.coldefs + uadd

    def print_create_table_statement(self):
        coldefs = self.coldefs
        if self.serial_extra_coldefs is not None:
            coldefs = self.serial_extra_coldefs + coldefs

        pcoldefs = ['          '+cd.lower() for cd in coldefs]
        pcoldefs = ',\n'.join(pcoldefs)
        printable = """
        CREATE TABLE 
            %s 
        (
%s
        )\n""" % (self.tablename, pcoldefs)
        sys.stdout.write(printable)

    def set_create_table_statement(self):
        coldefs = self.coldefs
        if self.serial_extra_coldefs is not None:
            coldefs = self.serial_extra_coldefs + coldefs
        coldefs = ',\n '.join(coldefs)
        tabledef = \
            'create table %s (%s)' % (self.tablename, coldefs)
        self.create_table_statement = tabledef
            


