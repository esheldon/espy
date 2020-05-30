"""
Package:
    pgnumpy
Description

    A class and a set of functions for interacting with a PostgreSql database.
    A C++ extension module allows returning results as a NumPy array.  Numpy
    arrays can also be written to tables.  
    
    The workhorse class is called PgNumpy

    This class has limited functionality compared to the full Python database
    api specification.  It can execute arbitrary queries and extract results
    into numpy arrays.  However, cursors are not yet supported.  For getting
    results, only the fetchall() command is available, as the goal is always to
    extract all rows into a single numpy structure rather than work row by row.
    
    More generic DB-API compliant packges like psycopg are more suitable when
    more flexible operations are needed.

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

import pgnumpy
import cpgnumpy

from pgnumpy import connect
from pgnumpy import PgNumpy
from pgnumpy import PgInput
from pgnumpy import ArrayWriter
from pgnumpy import ArrayStringifier
from pgnumpy import array2table

#from pgnumpy import tables
#from pgnumpy import table_exists
#from pgnumpy import describe
from pgnumpy import test
from pgnumpy import test_simple
#from pgnumpy import obliterate
#from pgnumpy import compare_arrays

# attempt to import the connect method from psycopg2
try:
    from psycopg2 import connect
except:
    pass
