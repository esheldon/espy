"""
Package:
    npypg pg_numpy 
Description
    A set of modules for interacting with a PostgreSql database.  A C++
    extension module allows returning results as a NumPy array.  Numpy
    arrays can also be written to tables.

    The query functions are not designed for generic database interactions,
    but primarily to efficiently get results into a numpy array.  More
    generic DB-API compliant packges like psycopg are more suitable when
    more general and flexible operations are needed.

some functions:
    query(query_string, connect_info=None, flengths={}) 
        Query the database and return the result as a numpy array

    array2table(array, tablename, connect_info=None)
        Write the array to a postgres table.  If the table does not yet exist
        it is created with column definitions based on the input array. If it
        does exist the data are appended as new rows in the table. 

Modules:
    PGConn: 
        A C++ extension for connecting to a PostgreSql database and executing a
        query. 
        
        The PGConn module defines the query() function described above, which is
        also exported to the primary namespace.   Most of the other functions
        in this package use the query() function.

Modification history:
    Created: 2008-04-07, Erin Sheldon, NYU
"""

import PGConn
from PGConn import query
import tools
from tools import array2table
from tools import tables
from tools import table_exists
from tools import describe
from tools import test
from tools import test_simple
from tools import obliterate
from tools import compare_arrays

# attempt to import the connect method from psycopg2
try:
    from psycopg2 import connect
except:
    pass
