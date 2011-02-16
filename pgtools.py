import os
import numpy

from sys import stdout, stderr

have_psycopg = False
try:
    import psycopg2
    have_psycopg = True
except:
    stderr.write('Could not import psycopg')

have_npypg = False
try:
    import npypg
    have_npypg = True
except:
    stderr.write('Could not import npypg')

dbnamedef = None
userdef = None
if 'PGDATABASE' in os.environ:
    dbnamedef=os.environ['PGDATABASE']
if 'PGUSER' in os.environ:
    userdef=os.environ['PGUSER']

pg2numpy_scalars = {
    21:   'i2',
    23:   'i4',
    20:   'i8',
    700:  'f4',
    701:  'f8',
    1042: 'S',       # length to be determined
    1043: 'S',      # varchar, possibly with max length
    25:   'S'} # text, don't know how long

pg2numpy_arrays = {
    1005: 'i2', 
    1007: 'i4',
    1016: 'i8',
    1009: 'S',   # text, len undefined, desc[3] is -1
    1014: 'S',      # fixed length character.  len in desc[3]
    1015: 'S',      # varchar.  desc[3] != -1 means max lenght is specified
    1021: 'f4',
    1022: 'f8'}

# Default max length for string fields. This is used if the field in the
# database does not specify a length and if the user did not sent a 
# max length
def_textmax = 10


def PgDescr2NumpyDescr(pg_descr, textmax=None, lendict={}):
    """
    Determine the numpy type descriptor for each of the input types
    For array types we will still have to deal with the array length. 
    """

    if lendict is None:
        lendict={'default':10}

    descr = []
    has_arrays=False
    for d in pg_descr:
        name=d[0]
        pgtype=d[1]

        if pgtype in pg2numpy_scalars:
            field_is_array=False
            npt = pg2numpy_scalars[pgtype]
        elif pgtype in pg2numpy_arrays:
            field_is_array=True
            has_arrays=True
            npt = pg2numpy_arrays[pgtype]
        else:
            raise ValueError('type '+pgtype+' not supported')


        if npt is 'S':
            # Get size.  We only deal with fixed length at this point
            sz = d[3]
            if sz != -1:
                npt = 'S'+str(sz)
            else:
                if name in lendict:
                    # See if the user sent something in lendict
                    npt = 'S'+str(lendict[name])
                elif textmax is not None:
                    # See if textmax is sent
                    npt = 'S'+str(textmax)
                else:
                    npt = 'S'+str(def_textmax)



        if field_is_array:
            # Will fill in the length later
            descr.append( (name,npt,-1) )
        else:
            descr.append( (name,npt) )

    return descr, has_arrays

def PgGetArrayDims(data):
    """
    Takes a list and determines the dimensions assuming all dims are the
    same
    """

    if isinstance(data, str):
        return [0]
    
    try:
        dlen = [len(data)]
    except:
        dlen = [0]

    if dlen[0] != 0:
        try:
            nextlen = PgGetArrayDims(data[0])
            if nextlen[0] != 0:
                dlen = dlen+nextlen
        except:
            pass
    return dlen

def PgGetDescrWithArrays(descr, data):
    newdescr=[]
    i=0
    for i in range(len(descr)):
        tdescr=descr[i]
        tdata=data[0][i]

        if len(tdescr) == 3:
            tdim = tuple( PgGetArrayDims(tdata) )
            if tdim[0] > 0:
                newdescr.append( (tdescr[0], tdescr[1], tdim) )
            else:
                newdescr.append(tdescr)
        else:
            newdescr.append(tdescr)

    return newdescr


def fetch_arrays(query_string,conninfo=None, user=userdef, database=dbnamedef, dtype=None,textmax=None,lendict={}):

    if have_npypg:
        res = npypg.query(query_string)
    elif have_psycopg:

        conn = psycopg2.connect(conninfo, user=user, database=database)

        cur = conn.cursor()

        cur.execute(query_string)

        # For queries that return nothing
        if cur.description is None or cur.rowcount == 0:
            return None

        # If the user sent dtype, we don't have to compute the dtype can 
        # set has_arrays to false
        has_arrays=False
        if dtype is None:
            dtype, has_arrays = PgDescr2NumpyDescr(cur.description, 
                                                   textmax=textmax, 
                                                   lendict=lendict)
    
        if has_arrays:
            #print 'Result contains arrays.  We must prefetch the data into'
            #print 'lists to determine the array sizes; this is less memory '
            #print 'efficient.  Suggestion: select a single row to get the dtype, '
            #print 'then get the full query sending dtype as a keyword'
            t = cur.fetchall()
            newdtype = PgGetDescrWithArrays(dtype, t)
            res = numpy.array(t, dtype=newdtype)
        else:
            res = numpy.fromiter(cur, dtype=dtype)

        conn.close()
    else:
        raise RuntimeError('Either psycopg or npypg must be available')

    return res


def test():
    sys.stdout.write('First check\n')
    tt = fetch('select * from test_sweep limit 1')
    dt=tt.dtype
    sys.stdout.write('Full get on dtype: %s\n' % dt)
    t = fetch('select * from test_sweep limit 10000', dtype=dt)
    return t
