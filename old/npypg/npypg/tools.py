import os
import sys
from sys import stdout,stderr
import numpy
import re
import tempfile
import npypg
from npypg import query as pgquery

# only unix for now.  This is because on Darwin
# the default tempdir is only readable by the owner
# of the file
tempfile.tempdir = '/tmp'

stringmatch = re.compile('^string.*')

adminconn_def = 'user=postgres'


def array2table(array, tablename, verbose=False, outdir=None,
                conninfo=None, keep=False, pad_strings=True,
                primary_key=None, unique=None, serial_extra=None):

    pgi = PgInput(verbose=verbose, conninfo=conninfo, pad_strings=pad_strings)
    pgi.array2table(array, tablename, outdir=outdir, keep=keep,
                    primary_key=primary_key, unique=unique, 
                    serial_extra=serial_extra)


def pgqueryln(query, conninfo=None):
    npypg.query(query, conninfo=conninfo)
    sys.stdout.write('\n')

def tables(conninfo=None):
    query="""
    select 
        tablename 
    from 
        pg_tables 
    where 
        schemaname = 'public'
    """
    tb=pgquery(query, flengths={'tablename':100}, 
               conninfo=conninfo)
    return tb['tablename']

def describe(table, conninfo=None):
    """
    Could do this with psycopg
    """
    q="select * from %s limit 1" % table
    res = npypg.query(q, conninfo=conninfo)
    dtypes=[]
    max_lens = {'name':0, 'type':0, 'shape':0}
    for d in res.dtype.descr:
        info = {}
        info['name'] = d[0]
        info['type'] = d[1][1:]
        if len(d) > 2:
            info['shape'] = d[2]
        else:
            info['shape'] = None

        for lentype in max_lens:
            try:
                if len(info[lentype]) > max_lens[lentype]:
                    max_lens[lentype] = len(info[lentype])
            except:
                # usually due to None
                pass

        dtypes.append(info)

    padding = 3

    forms = {}
    forms['name'] = '%-'+str(max_lens['name']+padding)+'s'
    forms['type'] = '%-'+str(max_lens['type']+padding)+'s'
    forms['shape'] = '%-'+str(max_lens['shape']+padding)+'s'

    for n in ['name','type','shape']:
        stdout.write(forms[n] % n)
    stdout.write('\n')

    #line = '-'*len(list(forms.keys()))

    for d in dtypes:
        for n in ['name','type','shape']:
            stdout.write(forms[n] % d[n])

        stdout.write('\n')


def table_exists(tablename, conninfo=None):
    tb=tables(conninfo=conninfo)
    w,=numpy.where(tb == tablename)
    if w.size > 0:
        return True
    else:
        return False

def obliterate(tablename, conninfo=None):
    pgqueryln('truncate %s' % tablename, conninfo=conninfo)
    pgqueryln('vacuum full %s' % tablename, conninfo=conninfo)
    pgqueryln('drop table %s' % tablename, conninfo=conninfo)


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

types={}
types['i2'] = {'postgres':'int2', 'numpy':'i2'}
types['int16'] = types['i2']

types['i4'] = {'postgres':'int4', 'numpy':'i4'}
types['int32'] = types['i4']

types['i8'] = {'postgres':'int8', 'numpy':'i8'}
types['int64'] = types['i8']

types['f4'] = {'postgres':'float4', 'numpy':'f4'}
types['float32'] = types['f4']

types['f8'] = {'postgres':'float8', 'numpy':'f8'}
types['float64'] = types['f8']


#string_store = 'varchar'
string_store = 'character'
#string_store = 'text'
types['string'] = {'postgres':string_store, 'numpy':'string'}

def npy2pg(npy_typename):
    if stringmatch.match(npy_typename):
        return types['string']['postgres']

    if npy_typename not in types:
        raise ValueError('Unsupported type: %s' % npy_typename)
    return types[npy_typename]['postgres']



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
    input = PgInput(tablename=None, dtype=None,
                    conninfo='user=postgres', verbose=False)

    oids and extended header not yet implemented

    example 1:
        pgi = PgInput()
        pgi.array2table(array, tablename, outdir=None)
        # see also the array2table convenience function

    example 2:  # copy in multiple structures to an existing table
        if array1.dtype != array2.dtype:
            raise ValueError('all array data types must be the same')
        pgi = PgInput('mydata', array1.dtype)
        pgi.open(fname)
        pgi.write_header()
        pgi.write_data(array1)
        pgi.write_data(array2)
        pgi.write_trailer()
        pgi.close()

        pgi.copy_input_to_table(fname)

    """
    def __init__(self, tablename=None, dtype=None, 
                 conninfo=None, 
                 pad_strings=True,
                 verbose=False):
        if dtype is not None:
            self.set_dtype(dtype)
        self.verbose=verbose
        self.tablename = tablename
        self.pad_strings = pad_strings
        self.conninfo=conninfo

    def array2table(self, array_in, tablename, outdir=None, keep=False,
                    primary_key=None, unique=None, serial_extra=None):
        array = array_in.view(numpy.ndarray)

        pref='pgcopy-'+tablename+'-'
        fname=tempfile.mktemp(prefix=pref,dir=outdir, suffix='.pgsql')
        fname=os.path.expanduser(fname)

        self.tablename = tablename
        self.serial_extra=serial_extra
        self.set_dtype(array.dtype)

        self.set_coldefs(primary_key=primary_key, unique=unique)
        self.set_create_table_statement()

        # does the table exist?  If not, create it
        if not table_exists(tablename, conninfo=self.conninfo):
            self.print_create_table_statement()
            pgqueryln(self.create_table_statement, conninfo=self.conninfo)

        # Noe write the file for input
        self.open(fname)
        self.write_header()
        self.write_data(array)
        self.write_trailer()
        self.close()

        t=os.popen('chmod a+rwx '+fname)
        t.read()
        t.close()

        self.copy_input_to_table(fname)

        os.remove(fname)

    def copy_input_to_table(self, fobj):
        """
        This requires admin privledges.

        If serial_extra not None, then we specify all the columns that will be
        copied not including the serial column
        """
        if isinstance(fobj, file):
            fname = file.name
        else:
            fname=fobj

        columns = ''
        if self.serial_extra is not None:
            # the serial extra columns were put at the beginning
            columns = '(%s)' % ','.join( self.names )

        q="""
        COPY BINARY 
            %s %s
        FROM
            '%s'
        """ % (self.tablename, columns, fname)
        if self.verbose:
            sys.stdout.write('%s\n' % q)
        conn=self.conninfo
        if conn is None:
            conn = adminconn_def 
        else:
            # libpq is smart and will take the last keyword-value pair when
            # there are duplicates
            conn = conn + ' ' + adminconn_def 
        
        conn = conn.strip()
        if self.verbose:
            sys.stdout.write("connection info = '%s'\n" % conn)
        pgqueryln(q, conninfo=conn)


    def set_dtype(self, dtype):
        # This can be list, dict, or already a dtype
        self.dtype = numpy.dtype(dtype)
        self.set_field_info()
        self.set_output_dtype()

        # we could implement these in the future
        self.flags = numpy.array([0], dtype='>i4')
        self.hdext_len = numpy.array([0], dtype='>i4')
        self.trailer = numpy.array([-1], dtype='>i2')
        
        # nothing sent for files yet
        self.fileobj_entered = False

    def open(self, fileobj=None):
        if fileobj is None:
            raise ValueError('need to implement temp files')
        elif isinstance(fileobj, file):
            self.fileobj_entered = True
            self.fobj = fileobj
        elif isinstance(fileobj, str):
            self.fileobj_entered = False
            self.fobj = open(fileobj, 'w')
        else:
            raise ValueError('fileobj entered must be None,file, or str')

    def write_header(self):
        self.write_sig()
        self.write_flags()
        self.write_hdext_len()

    def close(self):
        if not self.fileobj_entered:
            # If we opened the file, we should close it
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
        self.trailer.tofile(self.fobj)


    def write_data(self, array, mbperwrite=100):
        """
        This creates a temporary structure that contains field lengths 
        embedded, which matches the postgres output file structure
        """
    
        n=len(array)
        output = self.init_output(n)
        self.copy_to_output(array, output)
        if self.verbose:
            sys.stdout.write('Input dtype:\n')
            for el in array.dtype.descr:
                sys.stdout.write('%s\n' % str(el))
            sys.stdout.write('Output dtype:\n')
            for el in output.dtype.descr:
                sys.stdout.write('%s\n' % str(el))
        make_bigendian(output, each=True, verbose=self.verbose)
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
            if types['string']['postgres'] == 'character':
                sig[:] = [0,0,4,18]
            elif types['string']['postgres'] == 'varchar':
                # This is for varchar(n) fields
                sig[:] = [0,0,4,19]
            elif types['string']['postgres'] == 'text':
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
                w,=numpy.where(self.dims[i] > 0)
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

            if self.verbose:
                sys.stdout.write("field='"+dt.names[i]+"'\n")
                sys.stdout.write('\tnbytes= %s\n' % nbytes)
                sys.stdout.write('\tndim= %s\n' % ndim)
                sys.stdout.write('\tdims= %s\n' % self.dims[i])
                sys.stdout.write('\tnel= %s\n' % self.nel[i])
                sys.stdout.write('\ttotlen= %s\n' % self.totlen[i])
                sys.stdout.write('\tarrsig= %s\n' % self.arrsig[i])

        if self.verbose:
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
                if (val.dtype.byteorder != inarr[j][name].dtype.byteorder):
                    val = val.byteswap()

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
                if (val.dtype.byteorder != inarr[name][j].dtype.byteorder):
                    val = val.byteswap()

                if stringmatch.match(self.tname[index]):
                    outarr[oname][j]=\
                        self.copy_string(val,self.nbytes[index])
                else:
                    outarr[oname][j] = val



    def copy_field(self, inarr, outarr, index):
        name = self.names[index]
        if self.ndim[index] == 0:
            if self.verbose:
                sys.stdout.write('%s: scalar\n' % name)
            self.copy_scalar_field(inarr, outarr, index)
        else:
            if self.verbose:
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
        w,=numpy.where(dims > 0)
        
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
            if self.verbose:
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
            w,=numpy.where(arr1[n1].ravel() != arr2[n2].ravel())
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
