#include "cpgnumpy.h"

//Constructor and destructor
cPgNumpy::cPgNumpy() throw (const char*)
{
    import_array(); // Must be present for numpy.
    _set_defaults();
    _print_debug("Opening with no explicit connect info");
    open();
}

//Constructor and destructor
cPgNumpy::cPgNumpy(string conninfo) throw (const char*)
{
    import_array(); // Must be present for numpy.
    _set_defaults();
    _print_debug("Opening with no explicit connect info");
    open(conninfo.c_str());
}


//Constructor and destructor
cPgNumpy::cPgNumpy(const char* conninfo) throw (const char*)
{
    import_array(); // Must be present for numpy.
    _set_defaults();
    _print_debug("Opening with no explicit connect info");
    open(conninfo);
}



cPgNumpy::~cPgNumpy()
{
    Py_XDECREF(mFlenDict);
    close();
}


// Open a connection
void cPgNumpy::open() throw (const char*)
{
    open("");
}
void cPgNumpy::open(string conninfo) throw (const char*)
{
    open(conninfo.c_str());
}

void cPgNumpy::open(const char* conninfo) throw (const char*)
{
    close();
    if (conninfo != NULL)
        mConnectInfo = conninfo;

    mConn = PQconnectdb(mConnectInfo.c_str());

    if (PQstatus(mConn) != CONNECTION_OK) {
        PQfinish(mConn);
        stringstream err;
        err<<"Could not establish connection\n";
        err<<"connection info is '"<<mConnectInfo<<"'\n";
        throw err.str().c_str();
    }

}


// execute a query
void cPgNumpy::execute(string query_string) throw (const char*)
{
    execute(query_string.c_str());
}

void cPgNumpy::execute(const char* query_string) throw (const char*)
{
    // first clear any existing results
    clear();
    if (debug) {
        cout<<"Executing query: '"<<query_string<<"'"<<endl;
        fflush(stdout);
    }
    mResult = PQexecParams(mConn,
            query_string,
            0,    /* number of parameters (none) */
            NULL, /* see doc for the parameter info */
            NULL,
            NULL,
            NULL,
            mBinary); /* 0 for text, 1 for binary (network order) */

    // sets mResultStatus, mNtuples, mNfields
    // also can set error strings
    _check_status();
}

// Set this to use binary
void cPgNumpy::use_text() {
    mBinary = 0;
}
void cPgNumpy::set_decorators(bool decorate) {
    mUseArrayDecorators=decorate;
}

long long cPgNumpy::write(const char* filename) throw (const char*) {
    if (mBinary) {
        throw "You must execute the query in text retrieval mode to use write()";
    }

    long long nfields;
    long long ntuples;
    long long field;
    long long row;

    FILE *OFILEPTR=NULL;

    if (filename != NULL) {
        OFILEPTR = fopen(filename, "w");
        if (!OFILEPTR) {
            stringstream err;
            err<<"Could not open file for writing: "<<filename<<": ";
            err<<strerror(errno);
            throw err.str().c_str();
        }
    } else {
        OFILEPTR=stdout;
    }

    nfields = PQnfields(mResult);
    ntuples = PQntuples(mResult);

    // if we are not using array decorators, determine whic
    // are arrays in order to remove them later
    vector<int> isarray;
    isarray.assign(nfields,0);
    if (!mUseArrayDecorators) {
        for (long long i=0; i<nfields; i++) {
            int ftype = PQftype(mResult, i);

            string nptype;
            int nbytes=0;
            int this_isarray=0;
            _pg_type_info(ftype, nptype, nbytes, this_isarray);
            if (this_isarray) {
                isarray[i] = 1;
            }
        }
    }
    for (row=0; row<ntuples; row++) {
        for (field=0;field<nfields;field++)
        {
            const char* value = PQgetvalue(mResult, row, field);

            // should we remove { and }?
            if (isarray[field] && !mUseArrayDecorators) {
                int slen = strlen(value);
                for (int i=0; i<slen; i++) {
                    char c=value[i];
                    if (c != '{' && c != '}') {
                        fputc(c, OFILEPTR);
                    }
                }
            } else {
                fprintf(OFILEPTR,"%s",value); 
            }

            /* If last field print a new line, else a tab*/
            if (field == nfields-1)
                fprintf(OFILEPTR,"\n");
            else
                fprintf(OFILEPTR,",");
        }
    }

    if (filename != NULL) {
        fclose(OFILEPTR);
    }

    clear();

    return ntuples;
}



// size of fetch chunks when doing an execwrite
void cPgNumpy::set_fetch_count(long long fetchcount) {
    mFetchCount = fetchcount;
}

// Execute the query and use a cursor to retrieve the results in
// chunks, writing to a file or by default stdout.  Unlike write(),
// does *not* require use_text() to be set.
//
// This does not use the query execution infrastructure of PnNumpy,
// but rather is rolled specially for the cursor
//
// The return value is the number of rows returned
long long cPgNumpy::execwrite(const char* query, const char* filename) throw (const char*) {

    PGresult *res=NULL;

    //PQsetnonblocking(mConn,1);

    // Start a transaction block
    res = PQexec(mConn, "BEGIN");
    if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
        stringstream err;
        err<<"BEGIN command failed: "<<PQerrorMessage(mConn);
        PQclear(res);
        PQsetnonblocking(mConn,0);
        throw err.str().c_str();
    }
    /*
     * Should PQclear PGresult whenever it is no longer needed to avoid
     * memory leaks
     */
    PQclear(res);



    // set up the cursor
    string cursor_query = "DECLARE cur CURSOR FOR ";
    cursor_query += query;
    res = PQexec(mConn, cursor_query.c_str());
    if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
        stringstream err;
        err<<"DECLARE CURSOR failed: "<<PQerrorMessage(mConn);
        PQclear(res);
        PQsetnonblocking(mConn,0);
        throw err.str().c_str();
    }
    PQclear(res);




    long long nfields;
    long long ntuples;
    long long field;
    long long row;

    FILE *OFILEPTR=NULL;

    if (filename != NULL) {
        OFILEPTR = fopen(filename, "w");
        if (!OFILEPTR) {
            stringstream err;
            err<<"Could not open file for writing: "<<filename<<": ";
            err<<strerror(errno);
            PQsetnonblocking(mConn,0);
            throw err.str().c_str();
        }
    } else {
        OFILEPTR=stdout;
    }


    // true until we get the info
    bool getarrayinfo=true;
    vector<int> isarray;

    // grab the first set of results
    stringstream fcstream;
    fcstream << "FETCH "<<mFetchCount<<" in cur ";
    string fetch_command = fcstream.str();
    int rc = PQsendQuery(mConn, fetch_command.c_str());

    if (!rc) {
        stringstream err;
        err<<"Error fetching: "<<PQerrorMessage(mConn);
        throw err.str().c_str();
    }


    long long rowcount=0;
    while (1) {

        //rc = PQconsumeInput(mConn);
        res = PQgetResult(mConn);

        if (res == NULL) {
            // we need to send another fetch command
            rc = PQsendQuery(mConn, fetch_command.c_str());
            if (!rc) {
                stringstream err;
                err<<"Error fetching: "<<PQerrorMessage(mConn);
                throw err.str().c_str();
            }
        } else {

            nfields = PQnfields(res);
            ntuples = PQntuples(res);
            if (ntuples == 0) {
                // this means we are done
                PQclear(res);
                break;
            }


            rowcount += ntuples;

            fflush(stdout);

            // if we are not using array decorators, determine which
            // are arrays in order to remove them later
            isarray.assign(nfields,0);
            if (!mUseArrayDecorators && getarrayinfo) {

                for (long long i=0; i<nfields; i++) {
                    int ftype = PQftype(res, i);

                    string nptype;
                    int nbytes=0;
                    int this_isarray=0;
                    _pg_type_info(ftype, nptype, nbytes, this_isarray);
                    if (this_isarray) {
                        isarray[i] = 1;
                    }
                }
                getarrayinfo=false;
            }

            // now print results
            for (row=0; row<ntuples; row++) {
                for (field=0;field<nfields;field++) {
                    const char* value = PQgetvalue(res, row, field);

                    // should we remove { and }?
                    if (isarray[field] && !mUseArrayDecorators) {
                        int slen = strlen(value);
                        for (int i=0; i<slen; i++) {
                            char c=value[i];
                            if (c != '{' && c != '}') {
                                fputc(c, OFILEPTR);
                            }
                        }
                    } else {
                        fprintf(OFILEPTR,"%s",value); 
                    }

                    /* If last field print a new line, else a tab*/
                    if (field == nfields-1)
                        fprintf(OFILEPTR,"\n");
                    else
                        fprintf(OFILEPTR,",");
                }
            }

            PQclear(res);

        }


    }

    res = PQexec(mConn, "CLOSE cur");
    PQclear(res);
    /* end the transaction */
    res = PQexec(mConn, "END");
    PQclear(res);

    PQsetnonblocking(mConn,0);

    if (filename != NULL) {
        fclose(OFILEPTR);
    }

    return rowcount;
}



// Complex checking of the query status.  This is used to possibly set an
// error string for later raised exceptions

void cPgNumpy::_check_status() throw (const char*) {

    //   For successful queries, there are two options: either the query
    //   returns results or not.  If it does not return data, but
    //   successfully completed, then the status will be set to
    //   PGRES_COMMAND_OK.  If success and returns tuples, then it will be
    //   set to PGRES_TUPLES_OK 

    mResultStatus = PQresultStatus(mResult);
    mNtuples = (npy_intp) PQntuples(mResult);
    mNfields = (npy_intp) PQnfields(mResult);

    if (debug) {
        cout<<"    ntuple: "<<mNtuples<<"\n";
        cout<<"    nfields: "<<mNfields<<"\n";
    }
    // in this case no results: print what happened
    if (PGRES_COMMAND_OK == mResultStatus) {
        _print_debug("result status is PGRES_COMMAND_OK"); 
        // Success but no results
        cout<<PQcmdStatus(mResult)<<"\n";
        return;
    }
    // in 8.4 empty queries will return this instead of PGRES_COMMAND_OK
    if (PGRES_EMPTY_QUERY == mResultStatus) {
        _print_debug("result status is PGRES_EMPTY_QUERY"); 
        // Success but no results
        cout<<PQcmdStatus(mResult)<<"\n";
        return;
    }
    if (PGRES_TUPLES_OK != mResultStatus) {
        switch (mResultStatus) {
            // No results
            case PGRES_COPY_IN: 
                _print_debug("result status is PGRES_COPY_IN"); 
            return;
            case PGRES_COPY_OUT: 
                _print_debug("result status is PGRES_COPY_OUT"); 
                return;
            case PGRES_BAD_RESPONSE:
                // might be OK?
                _print_debug("result status is PGRES_BAD_RESPONSE"); 
                _print_debug(PQresultErrorMessage(mResult));
                break;
            case PGRES_NONFATAL_ERROR: 
                // this is OK, just print it
                _print_debug("result status is PGRES_NONFATAL_ERROR"); 
                _print_debug(PQresultErrorMessage(mResult));
                break;
            case PGRES_FATAL_ERROR:    
                _print_debug("result status is PGRES_FATAL_ERROR"); 
                throw PQresultErrorMessage(mResult);
            default: break;
        } 
    }

    // Result is empty, either an error or this query returns no data
    if ( PQnfields(mResult) == 0 && PQntuples(mResult) == 0 ) {
        // Query returns no data. Try to print a message to clarify
        cout<<PQcmdStatus(mResult);
        return;
    }

}

void cPgNumpy::set_field_lengths(PyObject* flenDict) throw (const char*)
{
    // make sure it's a dict
    if (!PyDict_Check(flenDict)) {
        throw "Input length dict must be a dictionary";
    }
    Py_XDECREF(mFlenDict);
    mFlenDict = PyDict_Copy(flenDict);
}

void cPgNumpy::clear_field_lengths() throw (const char*)
{
    // when we make a copy we should do a decref
    Py_XDECREF(mFlenDict);
}


PyObject* cPgNumpy::fetchall() throw (const char*) {

    PyObject* ret=NULL;
    npy_intp dims[1] = {0};
    dims[0] = mNtuples;

    if (mNtuples == 0) {
        // return empty array instead of None
        ret = PyArray_ZEROS(1, dims, NPY_INT, 0);
        clear();
        return ret;
    }

    // create the descriptor for our output data
    _make_descr();

    // pars: ndim, shape, descr, fortran_flag
    ret = PyArray_Zeros(1, dims, mDescr, NPY_FALSE);

    if (!ret) {
        // Still no decref of mDescr; it is stolen even on failure.
        throw "could not create output array";
        //preturn(ret);
    }

    // Creation of array was successfull
    // Note, we do NOT decref the mDescr since references are kept
    // In the array itself.  It "steals" the reference meaning it must take
    // care of it

    // only fill array if not empty
    _fill_array(ret);
    clear();
    return ret;

}

void cPgNumpy::_fill_array(PyObject* ret) {

    // Fill in the array
    PyArrayObject* retarray = (PyArrayObject* ) ret;
    char* data = retarray->data;

    for (npy_intp row=0; row<mNtuples; row++) {
        for (npy_intp field=0; field<mNfields; field++) {

            // For scalars, this is the actual length of this field, which
            // may be less or more than the bytes we have allocated.  For
            // arrays it the entire length of the array data including
            // dimension specifications, etc

            int flen   = PQgetlength(mResult, row, field);

            // Store the data
            _store_field(
                    PQgetvalue(mResult, row, field),
                    flen,
                    data,
                    mFieldInfo.nbytes[field], 
                    mFieldInfo.isarray[field], 
                    mFieldInfo.isstring[field]);

            data += mFieldInfo.nbytes[field]*mFieldInfo.nel[field];
        }
    }

}


void cPgNumpy::_store_field(
        char* input,    // Input postgresql data 
        int flen,       // length of this field in the input data
        char* output,   // Output array data
        int nbytes,     // number of available bytes in the output array
        int isarray,    // is this an array?
        int isstring)   // is this a string?
{
    if (isarray) {
        _fill_array(input, output, nbytes, isstring);
	} else {
		_fill_scalar(input, flen, output, nbytes, isstring);
	}
}

// Fill a scalar with the appropriate number of bytes.  Swap bytes if 
// necessary
void cPgNumpy::_fill_scalar(
        char* input,     // Input postgresql data 
        int32_t flen,  // Lenght of this field in the input data 
        char* output,    // Output array data
        int nbytes,      // number of bytes available in the output data 
        int isstring)    // Is this a string?
{
    // copy the smallest of the two
    int bytes2copy=0;
    if (flen <= nbytes)
        bytes2copy = flen;
    else
        bytes2copy = nbytes;

    if (debug2) {
        cout<<
            "        "<<
            "flen = "<<flen<<
            " nbytes="<<nbytes<<
            " bytes2copy="<<bytes2copy<<endl;
        fflush(stdout);
    }

    // Just copy over the bytes
    memcpy(output, input, bytes2copy);

    // Only swap if not string
    if (!isstring) {
		_tohost(output, nbytes);
	}
}


void cPgNumpy::_fill_array(
        char* input,    // Input postgresql data 
        char* output,   // Output array data
        int nbytes,     // Number of bytes available in output data field
        int isstring)   // Is this a string?
{
    
    // This is the most generic possible version of this program.  We can
    // probably speed it up later

    int32_t ndim = _tohost( *( (int32_t *) input) );
    input += 4;

    if (debug2) cout<<"    Ndim = "<<ndim<<endl;

    // Skip past the zero that is always there, and the array signature
    input += 2*4;

    int32_t n_elements=1;
    int32_t tdim=0;
    for (int32_t i=0; i<ndim; i++) {
        tdim = _tohost( *( (int32_t *) input) );
        n_elements *= tdim;
        // skip past this and the "one" that is always there
        input += 2*4;
    }

    if (debug2) cout<<"        Nel = "<<n_elements<<endl;

    
    int32_t flen=0;
    for (int32_t i=0; i< n_elements; i++) {

        // Read the size of this array element. This is always same for
        // numbers

        flen = *(int32_t *) input;
        flen = _tohost(flen);
        input += 4;


        // store the value
        _fill_scalar(input, flen, output, nbytes, isstring);

        // Skip to the next element
        input += flen;
        output += nbytes; 
    }

}

// convert from network to host byte ordering
void cPgNumpy::_tohost(char* data, int nbytes)
{
    // no need if this is already big endian!
    if (IS_BIG_ENDIAN) {
        return;
    }
    switch(nbytes) 
	{
        case 2:
            {
                npy_int16* tptr = (npy_int16 *) data;
                *tptr = ntohs( *tptr );
            }
            break;
        case 4:
            {
                int32_t* tptr = (int32_t *) data;
                *tptr = ntohl( *tptr );
            }
            break;
        case 8:
            {
                int64_t* tptr = (int64_t *) data;
                *tptr = ntohll( *tptr );
            }
            break;
        default:
            break;
    }

}

int32_t cPgNumpy::_tohost(int32_t num) {
    if (IS_BIG_ENDIAN) {
        return num;
    } else {
        return ntohl(num);
    }
}

// Make a desciption list, or format list, one element for each field
// in the output array
void cPgNumpy::_make_descr()
{

    // Make a list of tuples.  This is temporary, we will decref it
    PyObject* dlist=PyList_New(0);
    

    mFieldInfo.nbytes.assign(mNfields,0);
    mFieldInfo.isarray.assign(mNfields,0);
    mFieldInfo.nel.assign(mNfields,0);
    mFieldInfo.isstring.assign(mNfields,0);

    for (int i=0; i<mNfields; i++) {

        string fname=PQfname(mResult,i);
        if (debug) {
            cout<<endl<<"Working on field: "<<fname<<"\n";
        }

        // Determine the type, including if it is an array
        string nptype;
        int nbytes=0;
        int isarray=0;
        int ftype = PQftype(mResult, i);

        _pg_type_info(ftype, nptype, nbytes, isarray);


        int dsize;
        if (isarray)
            dsize=3;
        else
            dsize=2;


        // Tuple holding descr info for this field
        PyObject* dtup = PyTuple_New(dsize);

        // We need to describe the shape for arrays. max_flen is the maximum
        // size over all elements in this field for the first row
        int32_t max_flen=0;
        if (isarray) {
            PyObject* dims = _get_array_field_shape(
                    i, 
                    mFieldInfo.nel[i], 
                    max_flen);
            _print_debug("  Setting dims in descr tuple");
            PyTuple_SetItem(dtup, 2, dims);
        } else {
            _print_debug("  scalar");
            mFieldInfo.nel[i] = 1;
        }

        // We have a few options here.  There is the size of the actual
        // field, that is what we get from PQgetlength().  This is perfectly
        // good for character(n) fixed length fields.  For varchar(n) we
        // should try to find the max length.  
        if (nbytes == -1) {
            nbytes = _pg_string_nbytes(fname, i, ftype, isarray, max_flen);
        }

        mFieldInfo.nbytes[i] = nbytes;
        mFieldInfo.isarray[i] = isarray;

        // Set the length in the format
        if (nptype == "S") {
            stringstream sstr;
            sstr << nbytes;
            nptype=nptype+sstr.str();
            mFieldInfo.isstring[i] = 1;
        }
        if (debug) {
            cout<<"  Format = "<<nptype<<endl;
		}
        PyObject* ttype = PyString_FromString(nptype.c_str());
        PyObject* tname = PyString_FromString(fname.c_str());

        // Don't have to decref objects put into sequence.  They will get
        // decrefed when the sequence is decrefed
        PyTuple_SetItem(dtup, 0, tname);
        PyTuple_SetItem(dtup, 1, ttype);


        // Same goes for this tuple: no decref needed
        PyList_Append(dlist, dtup);
        
    }

    _print_debug("Calling PyArray_DescrConverter");
    if (!PyArray_DescrConverter(dlist, &mDescr)) {
        throw "Could not create descriptor\n";
    }
    // put error checking here

    // We decref this but not the Descr
    Py_XDECREF(dlist);

}


// Return a string containing the type of the field, the number of bytes
// if known, and if this is an array
void cPgNumpy::_pg_type_info(
        int pgtype, 
        string &nptype, 
        int& nbytes, 
        int& isarray)
{
    switch (pgtype) {
		// Scalar types
		case  18: nptype="S1"; // char
				  nbytes=1;
                  break;
        case  20: nptype="i8";
                  nbytes=8;
                  break;
        case  21: nptype="i2";
                  nbytes=2;
                  break;
        case  23: nptype="i4";
                  nbytes=4;
                  break;
        case  25: nptype="S"; // text, length unknown
                  nbytes=-1;
                  break;
        // This is an OID.  They are only 4 byte
        case  26: nptype="i4";
                  nbytes=4;
                  break;
        case 700: nptype="f4";
                  nbytes=4;
                  break;
        case 701: nptype="f8";
                  nbytes=8;
                  break;
        case 1042: nptype="S"; // Fixed length character(n)
                   nbytes=-1;
                   break;
        case 1043: nptype="S"; // varchar, possibly with max len TBD
                   nbytes=-1;
                   break;

                   // array types
        case 1002: nptype="S1";  //char
                   nbytes=1;
                   isarray=1;
                   break;
        case 1005: nptype="i2"; 
                   nbytes=2;
                   isarray=1;
                   break;
        case 1007: nptype="i4";
                   nbytes=4;
                   isarray=1;
                   break;
        case 1016: nptype="i8";
                   nbytes=8;
                   isarray=1;
                   break;
        case 1009: nptype="S"; // text, len undefined, desc[3] is -1
                   isarray=1;
                   nbytes=-1;
                   break;
        case 1014: nptype="S"; // fixed length character(n).
                   isarray=1;
                   nbytes=-1;
                   break;
        case 1015: nptype="S"; // varchar, possibly with max len TBD
                   isarray=1;
                   nbytes=-1;
                   break;
        case 1021: nptype="f4";
                   isarray=1;
                   nbytes=4;
                   break;
        case 1022: nptype="f8";
                   nbytes=8;
                   isarray=1;
                   break;

        default: 
                   if (debug) {
                       stringstream out;
                       out<<"Fell through type case with: "<<pgtype<<"\n";
                       _print_debug(out.str().c_str());
                   }
                   nptype="S";
                   nbytes=-1;
                   break;
    }

}


// Return the shape of the array and the maximum length of the fields
// in the first row
PyObject* cPgNumpy::_get_array_field_shape(
        int fieldnum, 
        npy_intp& totlen,
        int32_t& max_flen)
{

    PyObject* dims=NULL;
    char* mptr = PQgetvalue(mResult, 0, fieldnum); 

    // 32-bit integer, 4-bytes
    int32_t ndim = _tohost( *(int32_t *) mptr );

    // skip over ndim we just read
    mptr += 4;
    // Skip past the zero that is always there, and the array signature
    mptr += 2*4;

    totlen=1;
    if (ndim > 0) {
        dims = PyTuple_New(ndim);

        for (int32_t i=0;i<ndim;i++) {
            int32_t tdim = _tohost( *(int32_t *) mptr );
            PyTuple_SetItem(dims, i, PyLong_FromLong(tdim));

            totlen*= tdim;
            // skip past this and the "one" that is always there
            mptr += 2*4;
        }

        // Extract the largest of the fields.  Useful for guessing at
        // field lengthts of character arrays
        int32_t tflen=0;
        max_flen=0;
        for (int32_t i=0; i<totlen; i++) {
            tflen = _tohost( *(int32_t *) mptr );
            if (tflen > max_flen)
                max_flen=tflen;
            mptr += 4 + tflen;
        }
    }

    return(dims);

}

// This is for those cases where we can't just get it from the type
// definition
int32_t cPgNumpy::_pg_string_nbytes(
        string& fname, 
        int fnum,          // field number
        int ftype,         // postgres field type number
        int isarray, 
        int32_t max_flen)  // max_flen is the longest of the first row for arrays, zero otherwise
{
    int32_t nbytes=0;
    int32_t textmax=10;

    // First see if it was sent in the length dictionary
    int32_t flen = _long_from_dict(mFlenDict, fname.c_str());
    if (debug) cout<<"    Got flen="<<flen<<" from dictionary"<<endl;
    if (flen > 0) {
        // User the user specified size
        nbytes=flen;

    } else {

        // no user-specified size
        if (isarray) {
            flen = max_flen;
        } else {
            flen = PQgetlength(mResult, 0, fnum);
        }

        int fmod = PQfmod(mResult, fnum) - 4; // minus 4 for lenght field?

        if (debug) { 
			cout<<"    Got flen="<<flen<<" for character field"<<endl;
		}
        if (debug) { 
			cout<<"    Got fmod="<<fmod<<" for character field"<<endl;
		}
        if (ftype == 1042 || ftype == 1014 ) {
            // fixed length fields
            nbytes = flen;
        } else {
            // this is good for varchar(n) fields
            if (fmod >= flen) {
                _print_debug("      varchar(n) field found");
                nbytes=fmod;
            } else {

                // For scalars, we can run through all the rows and get the longest
                // might also be possible for arrays, need to look into it
                if (!isarray) {
                    if (debug) {
                        cout<<"Getting string size from all rows for column: '"<<fname<<"'\n";
                        fflush(stdout);
                    }
                    flen=0;
                    for (long long row=0;row<mNtuples; row++) {
                        int32_t t_nbytes = PQgetlength(mResult, row, fnum);
                        if (t_nbytes > flen) {
                            flen = t_nbytes;
                        }
                    }
                }

                if (flen > 0) {
                    nbytes=flen;
				} else {
					nbytes=textmax;
				}
                _print_debug("      text field found; size might be too small");
            }
        }
    } // length was not sent in the field length dictionary

    if (debug) {
		cout<<"  Using "<<nbytes<<" bytes for character field"<<endl;
	}
    return(nbytes);

}




long int cPgNumpy::_long_from_dict(
		PyObject* dict,
		const char* name)
{
	if (dict != NULL) {
		PyObject* item = PyDict_GetItemString(dict, name);
		if (item == NULL) {
			return(-9999);
		} else {
			return( PyLong_AsLong(item) );
		}
	} else {
		return(-9999);
	}
}







// Clear a result
void cPgNumpy::clear()
{
    if (mResult != NULL) {
        _print_debug("clearing existing result");
        PQclear(mResult);
    }
    mResult=NULL;

    mResultStatus = PGRES_FATAL_ERROR;
    mNtuples=0;
    mNfields=0;
}

void cPgNumpy::close()
{
    clear();
	if (PQstatus(mConn) == CONNECTION_OK) {
        _print_debug("closing connection");
        PQfinish(mConn);
        mConn = NULL;
    }
}

void cPgNumpy::_print_flush(const char *text)
{
    if (text) {
		cout<<text<<"\n";
        fflush(stdout);
	}
}

void cPgNumpy::_print_debug(const char *text)
{
    if (debug) {
		_print_flush(text);
	}
}

void cPgNumpy::_set_defaults() {
    // binary data retrieval
    mBinary=1;
    // for printig we would set mBinary=0 and
    // might turn off array decorators
    mUseArrayDecorators=true;
    mFetchCount = 1000;

    mConnectInfo="";
    mConn=NULL;
    mResult=NULL;

    mResultStatus = PGRES_FATAL_ERROR;

    mFlenDict=NULL;

    mNtuples=0;
    mNfields=0;

    mDescr=NULL;
}

