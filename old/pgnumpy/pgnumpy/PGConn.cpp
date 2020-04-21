#include <Python.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include "numpy/arrayobject.h"
#include "libpq-fe.h"

// For ntohl, etc
#include <netinet/in.h>


using namespace std;

#define debug 0
#define debug2 0


/* Macro to convert to host byte-ordering, 8-byte int or float */
#define ntohll(x) ((((npy_uint64)ntohl(x)) << 32) + ntohl(x >> 32))

// Hold some info about the fields in a result structure.  
typedef struct {
    // The number of bytes in each field.
    vector<npy_intp> nbytes;
    // The number of elements with the above length in bytes
    vector<npy_intp> nel;
    // Is this an array?
    vector<int> isarray;
    // Is this a string?  We don't count char(1) here
    vector<int> isstring; 
} FIELD_INFO;


class PGConn {

    public:
        PGConn();

        ~PGConn();

        int Open();
        int Open(char* conninfo);

        int Execute(char* query_string);
        int Execute(string query_string);

        void Clear();
        void Close();

        // Python related methods
        int PyCheckStatus();
        PyObject* PyFetchAsArray();

        void PyPgTypeInfo(
                int pgtype, 
                string &nptype, 
                int& nbytes, 
                int& isarray);
        void PySetFieldLengths(PyObject* flenDict);
        PyObject* PyGetArrayFieldShape(
                int fieldnum, 
                npy_intp& totlen,
                npy_int32& max_flen);
        long int LongFromDict(
                PyObject* dict,
                const char* name);


        npy_int32 GetStringNbytes(
                string& fname, 
                int fnum, 
                int ftype,
                int isarray, 
                npy_int32 max_flen);
        int PyMakeArrayDescr();

        //
        // For filling in the data
        // 

        // convert from network to host byte ordering
        void Network2Host(char* data, int nbytes);
        // Fill a scalar with the appropriate number of bytes.  Swap bytes if 
        // necessary
        void FillScalar(
                char* input,     // Input postgresql data 
                npy_int32 flen,  // Lenght of this field in the input data 
                char* output,    // Output array data
                int nbytes,      // number of bytes available in the output 
                int isstring);   // Is this a string?
        void FillArray(
                char* input,    // Input postgresql data 
                char* output,   // Output array data
                int nbytes,     // Number of bytes available in output 
                int isstring);  // Is this a string?

        void StoreBinary(
                char* input,    // Input postgresql data 
                int flen,       // length of this field in the input data
                char* output,   // Output array data
                int nbytes,     // number of available bytes in the output
                int isarray,    // is this an array?
                int isstring);  // is this a string?


        void PrintFlush(const char *text);
        void PrintDebug(const char *text);


    protected:

        int mBinary;  // always True

        FIELD_INFO mFieldInfo;
        string mConnectInfo;

        PGconn *mConn;
        PGresult* mResult;

        int mResultStatus;

        // User input field length dictionary
        PyObject* mFlenDict;

        // Number of tuples, fields
        npy_intp mNtuples;
        npy_intp mNfields;

        // description list
        PyArray_Descr* mDescr;

};

//Constructor and destructor
PGConn::PGConn()
{
    mConnectInfo="";
    mConn=NULL;
    mResult=NULL;

    mBinary=1; // binary data retrieval
    mResultStatus = PGRES_FATAL_ERROR;

    mFlenDict=NULL;
}
PGConn::~PGConn()
{
    Clear();
    Close();
}

// Open a connection
int PGConn::Open()
{
    return( Open("") );
}
int PGConn::Open(char* conninfo)
{
    if (conninfo != NULL)
        mConnectInfo = conninfo;

    mConn = PQconnectdb(mConnectInfo.c_str());

    if (PQstatus(mConn) != CONNECTION_OK) {
        PQfinish(mConn);
        PyErr_SetString(PyExc_RuntimeError,"Could not establish connection");
        return(0);  // sending NULL will trigger exception
    }

    return(1);
}

// Clear a result
void PGConn::Clear()
{
    if (mResult != NULL) {
        if (debug) {cout<<"Clearing result"<<endl;fflush(stdout);}
        PQclear(mResult);
        mResult=NULL;
    }
}

void PGConn::Close()
{
    //if (mConn != NULL)
	if (PQstatus(mConn) == CONNECTION_OK) {
        if (debug) {cout<<"Closing connection";fflush(stdout);}
        PQfinish(mConn);
        mConn = NULL;
        if (debug) {cout<<"...OK"<<endl;fflush(stdout);}
    }
}


// Execute a query
int PGConn::Execute(string query_string)
{
    return(  Execute(query_string.c_str())  );
}

int PGConn::Execute(char* query_string)
{
    Clear();
    if (debug) cout<<"Executing query: '"<<query_string<<"'"<<endl;
    mResult = PQexecParams(mConn,
            query_string,
            0,    /* number of parameters (none) */
            NULL, /* see doc for the parameter info */
            NULL,
            NULL,
            NULL,
            mBinary); /* 0 for text, 1 for binary (network order) */

    mResultStatus = PQresultStatus(mResult);
    mNtuples = (npy_intp) PQntuples(mResult);
    mNfields = (npy_intp) PQnfields(mResult);
    return(mResultStatus);
}

void PGConn::PrintFlush(const char *text)
{
    if (text) {
		cout<<text<<endl;fflush(stdout);
	}
}

void PGConn::PrintDebug(const char *text)
{
    if (debug) {
		PrintFlush(text);
	}
}




//
//
//  Python related routines
//
//

// Complex checking of the query status.  This is used to possibly set an
// error string for later raised exceptions
int PGConn::PyCheckStatus() {

    //   For successful queries, there are two options: either the query
    //   returns results or not.  If it does not return data, but successfully
    //   completed, then the status will be set to PGRES_COMMAND_OK.  If
    //   success and returns tuples, then it will be set to PGRES_TUPLES_OK 

    int mResultStatus = PQresultStatus(mResult);

    if (PGRES_COMMAND_OK == mResultStatus) {
        // Success but no results
        cout<<PQcmdStatus(mResult);
        return(mResultStatus);
    }

    if (PGRES_TUPLES_OK != mResultStatus) {
        switch (mResultStatus) {
            // No results
            case PGRES_COPY_IN:        return(mResultStatus);
            case PGRES_COPY_OUT:       return(mResultStatus);
            case PGRES_BAD_RESPONSE:
                PyErr_SetString(
                        PyExc_RuntimeError,
                        PQresultErrorMessage(mResult));
                return(mResultStatus);
            case PGRES_NONFATAL_ERROR: 
                PyErr_SetString(
                        PyExc_RuntimeError,
                        PQresultErrorMessage(mResult));
                return(mResultStatus);
            case PGRES_FATAL_ERROR:    
                PyErr_SetString(
                        PyExc_RuntimeError,
                        PQresultErrorMessage(mResult));
                return(mResultStatus);
            default: break;
        } 
    }

    // Result is empty, either an error or this query returns no data
    if ( PQnfields(mResult) == 0 && PQntuples(mResult) == 0 ) {
        // Query returns no data. Try to print a message to clarify
        cout<<PQcmdStatus(mResult);
        return(mResultStatus);
    }

    return(mResultStatus);

}

// Return a string containing the type of the field, the number of bytes
// if known, and if this is an array
void PGConn::PyPgTypeInfo(
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

        default: nptype="S";
				 nbytes=-1;
				 break;
    }

}


// Return the shape of the array and the maximum length of the fields
// in the first row
PyObject* PGConn::PyGetArrayFieldShape(
        int fieldnum, 
        npy_intp& totlen,
        npy_int32& max_flen)
{

    PyObject* dims=NULL;
    char* mptr = PQgetvalue(mResult, 0, fieldnum); 

    // 32-bit integer, 4-bytes
    npy_int32 ndim = 
        (npy_int32) ntohl( *(npy_int32 *) mptr );

    // skip over ndim we just read
    mptr += 4;
    // Skip past the zero that is always there, and the array signature
    mptr += 2*4;

    totlen=1;
    if (ndim > 0) {
        dims = PyTuple_New(ndim);

        for (npy_int32 i=0;i<ndim;i++) {
            npy_int32 tdim = 
                ntohl( *(npy_int32 *) mptr );
            PyTuple_SetItem(dims, i, PyLong_FromLong(tdim));

            totlen*= tdim;
            // skip past this and the "one" that is always there
            mptr += 2*4;
        }

        // Extract the largest of the fields.  Useful for guessing at
        // field lengthts of character arrays
        npy_int32 tflen=0;
        max_flen=0;
        for (npy_int32 i=0; i<totlen; i++) {
            tflen = ntohl(  *(npy_int32 *) mptr );
            if (tflen > max_flen)
                max_flen=tflen;
            mptr += 4 + tflen;
        }
    }

    return(dims);

}

long int PGConn::LongFromDict(
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


// This is for those cases where we can't just get it from the type
// definition
npy_int32 PGConn::GetStringNbytes(
        string& fname, 
        int fnum, 
        int ftype,
        int isarray, 
        npy_int32 max_flen)
{
    npy_int32 nbytes=0;
    npy_int32 textmax=10;

    // First see if it was sent in the length dictionary
    npy_int32 flen = LongFromDict(mFlenDict, fname.c_str());
    if (debug) cout<<"    Got flen="<<flen<<" from dictionary"<<endl;
    if (flen > 0) {
        // User the user specified size
        nbytes=flen;
    } else {
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
                PrintDebug("      varchar(n) field found");
                nbytes=fmod;
            } else {
                if (flen > 0) {
                    nbytes=flen;
				} else {
					nbytes=textmax;
				}
                PrintDebug("      text field found; size might be too small");
            }
        }
    } // length was not sent in the field length dictionary

    if (debug) {
		cout<<"  Using "<<nbytes<<" bytes for character field"<<endl;
	}
    return(nbytes);

}


// Make a desciption list, or format list, one element for each field
// in the output array
int PGConn::PyMakeArrayDescr()
{

    // Make a list of tuples.  This is temporary, we will decref it
    PyObject* dlist=PyList_New(0);
    

    mFieldInfo.nbytes.resize(mNfields,0);
    mFieldInfo.isarray.resize(mNfields,0);
    mFieldInfo.nel.resize(mNfields,0);
    mFieldInfo.isstring.resize(mNfields,0);

    for (int i=0; i<mNfields; i++) {

        string fname=PQfname(mResult,i);
        if (debug) cout<<endl<<"Working on field: "<<fname<<endl;

        // Determine the type, including if it is an array
        string nptype;
        int nbytes=0;
        int isarray=0;
        int ftype = PQftype(mResult, i);

        PyPgTypeInfo(ftype, nptype, nbytes, isarray);


        int dsize;
        if (isarray)
            dsize=3;
        else
            dsize=2;


        // Tuple holding descr info for this field
        PyObject* dtup = PyTuple_New(dsize);

        // We need to describe the shape for arrays max_flen is the maximum
        // size over all elements in this field for the first row
        npy_int32 max_flen=0;
        if (isarray) {
            PyObject* dims = PyGetArrayFieldShape(
                    i, 
                    mFieldInfo.nel[i], 
                    max_flen);
            PrintDebug("  Setting dims in descr tuple");
            PyTuple_SetItem(dtup, 2, dims);
        } else {
            mFieldInfo.nel[i] = 1;
        }

        // We have a few options here.  There is the size of the actual
        // field, that is what we get from PQgetlength().  This is perfectly
        // good for character(n) fixed length fields.  For varchar(n) we
        // should try to find the max length.  
        if (nbytes == -1) {
            nbytes = GetStringNbytes(fname, i, ftype, isarray, max_flen);
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

    int dstat = PyArray_DescrConverter(dlist, &mDescr);

    // We decref this but not the Descr
    Py_XDECREF(dlist);

    return(dstat);
}


// convert from network to host byte ordering
void PGConn::Network2Host(char* data, int nbytes)
{
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
                npy_int32* tptr = (npy_int32 *) data;
                *tptr = ntohl( *tptr );
            }
            break;
        case 8:
            {
                npy_int64* tptr = (npy_int64 *) data;
                *tptr = ntohll( *tptr );
            }
            break;
        default:
            break;
    }

}

// Fill a scalar with the appropriate number of bytes.  Swap bytes if 
// necessary
void PGConn::FillScalar(
        char* input,     // Input postgresql data 
        npy_int32 flen,  // Lenght of this field in the input data 
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
		Network2Host(output, nbytes);
	}
}



void PGConn::FillArray(
        char* input,    // Input postgresql data 
        char* output,   // Output array data
        int nbytes,     // Number of bytes available in output data field
        int isstring)   // Is this a string?
{
    
    // This is the most generic possible version of this program.  We can
    // probably speed it up later

    npy_int32 ndim = ntohl( *( (npy_int32 *) input) );
    input += 4;

    if (debug2) cout<<"    Ndim = "<<ndim<<endl;

    // Skip past the zero that is always there, and the array signature
    input += 2*4;

    npy_int32 n_elements=1;
    npy_int32 tdim=0;
    for (npy_int32 i=0; i<ndim; i++) {
        tdim = ntohl( *( (npy_int32 *) input) );
        n_elements *= tdim;
        // skip past this and the "one" that is always there
        input += 2*4;
    }

    if (debug2) cout<<"        Nel = "<<n_elements<<endl;

    
    npy_int32 flen=0;
    for (npy_int32 i=0; i< n_elements; i++) {

        // Read the size of this array element. This is always same for numbers
        flen = *(npy_int32 *) input;
        flen = ntohl(flen);
        input += 4;


        // store the number
        FillScalar(input, flen, output, nbytes, isstring);

        // Skip to the next element
        input += flen;
        output += nbytes; 
    }

}


void PGConn::StoreBinary(
        char* input,    // Input postgresql data 
        int flen,       // length of this field in the input data
        char* output,   // Output array data
        int nbytes,     // number of available bytes in the output array
        int isarray,    // is this an array?
        int isstring)   // is this a string?
{
    if (isarray) {
        FillArray(input, output, nbytes, isstring);
	} else {
		FillScalar(input, flen, output, nbytes, isstring);
	}
}



// create the output array
PyObject* PGConn::PyFetchAsArray()
{

    PyObject* ret;

    // Make the list of descr
    if (!PyMakeArrayDescr()) {
        return(NULL);  // sending NULL will trigger an exception
	}

    npy_intp dims[1] = {0};
    dims[0] = mNtuples;
    // pars: ndim, shape, descr, fortran_flag
    ret = PyArray_Zeros(1, dims, mDescr, NPY_FALSE);

    if (!ret) {
        // Still no decref of mDescr; it is stolen even on failure.
        return(ret);
    }

    // Creation of array was successfull
    // Note, we do NOT decref the mDescr since references are kept
    // In the array itself.  It "steals" the reference meaning it must take
    // care of it

    // Fill in the array
    PyArrayObject* retarray = (PyArrayObject* ) ret;
    char* data = retarray->data;

    for (npy_intp row=0; row<mNtuples; row++) {
        for (npy_intp field=0; field<mNfields; field++) {

            // For scalars, this is the actual length of this field, which may
            // be less or more than the bytes we have allocated.  For arrays
            // it the entire length of the array data including dimension
            // specifications, etc
            int flen   = PQgetlength(mResult, row, field);

            // Store the data
            StoreBinary(
                    PQgetvalue(mResult, row, field),
                    flen,
                    data,
                    mFieldInfo.nbytes[field], 
                    mFieldInfo.isarray[field], 
                    mFieldInfo.isstring[field]);

            data += mFieldInfo.nbytes[field]*mFieldInfo.nel[field];
        }
    }

    return(ret);
}

void PGConn::PySetFieldLengths(PyObject* flenDict)
{
    // Borrow a reference.  No need to decref.
    mFlenDict = flenDict;
}

// Parse the input and send the query.  Return any results.
static PyObject*
query(PyObject* ignored, PyObject* args, PyObject* kwds)
{

    PyObject* ret=NULL;         // The return value
    char* query_string=NULL;
    char* conninfo=NULL;   
    PyObject* flenDict=NULL;    // optional format strings in a dict
    static char *kwlist[] = {
        "query_string",
        "conninfo",
        "flengths",
        NULL};
    // s -> string
    // z -> string but may also be None, in which case it is copied here as 
    //    NULL
    // O! -> An object which must be of the specified type
    // To right side of | are optional args
    static char kwformat[] = "s|zO!:query";
   
    if (!PyArg_ParseTupleAndKeywords(
                args, kwds, kwformat, kwlist,
                &query_string, 
                &conninfo,
                &PyDict_Type, &flenDict)) {
        // returning NULL will raise an exception, the string of which
        // is set in the parser above.
        return(NULL);
    }

    if (debug) cout<<"Creating connection object"<<endl;
    PGConn conn;
    if (!conn.Open(conninfo))
        return(NULL);

    conn.PySetFieldLengths(flenDict);

    conn.Execute(query_string);

    int status = conn.PyCheckStatus();
    if (PGRES_TUPLES_OK != status) {
        if (PGRES_FATAL_ERROR == status) {
            // Return NULL to raise the exception
            ret=NULL;
        } else {
            Py_XINCREF(Py_None);
            ret = Py_None;
        }

    } else {
        ret = conn.PyFetchAsArray();
    }

    if (debug) fflush(stdout);

    return(ret);

}


static PyMethodDef NpPgsqlMethods[] = {
    {"query", 
    (PyCFunction)query, METH_VARARGS|METH_KEYWORDS,
    "Function: \n"
    "    pg_numpy.query()\n" 
    "\n"
    "    Execute a query and return the results.  If there are data\n"
    "    returned, they are returned as a NumPy array.  If there are\n"
    "    no results the return value is None.  This is pure C++ exported\n"
    "    from the PGConn module.\n"
    "\n"
    "    The query function is not designed for generic database\n"
    "    interaction, but primarily to efficiently get results into \n"
    "    a numpy array.  More generic DB-API compliant packages like \n"
    "    psycopg are more suitable when general and flexible operations \n"
    "    such as prepared statements are needed.\n"
    "\n"
    "    All basic numerical and character data types are supported.  \n"
    "    Types such as dates are returned as string types, but this \n"
    "    has not been fully tested.  Array columns are also supported.\n"
    "    with any dimensionality supported by both PostgreSql and NumPy\n"
    "\n"
    "Calling Sequence:\n"
    "    query(query_string, conninfo=None, flengths={})\n"
    "\n"
    "Inputs:\n"
    "    query_string:  A query string to be executed.  The query can be\n"
    "        any legal postgresql query.\n"
    "Optional Inputs:\n"
    "    conninfo:  Connection information string.  This is added to\n"
    "        or overrides any info defined in the user's environment such\n"
    "        as PGUSER or PGDATABASE.\n"
    "    flengths:  A dictionary specifying the length of columns in bytes.\n"
    "        This is useful for variable length character columns.  It is\n"
    "        ignored for numerical columns.  In the absence of this input\n"
    "        the size of variable length data is gotten from the first row\n"
    "        of the result, which may truncate later, longer data.\n"
    "\n"
    "Examples:\n"
    "    import numpypg as npg\n"
    "    a=npg.query('select x,y from mytable where x > 3')\n"
    "    print a['x']\n"
    "    prod = a['x']*a['y']\n"
    "\n"
    "    a=npg.query('select name, address from mytable',flengths={'name':20, 'address':100})\n"
    "\n"
    "Modification History:\n"
    "    Created, 2008-04-07, Erin Sheldon, NYU\n"},
    {NULL,NULL,0,NULL} // Sentinel
};

// Here interface is the name of the *module*, not the name of the 
// functions contained therein.
extern "C" PyMODINIT_FUNC
initPGConn(void)
{
    if (debug) {cout<<"Initialization"<<endl;}
    (void) Py_InitModule("PGConn", NpPgsqlMethods);
    import_array(); // Must be present for numpy.
}



