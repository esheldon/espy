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


void PrintSomething(const char *text)
{
    if (text)
        cout<<text<<endl;fflush(stdout);
}

void PrintDebug(const char *text)
{
    if (debug) PrintSomething(text);
}

// convert from network to host byte ordering
void Network2Host(char* data, int nbytes)
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
void FillScalar(
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

    if (debug) 
    {
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
    if (!isstring)
        Network2Host(output, nbytes);
}



void FillArray(
        char* input,    // Input postgresql data 
        char* output,   // Output array data
        int nbytes,     // Number of bytes available in output data field
        int isstring)   // Is this a string?
{
    
    // This is the most generic possible version of this program.  We can
    // probably speed it up later

    npy_int32 ndim = ntohl( *( (npy_int32 *) input) );
    input += 4;

    if (debug) cout<<"    Ndim = "<<ndim<<endl;

    // Skip past the zero that is always there, and the array signature
    input += 2*4;

    npy_int32 n_elements=1;
    npy_int32 tdim=0;
    for (npy_int32 i=0; i<ndim; i++)
    {
        tdim = ntohl( *( (npy_int32 *) input) );
        n_elements *= tdim;
        // skip past this and the "one" that is always there
        input += 2*4;
    }

    if (debug) cout<<"        Nel = "<<n_elements<<endl;

    
    npy_int32 flen=0;
    for (npy_int32 i=0; i< n_elements; i++)
    {

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


// Should we implement blobs?
void StoreBinary(
        char* input,    // Input postgresql data 
        int flen,       // length of this field in the input data
        char* output,   // Output array data
        int nbytes,     // number of available bytes in the output array
        int isarray,    // is this an array?
        int isstring)   // is this a string?
{
    if (isarray)
        FillArray(input, output, nbytes, isstring);
    else
        FillScalar(input, flen, output, nbytes, isstring);
}






// Return a string containing the type of the field, the number of bytes
// if known, and if this is an array
void PgType2NpTypeString(int pgtype, string &nptype, int& nbytes, int& isarray)
{
    switch (pgtype)
    {
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
    }

}


// Return the shape of the array and the maximum length of the fields
// in the first row
PyObject* GetArrayFieldShape(
        PGresult *res, 
        int fieldnum, 
        npy_intp& totlen,
        npy_int32& max_flen)
{

    PyObject* dims=NULL;
    char* mptr = PQgetvalue(res, 0, fieldnum); 

    // 32-bit integer, 4-bytes
    npy_int32 ndim = 
        (npy_int32) ntohl( *(npy_int32 *) mptr );

    // skip over ndim we just read
    mptr += 4;
    // Skip past the zero that is always there, and the array signature
    mptr += 2*4;

    totlen=1;
    if (ndim > 0)
    {
        dims = PyTuple_New(ndim);

        for (npy_int32 i=0;i<ndim;i++)
        {
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
        for (npy_int32 i=0; i<totlen; i++)
        {
            tflen = ntohl(  *(npy_int32 *) mptr );
            if (tflen > max_flen)
                max_flen=tflen;
            mptr += 4 + tflen;
        }
    }

    return(dims);

}

long int LongFromDict(
        PyObject* dict,
        const char* name)
{
    PyObject* item = PyDict_GetItemString(dict, name);
    if (item == NULL)
        return(-9999);
    else
    {
        return( PyLong_AsLong(item) );
    }
}

// Make a desciption list, or format list, one element for each field
// in the output array
PyObject* MakeDescrList(
        PGresult *res, 
        PyObject* flenDict,
        FIELD_INFO& finfo)
{

    int textmax=10;

    // Make a list of tuples
    PyObject* dlist=PyList_New(0);
    
    int nfields = PQnfields(res);

    finfo.nbytes.resize(nfields,0);
    finfo.isarray.resize(nfields,0);
    finfo.nel.resize(nfields,0);
    finfo.isstring.resize(nfields,0);

    for (int i=0; i<nfields; i++)
    {

        string fname=PQfname(res,i);
        if (debug) cout<<"Working on field: "<<fname<<endl;

        // Determine the type, including if it is an array
        string nptype;
        int nbytes=0;
        int isarray=0;
        int ftype = PQftype(res, i);

        PgType2NpTypeString(ftype, nptype, nbytes, isarray);


        int dsize;
        if (isarray)
            dsize=3;
        else
            dsize=2;


        // Tuple holding descr info for this field
        PyObject* dtup = PyTuple_New(dsize);
        npy_int32 max_flen=0;

        if (isarray) 
        {
            // We need to describe the shape for arrays
            PyObject* dims = GetArrayFieldShape(res, i, finfo.nel[i], max_flen);
            PrintDebug("Setting dims in descr tuple");
            PyTuple_SetItem(dtup, 2, dims);
        }
        else
        {
            finfo.nel[i] = 1;
        }


        // We have a few options here.  There is the size of the actual
        // field, that is what we get from PQgetlength().  This is perfectly
        // good for character(n) fixed length fields.  For varchar(n) we
        // should try to find the max length.  
        if (nbytes == -1)
        {

            npy_int32 flen=0;
            // First see if it was sent in the length dictionary
            flen = LongFromDict(flenDict, fname.c_str());
            if (debug) cout<<"  Got flen="<<flen<<" from dictionary"<<endl;
            if (flen > 0)
            {
                nbytes=flen;
            }
            else
            {

                if (isarray)
                {
                    flen = max_flen;
                }
                else
                {
                    flen = PQgetlength(res, 0, i);
                }

                int fmod   = PQfmod(res, i) - 4; // minus 4 for lenght field?

                if (debug) 
                    cout<<"  Got flen="<<flen<<" for character field"<<endl;
                if (debug) 
                    cout<<"  Got fmod="<<fmod<<" for character field"<<endl;
                if (ftype == 1042 || ftype == 1014 )
                {
                    // fixed length fields
                    nbytes = flen;
                }
                else
                {
                    if (fmod >= flen)
                    {
                        PrintDebug("    varchar(n) field found");
                        nbytes=fmod;
                    } 
                    else
                    {
                        if (textmax > flen)
                            nbytes=textmax;
                        else
                            nbytes=flen;
                        PrintDebug("    text field found");
                    }
                }
            } // length was not sent in the field length dictionary

            if (debug) cout<<"  Using "<<nbytes<<
                " bytes for character field"<<endl;
        }

        finfo.nbytes[i] = nbytes;
        finfo.isarray[i] = isarray;



        // Set the length in the format
        if (nptype == "S")
        {
            stringstream sstr;
            sstr << nbytes;
            nptype=nptype+sstr.str();
            finfo.isstring[i] = 1;
        }
        PrintDebug(nptype.c_str());
        PyObject* ttype = PyString_FromString(nptype.c_str());
        PyObject* tname = PyString_FromString(fname.c_str());

        // Don't have to decref objects put into sequence.  They will get
        // decrefed when the sequence is decrefed
        PyTuple_SetItem(dtup, 0, tname);
        PyTuple_SetItem(dtup, 1, ttype);


        // Same goes for this tuple: no decref needed
        PyList_Append(dlist, dtup);
        
    }


    return(dlist);

}

// create the output array
PyObject* CreateOutputArray(PGresult* res, PyObject* flenDict)
{

    PyObject* ret;
    FIELD_INFO finfo;

    npy_intp ntuples = (npy_intp) PQntuples(res);
    PyObject* descrlist = MakeDescrList(res, flenDict, finfo);

    PyArray_Descr* descr;
    int dstat = PyArray_DescrConverter(descrlist, &descr);
    if (!dstat)
    {
        PyErr_SetString(PyExc_ValueError,"Error making array descriptor");
        return(NULL);  // sending NULL will trigger an exception
    } else 
    {
        npy_intp dims[1] = {0};
        dims[0] = ntuples;
        // pars: ndim, shape, descr, fortran_flag
        ret = PyArray_Zeros(1, dims, descr, NPY_FALSE);
    }

    // Note, we do NOT decref the descr since references are kept
    Py_XDECREF(descrlist);

    // Fill in the array
    PyArrayObject* retarray = (PyArrayObject* ) ret;
    char* data = retarray->data;

    npy_intp nfields = finfo.nbytes.size();
    for (npy_intp row=0; row<ntuples; row++)
    {
        for (npy_intp field=0; field<nfields; field++)
        {

            // For scalars, this is the actual length of this field, which may
            // be less or more than the bytes we have allocated.  For arrays
            // it the entire length of the array data including dimension
            // specifications, etc
            int flen   = PQgetlength(res, row, field);

            // Store the data
            StoreBinary(
                    PQgetvalue(res, row, field),
                    flen,
                    data,
                    finfo.nbytes[field], 
                    finfo.isarray[field], 
                    finfo.isstring[field]);

            data += finfo.nbytes[field]*finfo.nel[field];
        }
    }

    return(ret);
}

/* Complex checking of the query status */
int CheckStatus(PGresult *res)
{
    /* 
       For successful queries, there are two options: either the query returns
       results or not.  If it does not return data, but successfully
       completed, then the status will be set to PGRES_COMMAND_OK.  If success
       and returns tuples, then it will be set to PGRES_TUPLES_OK 
       */

    int status = PQresultStatus(res);

    if (PQresultStatus(res) == PGRES_COMMAND_OK)
    {
        /* Success but no results */
        cout<<PQcmdStatus(res);
        return(status);
    }

    if (PQresultStatus(res) != PGRES_TUPLES_OK)
    {
        switch (status)
        {
            /* No results */
            case PGRES_COPY_IN:        return(status);
            case PGRES_COPY_OUT:       return(status);
            case PGRES_BAD_RESPONSE:
                PyErr_SetString(PyExc_RuntimeError,PQresultErrorMessage(res));
                return(status);
            case PGRES_NONFATAL_ERROR: 
                PyErr_SetString(PyExc_RuntimeError,PQresultErrorMessage(res));
                return(status);
            case PGRES_FATAL_ERROR:    
                PyErr_SetString(PyExc_RuntimeError,PQresultErrorMessage(res));
                return(status);
            default: break;
        } 
    }

    /* Result is empty, either an error or this query returns no data */
    if ( PQnfields(res) == 0 && PQntuples(res) == 0 ) 
    {
        /* Query returns no data. Try to print a message to clarify */
        cout<<PQcmdStatus(res);
        return(status);
    }

    return(status);

}


// Send the query and return the results
PyObject* SendQuery(
        char* query_string, 
        char* connect_info, 
        PyObject* flenDict)
{
    int binary=1;
    PGconn *conn=NULL;
    PGresult *res=NULL;
    PyObject* ret;

    /* Attempt to establish the connection */
    conn = PQconnectdb(connect_info);

    if (PQstatus(conn) != CONNECTION_OK)
    {
        PQfinish(conn);
        PyErr_SetString(PyExc_RuntimeError,"Could not establish connection");
        return(NULL);  // sending NULL will trigger exception
    }

    if (debug) cout<<"Executing query: '"<<query_string<<"'"<<endl;
    res = PQexecParams(conn,
            query_string,
            0,    /* number of parameters (none) */
            NULL, /* see doc for the parameter info */
            NULL,
            NULL,
            NULL,
            binary); /* 0 for text, 1 for binary (network order) */

    int status = CheckStatus(res);

    // Either the query returns nothing or some error occurred
    if (status != PGRES_TUPLES_OK)
    {
        PQclear(res);
        PQfinish(conn);
        if (status == PGRES_FATAL_ERROR)
        {
            // Return NULL to raise the exception
            return(NULL);
        }
        else
        {
            Py_XINCREF(Py_None);
            return(Py_None);
        }
    }


    if ( PQntuples(res) > 0 )
    {
        // We have some results.  Copy into the output array.
        ret = CreateOutputArray(res, flenDict);
    }
    else
    {
        Py_XINCREF(Py_None);
        ret = Py_None;
    }

    PQclear(res);
    PQfinish(conn);

    return(ret);
}




// Parse the input and send the query.  Return any results.
static PyObject*
query(PyObject* ignored, PyObject* args, PyObject* kwds)
{

    char* query_string=NULL;
    char* connect_info=NULL;   
    PyObject* flenDict;         // optional format strings in a dict
    static char *kwlist[] = {
        "query_string",
        "connect_info",
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
                &connect_info,
                &PyDict_Type, &flenDict))
    {
        // returning NULL will raise an exception, the string of which
        // is set in the parser above.
        return(NULL);
    }

    // If not set to "" we get a Bus error from PQconnectdb
    // Note, if set to "", the environment variables are still examined
    // such as PGUSER and PGDATABASE.
    if (connect_info == NULL)
        connect_info = "";

    PyObject* ret = SendQuery(query_string, connect_info, flenDict);
    return(ret);
}


static PyMethodDef NpPgsqlMethods[] = {
    {"query", 
    (PyCFunction)query, METH_VARARGS|METH_KEYWORDS,
    "Module: readfields.\n" 
    "    Efficiently read fields from a file into a numpy array.  Allows\n" 
    "    reading a subset of the rows and fields. Currently only supports\n"
    "    binary files.\n"
    "\n"
    "Calling Sequence:\n"
    "    readfields(filename_or_obj, dtype, nrows, rows=, fields=)\n\n"
    "Inputs:\n"
    "    filename or file object: Filename in string form or opened \n"
    "        file object\n"
    "    dtype: Numpy dtype in list of tuples form.  See example below.\n" 
    "    nrows: The number of rows in the file. This should be relative\n" 
    "        to current position if a file object is entered\n"
    "Optional Inputs:\n"
    "    rows: An ordered, unique set of rows to read.  Must be convertible\n"
    "        to a numpy array\n"
    "    fields: A unique sequence of fields to read. Order does not matter.\n"
    "Example:\n"
    "    from readfields import readfields\n"
    "    dtype=[('id', 'i8'), ('x', 'f4'), ('y', 'f4'), ('flux','5f4')]\n"
    "    nrows=12235221\n"
    "    result = readfields('data.bin', dtype, nrows, rows=[35,10002,532251], fields=['id', 'flux'])\n"
    "MODIFICATION HISTORY:\n"
    "    Created, 2007-05-25, Erin Sheldon, NYU\n"},
    {NULL,NULL,0,NULL} // Sentinel
};

// Here interface is the name of the *module*, not the name of the 
// functions contained therein.
extern "C" PyMODINIT_FUNC
initinterface(void)
{
    if (debug) {cout<<"Initialization"<<endl;}
    (void) Py_InitModule("interface", NpPgsqlMethods);
    import_array(); // Must be present for numpy.
}



