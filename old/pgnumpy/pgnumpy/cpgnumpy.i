%module cpgnumpy
%{
#include "cpgnumpy.h"
%}

%typemap(throws) const char * %{
    PyErr_SetString(PyExc_RuntimeError, $1);
    SWIG_fail;
%}

class cPgNumpy {

    public:
        cPgNumpy() throw (const char*);
        cPgNumpy(string conninfo) throw (const char*);
        cPgNumpy(const char* conninfo) throw (const char*);
        ~cPgNumpy();

        void open() throw (const char*);
        void open(string conninfo) throw (const char*);
        void open(const char* conninfo) throw (const char*);

        // return results as text.  Only useful if using
        // the write() method to write to a file or stdout
        void use_text();
        // should we use array decorators for array writing?
        void set_decorators(bool decorate);
        // for when doing an execwrite
        void set_fetch_count(long long fetchcount);

        void execute(string query_string) throw (const char*);
        void execute(const char* query_string) throw (const char*);

        // Write the results to a file or stdout. Requires running
        // use_text() to get results in text format
        long long write(const char* filename=NULL) throw (const char*);

        // Execute the query and use a cursor to retrieve the results in
        // chunks, writing to a file or by default stdout.  Unlike write(),
        // does *not* require use_text() to be set.
        //
        // This does not use the query execution infrastructure of PnNumpy,
        // but rather is rolled specially for the cursor
        //
        // The return value is the number of rows returned
        long long execwrite(
                       const char* query,
                       const char* filename=NULL) throw (const char*);

        long long ntuples() {
            return mNtuples;
        }
        long long nfields() {
            return mNfields;
        }
        int status() {
            return mResultStatus;
        }

        void set_field_lengths(PyObject* flenDict) throw (const char*);
        void clear_field_lengths() throw (const char*);

        // fetch results as an array
        PyObject* fetchall() throw (const char*);

        void clear();
        void close();
};
