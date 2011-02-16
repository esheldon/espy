#include <Python.h>
#include <iostream>
#include "numpy/arrayobject.h"

#include <setjmp.h>
#include <signal.h>

#ifndef _tarray_h
#define _tarray_h

class tarray {
	public:
		tarray() throw (const char *);
		PyObject* makearray(unsigned long long size) throw (const char *);

		int getint() {
			return 25;
		};

        // Test signal handler
        void testsig();

        // This must be defined fully or linking will fail.
		~tarray() {};
};

#endif
