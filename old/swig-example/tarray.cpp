#include "tarray.hpp"

tarray::tarray() throw (const char *) {
	// DONT FORGET THIS!!!!
	import_array();
}

PyObject* tarray::makearray(unsigned long long size)  throw (const char *) {

	PyArrayObject* outarr;
	int ndim=1;
	npy_intp intpsize = size;


	std::cout<<"Creating float64 output array with size: "<<size<<"\n";
	std::cout<<"ndim: "<<ndim<<"\n";
	std::cout<<"Dims value is: "<<intpsize<<"\n";
	std::cout<<"NPY_FLOAT64: "<<NPY_FLOAT64<<"\n";
	std::cout<<"NPY_DOUBLE: "<<NPY_DOUBLE<<"\n";
	std::cout<<"NPY_FALSE: "<<NPY_FALSE<<"\n";
	outarr = (PyArrayObject* ) PyArray_ZEROS(
			ndim, 
			&intpsize,
			NPY_FLOAT64,
			NPY_FALSE);

	std::cout<<"Done\n";
	if (outarr ==NULL) {
		throw "Could not allocate array";
	}
	

	// Make a pointer to the data area
	//mData = mReturnObject->data;

	std::cout<<"Returning data\n";
	return (PyObject* ) outarr;

}


/* Interrupt handler */

static jmp_buf jbuf;

/* ARGSUSED */
static void
onintr(int sig)
{
	longjmp(jbuf, 1);
}



void tarray::testsig() {

	PyOS_sighandler_t old_inthandler;

	old_inthandler = PyOS_setsig(SIGINT, onintr);

	if (setjmp(jbuf)) {
#ifdef HAVE_SIGRELSE
		/* This seems necessary on SunOS 4.1 (Rasmus Hahn) */
		sigrelse(SIGINT);
#endif
		PyOS_setsig(SIGINT, old_inthandler);
        printf("Interrupt found");
        return;
    }

    // just sleep
    sleep(1000);
	PyOS_setsig(SIGINT, old_inthandler);
}
