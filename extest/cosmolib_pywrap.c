#include <Python.h>
#include "cosmolib.h"
#include <numpy/arrayobject.h> 

struct PyCosmoObject {
  PyObject_HEAD
  struct cosmo* cosmo;
};





static void
PyCosmoObject_dealloc(struct PyCosmoObject* self)
{
    if (self->cosmo != NULL) {
       free(self->cosmo);
    } 
    self->ob_type->tp_free((PyObject*)self);
}

/* this never seems to get called */
/*
static PyObject *
PyCosmoObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    struct PyCosmoObject* self;

    self = (struct PyCosmoObject *)type->tp_alloc(type, 0);
    self->cosmo=NULL;
    return (PyObject *)self;
    if (self != NULL) {
        printf("initializing cosmo\n");fflush(stdout);
        self->cosmo = cosmo_new();
        if (self->cosmo == NULL) {
            printf("failed to allocate Cosmo");
            Py_DECREF(self);
            return NULL;
        }
    } else {
        printf("failed to allocate PyCosmoObject");
    }

    return (PyObject *)self;
}
*/

static int
PyCosmoObject_init(struct PyCosmoObject* self, PyObject *args, PyObject *kwds)
{
    double DH;
    int flat;
    double omega_m, omega_l, omega_k;

    if (self->cosmo != NULL) {
        free(self->cosmo);
    }

    if (!PyArg_ParseTuple(args, 
                          (char*)"diddd", 
                          &DH, &flat, &omega_m, &omega_l, &omega_k)) {
        printf("failed to Parse init");
        return -1;
    }

    self->cosmo = cosmo_new(DH, flat, omega_m, omega_l, omega_k);
    if (self->cosmo == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate struct cosmo");
        return -1;
    }
    return 0;
}

static PyObject *
PyCosmoObject_repr(struct PyCosmoObject* self) {
    char repr[255];
    if (self->cosmo != NULL) {
        sprintf(repr, "flat:    %d\n"
                      "DH:      %f\n"
                      "omega_m: %f\n" 
                      "omega_l: %f\n" 
                      "omega_k: %f", 
                      self->cosmo->flat, 
                      self->cosmo->DH, 
                      self->cosmo->omega_m, 
                      self->cosmo->omega_l, 
                      self->cosmo->omega_k);
        return PyString_FromString(repr);
    }  else {
        return PyString_FromString("");
    }
}


// helper to get our arrays
PyObject* pyarray_getdouble(PyObject* op) {
    PyObject* darr;
    darr=PyArray_ContiguousFromAny(op, NPY_DOUBLE, 0, 0);
    return darr;
}



/*
   The wrapper methods and vectorizations.

   For the array inputs, the caller is responsible for making sure the input is
   an array, contiguous, of the right type
 */


static PyObject*
PyCosmoObject_ez_inverse(struct PyCosmoObject* self, PyObject* args) {
    double z;
    double ezinv;

    if (!PyArg_ParseTuple(args, (char*)"d", &z)) {
        return NULL;
    }

    ezinv = ez_inverse(self->cosmo, z);
    return PyFloat_FromDouble(ezinv);
}

static PyObject*
PyCosmoObject_ez_inverse_integral(struct PyCosmoObject* self, PyObject* args) {
    double zmin, zmax;
    double ezinv_int;

    if (!PyArg_ParseTuple(args, (char*)"dd", &zmin, &zmax)) {
        return NULL;
    }

    ezinv_int = ez_inverse_integral(self->cosmo, zmin, zmax);
    return PyFloat_FromDouble(ezinv_int);
}





// comoving distance and vectorizations
static PyObject*
PyCosmoObject_cdist(struct PyCosmoObject* self, PyObject* args) {
    double zmin, zmax;
    double d;

    if (!PyArg_ParseTuple(args, (char*)"dd", &zmin, &zmax)) {
        return NULL;
    }

    d = cdist(self->cosmo, zmin, zmax);
    return PyFloat_FromDouble(d);

}


static PyObject*
PyCosmoObject_cdist_vec1(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zminObj=NULL, *resObj=NULL;;
    double *zmin, zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"Od", &zminObj, &zmax)) {
        return NULL;
    }

    n = PyArray_SIZE(zminObj);
    zmin = (double* )PyArray_DATA(zminObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = self->cosmo->DH*ez_inverse_integral(self->cosmo, zmin[i], zmax); 
    }

    return resObj;

}

static PyObject*
PyCosmoObject_cdist_vec2(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zmaxObj=NULL, *resObj=NULL;;
    double zmin, *zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"dO", &zmin, &zmaxObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zmaxObj);
    zmax = (double* )PyArray_DATA(zmaxObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = self->cosmo->DH*ez_inverse_integral(self->cosmo, zmin, zmax[i]); 
    }

    return resObj;
}

static PyObject*
PyCosmoObject_cdist_2vec(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zmaxObj, *zminObj=NULL, *resObj=NULL;
    double *zmin, *zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"OO", &zminObj, &zmaxObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zminObj);
    zmin = (double* )PyArray_DATA(zminObj);
    zmax = (double* )PyArray_DATA(zmaxObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = self->cosmo->DH*ez_inverse_integral(self->cosmo, zmin[i], zmax[i]); 
    }

    return resObj;
}

// transverse comoving distance and vectorizations
static PyObject*
PyCosmoObject_tcdist(struct PyCosmoObject* self, PyObject* args) {
    double zmin, zmax;
    double d;

    if (!PyArg_ParseTuple(args, (char*)"dd", &zmin, &zmax)) {
        return NULL;
    }

    d = tcdist(self->cosmo, zmin, zmax);
    return PyFloat_FromDouble(d);

}

static PyObject*
PyCosmoObject_tcdist_vec1(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zminObj=NULL, *resObj=NULL;;
    double *zmin, zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"Od", &zminObj, &zmax)) {
        return NULL;
    }

    n = PyArray_SIZE(zminObj);
    zmin = (double* )PyArray_DATA(zminObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = tcdist(self->cosmo, zmin[i], zmax); 
    }

    return resObj;

}

static PyObject*
PyCosmoObject_tcdist_vec2(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zmaxObj=NULL, *resObj=NULL;;
    double zmin, *zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"dO", &zmin, &zmaxObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zmaxObj);
    zmax = (double* )PyArray_DATA(zmaxObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = tcdist(self->cosmo, zmin, zmax[i]); 
    }

    return resObj;
}

static PyObject*
PyCosmoObject_tcdist_2vec(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zmaxObj, *zminObj=NULL, *resObj=NULL;
    double *zmin, *zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"OO", &zminObj, &zmaxObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zminObj);
    zmin = (double* )PyArray_DATA(zminObj);
    zmax = (double* )PyArray_DATA(zmaxObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = tcdist(self->cosmo, zmin[i], zmax[i]); 
    }

    return resObj;
}


// Angular diameter distance
static PyObject*
PyCosmoObject_angdist(struct PyCosmoObject* self, PyObject* args) {
    double zmin, zmax;
    double d;

    if (!PyArg_ParseTuple(args, (char*)"dd", &zmin, &zmax)) {
        return NULL;
    }

    d = angdist(self->cosmo, zmin, zmax);
    return PyFloat_FromDouble(d);

}

static PyObject*
PyCosmoObject_angdist_vec1(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zminObj=NULL, *resObj=NULL;;
    double *zmin, zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"Od", &zminObj, &zmax)) {
        return NULL;
    }

    n = PyArray_SIZE(zminObj);
    zmin = (double* )PyArray_DATA(zminObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = angdist(self->cosmo, zmin[i], zmax); 
    }

    return resObj;

}

static PyObject*
PyCosmoObject_angdist_vec2(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zmaxObj=NULL, *resObj=NULL;;
    double zmin, *zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"dO", &zmin, &zmaxObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zmaxObj);
    zmax = (double* )PyArray_DATA(zmaxObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = angdist(self->cosmo, zmin, zmax[i]); 
    }

    return resObj;
}

static PyObject*
PyCosmoObject_angdist_2vec(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zmaxObj, *zminObj=NULL, *resObj=NULL;
    double *zmin, *zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"OO", &zminObj, &zmaxObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zminObj);
    zmin = (double* )PyArray_DATA(zminObj);
    zmax = (double* )PyArray_DATA(zmaxObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = angdist(self->cosmo, zmin[i], zmax[i]); 
    }

    return resObj;
}


// luminosity distance
static PyObject*
PyCosmoObject_lumdist(struct PyCosmoObject* self, PyObject* args) {
    double zmin, zmax;
    double d;

    if (!PyArg_ParseTuple(args, (char*)"dd", &zmin, &zmax)) {
        return NULL;
    }

    d = lumdist(self->cosmo, zmin, zmax);
    return PyFloat_FromDouble(d);

}

static PyObject*
PyCosmoObject_lumdist_vec1(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zminObj=NULL, *resObj=NULL;;
    double *zmin, zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"Od", &zminObj, &zmax)) {
        return NULL;
    }

    n = PyArray_SIZE(zminObj);
    zmin = (double* )PyArray_DATA(zminObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = lumdist(self->cosmo, zmin[i], zmax); 
    }

    return resObj;

}

static PyObject*
PyCosmoObject_lumdist_vec2(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zmaxObj=NULL, *resObj=NULL;;
    double zmin, *zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"dO", &zmin, &zmaxObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zmaxObj);
    zmax = (double* )PyArray_DATA(zmaxObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = lumdist(self->cosmo, zmin, zmax[i]); 
    }

    return resObj;
}

static PyObject*
PyCosmoObject_lumdist_2vec(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zmaxObj, *zminObj=NULL, *resObj=NULL;
    double *zmin, *zmax, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"OO", &zminObj, &zmaxObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zminObj);
    zmin = (double* )PyArray_DATA(zminObj);
    zmax = (double* )PyArray_DATA(zmaxObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = lumdist(self->cosmo, zmin[i], zmax[i]); 
    }

    return resObj;
}

// Comoving volume element and vectorization
static PyObject*
PyCosmoObject_dV(struct PyCosmoObject* self, PyObject* args) {
    double z;
    double dv;

    if (!PyArg_ParseTuple(args, (char*)"d", &z)) {
        return NULL;
    }

    dv = dV(self->cosmo, z);
    return PyFloat_FromDouble(dv);

}

static PyObject*
PyCosmoObject_dV_vec(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zObj=NULL, *resObj=NULL;;
    double *z, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"O", &zObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zObj);
    z = (double* )PyArray_DATA(zObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = dV(self->cosmo, z[i]); 
    }

    return resObj;

}

// Comoving volume between zmin and zmax
static PyObject*
PyCosmoObject_V(struct PyCosmoObject* self, PyObject* args) {
    double zmin, zmax;
    double v;

    if (!PyArg_ParseTuple(args, (char*)"dd", &zmin, &zmax)) {
        return NULL;
    }

    v = volume(self->cosmo, zmin, zmax);
    return PyFloat_FromDouble(v);

}



// Inverse critical density
static PyObject*
PyCosmoObject_scinv(struct PyCosmoObject* self, PyObject* args) {
    double zl, zs;
    double d;

    if (!PyArg_ParseTuple(args, (char*)"dd", &zl, &zs)) {
        return NULL;
    }

    d = scinv(self->cosmo, zl, zs);
    return PyFloat_FromDouble(d);

}

static PyObject*
PyCosmoObject_scinv_vec1(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zlObj=NULL, *resObj=NULL;;
    double *zl, zs, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"Od", &zlObj, &zs)) {
        return NULL;
    }

    n = PyArray_SIZE(zlObj);
    zl = (double* )PyArray_DATA(zlObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = scinv(self->cosmo, zl[i], zs); 
    }

    return resObj;

}

static PyObject*
PyCosmoObject_scinv_vec2(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zsObj=NULL, *resObj=NULL;;
    double zl, *zs, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"dO", &zl, &zsObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zsObj);
    zs = (double* )PyArray_DATA(zsObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = scinv(self->cosmo, zl, zs[i]); 
    }

    return resObj;
}

static PyObject*
PyCosmoObject_scinv_2vec(struct PyCosmoObject* self, PyObject* args) {
    PyObject* zsObj, *zlObj=NULL, *resObj=NULL;
    double *zl, *zs, *res;
    npy_intp n, i;

    if (!PyArg_ParseTuple(args, (char*)"OO", &zlObj, &zsObj)) {
        return NULL;
    }

    n = PyArray_SIZE(zlObj);
    zl = (double* )PyArray_DATA(zlObj);
    zs = (double* )PyArray_DATA(zsObj);

    resObj = PyArray_ZEROS(1, &n, NPY_FLOAT64, 0);
    res = (double* )PyArray_DATA(resObj);

    for (i=0; i<n; i++) {
        res[i] = scinv(self->cosmo, zl[i], zs[i]); 
    }

    return resObj;
}




static PyMethodDef PyCosmoObject_methods[] = {
    {"ez_inverse", (PyCFunction)PyCosmoObject_ez_inverse, METH_VARARGS, "ez_inverse(z)\n\nGet 1/E(z)"},
    {"ez_inverse_integral", (PyCFunction)PyCosmoObject_ez_inverse_integral, METH_VARARGS, "ez_inverse_integral(zmin, zmax)\n\nGet integral of 1/E(z) from zmin to zmax"},
    {"cdist", (PyCFunction)PyCosmoObject_cdist, METH_VARARGS, "cdist(zmin,zmax)\n\nComoving distance between zmin and zmax"},
    {"cdist_vec1", (PyCFunction)PyCosmoObject_cdist_vec1, METH_VARARGS, "cdist_vec1(zmin,zmax)\n\nComoving distance between zmin(array) and zmax"},
    {"cdist_vec2", (PyCFunction)PyCosmoObject_cdist_vec2, METH_VARARGS, "cdist_vec2(zmin,zmax)\n\nComoving distance between zmin and zmax(array)"},
    {"cdist_2vec", (PyCFunction)PyCosmoObject_cdist_2vec, METH_VARARGS, "cdist_2vec(zmin,zmax)\n\nComoving distance between zmin and zmax both arrays"},
    {"tcdist", (PyCFunction)PyCosmoObject_tcdist, METH_VARARGS, "tcdist(zmin,zmax)\n\nTransverse comoving distance between zmin and zmax"},
    {"tcdist_vec1", (PyCFunction)PyCosmoObject_tcdist_vec1, METH_VARARGS, "tcdist_vec1(zmin,zmax)\n\nTransverse Comoving distance between zmin(array) and zmax"},
    {"tcdist_vec2", (PyCFunction)PyCosmoObject_tcdist_vec2, METH_VARARGS, "tcdist_vec2(zmin,zmax)\n\nTransverse Comoving distance between zmin and zmax(array)"},
    {"tcdist_2vec", (PyCFunction)PyCosmoObject_tcdist_2vec, METH_VARARGS, "tcdist_2vec(zmin,zmax)\n\nTransverse Comoving distance between zmin and zmax both arrays"},
    {"angdist", (PyCFunction)PyCosmoObject_angdist, METH_VARARGS, "angdist(zmin,zmax)\n\nAngular diameter distance distance between zmin and zmax"},
    {"angdist_vec1", (PyCFunction)PyCosmoObject_angdist_vec1, METH_VARARGS, "angdist_vec1(zmin,zmax)\n\nAngular diameter distance distance between zmin(array) and zmax"},
    {"angdist_vec2", (PyCFunction)PyCosmoObject_angdist_vec2, METH_VARARGS, "angdist_vec2(zmin,zmax)\n\nAngular diameter distance distance between zmin and zmax(array)"},
    {"angdist_2vec", (PyCFunction)PyCosmoObject_angdist_2vec, METH_VARARGS, "angdist_2vec(zmin,zmax)\n\nAngular diameter distance distance between zmin and zmax both arrays"},
    {"lumdist", (PyCFunction)PyCosmoObject_lumdist, METH_VARARGS, "lumdist(zmin,zmax)\n\nLuminosity distance distance between zmin and zmax"},
    {"lumdist_vec1", (PyCFunction)PyCosmoObject_lumdist_vec1, METH_VARARGS, "lumdist_vec1(zmin,zmax)\n\nLuminosity distance distance between zmin(array) and zmax"},
    {"lumdist_vec2", (PyCFunction)PyCosmoObject_lumdist_vec2, METH_VARARGS, "lumdist_vec2(zmin,zmax)\n\nLuminosity distance distance between zmin and zmax(array)"},
    {"lumdist_2vec", (PyCFunction)PyCosmoObject_lumdist_2vec, METH_VARARGS, "lumdist_2vec(zmin,zmax)\n\nLuminosity distance distance between zmin and zmax both arrays"},
    {"dV", (PyCFunction)PyCosmoObject_dV, METH_VARARGS, "dV(z)\n\nComoving volume element at redshift z"},
    {"dV_vec", (PyCFunction)PyCosmoObject_dV_vec, METH_VARARGS, "dV(z)\n\nComoving volume element at redshift z(array)"},
    {"V", (PyCFunction)PyCosmoObject_V, METH_VARARGS, "volume(z)\n\nComoving volume between zmin and zmax"},
    {"scinv", (PyCFunction)PyCosmoObject_scinv, METH_VARARGS, "scinv(zl,zs)\n\nInverse critical density distance between zl and zs"},
    {"scinv_vec1", (PyCFunction)PyCosmoObject_scinv_vec1, METH_VARARGS, "scinv_vec1(zl,zs)\n\nInverse critical density distance between zl(array) and zs"},
    {"scinv_vec2", (PyCFunction)PyCosmoObject_scinv_vec2, METH_VARARGS, "scinv_vec2(zl,zs)\n\nInverse critical density distance between zl and zs(array)"},
    {"scinv_2vec", (PyCFunction)PyCosmoObject_scinv_2vec, METH_VARARGS, "scinv_2vec(zl,zs)\n\nInverse critical density distance between zl and zs both arrays"},

    {NULL}  /* Sentinel */
};





static PyTypeObject PyCosmoType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "cosmolib.cosmo",             /*tp_name*/
    sizeof(struct PyCosmoObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PyCosmoObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    //0,                         /*tp_repr*/
    (reprfunc)PyCosmoObject_repr,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Cosmology Structure Pointer",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    PyCosmoObject_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    //0,     /* tp_init */
    (initproc)PyCosmoObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    //PyCosmoObject_new,                 /* tp_new */
    PyType_GenericNew,                 /* tp_new */
};


static PyMethodDef cosmotype_methods[] = {
    {NULL}  /* Sentinel */
};


#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initcosmolib(void) 
{
    PyObject* m;

    PyCosmoType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PyCosmoType) < 0)
        return;

    m = Py_InitModule3("cosmolib", cosmotype_methods, "Define cosmo type and methods.");

    Py_INCREF(&PyCosmoType);
    PyModule_AddObject(m, "cosmo", (PyObject *)&PyCosmoType);

    import_array();
}
