#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "kmpp.h"
#include "log.h"
#include "mat.h"
#include "nsc.h"

int prj_scan_input(size_t n, size_t d, PyObject * py_data, double * data)
{
  size_t i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < d; ++j) {
      data[(d * i) + j] = PyFloat_AsDouble(PyList_GetItem(py_data, (i * d) + j));
    }
  }

  return 0;
}

static PyObject * prj_main(PyObject *self, PyObject *args)
{
    (void)self;

    size_t d = 0, k = 0, n = 0, m = 0;
    PyObject * py_data = NULL;
    PyObject * py_labels_tuple = NULL;
    double * data = NULL;
    size_t * kmpp_labels = NULL;
    size_t * nsc_labels = NULL;

    if(!PyArg_ParseTuple(args, "IIIIO:prj_main", &k, &n, &d, &m, &py_data)) {
	return NULL;
    }

    log_init();

    log_emit("starting");

    if ((data = mat_allocate(n, d)) == NULL)
      {
	return NULL;
      }

    if (prj_scan_input(n, d, py_data, data) < 0)
      {
	perror("prj_main: failed to scan input data");
	return NULL;
      }
    if ((nsc_labels = (size_t *) malloc (sizeof(size_t) * n)) == NULL) {
      return NULL;
    }
    memset(nsc_labels, 0, sizeof(size_t) * n);

    log_emit("running nsc");

    k = normalized_spectral_clustering(k, n, d, data, 300, &nsc_labels);
    if (k == 0) {
      perror("prj_main: normalized_spectral_clustering failed");
      return NULL;
    }

    log_emit("running kmpp");

    kmpp_labels = kmpp(k, n, d, data, 300);

    log_emit("returning labels");

    py_labels_tuple = PyTuple_New(2);
    PyTuple_SET_ITEM(py_labels_tuple, 0,
		     PyByteArray_FromStringAndSize((const char *)kmpp_labels, sizeof(size_t) * n));
    PyTuple_SET_ITEM(py_labels_tuple, 1,
		     PyByteArray_FromStringAndSize((const char *)nsc_labels, sizeof(size_t) * n));
      //    Py_RETURN_NONE;
    return py_labels_tuple;
}

static PyMethodDef prj_methods[] = {
    {"prj_main", (PyCFunction) prj_main, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "prj_lib", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    prj_methods /* the PyMethodDef array from before containing the methods of the extension */
};


PyMODINIT_FUNC PyInit_prj_lib(void)
{
    return PyModule_Create(&moduledef);
}
