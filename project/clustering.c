#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "kmpp.h"
#include "mat.h"
#include "nsc.h"
#include "jaccard.h"

void scan_input(size_t n, size_t d, PyObject * py_data, double * data, PyObject * py_labels, size_t * labels)
{
  size_t i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < d; ++j) {
      data[(d * i) + j] = PyFloat_AsDouble(PyList_GetItem(py_data, (i * d) + j));
    }
    labels[i] = PyLong_AsSize_t(PyList_GetItem(py_labels, i));
  }
}

static PyObject * nsc_and_kmpp(PyObject *self, PyObject *args)
{
  bool err = false;

  /* input parameters */
  size_t d = 0, k = 0, n = 0, m = 0;
  PyObject * py_data = NULL;
  PyObject * py_labels = NULL;
  PyObject * py_labels_tuple = NULL;
  double * data = NULL;
  size_t * orig_labels = NULL;

  /* generated values */
  size_t * kmpp_labels = NULL;
  size_t * nsc_labels = NULL;

  (void)self;

  if(!PyArg_ParseTuple(args, "IIIIOO:nsc_and_kmpp", &k, &n, &d, &m, &py_data, &py_labels)) {
    err = true; goto cleanup;
  }

  if ((data = mat_allocate(n, d)) == NULL) {
    err = true; goto cleanup;
  }

  if ((orig_labels = (size_t *)malloc(sizeof(size_t) * n)) == NULL) {
    err = true; goto cleanup;
    }
  memset(orig_labels, 0, sizeof(size_t) * n);

  scan_input(n, d, py_data, data, py_labels, orig_labels);

  if ((nsc_labels = (size_t *)malloc(sizeof(size_t) * n)) == NULL) {
    err = true; goto cleanup;
  }
  memset(nsc_labels, 0, sizeof(size_t) * n);

  k = normalized_spectral_clustering(k, n, d, data, 300, &nsc_labels);
  if (k == 0) {
    err = true; goto cleanup;
  }

  kmpp_labels = kmpp(k, n, d, data, 300);

  py_labels_tuple = PyTuple_New(5);
  PyTuple_SET_ITEM(py_labels_tuple, 0,
		   PyLong_FromSize_t(k));
  PyTuple_SET_ITEM(py_labels_tuple, 1,
		   PyByteArray_FromStringAndSize((const char *)kmpp_labels, sizeof(size_t) * n));
  PyTuple_SET_ITEM(py_labels_tuple, 2,
		   PyFloat_FromDouble(jaccard_measure(n, kmpp_labels, orig_labels)));
  PyTuple_SET_ITEM(py_labels_tuple, 3,
		   PyByteArray_FromStringAndSize((const char *)nsc_labels, sizeof(size_t) * n));
  PyTuple_SET_ITEM(py_labels_tuple, 4,
		   PyFloat_FromDouble(jaccard_measure(n, nsc_labels, orig_labels)));
 cleanup:
  if (nsc_labels != NULL) {
    free(nsc_labels);
  }

  if (kmpp_labels != NULL) {
    free(kmpp_labels);
  }

  if (orig_labels != NULL) {
    free(orig_labels);
  }

  if (data != NULL) {
    free(data);
  }

  if (err) {
    if (py_labels_tuple != NULL) {
      Py_DECREF(py_labels_tuple);
      py_labels_tuple = NULL;
    }
  }

  return py_labels_tuple;
}

static PyObject * clustering_nsc(PyObject *self, PyObject *args)
{
  bool err = false;

  /* input parameters */
  size_t d = 0, k = 0, n = 0, m = 0;
  PyObject * py_data = NULL;
  PyObject * py_labels = NULL;
  PyObject * py_labels_tuple = NULL;
  double * data = NULL;
  size_t * orig_labels = NULL;

  /* generated values */
  size_t * nsc_labels = NULL;

  (void)self;

  if(!PyArg_ParseTuple(args, "IIIIOO:clustering_nsc", &k, &n, &d, &m, &py_data, &py_labels)) {
    err = true; goto cleanup;
  }

  if ((data = mat_allocate(n, d)) == NULL) {
    err = true; goto cleanup;
  }

  if ((orig_labels = (size_t *)malloc(sizeof(size_t) * n)) == NULL) {
    err = true; goto cleanup;
  }
  memset(orig_labels, 0, sizeof(size_t) * n);

  scan_input(n, d, py_data, data, py_labels, orig_labels);

  if ((nsc_labels = (size_t *)malloc(sizeof(size_t) * n)) == NULL) {
    err = true; goto cleanup;
  }
  memset(nsc_labels, 0, sizeof(size_t) * n);

  k = normalized_spectral_clustering(k, n, d, data, m, &nsc_labels);
  if (k == 0) {
    err = true; goto cleanup;
  }

  py_labels_tuple = PyTuple_New(3);
  PyTuple_SET_ITEM(py_labels_tuple, 0,
		   PyLong_FromSize_t(k));
  PyTuple_SET_ITEM(py_labels_tuple, 1,
		   PyByteArray_FromStringAndSize((const char *)nsc_labels, sizeof(size_t) * n));
  PyTuple_SET_ITEM(py_labels_tuple, 2,
		   PyFloat_FromDouble(jaccard_measure(n, nsc_labels, orig_labels)));

 cleanup:
  if (nsc_labels != NULL) {
    free(nsc_labels);
  }

  if (orig_labels != NULL) {
    free(orig_labels);
  }

  if (data != NULL) {
    free(data);
  }

  if (err) {
    if (py_labels_tuple != NULL) {
      Py_DECREF(py_labels_tuple);
      py_labels_tuple = NULL;
    }
  }

  return py_labels_tuple;
}

static PyObject * clustering_kmpp(PyObject *self, PyObject *args)
{
  bool err = false;

  /* input parameters */
  size_t d = 0, k = 0, n = 0, m = 0;
  PyObject * py_data = NULL;
  PyObject * py_labels = NULL;
  PyObject * py_labels_tuple = NULL;
  double * data = NULL;
  size_t * orig_labels = NULL;

  /* generated values */
  size_t * labels = NULL;

  (void)self;

  if(!PyArg_ParseTuple(args, "IIIIOO:clustering_kmpp", &k, &n, &d, &m, &py_data, &py_labels)) {
    err = true; goto cleanup;
  }

  if ((data = mat_allocate(n, d)) == NULL) {
    err = true; goto cleanup;
  }

  if ((orig_labels = (size_t *)malloc(sizeof(size_t) * n)) == NULL) {
    err = true; goto cleanup;
  }
  memset(orig_labels, 0, sizeof(size_t) * n);

  scan_input(n, d, py_data, data, py_labels, orig_labels);

  labels = kmpp(k, n, d, data, m);

  if (labels == NULL) {
    err = true; goto cleanup;
  }

  py_labels_tuple = PyTuple_New(2);
  PyTuple_SET_ITEM(py_labels_tuple, 0,
		   PyByteArray_FromStringAndSize((const char *)labels, sizeof(size_t) * n));
  PyTuple_SET_ITEM(py_labels_tuple, 1,
		   PyFloat_FromDouble(jaccard_measure(n, labels, orig_labels)));

 cleanup:
  if (labels != NULL) {
    free(labels);
  }

  if (orig_labels != NULL) {
    free(orig_labels);
  }

  if (data != NULL) {
    free(data);
  }

  if (err) {
    if (py_labels_tuple != NULL) {
      Py_DECREF(py_labels_tuple);
      py_labels_tuple = NULL;
    }
  }

  return py_labels_tuple;
}

static PyMethodDef clustering_methods[] = {
    {"nsc_and_kmpp", (PyCFunction) nsc_and_kmpp, METH_VARARGS, ""},
    {"nsc", (PyCFunction) clustering_nsc, METH_VARARGS, ""},
    {"kmpp", (PyCFunction) clustering_kmpp, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "clustering", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    clustering_methods /* the PyMethodDef array from before containing the methods of the extension */
};


PyMODINIT_FUNC PyInit_clustering(void)
{
    return PyModule_Create(&moduledef);
}
