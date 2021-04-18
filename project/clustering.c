/*
 * clustering.c
 * Python C-API extension for clustering real valued matrices
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "kmpp.h"
#include "mat.h"
#include "nsc.h"
#include "jaccard.h"

enum clustering_algorithm {
			   ALGORITHM_NORMALIZED_SPECTRAL_CLUSTERING,
			   ALGORITHM_KMEANS_PP,
};

/*
 * Fills data and labels with entries from the given respective Python lists
 */
static void scan_py_input(size_t n, size_t d, PyObject * py_data, double * data, PyObject * py_labels, size_t * labels)
{
  size_t i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < d; ++j) {
      data[(d * i) + j] = PyFloat_AsDouble(PyList_GetItem(py_data, (i * d) + j));
    }
    labels[i] = PyLong_AsSize_t(PyList_GetItem(py_labels, i));
  }
}

/*
 * Returns a Python 3-tuple of (K, Labels, Jaccard-Measure) where:
 * K is the number of given or inferred clusters;
 * Labels is an array of resulting labels for the given data points
 * computed with the given algorithm, packed as a Python bytes-array object; and
 * Jaccard-Measure is the Jaccard similiarity measure for the given labels and
 * the resulting labels.
 */
static PyObject * clustering_gen(PyObject * args, enum clustering_algorithm alg)
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

  if(!PyArg_ParseTuple(args, "IIIIOO:clustering_gen", &k, &n, &d, &m, &py_data, &py_labels)) {
    err = true; goto cleanup;
  }

  if ((data = mat_allocate(n, d)) == NULL) {
    err = true; goto cleanup;
  }

  if ((orig_labels = (size_t *)malloc(sizeof(size_t) * n)) == NULL) {
    err = true; goto cleanup;
  }
  memset(orig_labels, 0, sizeof(size_t) * n);

  scan_py_input(n, d, py_data, data, py_labels, orig_labels);

  switch (alg) {

  case ALGORITHM_NORMALIZED_SPECTRAL_CLUSTERING:
    k = normalized_spectral_clustering(k, n, d, data, m, &labels);
    break;

  case ALGORITHM_KMEANS_PP:
    k = k_means_pp(k, n, d, data, m, &labels);
    break;

  default:
    err = true; goto cleanup;
    break;
  }

  if (k == 0) {
    err = true; goto cleanup;
  }

  py_labels_tuple = PyTuple_New(3);
  PyTuple_SET_ITEM(py_labels_tuple, 0,
		   PyLong_FromSize_t(k));
  PyTuple_SET_ITEM(py_labels_tuple, 1,
		   PyByteArray_FromStringAndSize((const char *)labels, sizeof(size_t) * n));
  PyTuple_SET_ITEM(py_labels_tuple, 2,
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

static PyObject * clustering_nsc(PyObject *self, PyObject *args)
{
  (void)self;
  return clustering_gen(args, ALGORITHM_NORMALIZED_SPECTRAL_CLUSTERING);
}

static PyObject * clustering_kmpp(PyObject *self, PyObject *args)
{
  (void)self;
  return clustering_gen(args, ALGORITHM_KMEANS_PP);
}

static PyMethodDef clustering_methods[] = {
    {"nsc", (PyCFunction) clustering_nsc, METH_VARARGS, ""},
    {"kmpp", (PyCFunction) clustering_kmpp, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "clustering",
    NULL,
    -1,
    clustering_methods,
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC PyInit_clustering(void)
{
    return PyModule_Create(&moduledef);
}
