#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct km_ctx_s {
    size_t   dimension; /* dimension                                           */
    size_t   nclusters; /* # of clusters                                       */
    size_t   nobserves; /* total # of observations                             */
    size_t   max_iters; /* MAX_ITER                                            */
    size_t * ob_clusts; /* array of size N, observation index -> cluster index */
    size_t * cardinals; /* cardinalities of clusters, array of size K          */
    double * mean_vals; /* mean values, array of K * d                         */
    double * data_vals; /* input values, array of N * d                        */
};

void km_dump(struct km_ctx_s * ctx)
{
    size_t i, j;
    for (i = 0; i < ctx->nclusters; ++i) {
        for (j = 0; j < ctx->dimension - 1; ++j) {
            printf("%.14f,", *(ctx->mean_vals + (ctx->dimension * i) + j));
        }
        printf("%.14f\n", *(ctx->mean_vals + (ctx->dimension * i) + j));
    }
}

void km_destroy(struct km_ctx_s * ctx)
{
    if (ctx != NULL) {
        if (ctx->cardinals != NULL) {
            free(ctx->cardinals);
        }
        if (ctx->mean_vals != NULL) {
            free(ctx->mean_vals);
        }
        if (ctx->data_vals != NULL) {
            free(ctx->data_vals);
        }
        free(ctx);
    }
}

struct km_ctx_s * km_create(size_t d, size_t k, size_t n, size_t m)
{
    struct km_ctx_s * ctx = (struct km_ctx_s *) malloc (sizeof(struct km_ctx_s));
    if (ctx == NULL) {
        perror("km_create: allocate context");
        return NULL;
    }
    memset(ctx, 0, sizeof(struct km_ctx_s));

    ctx->dimension = d;
    ctx->nclusters = k;
    ctx->nobserves = n;
    ctx->max_iters = m;

    ctx->cardinals = (size_t *) malloc (ctx->nclusters * sizeof(size_t));
    if (ctx->cardinals == NULL) {
        perror("km_create: allocate cardinalities vector");
        km_destroy(ctx);
        return NULL;
    }
    memset(ctx->cardinals, 0, sizeof(size_t) * ctx->nclusters);

    ctx->mean_vals = (double *) malloc (ctx->nclusters * ctx->dimension * sizeof(double));
    if (ctx->mean_vals == NULL) {
        perror("km_create: allocate centroids matrix");
        km_destroy(ctx);
        return NULL;
    }
    memset(ctx->mean_vals, 0, sizeof(double) * ctx->nclusters * ctx->dimension);

    ctx->data_vals = (double *) malloc (ctx->nobserves * ctx->dimension * sizeof(double));
    if (ctx->data_vals == NULL) {
        perror("km_create: allocate observation matrix");
        km_destroy(ctx);
        return NULL;
    }
    memset(ctx->data_vals, 0, sizeof(double) * ctx->nobserves * ctx->dimension);

    return ctx;
}

double km_distance_squared(struct km_ctx_s * ctx, size_t cluster_id, double * w)
{
    double d = 0, a = 0;
    size_t j;

    for (j = 0; j < ctx->dimension; ++j) {
        a = (*(ctx->mean_vals + (ctx->dimension * cluster_id) + j)) - (*(w + j));
        d += a * a;
    }

    return d;
}

size_t km_cluster(struct km_ctx_s * ctx, double * w)
{
    size_t i, s = 0;
    double last, curr = 0;
    last = km_distance_squared(ctx, 0, w);
    for (i = 1; i < ctx->nclusters; ++i) {
        curr = km_distance_squared(ctx, i, w);
        if (curr < last) {
            s = i;
            last = curr;
        }
    }
    return s;
}

int prj_scan_input(size_t n, size_t d, PyObject * py_data, double * data, PyObject * py_labels, size_t * labels)
{
  size_t i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < d; ++j) {
      data[(d * i) + j] = PyFloat_AsDouble(PyList_GetItem(py_data, (i * d) + j));
    }
    labels[i] = PyLong_AsSize_t(PyList_GetItem(labels, i));
  }
}

int km_iterate(struct km_ctx_s * ctx)
{
    size_t i, j, s = 0;
    double * new_means = NULL, * temp = NULL;
    size_t * new_cards = NULL;

    new_means = (double *) malloc (ctx->nclusters * ctx->dimension * sizeof(double));
    if (new_means == NULL) {
        perror("km_iterate allocate new centroids matrix");
        return -1;
    }
    memset(new_means, 0, sizeof(double) * ctx->nclusters * ctx->dimension);

    new_cards = (size_t *) malloc (ctx->nclusters * sizeof(size_t));
    if (new_cards == NULL) {
        perror("km_iterate allocate new cardinalities vector");
        return -1;
    }
    memset(new_cards, 0, sizeof(size_t) * ctx->nclusters);

    /* iterate over all observations */
    for (i = 0; i < ctx->nobserves; ++i) {
        /* calculate cluster index for the current observation */
        s = km_cluster(ctx, ctx->data_vals + (i * ctx->dimension));

        /* iterate over the observation's entries */
        for (j = 0; j < ctx->dimension; ++j) {
            /* new_means[s][j] = (new_means[s][j] * new_cards[s] + observation[j]) / (new_cards[s] + 1) */
            *(new_means + (s * ctx->dimension) + j) *= new_cards[s];
            *(new_means + (s * ctx->dimension) + j) += *(ctx->data_vals + (i * ctx->dimension) + j);
            *(new_means + (s * ctx->dimension) + j) /= new_cards[s] + 1;
        }
        /* new_cards[s] += 1 */
        new_cards[s]++;
    }

    /* compare new_means to the current ones */
    for (i = 0; i < ctx->nclusters; ++i) {
        for (j = 0; j < ctx->dimension; ++j) {
            if ( *(new_means + (i * ctx->dimension) + j) != *(ctx->mean_vals + (i * ctx->dimension) + j) ) {
                temp = ctx->mean_vals;
                ctx->mean_vals = new_means;
                free(temp);
                free(new_cards);
                return 0;
            }
        }
    }
    free(new_means);
    free(new_cards);
    return 1;
}

int km_converge(struct km_ctx_s * ctx)
{
    int r = -1;
    size_t iter;

    for (iter = 0; iter < ctx->max_iters; ++iter) {
        if ((r = km_iterate(ctx)) == 1) {
            return 1;
        } else if (r < 0) {
            perror("km_converge: km_iterate failed");
            return  -1;
        }
    }

    return 0;
}

double distance_squared(double * v1, double * v2, size_t d)
{
  double r = 0, a = 0;
  size_t j;
  
  for (j = 0; j < d; ++j) {
    a = v1[j] - v2[j];
    r += a * a;
  }
  
  return r;
}

double col_norm_squared(size_t col, size_t n, double * m)
{
  double r = 0, a = 0;
  size_t i;
  
  for (i = 0; i < n; ++i) {
    a = v[(i * n) + col];
    r += a * a;
  }
  
  return r;
}

double * allocate_vector(size_t d)
{
  double * result = NULL;

  if ((result = (double *) malloc (sizeof(double) * d)) == NULL)
    {
      perror("out of memory");
      return NULL
    }

  memset(result, 0, sizeof(double) * d);

  return result;
}

double * allocate_matrix(size_t rows, size_t cols)
{
  return allocate_vector(rows * cols);
}

double * weighted_adjacency_matrix(size_t n, size_t d, double * x)
{
  double * w = NULL;
  size_t i, j;

  if ((w = allocate_matrix(n, n)) == NULL)
    {
      perror("failed to allocate weighted adjacency matrix");
      return NULL;
    }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      w[(i * n) + j] = (i != j ? 0 : exp(-(distance_squared(x[i*n], x[j*n]) / 2)));
    }
  }

  return w;
}

double * normalized_graph_laplacian(size_t n, double * x)
{
  double * l = NULL, *w = NULL, * d_diag = NULL;
  size_t i, j;

  if ((w = weighted_adjacency_matrix(n, d, x)) == NULL)
    {
      perror("failed to create weighted adjacency matrix");
      goto cleanup;
    }

  if ((l = allocate_matrix(n, n)) == NULL)
    {
      perror("failed to allocate normalized graph laplacian");
      goto cleanup;

    }

  if ((d_diag = allocate_vector(n)) == NULL)
    {
      perror("failed to allocate d_diag");
      goto cleanup;
    }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      d_diag[i] += w[(i * n) + j];
    }
    d_diag[i] = 1 / sqrt(d_diag[i]);
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      l[(i * n) + j] = (i == j) - (w[(i * n) + j] * d_diag[j] * d_diag[i]);
    }
  }
  
  return l;
  
cleanup:
  if (w != NULL) {
    free(w);
  }
  
  if (l != NULL) {
    free(w);
  }
  
  if (d != NULL) {
    free(w);
  }

  return NULL;
}

double * identity_matrix(size_t n)
{
  double * r = NULL;
  if ((r = allocate_matrix(n, n)) == NULL)
    {
      return NULL;
    }

  for (i = 0; i < n; i++) {
    r[(i * n) + i] = 1;
  }

  return r;
}

int modified_gram_schmidt(size_t n, double * a, double * q, double * r)
{
  double rij, * u = a;
  size_t i, j, k;
  for (i = 0; i < n; i++) {
    r[(i * n) + i] = col_norm_squared(i, n, a);
    for (j = 0; j < n; j++) {
      q[(j * n) + i] = u[(j * n) + i] / r[(i * n) + i];
    }
    for (j = 0; j < n; j++) {
      rij = 0;
      for (k = 0; k < n; k++) {
	rij += q[(i * n) + k] * u[(k * n) + j];
      }
      r[(i * n) + j] = rij;
      for (k = 0; k < n; k++) {
	u[(k * n) + j] = u[(k * n) + j] - (rij * q[(k * n) + i]);
      }
    }
  }
}

double * multiply_matrix(size_t n, double * m1, double * m2)
{
  double * r = NULL, rij;
  size_t i, j, k;
  if ((r = allocate_matrix(n, n)) == NULL)
    {
      return NULL;
    }

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      rij = 0;
      for (k = 0; k < n; ++k) {
	rij += m1[(i * n) + k] * m2[(k * n) + j];
      }
      r[(i * n) + j] = rij;
    }
  }

  return r;
}

int qr_iteration(size_t n, double * a, double * q_out, double * a_out)
{
  double * q = NULL, * r = NULL;
  size_t i;
  
  if ((q = allocate_matrix(n, n)) == NULL)
    {
      return -1;
    }
  
  if ((r = allocate_matrix(n, n)) == NULL)
    {
      free(q);
      return -1;
    }

  for (i = 0; i < n; i++) {
    modified_gram_schmidt(n, a, q, r);
    multiply_matrix(n, q, r, a_out);
    multiply_matrix(n, q_out, q, q_out);
    /* TODO - check for early convergance and bail */
  }
  
  return 0;
}

double * normalize_cols(size_t n, size_t d, double * m)
{
  double mik = 0, s = 0;
  double * n = NULL;

  if ((n = allocate_matrix(n, d)) == NULL) {
    perror("failed to allocate normalized matrix");
    return -1;
  }
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < d; ++j) {
      s = 0;
      for (k = 0; k < d; ++k) {
	mik = m[(i * n) + k];
	s = mik * mik;
      }
      n[(i * d) + j] = m[(i * n) + j] / sqrt(s);
    }
  }

  return n;
}

int normalized_spectral_clustering(size_t n, size_t d, size_t * k_inout, double * x, size_t * y_out)
{
  double * l = NULL, * a = NULL, * q = NULL, t = NULL;
  
  l = normalized_graph_laplacian(n, w);
  if (l == NULL) {
    perror("failed to create normalized graph laplacian");
    goto cleanup;
  }

  if ((a = allocate_matrix(n, n)) = NULL) {
    perror("failed to allocate A matrix for QR iteration");
    goto cleanup;
  }

  if ((q = allocate_matrix(n, n)) = NULL) {
    perror("failed to allocate Q matrix for QR iteration");
    goto cleanup;
  }

  if (qr_iteration(n, l, q, a) != 0) {
    perror("QR iteration failed");
    goto cleanup;
  }

  if ((k = eigengap_heuristic(n, a)) < 0) {
    perror("eigengap heuristic failed");
    goto cleanup;
  }

  if ((t = normalize_cols(n, k, q)) == NULL) {
    perror("renormalizing failed");
    goto cleanup;
  }

  mat_dump(n, k, t);

  //  ... now do kmeans
}

static PyObject * prj_main(PyObject *self, PyObject *args)
{
    (void)self;

    size_t d = 0, k = 0, n = 0, m = 0;
    struct km_ctx_s * ctx = NULL;
    PyObject * py_data = NULL;
    PyObject * py_labels = NULL;
    double * data = NULL;
    double * labels = NULL;

    if(!PyArg_ParseTuple(args, "IIIIOO:km_fit", &k, &n, &d, &m, &py_data, &py_labels)) {
        return NULL; 
    }

    if ((data = allocate_matrix(n, d)) == NULL)
      {
	return NULL;
      }

    if ((labels = allocate_vector(n)) == NULL)
      {
	return NULL;
      }

    if (prj_scan_input(n, d, py_data, data, py_labels, labels) < 0)
      {
	perror("prj_main: failed to scan input data");
	return NULL;
      }

    output_data(n, d, data, labels);

    if (nsc_fit(n, d, &k, data, labels) < 0)
      {
	perror("prj_main: NSC fitting failed");
	return NULL;
      }

    output_k(k);
    output_labels(labels);

    if (km_fit(n, d, k, data, labels) < 0)
      {
	perror("prj_main: KM fitting failed");
	return NULL;
      }

    output_labels(labels);

    Py_RETURN_NONE;
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
