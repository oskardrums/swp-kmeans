#include "mat.h"
#include "log.h"
#include "kmpp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/*
 * O(n^2 * d)
 */
double * weighted_adjacency_matrix(size_t n, size_t d, const double * x)
{
  double * w = NULL, wij = 0;
  size_t i, j;

  if ((w = mat_allocate(n, n)) == NULL)
    {
      perror("failed to allocate weighted adjacency matrix");
      return NULL;
    }

  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      wij = exp(-sqrt(vec_distance_squared(d, &(x[i*d]), &(x[j*d])) / 2));
      w[(i * n) + j] = wij;
      w[(j * n) + i] = wij;
    }
  }

  return w;
}

/*
 * O(n^2 * d)
 */
double * normalized_graph_laplacian(size_t n, size_t d, const double * x)
{
  double * l = NULL, *w = NULL, * d_diag = NULL;
  size_t i, j;
  bool err = false;

  if ((w = weighted_adjacency_matrix(n, d, x)) == NULL)
    {
      perror("failed to create weighted adjacency matrix");
      err = true;
      goto cleanup;
    }

  if ((l = mat_allocate(n, n)) == NULL)
    {
      perror("failed to allocate normalized graph laplacian");
      err = true;
      goto cleanup;

    }

  if ((d_diag = vec_allocate(n)) == NULL)
    {
      perror("failed to allocate d_diag");
      err = true;
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

cleanup:

  if (w != NULL) {
    free(w);
  }

  if (d_diag != NULL) {
    free(d_diag);
  }

  if (err) {
    if (l != NULL) {
      free(l);
      l = NULL;
    }
  }

  return l;
}

/*
 * O(n^3)
 */
void modified_gram_schmidt(size_t n, double * u, double * q, double * r)
{
  double rij;
  size_t i, j, k;
  for (i = 0; i < n; i++) {
    r[(i * n) + i] = sqrt(mat_col_norm_squared(i, n, n, u));
    for (j = 0; j < n; j++) {
      q[(j * n) + i] = u[(j * n) + i] / r[(i * n) + i];
    }
    for (j = i + 1; j < n; j++) {
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

/*
 * O(n^4)
 */
double * qr_iteration(size_t n, double * a_out)
{
  bool err = false, converged = true;
  double * q = NULL, * q_out = NULL, * q_out_times_q = NULL;
  double * r = NULL;
  double * temp = NULL;
  size_t i;

  if ((q_out = mat_identity(n)) == NULL) {
    err = true;
    goto cleanup;
  }

  if ((q = mat_allocate(n, n)) == NULL) {
    err = true;
    goto cleanup;
  }

  if ((q_out_times_q = mat_allocate(n, n)) == NULL) {
    err = true;
    goto cleanup;
  }

  if ((r = mat_allocate(n, n)) == NULL) {
    err = true;
    goto cleanup;
  }

  for (i = 0; i < n; i++) {
    log_emit("running modified_gram_schmidt");
    modified_gram_schmidt(n, a_out, q, r);
    mat_multiply(n, n, r, n, n, q, a_out);

    mat_multiply(n, n, q_out, n, n, q, q_out_times_q);

    if (mat_abs_equals(n, n, q_out, q_out_times_q)) {
      log_emit("converged");
      converged = true;
    }

    temp = q_out;
    q_out = q_out_times_q;
    q_out_times_q = temp;

    if (converged) {
      break;
    }

  }

 cleanup:
  free(r);
  free(q);
  free(q_out_times_q);

  if (err) {
    free(q_out);
    q_out = NULL;
  }

  return q_out;
}

/*
 * O(n)
 */
size_t eigengap_heuristic(size_t n, size_t * sorted_eigenvalues, const double * a)
{
  double gap = 0, max_gap = 0;
  size_t max_ind = 0, i, j = 0, j_prev = 0;

  for (i = 1; i < ((n & 1) ? ((n + 1) / 2) : (n / 2)); ++i) {
    j = sorted_eigenvalues[i];
    j_prev = sorted_eigenvalues[i - 1];

    gap = a[(j * n) + j] - a[((j_prev) * n) + (j_prev)];

    if (gap > max_gap) {
      max_gap = gap;
      max_ind = i;
    }
  }

  return max_ind;
}

struct nsc_compare_diag_arg {
  size_t n;
  double * mat;
};

struct nsc_compare_diag_arg comp_arg;

int nsc_compare_diag(const void * i_ptr, const void * j_ptr)
{
  size_t i = *(size_t*)i_ptr;
  size_t j = *(size_t*)j_ptr;
  size_t n = comp_arg.n;
  double mat_ii = comp_arg.mat[(i * n) + i];
  double mat_jj = comp_arg.mat[(j * n) + j];
  if (mat_ii > mat_jj) {
    return 1;
  } else if (mat_jj > mat_ii) {
    return -1;
  } else {
    return 0;
  }
}

size_t * nsc_sort_eigenvalues(size_t n, double * a)
{
  size_t i;
  size_t * eigenvalues = NULL;

  eigenvalues = (size_t *)malloc(sizeof(size_t) * n);
  if (eigenvalues == NULL) {
    return NULL;
  }

  for (i = 0; i < n; ++i) {
    eigenvalues[i] = i;
  }

  comp_arg.n = n;
  comp_arg.mat = a;
  qsort(eigenvalues, n, sizeof(size_t), nsc_compare_diag);

  return eigenvalues;
}

/*
 * O(n^4)
 */
size_t normalized_spectral_clustering(size_t k, size_t n, size_t d, const double * x, size_t m, size_t ** y_out)
{
  bool err = false;
  double * l = NULL, * q = NULL, * t = NULL;
  size_t * sorted_eigenvalues = NULL;

  log_emit("running normalized_graph_laplacian");
  log_emit("x");
  mat_dump(n, d, x);

  l = normalized_graph_laplacian(n, d, x);
  if (l == NULL) {
    err = true;
    perror("failed to create normalized graph laplacian");
    goto cleanup;
  }

  log_emit("l");
  mat_dump(n, n, l);

  log_emit("running qr_iteration");
  if ((q = qr_iteration(n, l)) == NULL) {
    err = true;
    perror("QR iteration failed");
    goto cleanup;
  }
  log_emit("q");
  mat_dump(n, n, q);
  log_emit("a");
  mat_dump(n, n, l);

  sorted_eigenvalues = nsc_sort_eigenvalues(n, l);

  if (k == 0) {
    log_emit("running eigengap_heuristic");
    if ((k = eigengap_heuristic(n, sorted_eigenvalues, l)) == 0) {
      err = true;
      perror("eigengap heuristic failed");
      goto cleanup;
    }
  }

  log_emit("running mat_copy_cols");
  if ((t = mat_copy_cols(n, n, q, k, sorted_eigenvalues)) == NULL) {
    err = true;
    goto cleanup;
  }

  log_emit("t");
  mat_dump(n, k, t);

  log_emit("running mat_normalize_rows");
  mat_normalize_rows(n, k, t);

  log_emit("normalized t");
  mat_dump(n, k, t);

  log_emit("running kmpp");
  *y_out = kmpp(k, n, k, t, m);
  if (*y_out == NULL) {
    err = true;
    goto cleanup;
  }

cleanup:
    log_emit("running kmpp");
  free(t);
    log_emit("running kmpp");
  free(q);
    log_emit("running kmpp");
  free(l);
    log_emit("running kmpp");
  free(sorted_eigenvalues);

  if (err) {
    k = 0;
    *y_out = NULL;
  }

  return k;
}
