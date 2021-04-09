#include "mat.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

double * weighted_adjacency_matrix(size_t n, size_t d, double * x)
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
      wij = exp(vec_distance_squared(d, &(x[i*d]), &(x[j*d])) / -2);
      w[(i * n) + j] = wij;
      w[(j * n) + i] = wij;
    }
  }

  return w;
}

double * normalized_graph_laplacian(size_t n, size_t d, double * x)
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

int modified_gram_schmidt(size_t n, double * a, double * q, double * r)
{
  double rij, * u = a;
  size_t i, j, k;
  for (i = 0; i < n; i++) {
    r[(i * n) + i] = mat_col_norm_squared(i, n, n, a);
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
  return 0;
}

int qr_iteration(size_t n, double * a, double * q_out, double * a_out)
{
  double * q = NULL, * r = NULL;
  size_t i;

  if ((q = mat_allocate(n, n)) == NULL) {
    return -1;
  }

  if ((r = mat_allocate(n, n)) == NULL) {
    free(q);
    return -1;
  }

  for (i = 0; i < n; i++) {
    modified_gram_schmidt(n, a, q, r);
    mat_multiply(n, n, r, n, n, q, a_out);
    mat_multiply(n, n, q_out, n, n, q, q_out);
    /* TODO - check for early convergance and bail */
  }

  free(r);
  free(q);

  return 0;
}

size_t eigengap_heuristic(size_t n, double * a)
{
  double gap = 0, max_gap = 0;
  size_t max_ind = 0, i;

  for (i = 1; i < n; ++i) {
    gap = a[((i-1) * n) + (i-1)] - a[(i * n) + i];

    if (gap > max_gap) {
      max_gap = gap;
      max_ind = i;
    }
  }

  return max_ind;
}

int normalized_spectral_clustering(size_t n, size_t d, size_t * k_inout, double * x, size_t * y_out)
{
  double * l = NULL, * a = NULL, * q = NULL, * t = NULL;
  size_t k = 0;
  (void)k_inout;
  (void)y_out;

  l = normalized_graph_laplacian(n, d, x);
  if (l == NULL) {
    perror("failed to create normalized graph laplacian");
    goto cleanup;
  }

  if ((a = mat_allocate(n, n)) == NULL) {
    perror("failed to allocate A matrix for QR iteration");
    goto cleanup;
  }

  if ((q = mat_identity(n)) == NULL) {
    perror("failed to allocate Q matrix for QR iteration");
    goto cleanup;
  }

  if (qr_iteration(n, l, q, a) != 0) {
    perror("QR iteration failed");
    goto cleanup;
  }

  free(l);

  if ((k = eigengap_heuristic(n, a)) == 0) {
    perror("eigengap heuristic failed");
    goto cleanup;
  }

  free(a);

  if ((t = mat_allocate(n, k)) == NULL) {
    perror("failed to allocate T matrix");
    goto cleanup;
  }

  mat_trim_and_normalize_cols(n, n, q, k, t);

  free(q);

  mat_dump(n, k, t);

  free(t);

cleanup:
  return 0;
  //  ... now do kmeans
}

int main()
{
  double data[] =
    {
     1, 2,
     3, 4,
     5, 6,
     7, 8,
     9, 0,
     1, 3,
     2, 4,
     3, 5,
     4, 6,
     5, 7,
    };

  return normalized_spectral_clustering(10, 2, NULL, data, NULL);
}
