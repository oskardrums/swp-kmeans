/*
 * Normalized Spectral Clustering algorithm
 */
#include "nsc.h"
#include "mat.h"
#include "kmpp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>


static void nsc_qr_print_eta(size_t cur_iteration,
			     size_t max_iterations,
			     time_t start_time,
			     time_t cur_time)
{
  double eta = 0, time_delta = 0;
  if (cur_iteration > 0) {
    time_delta += cur_time - start_time;
    eta = time_delta * ((((double)max_iterations / cur_iteration)) - 1);
    printf("QR iteration %lu/%lu\tTime elapsed %.2fs\tETA %.2fs\tExpected total %.2fs\n",
	   cur_iteration,
	   max_iterations,
	   time_delta,
	   eta,
	   time_delta + eta);
  }
}

/*
 * Creates the Weighted Adjacency Matrix W for
 * the set of n d-dimensional samples in mat
 *
 * O(num_rows^2 * num_cols)
 */
static double * nsc_weighted_adjacency_matrix(size_t num_rows,
					      size_t num_cols,
					      const double * mat)
{
  double * w = NULL, wij = 0;
  size_t i, j;

  if ((w = mat_allocate(num_rows, num_rows)) == NULL)
    {
      perror("failed to allocate weighted adjacency matrix");
      return NULL;
    }

  for (i = 0; i < num_rows; i++) {
    for (j = i + 1; j < num_rows; j++) {
      wij = exp(-sqrt(vec_distance_squared(num_cols, &(mat[i*num_cols]), &(mat[j*num_cols])) / 2));
      w[(i * num_rows) + j] = wij;
      w[(j * num_rows) + i] = wij;
    }
  }

  return w;
}

/*
 * Creates the Normalized Graph Laplacian L for
 * the set of n d-dimensional samples in x
 *
 * O(num_rows^2 * num_cols)
 */
static double * nsc_normalized_graph_laplacian(size_t num_rows,
					       size_t num_cols,
					       const double * mat)
{
  double * l = NULL, *w = NULL, * d_diag = NULL;
  size_t i, j;
  bool err = false;

  if ((w = nsc_weighted_adjacency_matrix(num_rows, num_cols, mat)) == NULL)
    {
      err = true;
      goto cleanup;
    }

  if ((l = mat_allocate(num_rows, num_rows)) == NULL)
    {
      err = true;
      goto cleanup;

    }

  if ((d_diag = vec_allocate(num_rows)) == NULL)
    {
      err = true;
      goto cleanup;
    }

  for (i = 0; i < num_rows; i++) {
    for (j = 0; j < num_rows; j++) {
      d_diag[i] += w[(i * num_rows) + j];
    }
    d_diag[i] = 1 / sqrt(d_diag[i]);
  }

  for (i = 0; i < num_rows; i++) {
    for (j = 0; j < num_rows; j++) {
      l[(i * num_rows) + j] = (i == j) - (w[(i * num_rows) + j] * d_diag[j] * d_diag[i]);
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
 * Given an n x n matrix u, computes the decomposition of into q and r
 * where q is orthogonal and r is upper triangular.
 * The input matrix u may be altered in hte process in an undefined manner.
 *
 * O(n^3)
 */
static void nsc_modified_gram_schmidt(size_t n, double * u, double * q, double * r)
{
  double rij, l2norm = 0;
  size_t col, row, col2;
  for (col = 0; col < n; ++col) {
    l2norm = mat_col_norm_squared(col, n, n, u);

    if (l2norm > 0) {
      l2norm = sqrt(l2norm);

      for (row = 0; row < n; ++row) {
	q[(row * n) + col] = u[(row * n) + col] / l2norm;
      }
    }

    r[(col * n) + col] = l2norm;

    for (col2 = col + 1; col2 < n; ++col2) {
      rij = 0;
      for (row = 0; row < n; ++row) {
	rij += q[(row * n) + col] * u[(row * n) + col2];
      }

      r[(col * n) + col2] = rij;

      for (row = 0; row < n; ++row) {
	u[(row * n) + col2] -= (rij * q[(row * n) + col]);
      }
    }
  }
}

/*
 * Computes and returns an orthogonal matrix q composed of eigenvectors of
 * the input matrix a.
 * a itself is altered so that the diagonal of the ouput a is made up of eigenvalues
 * of the input a, in matching order to the eigenvectors in q.
 *
 * O(n^4)
 */
static double * nsc_qr_iteration(size_t n, double * a)
{
  bool err = false, converged = false;
  double * q = NULL, * q_out = NULL, * q_out_times_q = NULL;
  double * r = NULL;
  double * temp = NULL;
  size_t i = 0;
  time_t base_ts = 0, ts = 0, last_ts = 0;

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

  base_ts = time(NULL);
  last_ts = base_ts;

  for (i = 0; i < n; i++) {
    ts = time(NULL);

    if (ts - last_ts > 1) {
      nsc_qr_print_eta(i, n, base_ts, ts);
      last_ts = ts;
    }

    nsc_modified_gram_schmidt(n, a, q, r);

    mat_upper_triangular_multiply(n, n, r, n, n, q, a);

    mat_multiply(n, n, q_out, n, n, q, q_out_times_q);

    if (mat_abs_equals(n, n, q_out, q_out_times_q)) {
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
  if (r != NULL) {
    free(r);
  }

  if (q != NULL) {
    free(q);
  }

  if (q_out_times_q != NULL) {
    free(q_out_times_q);
  }

  if (err) {
    free(q_out);
    q_out = NULL;
  }

  return q_out;
}

/*
 * Returns the result of the Eigengap Heuristic given a ordered set of eigen
 * values stored in the main diagonal of eigenvalues_mat, where the ordering
 * of the eigen values is given by the index array eigenvalues_order.
 *
 * O(num_eigenvalues)
 */
static size_t nsc_eigengap_heuristic(size_t num_eigenvalues,
				     size_t * eigenvalues_order,
				     const double * eigenvalues_mat)
{
  double gap = 0, max_gap = 0;
  size_t max_ind = 0, i = 0, j = 0, j_prev = 0;
  size_t num_eigenvalues_half =
    (num_eigenvalues & 1) ?
    ((num_eigenvalues + 1) / 2) :
    (num_eigenvalues / 2);

  for (i = 1; i < num_eigenvalues_half; ++i) {
    j = eigenvalues_order[i];
    j_prev = eigenvalues_order[i - 1];

    gap =
      eigenvalues_mat[(j * num_eigenvalues) + j] -
      eigenvalues_mat[((j_prev) * num_eigenvalues) + (j_prev)];

    if (gap > max_gap) {
      max_gap = gap;
      max_ind = i;
    }
  }

  return max_ind;
}

/* Layout for the static variable nsc_compare_diag_arg used by nsc_compare_diag */
struct nsc_compare_diag_arg {
  size_t n;
  double * mat;
};

static struct nsc_compare_diag_arg comp_arg;

/*
 * Comparing routine for sorting an array of eigen value indices based on
 * the corresponding eigen values.
 */
static int nsc_compare_diag(const void * i_ptr, const void * j_ptr)
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

/*
 * Returns an array of indices [i_0, i_2, ..., i_n-1] ordered by the
 * corresponding entries in eigenvalues_mat main diagonal.
 *
 * O(num_eigenvalues * log(num_eigenvalues))
 */
static size_t * nsc_sort_eigenvalues(size_t num_eigenvalues, double * eigenvalues_mat)
{
  size_t i;
  size_t * eigenvalues = NULL;

  eigenvalues = (size_t *)malloc(sizeof(size_t) * num_eigenvalues);
  if (eigenvalues == NULL) {
    return NULL;
  }

  for (i = 0; i < num_eigenvalues; ++i) {
    eigenvalues[i] = i;
  }

  comp_arg.n = num_eigenvalues;
  comp_arg.mat = eigenvalues_mat;

  qsort(eigenvalues, num_eigenvalues, sizeof(size_t), nsc_compare_diag);

  comp_arg.n = 0;
  comp_arg.mat = NULL;

  return eigenvalues;
}

/*
 * Clusters the rows of mat into num_clusters clusters
 * and stores the resulting labels for each row in an array pointed to by
 * labels. If num_clusters is 0, computes a viable num_clusters parameter.
 * Returns the number of clusters or 0 in case of error.
 *
 * O(num_rows^4) when num_rows >> num_cols
 */
size_t normalized_spectral_clustering(size_t num_clusters,
				      size_t num_rows,
				      size_t num_cols,
				      const double * mat,
				      size_t max_iters,
				      size_t ** labels)
{
  bool err = false;
  double * l = NULL, * q = NULL, * t = NULL;
  size_t * sorted_eigenvalues = NULL;

  l = nsc_normalized_graph_laplacian(num_rows, num_cols, mat);
  if (l == NULL) {
    err = true; goto cleanup;
  }

  if ((q = nsc_qr_iteration(num_rows, l)) == NULL) {
    err = true; goto cleanup;
  }

  sorted_eigenvalues = nsc_sort_eigenvalues(num_rows, l);

  if (num_clusters == 0) {
    if ((num_clusters = nsc_eigengap_heuristic(num_rows, sorted_eigenvalues, l)) == 0) {
      err = true; goto cleanup;
    }
  }

  if ((t = mat_copy_cols(num_rows, num_rows, q, num_clusters, sorted_eigenvalues)) == NULL) {
    err = true;
    goto cleanup;
  }

  mat_normalize_rows(num_rows, num_clusters, t);

  num_clusters = k_means_pp(num_clusters, num_rows, num_clusters, t, max_iters, labels);
  if (num_clusters == 0) {
    err = true; goto cleanup;
  }

cleanup:
  if (t != NULL) free(t);

  if (q != NULL) free(q);

  if (l != NULL) free(l);

  if (sorted_eigenvalues != NULL) {
    free(sorted_eigenvalues);
  }

  if (err) {
    if (*labels != NULL) {
      free(*labels);
      *labels = NULL;
    }
    return 0;
  }

  return num_clusters;
}
