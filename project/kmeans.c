#include "kmpp.h"
#include "mat.h"
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

static size_t kmpp_cluster(double * centroids, size_t num_clusters, size_t num_cols, const double * w)
{
    size_t i, s = 0;
    double last, curr = 0;
    last = vec_distance_squared(num_cols, &(centroids[(0 * num_cols)]), w);

    for (i = 1; i < num_clusters; ++i) {
      curr = vec_distance_squared(num_cols, &(centroids[(i * num_cols)]), w);
      if (curr < last) {
	s = i;
	last = curr;
      }
    }
    return s;
}

static size_t * kmpp_converge(double * centroids, size_t num_clusters, size_t num_rows, size_t num_cols, const double * mat, size_t max_iters)
{
  bool err = false;
  size_t iter, i, j;
  size_t * cardinals = NULL, * clusters = NULL;
  double * new_centroids = NULL, * centroids_a = NULL, * centroids_b = NULL;

  cardinals = (size_t *)malloc(sizeof(size_t) * num_clusters);
  if (cardinals == NULL) {
    err = true;
    goto cleanup;
  }

  clusters = (size_t *)malloc(sizeof(size_t) * num_rows);
  if (clusters == NULL) {
    err = true;
    goto cleanup;
  }
  memset(clusters, 0, sizeof(size_t) * num_rows);

  centroids_a = centroids;

  centroids_b = (double *)malloc(sizeof(double) * num_clusters * num_cols);
  if (centroids_b == NULL) {
    err = true;
    goto cleanup;
  }
  memset(centroids_b, 0, sizeof(double) * num_clusters * num_cols);

  for (iter = 0; iter < max_iters; ++iter) {
    if (iter & 1) {
      centroids = centroids_b;
      new_centroids = centroids_a;
    } else {
      centroids = centroids_a;
      new_centroids = centroids_b;
    }

    memset(cardinals, 0, sizeof(size_t) * num_clusters);
    memset(new_centroids, 0, sizeof(size_t) * num_clusters * num_cols);

    for (i = 0; i < num_rows; ++i) {
      clusters[i] = kmpp_cluster(centroids, num_clusters, num_cols, &(mat[i * num_cols]));
      for (j = 0; j < num_cols; ++j) {
	new_centroids[(clusters[i] * num_cols) + j] *= cardinals[clusters[i]];
	new_centroids[(clusters[i] * num_cols) + j] += mat[(i * num_cols) + j];
	new_centroids[(clusters[i] * num_cols) + j] /= ++(cardinals[clusters[i]]);
      }
    }

    if (mat_equals(num_clusters, num_cols, centroids, new_centroids)) {
      break;
    }
  }

cleanup:
  if (err) {
    if (clusters != NULL) {
      free(clusters);
      clusters = NULL;
    }
  }

  if (cardinals != NULL) {
    free(cardinals);
  }

  if (centroids_b != NULL) {
      free(centroids_b);
  }

  return clusters;
}

static size_t kmpp_random_size_t(size_t n)
{
  return rand() % n;
}

static double kmpp_random_double(double r)
{
  return ((double)rand() * r /(double)RAND_MAX);
}

static double kmpp_min_dist_squared(const double * centroids, size_t num_centroids, const double * vec, size_t num_cols)
{
  double min_dist_squared = 0, cur_dist_squared = 0;
  size_t z = 0;

  min_dist_squared = vec_distance_squared(num_cols, vec, &(centroids[0 * num_cols]));

  for (z = 1; z < num_centroids; ++z) {
    cur_dist_squared = vec_distance_squared(num_cols, vec, &(centroids[z * num_cols]));
    if (cur_dist_squared < min_dist_squared) {
      min_dist_squared = cur_dist_squared;
    }
  }
  return min_dist_squared;
}

static double kmpp_half(size_t n)
{
  return (n & 1) ? ((n - 1) / 2) : (n / 2);
}

static size_t kmpp_binary_search(size_t n, double * cdf, double r)
{
  size_t i = 0, head = 0, tail = n;
  while (head < tail) {
    i = kmpp_half(head + tail);
    if (cdf[i] >= r) {
      if (cdf[i - 1] < r) {
	break;
      } else {
	tail = i;
      }
    } else {
      head = i;
    }
  }
  return i;
}

static size_t kmpp_random_size_t_with_cdf(size_t n, double * cdf)
{
  double random_double = 0;

  random_double = kmpp_random_double(cdf[n - 1]);

  return kmpp_binary_search(n, cdf, random_double);
}

static double * kmpp_initial_centroids(size_t num_clusters, size_t num_rows, size_t num_cols, const double * mat)
{
  double * cdf = NULL, * centroids = NULL;
  size_t j, i, k = 0;

  cdf = vec_allocate(num_rows);
  if (cdf == NULL) {
    return NULL;
  }

  centroids = mat_allocate(num_clusters, num_cols);
  if (centroids == NULL) {
    return NULL;
  }

  memcpy(&centroids[(0 * num_cols)], &(mat[(kmpp_random_size_t(num_rows) * num_cols)]), sizeof(double) * num_cols);

  for (j = 1; j < num_clusters; ++j) {
    cdf[0] = kmpp_min_dist_squared(centroids, j, mat, num_cols);
    for (i = 1; i < num_rows; ++i) {
      cdf[i] = cdf[i - 1] + kmpp_min_dist_squared(centroids, j, &(mat[(i * num_cols)]), num_cols);
    }
    k = kmpp_random_size_t_with_cdf(num_rows, cdf);
    memcpy(&centroids[(j * num_cols)], &(mat[k * num_cols]), sizeof(double) * num_cols);
  }

  free(cdf);

  return centroids;
}

size_t * kmpp(size_t num_clusters, size_t num_rows, size_t num_cols, const double * mat, size_t max_iters)
{
  double * initial_centroids = NULL;
  srand(time(NULL));

  if ((initial_centroids = kmpp_initial_centroids(num_clusters, num_rows, num_cols, mat)) == NULL) {
    return NULL;
  }

  return kmpp_converge(initial_centroids, num_clusters, num_rows, num_cols, mat, max_iters);
}
