#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void mat_dump(size_t n, size_t d, double * m)
{
  size_t i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < d; ++j) {
      printf("%.3f, ", m[(i * d) + j]);
    }
    printf("\n");
  }
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
    a = m[(i * n) + col];
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
      return NULL;
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
  double * w = NULL, wij = 0;
  size_t i, j;

  if ((w = allocate_matrix(n, n)) == NULL)
    {
      perror("failed to allocate weighted adjacency matrix");
      return NULL;
    }

  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      wij = exp(distance_squared(&(x[i*d]), &(x[j*d]), d) / -2);
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

  if ((w = weighted_adjacency_matrix(n, d, x)) == NULL)
    {
      perror("failed to create weighted adjacency matrix");
      goto cleanup;
    }

  mat_dump(n, n, w);

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
    free(l);
  }
  
  if (d_diag != NULL) {
    free(d_diag);
  }

  return NULL;
}

double * identity_matrix(size_t n)
{
  double * r = NULL;
  size_t i;
  
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

int multiply_matrix(size_t n, double * m1, double * m2, double * out)
{
  double rij;
  size_t i, j, k;
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      rij = 0;
      for (k = 0; k < n; ++k) {
	rij += m1[(i * n) + k] * m2[(k * n) + j];
      }
      out[(i * n) + j] = rij;
    }
  }

  return 0;
}

int qr_iteration(size_t n, double * a, double * q_out, double * a_out)
{
  double * q = NULL, * r = NULL;
  size_t i;
  
  if ((q = allocate_matrix(n, n)) == NULL)
    {
      return -1;
    }
  printf("q:\n");
  mat_dump(n, n, q);
  
  if ((r = allocate_matrix(n, n)) == NULL)
    {
      free(q);
      return -1;
    }

  for (i = 0; i < n; i++) {
    modified_gram_schmidt(n, a, q, r);
    printf("q:\n");
    mat_dump(n, n, q);
    printf("r:\n");
    mat_dump(n, n, r);
    multiply_matrix(n, r, q, a_out);
    printf("q:\n");
    mat_dump(n, n, q);
    printf("q_out:\n");
    mat_dump(n, n, q_out);
    multiply_matrix(n, q_out, q, q_out);
    printf("q_out:\n");
    mat_dump(n, n, q_out);
    printf("a_out:\n");
    mat_dump(n, n, a_out);
    /* TODO - check for early convergance and bail */
  }
  
  return 0;
}

int normalize_cols(size_t n, size_t d, double * m, double * out)
{
  double mij = 0, s = 0;
  size_t i, j;
  
  for (i = 0; i < n; ++i) {
    s = 0;
    
    for (j = 0; j < d; ++j) {
	mij = m[(i * n) + j];
	s += mij * mij;
    }
    
    for (j = 0; j < d; ++j) {
      out[(i * d) + j] = m[(i * n) + j] / sqrt(s);			     
    }
    
  }

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

  mat_dump(n, n, l);

  if ((a = allocate_matrix(n, n)) == NULL) {
    perror("failed to allocate A matrix for QR iteration");
    goto cleanup;
  }

  if ((q = identity_matrix(n)) == NULL) {
    perror("failed to allocate Q matrix for QR iteration");
    goto cleanup;
  }

  if (qr_iteration(n, l, q, a) != 0) {
    perror("QR iteration failed");
    goto cleanup;
  }

  printf("q\n");
  mat_dump(n, n, q);


  if ((k = eigengap_heuristic(n, a)) == 0) {
    perror("eigengap heuristic failed");
    goto cleanup;
  }
  printf("k=%lu\n", k);

  if ((t = allocate_matrix(n, k)) == NULL) {
    perror("failed to allocate T matrix");
    goto cleanup;
  }

  if (normalize_cols(n, k, q, t) != 0) {
    perror("renormalizing failed");
    goto cleanup;
  }

  mat_dump(n, k, t);

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
