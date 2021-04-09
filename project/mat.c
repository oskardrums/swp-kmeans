#include "mat.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void mat_dump(size_t num_rows, size_t num_cols, double * mat)
{
  size_t row, col;
  double mat_entry;
  
  for (row = 0; row < num_rows; ++row) {
    for (col = 0; col < num_cols - 1; ++col) {
      mat_entry = mat[(row * num_cols) + col];
      printf("%s%.3f, ", (mat_entry > 0) ? "+" : "", mat_entry);
    }
    printf("%.3f\n", mat[(row * num_cols) + col]);
  }
}

double vec_distance_squared(size_t d, double * v1, double * v2)
{
  double r = 0, a = 0;
  size_t j;
  
  for (j = 0; j < d; ++j) {
    a = v1[j] - v2[j];
    r += a * a;
  }
  
  return r;
}

double mat_col_norm_squared(size_t col, size_t num_rows, size_t num_cols, double * mat)
{
  double result = 0, temp = 0;
  size_t row;
  
  for (row = 0; row < num_rows; ++row) {
    temp = mat[(row * num_cols) + col];
    result += temp * temp;
  }
  
  return result;
}

double * vec_allocate(size_t d)
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

double * mat_allocate(size_t rows, size_t cols)
{
  return vec_allocate(rows * cols);
}

double * mat_identity(size_t n)
{
  double * r = NULL;
  size_t i;
  
  if ((r = mat_allocate(n, n)) == NULL) {
    perror("mat_identity: allocation failed");
    return NULL;
  }

  for (i = 0; i < n; i++) {
    r[(i * n) + i] = 1;
  }

  return r;
}

/*
 * num_cols1 == num_rows2
 */
void mat_multiply(size_t num_rows1, size_t num_cols1, double * mat1,
		  size_t num_rows2, size_t num_cols2, double * mat2,
		  double * mat_out)
{
  double mat_entry;
  size_t row, col, k;

  (void)num_rows2;
  
  for (row = 0; row < num_rows1; ++row) {
    for (col = 0; col < num_cols2; ++col) {
      mat_entry = 0;
      for (k = 0; k < num_cols1; ++k) {
	mat_entry += mat1[(row * num_cols1) + k] * mat2[(k * num_cols2) + col];
      }
      mat_out[(row * num_cols2) + col] = mat_entry;
    }
  }
}


/*
 * num_cols >= num_out_cols
 */
void mat_trim_and_normalize_cols(size_t num_rows, size_t num_cols, double * mat,
				 size_t num_out_cols, double * mat_out)
{
  double mat_entry = 0, row_sum = 0;
  size_t row, col;
  
  for (row = 0; row < num_rows; ++row) {
    row_sum = 0;
    
    for (col = 0; col < num_out_cols; ++col) {
	mat_entry = mat[(row * num_cols) + col];
	row_sum += mat_entry * mat_entry;
    }
    
    for (col = 0; col < num_out_cols; ++col) {
      mat_out[(row * num_out_cols) + col] = mat[(row * num_cols) + col] / sqrt(row_sum);
    }
  }
}

