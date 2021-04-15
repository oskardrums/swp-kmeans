#include "mat.h"
#include "log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void mat_dump(size_t num_rows, size_t num_cols, const double * mat)
{
  size_t row, col;
  double mat_entry;

  for (row = 0; row < num_rows; ++row) {
    for (col = 0; col < num_cols - 1; ++col) {
      mat_entry = mat[(row * num_cols) + col];
      printf("%s%.4f, ", (mat_entry > 0) ? "+" : "", mat_entry);
    }
    mat_entry = mat[(row * num_cols) + col];
    printf("%s%.4f\n", (mat_entry > 0) ? "+" : "", mat_entry);
  }
}

double vec_distance_squared(size_t d, const double * v1, const double * v2)
{
  double r = 0, a = 0;
  size_t j;

  for (j = 0; j < d; ++j) {
    a = v1[j] - v2[j];
    r += a * a;
  }

  return r;
}

/*
 * O(num_rows)
 */
double mat_col_norm_squared(size_t col, size_t num_rows, size_t num_cols, const double * mat)
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

  printf("allocating %lu doubles\n", d);

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
 * O(num_rows1 * num_cols2 * num_cols1)
 */
void mat_multiply(size_t num_rows1, size_t num_cols1, const double * mat1,
		  size_t num_rows2, size_t num_cols2, const double * mat2,
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
 * O(num_rows * num_out_cols)
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

void mat_normalize_rows(size_t num_rows, size_t num_cols, double * mat)
{
  double mat_entry = 0, row_squares_sum = 0;
  size_t row, col;

  for (row = 0; row < num_rows; ++row) {
    row_squares_sum = 0;

    for (col = 0; col < num_cols; ++col) {
	mat_entry = mat[(row * num_cols) + col];
	row_squares_sum += mat_entry * mat_entry;
    }

    if (row_squares_sum > 0) {
      for (col = 0; col < num_cols; ++col) {
	mat_entry = mat[(row * num_cols) + col];
	mat[(row * num_cols) + col] = mat_entry / sqrt(row_squares_sum);
      }
    }
  }
}

/*
 * num_cols >= num_out_cols
 * O(num_rows * num_out_cols)
 */
double * mat_copy_cols(size_t num_rows, size_t num_cols, const double * mat,
		       size_t num_out_cols, size_t * indices)
{
  size_t row, col, i;
  double * mat_out = NULL;

  if ((mat_out = mat_allocate(num_rows, num_out_cols)) == NULL) {
    return NULL;
  }

  for (i = 0; i < num_out_cols; ++i) {
    col = indices[i];
    printf("copying col=%lu\n", col);
    for (row = 0; row < num_rows; ++row) {
      printf("row=%lu\n", row);
      mat_out[(row * num_out_cols) + i] = mat[(row * num_cols) + col];
    }
  }

  return mat_out;
}

bool mat_abs_equals(size_t num_rows, size_t num_cols, const double * mat1, const double * mat2)
{
  size_t i, j;
  double delta = 0, mat1ij = 0, mat2ij = 0;

  for (i = 0; i < num_rows; ++i) {
    for (j = 0; j < num_cols; ++j) {
      mat1ij = mat1[(i * num_cols) + j];
      mat2ij = mat2[(i * num_cols) + j];
      delta = fabs(fabs(mat1ij) - fabs(mat2ij));

#define EPSILON 0.0001
      if (delta > EPSILON) {
	return false;
      }
    }
  }

  return true;

}

bool mat_equals(size_t num_rows, size_t num_cols, const double * mat1, const double * mat2)
{
  size_t i, j;
  double delta = 0, mat1ij = 0, mat2ij = 0;

  for (i = 0; i < num_rows; ++i) {
    for (j = 0; j < num_cols; ++j) {
      mat1ij = mat1[(i * num_cols) + j];
      mat2ij = mat2[(i * num_cols) + j];
      delta = fabs(mat1ij - mat2ij);

      if (delta > EPSILON) {
	return false;
      }
    }
  }

  return true;
}
