#ifndef __MAT_H
#define __MAT_H

#include <stddef.h>
#include <stdbool.h>

void mat_dump(size_t num_rows, size_t num_cols, const double * mat);

double mat_col_norm_squared(size_t col, size_t num_rows, size_t num_cols, const double * mat);

double vec_distance_squared(size_t, const double *, const double *);

double * vec_allocate(size_t);

double * mat_allocate(size_t num_rows, size_t num_cols);

double * mat_identity(size_t);

void mat_multiply(size_t num_rows1, size_t num_cols1, const double * mat1,
		  size_t num_rows2, size_t num_cols2, const double * mat2,
		  double * mat_out);

void mat_trim_and_normalize_cols(size_t num_rows, size_t num_cols, double * mat,
				 size_t num_out_cols, double * mat_out);

double * mat_copy_cols(size_t num_rows, size_t num_cols, const double * mat,
		       size_t num_out_cols, size_t * indices);

void mat_normalize_rows(size_t num_rows, size_t num_cols, double * mat);

bool mat_abs_equals(size_t num_rows, size_t num_cols, const double * mat1, const double * mat2);

bool mat_equals(size_t num_rows, size_t num_cols, const double * mat1, const double * mat2);

#endif /* __MAT_H */
