#ifndef __MAT_H
#define __MAT_H

#include <stddef.h>

void mat_dump(size_t num_rows, size_t num_cols, double * mat);

double mat_col_norm_squared(size_t col, size_t num_rows, size_t num_cols, double * mat);

double vec_distance_squared(size_t, double *, double *);

double * vec_allocate(size_t);

double * mat_allocate(size_t num_rows, size_t num_cols);

double * mat_identity(size_t);

void mat_multiply(size_t num_rows1, size_t num_cols1, double * mat1,
		  size_t num_rows2, size_t num_cols2, double * mat2,
		  double * mat_out);

void mat_trim_and_normalize_cols(size_t num_rows, size_t num_cols, double * mat,
				 size_t num_out_cols, double * mat_out);

#endif /* __MAT_H */
