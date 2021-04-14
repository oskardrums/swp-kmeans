#ifndef __KMPP_H
#define __KMPP_H

#include <stdlib.h>

size_t * kmpp(size_t num_clusters, size_t num_rows, size_t num_cols, const double * mat, size_t max_iters);

#endif
