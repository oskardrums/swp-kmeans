/* header file for Normalized Spectral Clustering */
#ifndef __NSC_H
#define __NSC_H

#include <stdlib.h>

size_t normalized_spectral_clustering(size_t num_clusters,
				      size_t num_rows,
				      size_t num_cols,
				      const double * mat,
				      size_t max_iters,
				      size_t ** labels);

#endif
