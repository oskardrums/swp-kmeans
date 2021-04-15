#include "jaccard.h"
#include <stdbool.h>

double jaccard_measure(size_t n, const size_t * v1, const size_t * v2)
{
  size_t i, j;
  size_t inter_size = 0, union_size = 0;
  bool cond1 = false, cond2 = false;

  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {

      cond1 = (v1[i] == v1[j]);
      cond2 = (v2[i] == v2[j]);

      if (cond1 && cond2) {
	++inter_size;
      }

      if (cond1 || cond2) {
	++union_size;
      }
    }
  }
  return ((double)inter_size / union_size);
}
