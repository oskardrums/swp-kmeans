/*
 * Jaccard similiarity measure calculation
 */
#include "jaccard.h"
#include <stdbool.h>

/*
 * Returns the Jaccard measure of the two label sets.
 *
 * O(num_labels^2)
 */
double jaccard_measure(size_t num_labels, const size_t * labels1, const size_t * labels2)
{
  size_t i, j;
  size_t inter_size = 0, union_size = 0;
  bool cond1 = false, cond2 = false;

  for (i = 0; i < num_labels; ++i) {
    for (j = i + 1; j < num_labels; ++j) {

      cond1 = (labels1[i] == labels1[j]);
      cond2 = (labels2[i] == labels2[j]);

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
