#include "log.h"

struct timespec ts;
struct timespec base_ts;

void log_init()
{
  clock_gettime(CLOCK_MONOTONIC, &base_ts);
}
