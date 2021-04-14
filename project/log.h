#ifndef __LOG_H
#define __LOG_H

#include <time.h>
extern struct timespec ts;
extern struct timespec base_ts;

#define log_emit(msg)							\
  do									\
    {									\
      clock_gettime(CLOCK_MONOTONIC, &ts);				\
      printf("[%lu] %s:%s:%u - %s\n",					\
	     ts.tv_sec - base_ts.tv_sec,				\
	     __FILE__, __FUNCTION__, __LINE__,				\
	     msg);							\
    }									\
  while (0)

void log_init();

#endif
