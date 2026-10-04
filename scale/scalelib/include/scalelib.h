#ifndef SCALE_H
#define SCALE_H

#include "scale_log.h"
#include "scale_openmp.h"

#ifdef _OPENACC

#ifndef LSIZE
#define LSIZE 1
#endif

#else

#ifndef LSIZE
#ifdef SINGLE
#define LSIZE 16
#else
#define LSIZE 8
#endif
#endif

#endif

#endif
