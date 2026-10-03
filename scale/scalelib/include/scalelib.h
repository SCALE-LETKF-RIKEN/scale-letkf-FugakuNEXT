#ifndef SCALE_H
#define SCALE_H

#include "scale_log.h"
#include "scale_openmp.h"

#ifdef _OPENACC

#ifndef LSIZE
#define LSIZE 1
#endif

! upper limit of KA for the column work arrays of fixed size in the kernels,
! which can be placed in the shared memory.
! KA larger than this is rejected at the setup; rebuild with -DKAMAX=<n> (n >= KA).
! An odd number avoids the bank conflicts when the vector lanes access
! the different columns of a work array (KAMAX,n) at the same k.
#ifndef KAMAX
#define KAMAX 129
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
