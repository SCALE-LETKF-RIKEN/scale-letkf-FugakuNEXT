#!/bin/sh 
#
#PJM -g <--GROUP--> 
#PJM -x PJM_LLIO_GFSCACHE=/vol0003:/vol0004
#PJM -L "rscgrp=small"
#PJM -L "node=<--NNODES-->"
#PJM -L "elapse=00:05:00"
#PJM --mpi "max-proc-per-node=<--PPN-->"
#PJM -j
#PJM -s

export OMP_NUM_THREADS=<--NTHREADS--> 
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
#export FLIB_BARRIER=HARD

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-c-4.9.0-g462kcd2ivou7ewax6wddywoyrbz2oib/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-fortran-4.6.0-mmdtg5243y4mwqsl3gcu3m2kh27raq5n/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/parallel-netcdf-1.12.3-avpnzm4pwv2tuu2mv73lacb4vhcwlnds/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/hdf5-1.12.2-kb4msz2kuwzsmqsshhpryqebui6tqcfs/lib:$LD_LIBRARY_PATH

mpiexec -n <--NPROCS--> -std-proc log/NOUT_scale-rm_ens bin/scale-rm_ens conf/scale-rm_ens_20000101002500.conf
mpiexec -n <--NPROCS--> -std-proc log/NOUT_letkf bin/letkf conf/letkf_20000101003000.conf

