#!/bin/sh 
#
#PJM -g <--GROUP--> 
#PJM -x PJM_LLIO_GFSCACHE=/vol0003:/vol0004
#PJM -L "rscgrp=small"
#PJM -L "node=<--NNODES-->"
#PJM -L "elapse=00:10:00"
#PJM --mpi "max-proc-per-node=<--PPN-->"
#PJM -j
#PJM -s

export OMP_NUM_THREADS=<--NTHREADS--> 
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
#export FLIB_BARRIER=HARD

export LD_LIBRARY_PATH=/lib64:/opt/FJSVxtclanga/tcsds-mpi-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:`cat /home/apps/oss/scale/llio.list | sed 's:\(.*/lib\)/.*:\1:' | uniq | sed -z 's/\n/:/g'`

###pp### mpiexec -n <--NPROCS_PP--> -std-proc log/NOUT_scale-rm_pp_ens bin/scale-rm_pp_ens conf/scale-rm_pp_ens_20210730060000.conf
mpiexec -n <--NPROCS--> -std-proc log/NOUT_scale-rm_init_ens bin/scale-rm_init_ens conf/scale-rm_init_ens_20210730060000.conf

echo "SCALE-RM   run starting ... `date`"
startms=$(date +'%s.%3N')
mpiexec -n <--NPROCS--> -std-proc log/NOUT_scale-rm_ens bin/scale-rm_ens conf/scale-rm_ens_20210730060000.conf
endms=$(date +'%s.%3N')
elapse=$(echo "$endms $startms" | awk '{printf "%.3f\n", $1 - $2}')
echo "SCALE-RM   run ending ... elapse (`echo $elapse` sec)"

echo "LETKF      run starting ... `date`"
startms=$(date +'%s.%3N')
mpiexec -n <--NPROCS--> -std-proc log/NOUT_letkf bin/letkf conf/letkf_20210730060030.conf
endms=$(date +'%s.%3N')
elapse=$(echo "$endms $startms" | awk '{printf "%.3f\n", $1 - $2}')
echo "LETKF      run ending ... elapse (`echo $elapse` sec)"
