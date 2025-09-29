#!/bin/sh 
#
#PJM -g hp150019
#PJM -x PJM_LLIO_GFSCACHE=/vol0003:/vol0004
#PJM -L "rscgrp=int"
#PJM -L "node=12"
#PJM -L "elapse=00:20:00"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -s
#
export PARALLEL=12
export OMP_NUM_THREADS=${PARALLEL}
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-c-4.9.0-g462kcd2ivou7ewax6wddywoyrbz2oib/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-fortran-4.6.0-mmdtg5243y4mwqsl3gcu3m2kh27raq5n/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/parallel-netcdf-1.12.3-avpnzm4pwv2tuu2mv73lacb4vhcwlnds/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/hdf5-1.12.2-kb4msz2kuwzsmqsshhpryqebui6tqcfs/lib:$LD_LIBRARY_PATH

echo "scale-rm_init_ens"
mpiexec -std-proc log/scale_init/NOUT -n 48 ./scale-rm_init_ens config/scale-rm_init_ens_20220101000000.conf
echo "scale-rm_ens"
mpiexec -std-proc log/scale/NOUT -n 48 ./scale-rm_ens config/scale-rm_ens_20220101000000.conf
echo "copy restart files"
for mem in $(seq -f %04g 1 5) mean;do
  for pe in $(seq -f %06g 0 7);do
    cp ${mem}/gues/init_20220101-060000.000.pe${pe}.nc ${mem}/anal/
  done
done
echo "letkf"
mpiexec -std-proc log/letkf/NOUT -n 48 ./letkf config/letkf_20220101060000.conf
echo "done."
