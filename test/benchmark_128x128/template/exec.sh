#!/bin/bash
#
#SBATCH --job-name=scale-letkf
#SBATCH --partition=qc-gh200
#SBATCH -t 24:00:00
#
module purge
module load system/qc-gh200
module load nvhpc/25.9

. /lvs0/rccs-sdt/tyamaura/lib_qc-gh200/setup-env.sh
export OMP_NUM_THREADS=1

###pp### mpiexec -n <--NPROCS_PP--> bin/scale-rm_pp_ens conf/scale-rm_pp_ens_20210730060000.conf
mpiexec -n <--NPROCS--> bin/scale-rm_init_ens conf/scale-rm_init_ens_20210730060000.conf

echo "SCALE-RM   run starting ... `date`"
startms=$(date +'%s.%3N')
mpiexec -n <--NPROCS--> bin/scale-rm_ens conf/scale-rm_ens_20210730060000.conf
endms=$(date +'%s.%3N')
elapse=$(echo "$endms $startms" | awk '{printf "%.3f\n", $1 - $2}')
echo "SCALE-RM   run ending ... elapse (`echo $elapse` sec)"

echo "LETKF      run starting ... `date`"
startms=$(date +'%s.%3N')
mpiexec -n <--NPROCS--> --oversubscribe bin/letkf conf/letkf_20210730060030.conf
endms=$(date +'%s.%3N')
elapse=$(echo "$endms $startms" | awk '{printf "%.3f\n", $1 - $2}')
echo "LETKF      run ending ... elapse (`echo $elapse` sec)"
