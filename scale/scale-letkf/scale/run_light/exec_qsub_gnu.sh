#!/bin/sh
#PBS -N test_dacycle
#PBS -q m
#PBS -l nodes=12:ppn=4
#PBS -l walltime=01:30:00
#

cd ${PBS_O_WORKDIR}

export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

source /etc/profile.d/modules.sh
module unload mpt/2.12
module unload intelcompiler/16.0.1.150
module unload intelmpi/5.1.2.150
module unload hdf5/1.8.16-intel
module unload netcdf4/4.3.3.1-intel
module unload netcdf4/fortran-4.4.2-intel
module load gcc/5.3.0
module load openmpi/2.1.6-gcc
module load hdf5/1.10.8
module load netcdf4/4.6.2
module load netcdf4/fortran-4.4.5
module load lapack/3.6.0

export OMP_NUM_THREADS=1
export KMP_AFFINITY=compact

export LD_LIBRARY_PATH="/home/seiya/lib:$LD_LIBRARY_PATH"

ulimit -s unlimited
umask 0007

#mpicommand="mpirun --mca btl openib,sm,self --bind-to core"
mpicommand="mpirun --mca btl tcp,vader,self --bind-to core"

echo "scale-rm_init_ens" 
 $mpicommand ./scale-rm_init_ens config/scale-rm_init_ens_20220101000000.conf 
echo "scale-rm_ens"
 $mpicommand ./scale-rm_ens config/scale-rm_ens_20220101000000.conf 
echo "copy restart files"
for mem in $(seq -f %04g 1 5) mean;do
  for pe in $(seq -f %06g 0 7);do
    cp ${mem}/gues/init_20220101-060000.000.pe${pe}.nc ${mem}/anal/
  done
done
echo "letkf"
 $mpicommand ./letkf config/letkf_20220101060000.conf >> run_progress
echo "done."
