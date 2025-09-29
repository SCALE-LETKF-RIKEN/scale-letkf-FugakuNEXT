#!/bin/sh
#PBS -N test_fcst
#PBS -q m
#PBS -l nodes=12:ppn=4
#PBS -l walltime=01:30:00
#

cd ${PBS_O_WORKDIR}

export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

source /etc/profile.d/modules.sh
module unload mpt/2.12
module load intelmpi/5.1.2.150
module load hdf5/1.8.16-intel
module load netcdf4/4.3.3.1-intel
module load netcdf4/fortran-4.4.2-intel

export OMP_NUM_THREADS=1
export KMP_AFFINITY=compact

export LD_LIBRARY_PATH="/home/seiya/lib:$LD_LIBRARY_PATH"

ulimit -s unlimited
umask 0007

echo "scale-rm_init_ens"
impijob ./scale-rm_init_ens config/scale-rm_init_ens_20220101000000.conf 
echo "scale-rm_ens"
impijob ./scale-rm_ens config/scale-rm_ens_20220101000000.conf 
echo "copy restart files"
for mem in $(seq -f %04g 1 5) mean;do
  for pe in $(seq -f %06g 0 7);do
    cp ${mem}/gues/init_20220101-060000.000.pe${pe}.nc ${mem}/anal/
  done
done
echo "letkf"
impijob ./letkf config/letkf_20220101060000.conf >> run_progress
echo "done."
