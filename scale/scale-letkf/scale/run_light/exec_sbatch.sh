#!/bin/bash
####SBATCH -A MST110250 
####SBATCH -p t3atm_NRQ
#SBATCH --job-name="test_dacycle"
#SBATCH --exclusive
#SBATCH -o test_dacycle
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=16
#SBATCH -n 48
#SBATCH -N 3

ulimit -s unlimited
umask 0007

module purge
module load  compiler/intel/2020u4 IntelMPI/2020 netcdf/4.7.4

echo "scale-rm_init_ens"
mpirun ./scale-rm_init_ens config/scale-rm_init_ens_20220101000000.conf 
echo "scale-rm_ens"
mpirun ./scale-rm_ens config/scale-rm_ens_20220101000000.conf 
for mem in $(seq -f %04g 1 5) mean;do
  for pe in $(seq -f %06g 0 7);do
    cp ${mem}/gues/init_20220101-060000.000.pe${pe}.nc ${mem}/anal/
  done
done
echo "letkf"
mpirun ./letkf config/letkf_20220101060000.conf 
echo "done."
