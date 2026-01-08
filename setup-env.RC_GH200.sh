#!/bin/bash

module purge
module load system/qc-gh200
module load nvhpc/25.9

SYS_LIB="/lvs0/rccs-nghpcadu/CX_input/SCALE-LETKF/lib_qc-gh200"

SCALE_HDF="$SYS_LIB/hdf5/1.10.10"
SCALE_NETCDF_C="$SYS_LIB/netcdf-c/4.9.2"
SCALE_NETCDF_F="$SYS_LIB/netcdf-fortran/4.6.1"
SCALE_LAPACK="$SYS_LIB/lapack/3.10.1"
SCALE_SZIP="$SYS_LIB/szip/2.1.1"

export SCALE_SYS="Linux64-nvidia"
export SCALE_NETCDF_INCLUDE="-I${SCALE_NETCDF_C}/include -I${SCALE_NETCDF_F}/include"
export SCALE_NETCDF_LIBS="-L${SCALE_NETCDF_C}/lib -L${SCALE_NETCDF_F}/lib -L${SCALE_HDF}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5"
export SCALE_MATHLIB_LIBS="-L${SCALE_LAPACK}/lib -llapack -lrefblas"
export SCALE_ENABLE_MATHLIB=T
export SCALE_ENABLE_OPENMP=F
export SCALE_ENABLE_OPENACC=T

export PATH="${SCALE_NETCDF_C}/bin:${SCALE_NETCDF_F}/bin:$PATH"
export LD_LIBRARY_PATH="${SCALE_HDF}/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_SZIP}/lib:${SCALE_LAPACK}/lib:$LD_LIBRARY_PATH"
