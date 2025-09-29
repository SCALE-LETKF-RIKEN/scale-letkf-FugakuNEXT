#!/bin/bash

. /vol0004/apps/oss/spack/share/spack/setup-env.sh

NC_HASH=`spack find -lx netcdf-c%fj@4.10.0 | grep netcdf-c | awk '{print $1}'`
NF_HASH=`spack find -lx netcdf-fortran%fj@4.10.0 | grep netcdf-fortran | awk '{print $1}'`
HDF_HASH=`spack find -l --deps /${NC_HASH} | grep hdf5 | awk '{print $1}'`
SCALE_HDF=`spack location --install-dir /${HDF_HASH}`
SCALE_NETCDF_C=`spack location --install-dir /${NC_HASH}`
SCALE_NETCDF_F=`spack location --install-dir /${NF_HASH}`
SCALE_PNETCDF=`spack location --install-dir parallel-netcdf%fj@4.10.0`

export SCALE_DB="`pwd`/scale_database"
export SCALE_NETCDF_INCLUDE="-I${SCALE_NETCDF_C}/include -I${SCALE_NETCDF_F}/include -I${SCALE_PNETCDF}/include"
export SCALE_NETCDF_LIBS="-L${SCALE_NETCDF_C}/lib -L${SCALE_NETCDF_F}/lib -L${SCALE_HDF}/lib -L${SCALE_PNETCDF}/lib -lpnetcdf -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lfjprofmpi -lmpi_cxx"
