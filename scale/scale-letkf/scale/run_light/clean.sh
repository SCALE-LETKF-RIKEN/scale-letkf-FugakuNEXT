#!/bin/sh

### prepare symbolic links to binary and data files

cd "$(dirname "$0")"
mydir=$(pwd)
myname="$(basename "$0")"

### binary
rm ./scale-rm_ens ./scale-rm_init_ens ./letkf  

### topo and landuse
rm -r ./const 

### SCALE data
rm ./dat

### boundary orig
rm mean/bdyorg*.grd
rm mean/gradsbdy.conf

### netcdf files
for mem in mean $(seq -f %04g 1 5) ;do
  rm -r ./$mem/anal
  rm -r ./$mem/gues
  rm -r ./$mem/hist
done
rm ./mean/boundary*.nc

### observation
rm -r ./obs

### log 
rm -r ./log 
