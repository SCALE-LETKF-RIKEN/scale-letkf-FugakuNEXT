#!/bin/sh

### prepare symbolic links to binary and data files

cd "$(dirname "$0")"
mydir=$(pwd)
myname="$(basename "$0")"

######
LETKFDIR=${mydir}/..
SCALEDIR=${mydir}/../../..
DATADIR="/share/hp150019/scale_database/scale-letkf-test-suite"
#DATADIR="${HOME}/scale_database/scale-letkf-test-suite"
copy="cp" ### tentative for Fugaku
#copy="ln -s"
######

### binary

ln -s $LETKFDIR/ensmodel/scale-rm_ens . 
ln -s $LETKFDIR/ensmodel/scale-rm_init_ens . 
ln -s $LETKFDIR/letkf/letkf . 

### topo and landuse
for item in topo landuse ; do
  mkdir -p const/$item
  for pe in $(seq -f %06g 0 7) ; do
    $copy $DATADIR/exp/18km_Japan/init/const/$item/${item}.pe${pe}.nc ./const/$item/
  done
done

### SCALE data
rm -rf ./dat
ln -s $SCALEDIR/data ./dat

### boundary orig
timef="2022-01-01 00:00:00"
time=$(date -ud "$timef" +%Y%m%d%H%M%S)
for it in $(seq 0 2); do
  atime=$(date -ud "$((it*21600)) second $timef" +%Y%m%d%H%M%S)
  $copy $DATADIR/exp/18km_Japan/bdy/${time}/mean/atm_${atime}.grd  mean/bdyorg_atm_${time:0:8}-${time:8:6}.000_0000${it}.grd
  $copy $DATADIR/exp/18km_Japan/bdy/${time}/mean/land_${atime}.grd mean/bdyorg_lnd_${time:0:8}-${time:8:6}.000_0000${it}.grd
  $copy $DATADIR/exp/18km_Japan/bdy/${time}/mean/sfc_${atime}.grd  mean/bdyorg_sfc_${time:0:8}-${time:8:6}.000_0000${it}.grd
done

cp mean/gradsbdy.conf.temp mean/gradsbdy.conf
sed -i -e "s#<--bdyorg_atm-->#${mydir}/mean/bdyorg_atm_${time:0:8}-${time:8:6}.000#g" mean/gradsbdy.conf
sed -i -e "s#<--bdyorg_sfc-->#${mydir}/mean/bdyorg_sfc_${time:0:8}-${time:8:6}.000#g" mean/gradsbdy.conf
sed -i -e "s#<--bdyorg_lnd-->#${mydir}/mean/bdyorg_lnd_${time:0:8}-${time:8:6}.000#g" mean/gradsbdy.conf

### init (prepared)
for mem in mean $(seq -f %04g 1 5) ;do
  mkdir -p $mem/gues
  mkdir -p $mem/anal
  mkdir -p $mem/hist
  for pe in $(seq -f %06g 0 7) ; do
    $copy $DATADIR/exp/18km_Japan/init/$time/anal/${mem}/init_${time:0:8}-${time:8:6}.000.pe${pe}.nc ./$mem/anal/
  done
done

### observation
 
mkdir -p ./obs
atime=$(date -ud "21600 second $timef" +%Y%m%d%H%M%S)
$copy ${DATADIR}/obs/prepbufr_Japan/obs_${atime}.dat ./obs/

### log 
mkdir -p log/scale_init log/scale log/letkf
