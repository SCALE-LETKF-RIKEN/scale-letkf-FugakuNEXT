#!/bin/sh

WKDIR=$(cd $(dirname $0) ; pwd)

GROUP="ra010016" # your group

MEMBER=2 ### up to 50
PRC_NUM_X=10
PRC_NUM_Y=10

IMAXG=1280
JMAXG=1280
IMAX=$((IMAXG / PRC_NUM_X))
JMAX=$((JMAXG / PRC_NUM_Y))

MAX_PROCS_PER_NODE=48 ### Fugaku : fixed
NTHREADS=12 ### Fugaku fixed
PPN=$((MAX_PROCS_PER_NODE / NTHREADS))
NPROCS_PER_MEM=$((PRC_NUM_X * PRC_NUM_Y))

NPROCS=$((NPROCS_PER_MEM * (MEMBER + 1)))
NNODES=$(( (NPROCS-1) / PPN + 1 ))

mkdir -p $WKDIR/bin
mkdir -p $WKDIR/log
mkdir -p $WKDIR/conf

#cp $TOPDIR/scale/scale-letkf/scale/ensmodel/scale-rm*_ens $WKDIR/bin/
#cp $TOPDIR/scale/scale-letkf/scale/letkf/letkf $WKDIR/bin/

mkdir -p $WKDIR/result/SC23/const/topo
mkdir -p $WKDIR/result/SC23/const/landuse

if [ -f $WKDIR/result/SC23/const/topo/topo.pe000000.nc ] ;then
  if [ $(/bin/ls -x1 $WKDIR/result/SC23/const/topo/ | wc -l | cut -f 1 -d " ") == $NPROCS_PER_MEM ] ;then
    RUN_PP=0
  else
    echo "The data with different NPROCS_PER_MEM exists in the output directory !"
    echo "Rename $WKDIR/result/SC23 to store the existing data, or remove topography files."
    exit 1
  fi
else
  RUN_PP=1
  mkdir -p $WKDIR/result/SC23/const/topo
  mkdir -p $WKDIR/result/SC23/const/landuse
fi

mkdir -p $WKDIR/result/SC23/20210730060030/obs
for mem in mean $(seq -f %04g 1 $MEMBER) ;do 
  mkdir -p $WKDIR/result/SC23/20210730060000/anal/$mem
  mkdir -p $WKDIR/result/SC23/20210730060000/bdy/$mem
  mkdir -p $WKDIR/result/SC23/20210730060030/anal/$mem
done

if [ $RUN_PP == 1 ]; then
mem=0001
mkdir -p conf/$mem
cat template/pp.conf | \
  sed -e "s#<--WKDIR-->#$WKDIR#g" | \
  sed -e "s#<--PRC_NUM_X-->#$PRC_NUM_X#g" | \
  sed -e "s#<--PRC_NUM_Y-->#$PRC_NUM_Y#g" | \
  sed -e "s#<--IMAX-->#$IMAX#g" | \
  sed -e "s#<--JMAX-->#$JMAX#g" | \
  sed -e "s#<--MEM-->#$mem#g" \
> conf/$mem/pp.d01_20210730060000.conf
cat template/scale-rm_pp_ens.conf | \
sed -e "s/<--MEMBER-->/$MEMBER/g" | \
sed -e "s/<--PPN-->/$PPN/g" | \
sed -e "s/<--NPROCS_PER_MEM-->/$NPROCS_PER_MEM/g" \
> conf/scale-rm_pp_ens_20210730060000.conf
fi

for mem in mean $(seq -f %04g 1 $MEMBER);do
  mkdir -p conf/$mem
  cat template/init.conf | \
  sed -e "s#<--WKDIR-->#$WKDIR#g" | \
  sed -e "s#<--PRC_NUM_X-->#$PRC_NUM_X#g" | \
  sed -e "s#<--PRC_NUM_Y-->#$PRC_NUM_Y#g" | \
  sed -e "s#<--IMAX-->#$IMAX#g" | \
  sed -e "s#<--JMAX-->#$JMAX#g" | \
  sed -e "s#<--MEM-->#$mem#g" \
> conf/$mem/init.d01_20210730060000.conf
  cat template/run.conf | \
  sed -e "s#<--WKDIR-->#$WKDIR#g" | \
  sed -e "s#<--PRC_NUM_X-->#$PRC_NUM_X#g" | \
  sed -e "s#<--PRC_NUM_Y-->#$PRC_NUM_Y#g" | \
  sed -e "s#<--IMAX-->#$IMAX#g" | \
  sed -e "s#<--JMAX-->#$JMAX#g" | \
  sed -e "s#<--MEM-->#$mem#g" \
> conf/$mem/run.d01_20210730060000.conf
done

cat template/scale-rm_init_ens.conf | \
sed -e "s/<--MEMBER-->/$MEMBER/g" | \
sed -e "s/<--MEMBER_RUN-->/$((MEMBER+1))/g" | \
sed -e "s/<--PPN-->/$PPN/g" | \
sed -e "s/<--NPROCS_PER_MEM-->/$NPROCS_PER_MEM/g" \
> conf/scale-rm_init_ens_20210730060000.conf

cat template/scale-rm_ens.conf | \
sed -e "s/<--MEMBER-->/$MEMBER/g" | \
sed -e "s/<--MEMBER_RUN-->/$((MEMBER+1))/g" | \
sed -e "s/<--PPN-->/$PPN/g" | \
sed -e "s/<--NPROCS_PER_MEM-->/$NPROCS_PER_MEM/g" \
> conf/scale-rm_ens_20210730060000.conf

cat template/letkf.conf | \
sed -e "s#<--WKDIR-->#$WKDIR#g" | \
sed -e "s/<--MEMBER-->/$MEMBER/g" | \
sed -e "s/<--MEMBER_RUN-->/$((MEMBER+1))/g" | \
sed -e "s/<--PPN-->/$PPN/g" | \
sed -e "s#<--PRC_NUM_X-->#$PRC_NUM_X#g" | \
sed -e "s#<--PRC_NUM_Y-->#$PRC_NUM_Y#g" | \
sed -e "s#<--IMAX-->#$IMAX#g" | \
sed -e "s#<--JMAX-->#$JMAX#g" | \
sed -e "s/<--NPROCS_PER_MEM-->/$NPROCS_PER_MEM/g" \
> conf/letkf_20210730060030.conf

cat template/exec.sh | \
sed -e "s/<--GROUP-->/$GROUP/g" | \
sed -e "s/<--NNODES-->/$NNODES/g" | \
sed -e "s/<--PPN-->/$PPN/g" | \
sed -e "s/<--NTHREADS-->/$NTHREADS/g" | \
sed -e "s/<--NPROCS-->/$NPROCS/g" \
> exec.sh  

if [ $RUN_PP == 1 ] ;then
  sed -i -e "s/###pp###//g" exec.sh
  sed -i -e "s/<--NPROCS_PP-->/$NPROCS_PER_MEM/g" exec.sh
fi
 
