#!/bin/sh

WKDIR=$(cd $(dirname $0) ;pwd)
TOPDIR=$(cd $WKDIR/../.. ; pwd)

NPROCS_PER_MEM=4 ### fixed
MAX_PROCS_PER_NODE=48 ### fixed

GROUP="ra010016" # your group

MEMBER=10
NTHREADS=12
PPN=$((MAX_PROCS_PER_NODE / NTHREADS))

NPROCS=$((NPROCS_PER_MEM * (MEMBER + 1)))
NNODES=$(( (NPROCS-1) / PPN + 1 ))

mkdir -p $WKDIR/bin
mkdir -p $WKDIR/log
mkdir -p $WKDIR/conf

for mem in mean $(seq -f %04g 1 $MEMBER);do
  mkdir -p $TOPDIR/result/supercell/20000101003000/anal/$mem
done
mkdir -p $TOPDIR/result/supercell/20000101003000/obs

cp $TOPDIR/scale/scale-letkf/scale/ensmodel/scale-rm*_ens $WKDIR/bin/
cp $TOPDIR/scale/scale-letkf/scale/letkf/letkf $WKDIR/bin/

for mem in mean $(seq -f %04g 1 $MEMBER);do
  mkdir -p conf/$mem
  cat template/run.conf | \
  sed -e "s#<--TOPDIR-->#$TOPDIR#g" | \
  sed -e "s#<--MEM-->#$mem#g" \
> conf/$mem/run.d01_20000101002500.conf
done

cat template/scale-rm_ens.conf | \
sed -e "s/<--MEMBER-->/$MEMBER/g" | \
sed -e "s/<--MEMBER_RUN-->/$((MEMBER+1))/g" | \
sed -e "s/<--PPN-->/$PPN/g" \
> conf/scale-rm_ens_20000101002500.conf

cat template/letkf.conf | \
sed -e "s#<--TOPDIR-->#$TOPDIR#g" | \
sed -e "s/<--MEMBER-->/$MEMBER/g" | \
sed -e "s/<--PPN-->/$PPN/g" \
> conf/letkf_20000101003000.conf

cat template/exec.sh | \
sed -e "s/<--GROUP-->/$GROUP/g" | \
sed -e "s/<--NNODES-->/$NNODES/g" | \
sed -e "s/<--PPN-->/$PPN/g" | \
sed -e "s/<--NTHREADS-->/$NTHREADS/g" | \
sed -e "s/<--NPROCS-->/$NPROCS/g" \
> exec.sh   
