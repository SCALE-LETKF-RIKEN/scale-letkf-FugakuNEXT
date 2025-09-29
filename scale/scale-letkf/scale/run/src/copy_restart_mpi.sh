#!/bin/bash

## 
df -t lliofs &> /dev/null
###

. ./config.main || exit $?
. ./src/func_datetime.sh || exit $?
. ./src/func_util.sh || exit $?
. ./node/distr || exit $?

DIR_IN=$1
DIR_OUT=$2
ATIME=$3

#
GLOBAL_RANK=${PMIX_RANK}
irank=$((GLOBAL_RANK+1))
imem=${proc2group[$irank]}
ipe=${proc2grpproc[$irank]}
ipe=$((ipe-1))

if [ ! -z "$imem" ];then
if [ $imem -le $MEMBER ]; then
  mem=$(printf %04d $imem)
else
  mem="mean" 
fi
fi

ifile=$DIR_IN/$mem/init_$(datetime_scale $ATIME).pe$(printf %06d $ipe).nc 
ofile=$DIR_OUT/$mem/init_$(datetime_scale $ATIME).pe$(printf %06d $ipe).nc 
if [ -f $ifile ]; then
#  echo "===CHECK=== $GLOBAL_RANK $imem $ipe copy file from $ifile to $ofile" > check.$GLOBAL_RANK
  mkdir -p $(dirname $ofile) 
  cp $ifile $ofile 
#else
#  echo "===CHECK=== $GLOBAL_RANK $imem $ipe no file $ifile > check.$GLOBAL_RANK"
fi
