#!/bin/bash
#===============================================================================
#
#  Nature run 
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_run.sh [..]
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"

#===============================================================================
# Configuration

. ./config.main.${PRESET} || exit $?
. ./config.obsmake || exit $?

. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_common_static.sh || exit $?
. src/func_fcst_static.sh || exit $?

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $myname $@"

setting "$@" || exit $?

echo
print_setting || exit $?
echo

#===============================================================================
# Create and clean the temporary directory

echo "[$(datetime_now)] Create and clean the temporary directory"
 
safe_init_tmpdir $TMP || exit $?

#===============================================================================
# Determine the distibution schemes

echo "[$(datetime_now)] Determine the distibution schemes"

safe_init_tmpdir $NODEFILE_DIR || exit $?
#distribute_fcst "$MEMBERS" $CYCLE - $NODEFILE_DIR || exit $?

if ((CYCLE == 0)); then
  CYCLE=$cycle_auto
fi

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

cat $SCRP_DIR/config.main.${PRESET} | \
    sed -e "/\(^DIR=\| DIR=\)/c DIR=\"$DIR\"" \
    > $TMP/config.main

echo "SCRP_DIR=\"\$TMPROOT\"" >> $TMP/config.main
echo "NODEFILE_DIR=\"\$TMPROOT/node\"" >> $TMPS/config.main
echo "RUN_LEVEL=4" >> $TMP/config.main

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMP/config.main

safe_init_tmpdir $STAGING_DIR || exit $?
#staging_list_static || exit $?

#config_file_list $TMPS/config || exit $?

#-------------------------------------------------------------------------------
# Add shell scripts and node distribution files into the staging list

cp -r ${SCRP_DIR}/src $TMP/src
cp ${SCRP_DIR}/config.rc $TMP/config.rc

#===============================================================================
# Stage in

echo "[$(datetime_now)] Initialization (stage in)"

stage_in server || exit $?

mkdir -p $TMP/log
mkdir -p $TMP/obsin
mkdir -p $TMP/config
mkdir -p $OBS

ln -s $DIR/obs/obsmake $TMP/. 
cp $OBSIN $TMP/obsin/obsin.dat

time=$STIME
time_list=""
while ((time <= ETIME)); do
nslot=1
timet=$BASETIME
while (( timet < time )); do
  nslot=$((nslot + 1))
  timet=$(datetime $timet $FCSTOUT)
done
if (( timet > time )); then 
  echo "observation time does not match with history file !"
  exit 1 
fi

conf_file=$TMP/config/obsmake_${time}.conf
cat $SCRP_DIR/config.nml.obsmake | sed \
    -e "/!--PPN--/a PPN=$PPN,"  \
    -e "/!--PRC_DOMAINS--/a PRC_DOMAINS=$SCALE_NP,"  \
    -e "/!--OBS_IN_NAME--/a OBS_IN_NAME=\"$TMP/obsin/obsin.dat\","  \
    -e "/!--OBS_IN_FORMAT--/a OBS_IN_FORMAT=\"${OBS_IN_FORMAT}\","  \
    -e "/!--LETKF_TOPOGRAPHY_IN_BASENAME--/a LETKF_TOPOGRAPHY_IN_BASENAME=\"$OUTDIR/const/topo/topo\"," \
    -e "/!--HISTORY_IN_BASENAME--/a HISTORY_IN_BASENAME=\"$OUTDIR/nature/hist/history\"," \
    -e "/!--GUES_IN_BASENAME--/a GUES_IN_BASENAME=\"$OUTDIR/nature/init/init_$(datetime_scale $time)\"," \
    -e "/!--SLOT_START--/a SLOT_START=$nslot,"  \
    -e "/!--SLOT_END--/a SLOT_END=$nslot,"  \
    -e "/!--SLOT_BASE--/a SLOT_BASE=$nslot,"  \
    -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL=$FCSTOUT,"  \
> $conf_file

cat $SCRP_DIR/config.nml.scale | sed \
    -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
    -e "/!--TIME_DURATION--/a TIME_DURATION = ${LCYCLE}.D0," \
    -e "/!--RESTART_IN_BASENAME--/a RESTART_IN_BASENAME=\"$OUTDIR/nature/init/init_$(datetime_scale $time)\","  \
    -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME=\"$OUTDIR/nature/init/init\","  \
    -e "/!--TOPOGRAPHY_IN_BASENAME--/a TOPOGRAPHY_IN_BASENAME=\"$OUTDIR/const/topo/topo\"," \
    -e "/!--FILE_HISTORY_DEFAULT_BASENAME--/a FILE_HISTORY_DEFAULT_BASENAME=\"$OUTDIR/nature/hist/history\"," \
>> $conf_file

  time_list="$time_list $time"
  time=$(datetime $time $LCYCLE s)
done

# Overwrite NNODES
NNODES_USE=$(( SCALE_NP / PPN ))
echo "NNODES=$NNODES_USE" >> $TMP/config.main

#===============================================================================
# Creat a job script

jobscrp="$TMP/obsmake_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

# FUGAKU
if [ "$PRESET" = 'FUGAKU' ]; then

  if (( NNODES_USE > 384 )) ; then
    RSCGRP="large"
  else
    RSCGRP="small"
  fi
  TPROC=$((NNODES_USE*PPN))

  CVOLUME=$(realpath $(pwd) | cut -d "/" -f 2) # current volume (e.g., /vol0X0Y or /vol000X)
  NUM_VOLUME=${CVOLUME:4:1} # get number of current volume 

  if [ "$NUM_VOLUME" = "0" ] ; then
    VOLUMES="/"${CVOLUME}
  else
    VOLUMES="/vol000${NUM_VOLUME}"
  fi

  if [ $VOLUMES != "/vol0004" ] ;then
    VOLUMES="${VOLUMES}:/vol0004" # spack
  fi

cat > $jobscrp << EOF
#!/bin/sh 
#
#PJM -g ${GROUP} 
#PJM -x PJM_LLIO_GFSCACHE=${VOLUMES}
#PJM -L "rscgrp=${RSCGRP}"
#PJM -L "node=$(((TPROC+3)/4))"
#PJM -L "elapse=${TIME_LIMIT}"
#PJM --mpi "max-proc-per-node=${PPN}"
#PJM -j
#PJM -s
#
#
export PARALLEL=${THREADS}
export OMP_NUM_THREADS=\${PARALLEL}
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD

EOF

  if (( USE_LLIO_BIN == 1 )); then
      echo "llio_transfer obsmake" >> $jobscrp 
    echo "" >> $jobscrp
  fi

  if (( USE_SPACK > 0 )); then

cat << EOF >>  $jobscrp 
SPACK_FJVER=${SPACK_FJVER}
. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load netcdf-c%fj@\${SPACK_FJVER}
spack load netcdf-fortran%fj@\${SPACK_FJVER}
spack load parallel-netcdf%fj@\${SPACK_FJVER}

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:\$LD_LIBRARY_PATH

EOF

  else

    if [ -z "$SCALE_NETCDF_C" ] || [ -z "$SCALE_NETCDF_F" ] || [ -z "$SCALE_PNETCDF" ] || [ -z "$SCALE_HDF" ] ; then
      echo "[Error] Export SCALE environmental parameters (e.g., SCALE_NETCDF_C)"
      exit 1
    fi

cat << EOF >>  $jobscrp 
export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:\$LD_LIBRARY_PATH

EOF

  fi

cat << EOF >>  $jobscrp 
. ./src/func_util.sh
. ./config.main

for time in $time_list ; do 
    mpirunf - ./obsmake $TMP/config/obsmake_\${time}.conf log/obsmake || exit \$?
    cp $TMP/obsin/obsin.dat.out $OBS/${OBSNAME[1]}_\${time}.dat
done
EOF

  if (( USE_LLIO_BIN == 1 )); then
      echo "llio_transfer --purge obsmake" >> $jobscrp 
  fi

  echo "[$(datetime_now)] Run obsmake job on PJM"
  echo
  
  job_submit_PJM $jobscrp
  echo
  
  job_end_check_PJM $jobid
  res=$?

# qsub
elif [ "$PRESET" = 'Linux_torque' ]; then

if [ $NNODES_USE -lt 4 ] ; then
  RSCGRP=s
elif [ $NNODES_USE -le 16 ] ; then
  RSCGRP=m
elif [ $NNODES_USE -le 24 ] ; then
  RSCGRP=l
else
  echo "too many nodes required. " $NNODES_USE " > 16"
  exit 1
fi

cat > $jobscrp << EOF
#!/bin/sh
#PBS -N exec_obsmake
#PBS -q ${RSCGRP}
#PBS -l nodes=${NNODES_USE}:ppn=${PPN}
#PBS -l walltime=${TIME_LIMIT}
#
#

cd \${PBS_O_WORKDIR}

export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

EOF

if [ "$SCALE_SYS" == "Linux64-gnu-ompi" ] ; then

cat >> $jobscrp << EOF

source /etc/profile.d/modules.sh
module unload mpt/2.12
module unload intelcompiler/16.0.1.150
module unload intelmpi/5.1.2.150
module unload hdf5/1.8.16-intel
module unload netcdf4/4.3.3.1-intel
module unload netcdf4/fortran-4.4.2-intel
module load gcc/4.7.2
module load openmpi/2.0.4-gcc
module load hdf5/1.8.16
module load netcdf4/4.3.3.1
module load netcdf4/fortran-4.4.2
module load lapack/3.6.0

export OMP_NUM_THREADS=1
export KMP_AFFINITY=compact

ulimit -s unlimited

cd $TMP

. ./src/func_util.sh
. ./config.main

for time in $time_list ; do 
    mpirun --mca btl openib,sm,self --bind-to core ./obsmake $TMP/config/obsmake_\${time}.conf log/obsmake || exit \$?
    cp $TMP/obsin/obsin.dat.out $OBS/${OBSNAME[1]}_\${time}.dat
done


EOF

else

cat >> $jobscrp << EOF

source /etc/profile.d/modules.sh 
module unload mpt/2.12
module load intelmpi/5.1.2.150

export OMP_NUM_THREADS=${THREADS}
export KMP_AFFINITY=compact

export LD_LIBRARY_PATH="/home/seiya/lib:$LD_LIBRARY_PATH"

ulimit -s unlimited

cd $TMP

. ./src/func_util.sh
. ./config.main
 
for time in $time_list ; do 
    mpirunf - ./obsmake $TMP/config/obsmake_\${time}.conf log/obsmake || exit \$?
    cp $TMP/obsin/obsin.dat.out $OBS/${OBSNAME[1]}_\${time}.dat
done

EOF

fi 

  echo "[$(datetime_now)] Run obsmake job on PJM"
  echo

  job_submit_torque $jobscrp
  echo
  
  job_end_check_torque $jobid
  res=$?
# direct
elif [ "$PRESET" = 'Linux' ]; then

  echo "[$(datetime_now)] Run ${job} job on PJM"
  echo

  cd $TMPS
  for time in $(echo $time_list) ; do 
    mpirunf - ./obsmake $TMPS/config/obsmake_${time}.conf log/obsmake || exit $?
    cp $TMPS/obsin/obsin.dat.out $OBS/${OBSNAME[1]}_${time}.dat
  done

else
  echo "PRESET '$PRESET' is not supported."
  exit 1
fi

#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)"

stage_out server || exit $?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $TMP $jobid ${job}_job.sh 'o e'

archive_log


#===============================================================================

echo "[$(datetime_now)] Finish $myname $@"

exit $res
