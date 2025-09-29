#!/bin/bash

#
# SNO is executed by a bulk job
#

GRADS=F
PLEV=F

INPUT_FROM_SNOW=F
INPUT_SNOW_NP=16

ALLVAR=F
SINGLE_VAR=T

tint=21600 # [second]
tstart='2018-07-05 0:10:00'
tend=$tstart

. ./config.main
RUNDIR="${TMP}_sno"


SCALEDIR="$(cd "$(pwd)/../../.." && pwd)"  

TYPE=fcst
TYPE=anal
#TYPE=gues
TYPE=hist

## Which domain do you want to convert?
#DOM=2 

# Output file (X & Y process number) for each member
NP_OFILE_X=1
NP_OFILE_Y=1

if [ "$INPUT_FROM_SNOW" == "T" ] ; then
  NP_OFILE_X=1
  NP_OFILE_Y=1
fi

# Do not edit!
NP_OFILE=$((${NP_OFILE_X} * ${NP_OFILE_Y})) # Output file (process number) for each member

# Specify members that will be processed
#SNO_MEMBERS=50
#SNO_MEM_L="mean "$(seq -f %04g ${SNO_MEMBERS})

SNO_MEM_L="mean "


if [ "$ALLVAR" == "T" ] ; then
  # Convert all variables
  VARS=""
elif [ "$SINGLE_VAR" == "T" ] ; then
  # Convert single variable 
  NVAR='MSLP'
  VARS="$NVAR"
else
  # Specify variables that will be converted
  VARS='"W", "QV"' 
fi

TOPO=0 # Process topography file? # 1: Yes, 0: No
if (( TOPO > 0 )) ; then
  VARS='"topo"'
  SNO_MEM_L="mean"
fi


if [ "$GRADS" == "T" ] ; then
  NP_OFILE_X=1
  NP_OFILE_Y=1
  OUTPUT_GRADS=".true."
  OUTPUT_GRADSCTL=".true."
else
  OUTPUT_GRADS=".false."
  OUTPUT_GRADSCTL=".false."
fi



###############################

# Path for SNO binary
SNOBIN_ORG=${SCALEDIR}/bin/sno
SNOBIN=${RUNDIR}/sno
if [ ! -e ${SNOBIN_ORG} ] ; then
  echo "No SNO binary!"
  exit
fi

###############################


rm -rf ${RUNDIR}

mkdir -p ${RUNDIR}/conf
mkdir -p ${RUNDIR}/log

# copy binary 
cp ${SNOBIN_ORG} ${SNOBIN}

conf_bulk="${RUNDIR}/conf/bulk_sno.conf"

cnt=0

time="$tstart"
while (($(date -ud "$time" '+%s') <= $(date -ud "$tend" '+%s'))); do # time loop

  timef=$(date -ud "$time" '+%Y-%m-%d %H:%M:%S')

  YYYYs=$(date -ud "$timef" '+%Y')
  MMs=$(date -ud "$timef" '+%m')
  DDs=$(date -ud "$timef" '+%d')
  HHs=$(date -ud "$timef" '+%H')
  MNs=$(date -ud "$timef" '+%M')
  SSs=$(date -ud "$timef" '+%S')

  DTIME=${YYYYs}${MMs}${DDs}${HHs}${MNs}${SSs}
  echo $DTIME

  SCALE_TIME=$(date -ud "$timef" '+%Y%m%d-%H%M%S.000')

  for mem in  ${SNO_MEM_L} # member loop
  do 
    cnt=$((${cnt} + 1))
    echo $mem
  
    SNO_BASENAME_IN="${OUTDIR}/${DTIME}/${TYPE}/${mem}/history"
  
    SNO_BASENAME_OUT="history"
  
    if [ "$TYPE" != "fcst" ] && [ "$TYPE" != "hist" ]  ; then
      SNO_BASENAME_OUT="$TYPE"
      SNO_BASENAME_IN="${OUTDIR}/${DTIME}/${TYPE}/${mem}/init_${SCALE_TIME}"
    fi
  
    if [ "$ALLVAR" != "T" ] && [ "$SINGLE_VAR" == "T" ]; then
      SNO_BASENAME_OUT="$NVAR"
    fi
    
  
    if [ "$INPUT_FROM_SNOW" == "T" ] ; then
      SNO_BASENAME_IN=${OUTDIR}/${DTIME}/${TYPE}_sno_np$(printf %05d ${INPUT_SNOW_NP})/${mem}/${SNO_BASENAME_OUT}
    fi
  
    if [ "$GRADS" == "T" ] ; then
      SNO_OUTPUT_PATH=${OUTDIR}/${DTIME}/${TYPE}_sno_grads/${mem}
    else
      SNO_OUTPUT_PATH=${OUTDIR}/${DTIME}/${TYPE}_sno_np$(printf %05d ${NP_OFILE})/${mem}
    fi
  
    if (( TOPO > 0 )) ; then
      SNO_OUTPUT_PATH=${OUTDIR}/const/topo_sno_np$(printf %05d ${NP_OFILE})
      SNO_BASENAME_IN="${OUTDIR}/const/topo/topo"
      SNO_BASENAME_OUT="topo"
    fi
  
    if [ ! -e ${SNO_OUTPUT_PATH} ] ; then
      mkdir -p ${SNO_OUTPUT_PATH}
    fi
  
  
    conf="${RUNDIR}/conf/sno_${mem}_${DTIME}.conf"
  
cat << EOF >> $conf
&PARAM_IO
 IO_LOG_BASENAME = "log/LOG_${mem}_${DTIME}",
 IO_LOG_ALLNODE = .false.,
 IO_LOG_SUPPRESS = .true.,
 IO_LOG_NML_SUPPRESS = .true.,
/
&PARAM_SNO
 basename_in  = "${SNO_BASENAME_IN}",
 basename_out  = "${SNO_BASENAME_OUT}",
 dirpath_out = "${SNO_OUTPUT_PATH}",
 vars         = ${VARS},
 nprocs_x_out = ${NP_OFILE_X},
 nprocs_y_out = ${NP_OFILE_Y},
 output_gradsctl = ${OUTPUT_GRADSCTL},
 output_grads = ${OUTPUT_GRADS},
 debug = .true.,
/
EOF

    if [ "$PLEV" == "T" ] &&  [ "$INPUT_FROM_SNOW" != "T" ] ; then
  
cat << EOF >> $conf
&PARAM_SNOPLGIN_VGRIDOPE
 SNOPLGIN_vgridope_type        = 'PLEV', 
 SNOPLGIN_vgridope_lev_num     = 15,
 SNOPLGIN_vgridope_lev_data    = 1000.e+2, 950.e+2, 925.e+2, 900.e+2, 850.e+2, 800.e+2, 700.e+2, 600.e+2, 500.e+2, 400.e+2, 300.e+2, 200.e+2, 100.e+2, 70.e+2, 50.e+2, 
/
EOF

    fi
  
    ln -s ${conf} ${conf_bulk}.${cnt}

  done # member loop

  time=$(date -ud "${tint} second $time" '+%Y-%m-%d %H:%M:%S')
done # time loop


# Total SNO processes  
echo "Total number of nodes: $(( NP_OFILE * cnt / PPN ))"

# Get total SNO NODEs
if (( NP_OFILE < PPN )) ; then
  SNO_NODE=1
else
  SNO_NODE=$(( NP_OFILE / PPN ))
  if (( SNO_NODE*PPN < NP_OFILE )) ; then
    SNO_NODE=$(( SNO_NODE + 1 ))
  fi
fi

jobsh="${RUNDIR}/job_sno.sh"

if [ "$PRESET" = 'FUGAKU' ]; then

#  if (( SNO_NODE < 12 )) ; then
#    SNO_NODE=12
#  fi

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

cat << EOF >> $jobsh
#!/bin/sh 
#
#PJM -g ${GROUP} 
#PJM -x PJM_LLIO_GFSCACHE=${VOLUMES}
#PJM -L "rscgrp=small"
#PJM -L "node=${SNO_NODE}"
#PJM -L "elapse=00:30:00"
#PJM --mpi "max-proc-per-node=${PPN}"
#PJM -j
#PJM -s
#
#
export PARALLEL=${THREADS}
export OMP_NUM_THREADS=${THREADS}
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD

EOF

  if (( USE_SPACK > 0 )); then
cat << EOF >> $jobsh
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

cat << EOF >>  $jobsh

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:\$LD_LIBRARY_PATH

EOF

  fi

  if (( USE_LLIO_BIN == 1 )); then
    echo "llio_transfer ${SNOBIN}" >> $jobsh
    echo "" >> $jobsh
  fi

cat << EOF >> $jobsh
echo "[\$(date "+%Y/%m/%d %H:%M:%S")] Start SNO"

mpiexec -std-proc log/NOUT -n $((NP_OFILE)) ${SNOBIN} ${conf_bulk}.\${PJM_BULKNUM}

echo "[\$(date "+%Y/%m/%d %H:%M:%S")] End SNO"

EOF

  if (( USE_LLIO_BIN == 1 )); then
    echo "llio_transfer --purge ${SNOBIN}" >> $jobsh
    echo "" >> $jobsh
  fi

# qsub
elif [ "$PRESET" = 'Linux_torque' ]; then

if [ $SNO_NODE -lt 4 ] ; then
  RSCGRP=s
elif [ $SNO_NODE -le 16 ] ; then
  RSCGRP=m
elif [ $SNO_NODE -le 24 ] ; then
  RSCGRP=l
else
  echo "too many nodes required. " $SNO_NODE " > 24"
  exit 1
fi

cat > $jobscrp << EOF
#!/bin/sh
#PBS -N $job
#PBS -q $RSCGRP
#PBS -l nodes=${SNO_NODE}:ppn=${PPN}
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

EOF

else

cat >> $jobscrp << EOF

source /etc/profile.d/modules.sh 
module unload mpt/2.12
module load intelmpi/5.1.2.150

export OMP_NUM_THREADS=${THREADS}
export KMP_AFFINITY=compact

export LD_LIBRARY_PATH="/home/seiya/lib:$LD_LIBRARY_PATH"
EOF


fi 

cat >> $jobscrp << EOF
ulimit -s unlimited

echo "[\$(date "+%Y/%m/%d %H:%M:%S")] Start SNO"

mpirun --mca btl openib,sm,self --bind-to core ${SNOBIN} ${conf_bulk}.\${PJM_BULKNUM} log/NOUT || exit \$?

echo "[\$(date "+%Y/%m/%d %H:%M:%S")] End SNO"

EOF

  echo "[$(datetime_now)] Run ${job} job on PJM"
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

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" &> run_progress || exit $?

else
  echo "PRESET '$PRESET' is not supported."
  exit 1
fi

cd ${RUNDIR}
pjsub --bulk --sparam 1-${cnt} job_sno.sh 
cd - > /dev/null

