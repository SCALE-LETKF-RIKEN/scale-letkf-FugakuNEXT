#!/bin/bash
#===============================================================================
#
#  Common functions for 'cycle/fcst' jobs
#
#===============================================================================

staging_list_common_static () {
#-------------------------------------------------------------------------------
# Usage: staging_list_common_static JOBTYPE
#
#   JOBTYPE  Job type (cycle/fcst)
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBTYPE="$1"

#-------------------------------------------------------------------------------
# common variables

for d in $(seq $DOMNUM); do
  if ((PNETCDF == 1)); then
    mem_np_[$d]=1
  else
    mem_np_[$d]=${SCALE_NP[$d]}
  fi
done

if [ "$TOPO_FORMAT" != 'prep' ]; then
  mkdir -p $TMP/dat/topo/${TOPO_FORMAT}
  #ln -sf ${DATADIR}/topo/${TOPO_FORMAT}/Products $TMP/dat/topo/${TOPO_FORMAT}/Products
  echo "${DATADIR}/topo/${TOPO_FORMAT}/Products/|dat/topo/${TOPO_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  mkdir -p $TMP/dat/landuse/${LANDUSE_FORMAT}
  #ln -sf ${DATADIR}/landuse/${LANDUSE_FORMAT}/Products $TMP/dat/landuse/${LANDUSE_FORMAT}/Products
  echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products/|dat/landuse/${LANDUSE_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi

#-------------------------------------------------------------------------------
# observations
      

if [ "$JOBTYPE" = 'cycle' ]; then
  mkdir -p $TMP/obs

  time=$(datetime $STIME $LCYCLE s)
  while ((time <= $(datetime $ETIME $LCYCLE s))); do
    for iobs in $(seq $OBSNUM); do
#      if [ "${OBSNAME[$iobs]}" != '' ] && [ -e ${OBS}/${OBSNAME[$iobs]}_${time}.dat ]; then
      if [ "${OBSNAME[$iobs]}" != '' ] ; then
        mkdir -p ${OBS}
        echo "${OBS}/${OBSNAME[$iobs]}_${time}.dat|obs/${OBSNAME[$iobs]}_${time}.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
#        ln -sf ${OBS}/${OBSNAME[$iobs]}_${time}.dat $TMP/obs/${OBSNAME[$iobs]}_${time}.dat
      fi
    done
    time=$(datetime $time $LCYCLE s)
  done
fi

#-------------------------------------------------------------------------------
# create empty directories

#cat >> ${STAGING_DIR}/${STGINLIST} << EOF
#|mean/
#|log/
#EOF
#
#if [ "$JOBTYPE" = 'cycle' ]; then
#  cat >> ${STAGING_DIR}/${STGINLIST} << EOF
#|sprd/
#EOF
#fi

#-------------------------------------------------------------------------------
# time-invariant outputs

#-------------------
# stage-in
#-------------------

# domain catalogue
#-------------------
#if ((BDY_FORMAT == 1)); then
#  if [ -s "$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt" ]; then
#    pathin="$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt"
#    path="latlon_domain_catalogue.bdy.txt"
#    echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#  else
#    echo "[Error] Cannot find a lat/lon domain catalogue file at" >&2
#    echo "        '$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt'" >&2
#    exit 1
#  fi
#fi

#-------------------
# stage-out
#-------------------

# domain catalogue
#-------------------
if ((LOG_OPT <= 3 && BDY_FORMAT == 1)); then
  for d in $(seq $DOMNUM); do
    path="latlon_domain_catalogue.d$(printf $DOMAIN_FMT $d).txt"
    pathout="${OUTDIR[$d]}/const/log/latlon_domain_catalogue.txt"
    echo "${pathout}|${path}|1" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$((${SCALE_NP_S[$d]}+1))]}
  done
fi

#-------------------------------------------------------------------------------
}

staging_list_core_bdyorg () {

### Local variables
m=$1
q=$2

if [ "$job" == "cycle" ];then
  FLEN=$CYCLEFLEN
elif [ "$job" == "fcst" ];then
  FLEN=$FCSTLEN
else
  echo '$job = '$job' not supported.'
  exit 1
fi

time=$STIME
loop=0
time_bdy_start_prev=0

while ((time <= ETIME)); do
  loop=$((loop+1))

  # bdy (parent)
  #-------------------
  if ((BDY_FORMAT > 0 && q == 1)) ; then
    bdy_setting $time $FLEN $BDYCYCLE_INT "$BDYINT" "$PARENT_REF_TIME" "$BDY_SINGLE_FILE"
    if ((BDY_ROTATING == 1)) || [ ${bdy_times[1]} != $time_bdy_start_prev ] ; then
      time_bdy_start_prev=${bdy_times[1]}
      nbdy_max=0
    fi
    if ((nbdy > nbdy_max)); then
      for ibdy in $(seq $((nbdy_max+1)) $nbdy); do
        time_bdy=${bdy_times[$ibdy]}

        if ((BDY_FORMAT == 1)); then

          if ((BDY_ENS == 1)); then
              if ((m == mmean)); then
                mem_bdy="$BDY_MEAN"
              elif (( m == mmgue && EFSO_RUN == 1 )); then
                mem_bdy="$BDY_MEAN"
              else
                mem_bdy="${name_m[$m]}"
              fi
              for qb in $(seq $mem_np_bdy_); do
                pathin="${DATA_BDY_SCALE}/${time_bdy}/${BDY_SCALE_DIR}/${mem_bdy}${CONNECTOR_BDY}history$(scale_filename_bdy_sfx $((qb-1)))"
                path="${name_m[$m]}/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))$(scale_filename_bdy_sfx $((qb-1)))"
                if ((DISK_MODE_BDYDATA >= 1)); then
                  echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
                else
                  ln -sf $pathin $TMP/$path
                fi
              done
          else
            if ((m == 1)); then
              for qb in $(seq $mem_np_bdy_); do
                pathin="${DATA_BDY_SCALE}/${time_bdy}/${BDY_SCALE_DIR}/${BDY_MEAN}${CONNECTOR_BDY}history$(scale_filename_bdy_sfx $((qb-1)))"
                path="mean/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))$(scale_filename_bdy_sfx $((qb-1)))"
                if ((DISK_MODE_BDYDATA >= 1)); then
                  echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
                else
                  ln -sf $pathin $TMP/$path
                fi
              done
            fi
          fi

        elif ((BDY_FORMAT == 2 || BDY_FORMAT == 4)); then
          if ((BDY_FORMAT == 2)); then
            data_bdy_i="$DATA_BDY_WRF"
            filenum=1
            filename_prefix[1]='wrfout_'
            filename_suffix[1]=''
            filenamein_prefix[1]=''
            filenamein_suffix[1]=''
          elif ((BDY_FORMAT == 4)); then
            data_bdy_i="$DATA_BDY_GRADS"
            filenum=3
            filename_prefix[1]='atm_'
            filename_suffix[1]='.grd'
            filenamein_prefix[1]='atm_'
            filenamein_suffix[1]='.grd'
            filename_prefix[2]='sfc_'
            filename_suffix[2]='.grd'
            filenamein_prefix[2]='sfc_'
            filenamein_suffix[2]='.grd'
            filename_prefix[3]='land_'
            filename_suffix[3]='.grd'
            filenamein_prefix[3]='lnd_'
            filenamein_suffix[3]='.grd'
          fi

          if ((BDY_ENS == 1)); then
            if ((m == mmean)); then
              mem_bdy="$BDY_MEAN"
            else
              mem_bdy="${name_m[$m]}"
            fi
            for ifile in $(seq $filenum); do
              if ((BDY_ROTATING == 1)); then
                pathin="${data_bdy_i}/${time}/${mem_bdy}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
              else
                pathin="${data_bdy_i}/${mem_bdy}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
              fi
              path="${name_m[$m]}/bdyorg_${filenamein_prefix[$ifile]}$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))${filenamein_suffix[$ifile]}"
              if ((DISK_MODE_BDYDATA >= 1)); then
                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
              else
                ln -sf  $pathin $TMP/$path
              fi
            done
          else
            if ((m==mmean || (job == "fcst" && m==1) ));then
              for ifile in $(seq $filenum); do
                if ((BDY_ROTATING == 1)); then
                  pathin="${data_bdy_i}/${time}/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                else
                  pathin="${data_bdy_i}/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                fi
                path="mean/bdyorg_${filenamein_prefix[$ifile]}$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))${filenamein_suffix[$ifile]}"
                if ((DISK_MODE_BDYDATA >= 1)); then
                  echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
                else
                  ln -sf  $pathin $TMP/$path
                fi
              done
            fi
          fi

        fi # [ BDY_FORMAT == 2 || BDY_FORMAT == 4 ]
      done # [ ibdy in $(seq $((nbdy_max+1)) $nbdy) ]
      nbdy_max=$nbdy
    fi
  fi # [ BDY_FORMAT > 0 ]

  time=$(datetime $time $LCYCLE s)
done # [ ((time <= ETIME)) ]

}

#===============================================================================

config_file_scale_launcher () {
#-------------------------------------------------------------------------------
# Generate the launcher configuration files for scale_pp/scale_init/scale
#
# Usage: config_file_scale_launcher MODEL_NAME CONF_NAME
#
#   JOBTYPE     Job type (cycle/fcst)
#   MODEL_NAME  (scale-rm_pp/scale-rm_init/scale-rm)
#   CONF_NAME   (pp/init/run)
#   MEMBER_RUN  Number of members needed to run
#
# Other input variables:
#   $nitmax
#   $time
#   $SCRP_DIR
#   $MEMBER
#   $mtot
#   $PPN_APPAR
#   $mem_nodes
#   $DOMNUM
#   $PRC_DOMAINS_LIST
#   $STAGING_DIR
#-------------------------------------------------------------------------------

local JOBTYPE="$1"; shift
local MODEL_NAME="$1"; shift
local CONF_NAME="$1"; shift
local MEMBER_RUN="$1"

#-------------------------------------------------------------------------------

local it
local conf_file

if [ "$JOBTYPE" = 'cycle' ]; then
  CONF_FILES_SEQNUM='.false.'
else
  CONF_FILES_SEQNUM='.true.'
fi

DET_RUN_TF='.false.'
if ((DET_RUN == 1)); then
  DET_RUN_TF='.true.'
fi

EFSO_RUN_TF='.false.'
if (( EFSO_RUN == 1 )); then
  EFSO_RUN_TF='.true.'
fi

conf_file="${MODEL_NAME}_${time}.conf"

echo "  $conf_file"
cat $SCRP_DIR/config.nml.ensmodel | \
    sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
        -e "/!--MEMBER_RUN--/a MEMBER_RUN = $MEMBER_RUN," \
        -e "/!--CONF_FILES--/a CONF_FILES = \"${CONF_NAME}.d<domain>_${time}.conf\"," \
        -e "/!--CONF_FILES_SEQNUM--/a CONF_FILES_SEQNUM = $CONF_FILES_SEQNUM," \
        -e "/!--DET_RUN--/a DET_RUN = $DET_RUN_TF," \
        -e "/!--EFSO_RUN--/a EFSO_RUN = $EFSO_RUN_TF," \
        -e "/!--PPN--/a PPN = $PPN_APPAR," \
        -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
        -e "/!--NUM_DOMAIN--/a NUM_DOMAIN = $DOMNUM," \
        -e "/!--PRC_DOMAINS--/a PRC_DOMAINS = $PRC_DOMAINS_LIST" \
    > $TMP/config/${conf_file}

#-------------------------------------------------------------------------------
}

#===============================================================================

config_file_save () {
#-------------------------------------------------------------------------------
# Save the runtime configuration files in $OUTDIR
#
# Usage: config_file_save [CONFIG_DIR]
#
#   CONFIG_DIR  Temporary directory of configuration files
#               '-': Use $TMPROOT
#
# Other input variables:
#   $TMPROOT
#-------------------------------------------------------------------------------

local CONFIG_DIR="${1:--}"

if [ "$CONFIG_DIR" = '-' ]; then
  CONFIG_DIR="$TMPROOT"
fi

#-------------------------------------------------------------------------------

mkdir -p ${OUTDIR[1]}/config

cp -fr $CONFIG_DIR/*.conf ${OUTDIR[1]}/config/

#-------------------------------------------------------------------------------
}

#===============================================================================
