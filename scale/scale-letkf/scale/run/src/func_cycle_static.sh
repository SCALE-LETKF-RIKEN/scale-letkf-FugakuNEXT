#!/bin/bash
#===============================================================================
#
#  Functions for 'cycle' jobs
#
#===============================================================================

staging_list_static () {
#-------------------------------------------------------------------------------
# Prepare all the staging list files
#
# Usage: staging_list_static
#-------------------------------------------------------------------------------

declare -a mem_np_

if ((PNETCDF == 1)); then
  local mem_np_=1
else
  local mem_np_=$mem_np
fi
if ((PNETCDF_BDY_SCALE == 1)); then
  local mem_np_bdy_=1
else
  local mem_np_bdy_=$((DATA_BDY_SCALE_PRC_NUM_X*DATA_BDY_SCALE_PRC_NUM_Y))
  if (( mem_np_bdy_ < 1 )) && (( BDY_FORMAT < 4 ))  && (( BDY_FORMAT > 0 )); then
    echo "[Error] $0: Specify DATA_BDY_SCALE_PRC_NUM_X/Y" >&2
    exit 1
  fi
fi

#-------------------------------------------------------------------------------
# common section

if ((DISK_MODE == 1)) ;then
  mkdir -p ${TMPROOT}/dat
  mkdir -p ${TMPROOT}/log
  staging_list_common_static cycle
fi

mtot=$(( MEMBER + 1 ))
mmean=$(( MEMBER + 1 ))
name_m[$mmean]='mean'
mkdir -p $TMP/mean
if (( DET_RUN == 1 )); then
  mtot=$(( mtot + 1 ))
  mmdet=$(( mmean + 1 ))
  name_m[$mmdet]='mdet'
  mkdir -p $TMP/mdet
else
  mmdet=$mmean
fi

if (( EFSO_RUN == 1 )); then
  mtot=$(( mtot + 1 ))
  mmgue=$(( mmdet + 1 ))
  name_m[$mmgue]='mgue'
  mkdir -p $TMP/mgue
fi


for m in $(seq $MEMBER); do
  name_m[$m]=$(printf $MEMBER_FMT $m)
  mkdir -p $TMP/${name_m[$m]}
done

totalnp=$((PPN*NNODES))
SCALE_NP_TOTAL=0
for d in `seq $DOMNUM`; do
  SCALE_NP_TOTAL=$((SCALE_NP_TOTAL+SCALE_NP[$d]))
done

repeat_mems=$((mtot*SCALE_NP_TOTAL/totalnp))
nitmax=$(( ( mtot - 1) * SCALE_NP_TOTAL / totalnp + 1 ))

if ((DISK_MODE >= 1)) ;then

#-------------------------------------------------------------------------------
# executable files

  # no link
cat >> ${STAGING_DIR}/${STGINLIST_SHARE} << EOF
${COMMON_DIR}/pdbash|pdbash
${COMMON_DIR}/datetime|datetime
${ENSMODEL_DIR}/scale-rm_pp_ens|scale-rm_pp_ens
${ENSMODEL_DIR}/scale-rm_init_ens|scale-rm_init_ens
${ENSMODEL_DIR}/scale-rm_ens|scale-rm_ens
EOF
 
if (( OBSOPE_RUN == 1 )) ;then
  echo "${OBSUTIL_DIR}/obsope|obsope" >> ${STAGING_DIR}/${STGINLIST_SHARE}
fi
if (( PAWR_DECODE == 1 )) ;then
  echo "${OBSUTIL_DIR}/dec_pawr|dec_pawr" >> ${STAGING_DIR}/${STGINLIST_SHARE}
fi
  echo "${LETKF_DIR}/letkf|letkf" >> ${STAGING_DIR}/${STGINLIST_SHARE} 

#-------------------------------------------------------------------------------
# database

cat >> ${STAGING_DIR}/${STGINLIST_CONSTDB} << EOF
${SCALEDIR}/data/rad|dat/rad
${SCALEDIR}/data/land|dat/land
${SCALEDIR}/data/urban|dat/urban
${SCALEDIR}/data/lightning|dat/lightning
EOF

if [ "${SOUNDING}" != "" ] ; then
  echo "${SOUNDING}|dat/${SOUNDING}" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi

else # DISK_MODE=0 : no staging 

#-------------------------------------------------------------------------------
# executable files

cp ${COMMON_DIR}/pdbash ${TMPROOT}/pdbash
cp ${COMMON_DIR}/datetime ${TMPROOT}/datetime
cp ${ENSMODEL_DIR}/scale-rm_pp_ens ${TMPROOT}/scale-rm_pp_ens
cp ${ENSMODEL_DIR}/scale-rm_init_ens ${TMPROOT}/scale-rm_init_ens
cp ${ENSMODEL_DIR}/scale-rm_ens ${TMPROOT}/scale-rm_ens
 
cp ${OBSUTIL_DIR}/obsope ${TMPROOT}/obsope
if (( PAWR_DECODE == 1 )) ;then
  cp ${OBSUTIL_DIR}/dec_pawr ${TMPROOT}/dec_pawr
fi
cp ${LETKF_DIR}/letkf ${TMPROOT}/letkf
cp ${LETKF_DIR}/efso ${TMPROOT}/efso

#-------------------------------------------------------------------------------
# database

mkdir -p ${TMPROOT}/dat
cp -r ${SCALEDIR}/data/rad ${TMPROOT}/dat/rad
cp -r ${SCALEDIR}/data/land ${TMPROOT}/dat/land
cp -r ${SCALEDIR}/data/urban ${TMPROOT}/dat/urban
cp -r ${SCALEDIR}/data/lightning ${TMPROOT}/dat/lightning

if [ "${SOUNDING}" != "" ] ; then
  cp ${SOUNDING} ${TMPROOT}/dat/
fi

fi ### DISK_MODE >= 1

#-------------------------------------------------------------------------------

ith=0
for m in $(seq $mtot) ; do
  for q in $(seq ${mem_np_}); do
      if ((DISK_MODE >= 1));then
        ith=$((ith+1))
        staging_list_core $m $q &
      else
        if ((BDY_FORMAT > 0 && q == 1)) ; then
          ith=$((ith+1))
          staging_list_core_bdyorg $m $q &
        fi
      fi
      if (( ith == SHELL_PROCS )) ; then 
         wait 
         ith=0
      fi
  done
done
wait

}

staging_list_core () {

### tentative 
if ((DOMNUM > 1)); then
  echo "not supported."
  exit 0 
fi

d=1
dom=".d$(printf $DOMAIN_FMT $d)" ###

### Local variables
m=$1
q=$2

### bdyorg
if ((BDY_FORMAT > 0 && q == 1)) ; then
  staging_list_core_bdyorg $m $q 
fi
 
time=$STIME

atime=$(datetime $time $LCYCLE s)
loop=0
while ((time <= ETIME)); do
  loop=$((loop+1))

  time=$(datetime $time $((LCYCLE * (loop-1) )) s)
  atime=$(datetime $time $((LCYCLE * loop )) s)

  sfx=$(scale_filename_sfx $((q-1)))
  tsfx="_"$(datetime_scale $time)$(scale_filename_sfx $((q-1)))
  atsfx="_"$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))

  sproc=$(((m-1)*mem_np+q))
  snode=${mem2node[${sproc}]}
  msnode=${mem2node[$(((mmean-1)*mem_np+q))]}


  #-------------------
  # stage-in
  #-------------------

  mkdir -p ${TMPROOT}/${name_m[$m]}

  # anal
  #-------------------
  if ((loop == 1 && MAKEINIT != 1)); then
      mkdir -p ${TMPROOT}/${name_m[$m]}
      pathin="${INDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init${tsfx}"
      path="${name_m[$m]}/anal${dom}${tsfx}"
      echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${snode}
  fi

  # topo
  #-------------------
  if ((loop == 1)) && [ "$TOPO_FORMAT" == 'prep' ] && [ "$TOPO_FORMAT" == 'none' ]; then
    if ((DISK_MODE == 3)); then
            pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo${sfx}"
            path="topo/topo${dom}${sfx}"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${snode}
    else
      if ((m == 1)); then
          pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo${sfx}"
          path="topo/topo${dom}${sfx}"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
      fi
    fi
  fi

    # topo (bdy_scale)
    #-------------------
    if ((loop == 1 && BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ] && [ "$TOPO_FORMAT" != 'none' ]  ; then
      if ((q == 1)) && ((m == 1)); then
        echo "|bdytopo/" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
        for qb in $(seq $mem_np_bdy_);do
          pathin="${DATA_TOPO_BDY_SCALE}/const/topo/topo$(scale_filename_bdy_sfx $((qb-1)))"
          path="bdytopo/bdytopo$(scale_filename_bdy_sfx $((qb-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
        done
        pathin="${DATA_TOPO_BDY_SCALE}/const/log/latlon_domain_catalogue.txt"
        path="bdytopo/latlon_domain_catalogue.txt"
        echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
      fi
    fi

  # landuse
  #-------------------
  if ((loop == 1)) && [ "$LANDUSE_FORMAT" = 'prep' ] && [ "$LANDUSE_FORMAT" = 'none' ]; then
    if ((DISK_MODE == 3)); then
            pathin="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse${sfx}"
            path="landuse/landuse${dom}${sfx}"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${snode}
    else
      if ((m == 1)); then
          pathin="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse${sfx}"
          path="landuse/landuse${dom}${sfx}"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
      fi
    fi
  fi

  # bdy (prepared)
  #-------------------
  if ((BDY_FORMAT == 0)); then
    if ((BDY_ENS == 0)); then
      if ((DISK_MODE == 3)); then
            pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary${sfx}"
            path="mean/bdy${tsfx}"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${snode}
          if ((USE_INIT_FROM_BDY == 1)); then
                pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy${sfx}"
                path="mean/init${dom}${tsfx}"
                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${snode}
          fi
      else
          pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary${sfx}"
          path="mean/bdy${tsfx}"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
        if ((USE_INIT_FROM_BDY == 1)); then
              pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy${sfx}"
              path="mean/init${dom}${tsfx}"
              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
        fi
      fi
    elif ((BDY_ENS == 1)); then
          pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary${sfx}"
          path="${name_m[$m]}/bdy${tsfx}"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${snode}
        if ((USE_INIT_FROM_BDY == 1)); then
              pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy${sfx}"
              path="${name_m[$m]}/init${dom}${tsfx}"
              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${snode}
        fi
    fi
  fi

  # additive inflation
  #-------------------
  if ((loop == 1 && ADDINFL == 1)); then
          pathin="${DATA_ADDINFL[$d]}/const/addi/${name_m[$m]}${CONNECTOR}init${sfx}"
          path="${name_m[$m]}/addi${dom}${sfx}"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${snode}
  fi

  #-------------------
  # stage-out
  #-------------------

  # anal (initial time)
  #-------------------
  if ((loop == 1 && MAKEINIT == 1)); then
          path="${name_m[$m]}/init${dom}${tsfx}"
          pathout="${OUTDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init${sfx}"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${snode}
  fi

  # topo
  #-------------------
  if ((loop == 1 )) && [ "$TOPO_FORMAT" != 'prep' ] && [ "$TOPO_FORMAT" != 'none' ]; then
    echo "|topo/" >> ${STAGING_DIR}/${STGINLIST}
    if ((m == 1)) && ((TOPOOUT_OPT <= 1)); then
      path="topo/topo${dom}${sfx}"
      pathout="${OUTDIR[$d]}/const/${CONNECTOR_TOPO}topo${sfx}"
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
    fi
  fi

  # landuse
  #-------------------
  if ((loop == 1)) && [ "$LANDUSE_FORMAT" != 'prep' ] && [ "$LANDUSE_FORMAT" != 'none' ]; then
    echo "|landuse/" >> ${STAGING_DIR}/${STGINLIST}
    if ((m == 1)) && ((LANDUSEOUT_OPT <= 1)); then
      path="landuse/landuse${dom}${sfx}"
      pathout="${OUTDIR[$d]}/const/${CONNECTOR_LANDUSE}landuse${sfx}"
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
    fi
  fi
  # bdy
  #-------------------
  if ((BDY_FORMAT != 0)); then
    if ((BDY_ENS == 1 && BDYOUT_OPT <= 1)); then
          path="${name_m[$m]}/bdy${tsfx}"
          pathout="${OUTDIR[1]}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary${sfx}"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${snode}
        if ((USE_INIT_FROM_BDY == 1)); then
              path="${name_m[$m]}/init${dom}${tsfx}"
              pathout="${OUTDIR[$d]}/${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy${sfx}"
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${snode}
        fi
    elif ((BDYOUT_OPT <= 2)) && ((m==mmean)); then
        path="mean/bdy${tsfx}"
        pathout="${OUTDIR[1]}/${time}/bdy/mean${CONNECTOR}boundary${sfx}"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$q]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${msnode}
      if ((USE_INIT_FROM_BDY == 1)); then
            path="mean/init${dom}${tsfx}"
            pathout="${OUTDIR[$d]}/${time}/bdy/mean${CONNECTOR}init_bdy${sfx}"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${msnode}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${msnode}
      fi
    fi
  fi

  # anal
  #-------------------
  if ((  ( (OUT_OPT <= 4 || (OUT_OPT <= 5 && loop % OUT_CYCLE_SKIP == 0) || atime > ETIME)  && m <= mtot ) || \
         ( OUT_OPT <= 7 && (m == mmean || (DET_RUN==1 && m==mmdet) ) ) )) ; then 
        path="${name_m[$m]}/anal${dom}${atsfx}"
        pathout="${OUTDIR[$d]}/${atime}/anal/${name_m[$m]}${CONNECTOR}init${sfx}"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${snode}
        if ((m == mmean && SPRD_OUT == 1)); then
          path="sprd/anal${dom}${atsfx}"
          pathout="${OUTDIR[$d]}/${atime}/anal/sprd${CONNECTOR}init${sfx}"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${snode}
        fi
  fi

  # gues
  #-------------------

 if  (( OUT_OPT <= 3 &&  m <= mtot )) || (( OUT_OPT <= 6 && (m == mmean || (DET_RUN==1 && m==mmdet) ) )) ; then 
        path="${name_m[$m]}/gues${dom}${atsfx}"
        pathout="${OUTDIR[$d]}/${atime}/gues/${name_m[$m]}${CONNECTOR}init${sfx}"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${snode}
        if ((m == mmean && SPRD_OUT == 1)); then
          path="sprd/gues${dom}${atsfx}"
          pathout="${OUTDIR[$d]}/${atime}/gues/sprd${CONNECTOR}init${sfx}"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${snode}
        fi
 fi

  # hist
  #-------------------

  if (( OUT_OPT <= 1 &&  m <= mtot )) || (( OUT_OPT <= 2 && (m == mmean || (DET_RUN==1 && m==mmdet) ) )) ; then 
        path="${name_m[$m]}/hist${dom}${tsfx}"
        pathout="${OUTDIR[$d]}/${time}/hist/${name_m[$m]}${CONNECTOR}history${sfx}"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${snode}
  fi

#    # diag
#    #-------------------
#    if (( q == 1 && m == mmean )); then
#      if ((RTPS_INFL_OUT == 1)); then
#        path="rtpsinfl.d01_$(datetime_scale $atime).nc"
#        pathout="${OUTDIR}/${atime}/diag/rtpsinfl.init.nc"
##        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${msnode}
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${msnode}
#      fi
#      if ((NOBS_OUT == 1)); then
#        path="nobs.d01_$(datetime_scale $atime).nc"
#        pathout="${OUTDIR}/${atime}/diag/nobs.init.nc"
##        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${msnode}
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${msnode}
#      fi
#    fi
# 
#    # obsgues
#    #-------------------
#    if (( q == 1 && m <= mtot )); then
#      if ((OBSOUT_OPT <= 2)); then
#          path="${name_m[$m]}/obsgues.d01_${atime}.dat"
#          pathout="${OUTDIR}/${atime}/obsgues/${name_m[$m]}.obsda.dat"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#      fi
#    fi

  # log
  #-------------------

#  log_nfmt="-${PROCESS_FMT}"
#  m_out=0
#  m_init_out=0
#  p_out=1
#  if ((LOG_TYPE == 1 && m == 1)); then
#    m_out=1
#    m_init_out=1
#  else
#    if (( m <= mtot)) ;then
#      m_out=1
#    fi 
#    if ((BDY_ENS == 1 && m <= mtot)) || (( DISK_MODE <= 2 && m == 1 )) ;then
#      m_init_out=1
#    fi
#  fi

#  if ((LOG_TYPE != 1 )) || (( LOG_TYPE == 1 && sproc == 1)) ; then
#    p_out=1
#    p=$sproc
#  fi

#  if ((BDY_FORMAT != 0 && LOG_OPT <= 2)); then
#    if (( m_out == 1));then
#     if (( q == 1 )); then
#        path="log/scale_init.${name_m[$m]}${dom}.LOG_${time}${SCALE_SFX_NONC_0}"
#        pathout="${OUTDIR[$d]}/${time}/log/scale_init/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#     fi
#    fi
#    if (( p_out == 1));then
#      if ((nitmax == 1)); then
#        path="log/scale-rm_init_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
#        pathout="${OUTDIR[1]}/${time}/log/scale_init/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#      else
#        for it in $(seq $((BDY_ENS == 1 ? nitmax : 1))); do
#          path="log/scale-rm_init_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
#          pathout="${OUTDIR[1]}/${time}/log/scale_init/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#        done
#      fi
#    fi
#  fi
#  if ((LOG_OPT <= 3)); then
#    for m in $mlist; do
#      if (( q == 1 )); then
#        path="log/scale.${name_m[$m]}${dom}.LOG_${time}${SCALE_SFX_NONC_0}"
#        pathout="${OUTDIR[$d]}/${time}/log/scale/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#        path="log/scale.${name_m[$m]}${dom}.monitor_${time}${SCALE_SFX_NONC_0}"
#        pathout="${OUTDIR[$d]}/${time}/log/scale/${name_m[$m]}_monitor${SCALE_SFX_NONC_0}"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#      fi
#    done
#    for p in $plist; do
#      if ((nitmax == 1)); then
#        path="log/scale-rm_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
#        pathout="${OUTDIR[1]}/${time}/log/scale/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#      else
#        for it in $(seq $nitmax); do
#          path="log/scale-rm_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
#          pathout="${OUTDIR[1]}/${time}/log/scale/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#        done
#      fi
#    done
#  fi
#  if ((LOG_OPT <= 4)); then
#    for p in $plist; do
#      path="log/letkf.NOUT_${atime}$(printf -- "${log_nfmt}" $((p-1)))"
#      pathout="${OUTDIR[1]}/${atime}/log/letkf/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
#      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${snode}
#    done
#  fi


  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done # [ ((time <= ETIME)) ]


#-------------------------------------------------------------------------------
}
#===============================================================================

config_file_list () {
#-------------------------------------------------------------------------------
# Prepare all runtime configuration files
#
# Usage: config_file_list [CONFIG_DIR]
#
#   CONFIG_DIR  Temporary directory of configuration files to be staged to $TMPROOT
#               '-': Do not use a temporary directory and stage;
#                    output configuration files directly to $TMPROOT
#
# Other input variables:
#   $TMPROOT
#   $STAGING_DIR
#-------------------------------------------------------------------------------

local CONFIG_DIR="${1:--}"

local stage_config=1
if [ "$CONFIG_DIR" = '-' ]; then
  CONFIG_DIR="$TMPROOT"
  stage_config=0
fi

#-------------------------------------------------------------------------------

echo
echo "Generate configration files..."

mkdir -p $CONFIG_DIR

FILE_AGGREGATE=".false"
if ((PNETCDF == 1)); then
  FILE_AGGREGATE=".true."
fi

PRC_DOMAINS_LIST=
for d in $(seq $DOMNUM); do
  PRC_DOMAINS_LIST="$PRC_DOMAINS_LIST${SCALE_NP[$d]}, "
done

if [ "$TOPO_FORMAT" = "GTOPO30" ] || [ "$TOPO_FORMAT" = "DEM50M" ] || [ "$LANDUSE_FORMAT" = "GLCCv2" ] || [ "$LANDUSE_FORMAT" = "LU100M" ] ; then

  mkdir -p $OUTDIR/const/topo
  mkdir -p $OUTDIR/const/landuse
  time=$STIME
  config_file_scale_launcher cycle scale-rm_pp_ens "f<member>/pp" 1
  OFFLINE_PARENT_BASENAME=

  if ((BDY_FORMAT == 1)); then
    if ((DISK_MODE >= 1)); then
      BDYTOPO=${TMP}/bdytopo/bdytopo
    else
      BDYTOPO=${DATA_TOPO_BDY_SCALE}/const/topo
    fi
  fi

#  if ((BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#    OFFLINE_PARENT_BASENAME="$COPYTOPO"
#  fi

  if [ "$TOPO_FORMAT" != 'prep' ] && [ "$TOPO_FORMAT" != 'none' ]; then
    CONVERT_TOPO='.true.'
  else
    CONVERT_TOPO='.false.'
  fi
  
  if [ "$LANDUSE_FORMAT" != 'prep' ] && [ "$LANDUSE_FORMAT" != 'none' ]; then
    CONVERT_LANDUSE='.true.'
  else
    CONVERT_LANDUSE='.false.'
  fi

  # assume Domain 1

  mkdir -p $OUTDIR/$time/log/scale_pp

  if ((DISK_MODE >= 1));then
    TOPO_PATH=${TMP}
    LANDUSE_PATH=${TMP}
    SRC_TOPO_PATH=${TMP}/dat/topo
    SRC_LANDUSE_PATH=${TMP}/dat/landuse
  else
    TOPO_PATH="${DATA_TOPO}/const"
    LANDUSE_PATH="${DATA_LANDUSE}/const"
    SRC_TOPO_PATH=${DATADIR}/topo
    SRC_LANDUSE_PATH=${DATADIR}/landuse
  fi

  conf_file_src=${TMP}/config.nml.scale_pp
  conf="$(cat $conf_file_src | \
           sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"$OUTDIR/$time/log/scale_pp/LOG\"," \
               -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
               -e "/!--TOPOGRAPHY_OUT_BASENAME--/a TOPOGRAPHY_OUT_BASENAME = \"${TOPO_PATH}/topo/topo\"," \
               -e "/!--LANDUSE_OUT_BASENAME--/a LANDUSE_OUT_BASENAME = \"${LANDUSE_PATH}/landuse/landuse\"," \
               -e "/!--CONVERT_TOPO--/a CONVERT_TOPO = $CONVERT_TOPO," \
               -e "/!--CONVERT_LANDUSE--/a CONVERT_LANDUSE = $CONVERT_LANDUSE," \
               -e "/!--CNVTOPO_name--/a CNVTOPO_name = \"$TOPO_FORMAT\"," \
               -e "/!--GTOPO30_IN_DIR--/a GTOPO30_IN_DIR = \"${SRC_TOPO_PATH}/GTOPO30/Products\"," \
               -e "/!--DEM50M_IN_DIR--/a DEM50M_IN_DIR = \"${SRC_TOPO_PATH}/DEM50M/Products\"," \
               -e "/!--CNVLANDUSE_name--/a CNVLANDUSE_name = '$LANDUSE_FORMAT'," \
               -e "/!--GLCCv2_IN_DIR--/a GLCCv2_IN_DIR = \"${SRC_LANDUSE_PATH}/GLCCv2/Products\"," \
               -e "/!--LU100M_IN_DIR--/a LU100M_IN_DIR = \"${SRC_LANDUSE_PATH}/LU100M/Products\"," \
               -e "/!--COPYTOPO_IN_BASENAME--/a COPYTOPO_IN_BASENAME = \"${BDYTOPO}\"," \
               -e "/!--OFFLINE_PARENT_BASENAME--/a OFFLINE_PARENT_BASENAME = \"${OFFLINE_PARENT_BASENAME}\"," \
               -e "/!--OFFLINE_PARENT_PRC_NUM_X--/a OFFLINE_PARENT_PRC_NUM_X = ${DATA_BDY_SCALE_PRC_NUM_X}," \
               -e "/!--OFFLINE_PARENT_PRC_NUM_Y--/a OFFLINE_PARENT_PRC_NUM_Y = ${DATA_BDY_SCALE_PRC_NUM_Y}," \
          )"
   mkdir -p $TMP/f$(printf $MEMBER_FMT 1)
   conf_file="$TMP/f$(printf $MEMBER_FMT 1)/pp.d01_${STIME}.conf"
   echo "$conf" > ${conf_file}

#   for q in $(seq ${SCALE_NP[1]}); do
#      pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
#      path="${TMPROOT}/topo/topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
#      ln -sf $pathin $path
#
#      pathin2="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
#      path2="${TMPROOT}/landuse/landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
#      ln -sf $pathin2 $path2
#   done
fi


mkdir -p ${OUTDIR[$d]}/score
mkdir -p ${OUTDIR[$d]}/obsdep

time=$STIME
atime=$(datetime $time $LCYCLE s)
time_bdy_start_prev=0
loop=0
while ((time <= ETIME)); do
  loop=$((loop+1))
  time_efso=$(datetime $atime -${EFSO_FCST_LENGTH} s)

#  for s in $(seq $nsteps); do
#    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then

#      if ((s == 1)); then
#        if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
#          continue
#        elif ((BDY_FORMAT == 0)); then
#          continue
#        elif ((LANDUSE_UPDATE != 1 && loop > 1)); then
#          continue
#        fi
#      fi
#      if ((s == 2)); then
#        if ((BDY_FORMAT == 0)); then
#          continue
#        fi
#      fi
#      if ((s == 4)); then
#        if ((OBSOPE_RUN == 0)); then
#          continue
#        fi
#      fi

#    fi
#  done

  obstime $time

  if ((BDY_FORMAT != 0)); then

  bdy_setting $time $CYCLEFLEN $BDYCYCLE_INT "$BDYINT" "$PARENT_REF_TIME" "$BDY_SINGLE_FILE"

    #---------------------------------------------------------------------------
    # scale_init (launcher)
    #---------------------------------------------------------------------------

    if ((BDY_ENS == 1)); then
      config_file_scale_launcher cycle scale-rm_init_ens "<member>/init" $mtot
    elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
      config_file_scale_launcher cycle scale-rm_init_ens "<member>/init" 1
    else # local run directory: run multiple members as needed
      config_file_scale_launcher cycle scale-rm_init_ens "<member>/init" $((repeat_mems <= mtot ? repeat_mems : mtot))
    fi

    #---------------------------------------------------------------------------
    # scale_init (each member)
    #---------------------------------------------------------------------------

#    if (((loop == 1 && MAKEINIT == 1) && ${bdy_times[1]} != time)); then
#      echo "[Error] $0: Unable to generate initial analyses (MAKEINIT) at this time" >&2
#      echo "        that does not fit to any boundary data." >&2
#      exit 1
#    fi

    if ((DISK_MODE >= 1));then
      TOPO_PATH=${TMP}
      LANDUSE_PATH=${TMP}
      RESTART_IN_PATH[$d]=${TMP}
      RESTART_OUT_PATH[$d]=${TMP}
      BOUNDARY_PATH[$d]=${TMP}
      CONSTDB_PATH=$TMPROOT_CONSTDB/dat
     else
      TOPO_PATH="${DATA_TOPO}/const"
      LANDUSE_PATH="${DATA_LANDUSE}/const"
      RESTART_IN_PATH[$d]=${INDIR[$d]}/$time/anal
      RESTART_OUT_PATH[$d]=${OUTDIR[$d]}/${time}/anal
      BOUNDARY_PATH[$d]=${OUTDIR[$d]}/$time/bdy

      CONSTDB_PATH=$SCALEDIR/scale-rm/test/data

      if [ $PRESET = 'FUGAKU' ] ; then

        CONSTDB_PATH=$TMPROOT_CONSTDB/dat

        if (( BDY_LLIO_TMP == 1 )) ; then
          if ((BDY_ENS ==1));then
            BOUNDARY_PATH[$d]=/local/$time/bdy
          else
            BOUNDARY_PATH[$d]=/share/$time/bdy
          fi
        fi

      fi

    fi

    if ((BDY_ROTATING == 1 || ${bdy_times[1]} != time_bdy_start_prev)); then
      time_bdy_start_prev=${bdy_times[1]}
    fi
 
    ith=0
    for m in $(seq $mtot); do
      ith=$((ith+1))
      config_file_init_core $m &
      if (( ith == SHELL_PROCS )) || ((m == mtot)) ; then 
         wait 
         ith=0
      fi
    done

  fi # [ BDY_FORMAT != 0 ]

  #-----------------------------------------------------------------------------
  # scale (launcher)
  #-----------------------------------------------------------------------------

  config_file_scale_launcher cycle scale-rm_ens "<member>/run" $mtot

  #-----------------------------------------------------------------------------
  # scale (each member)
  #-----------------------------------------------------------------------------

    if ((DISK_MODE >= 1));then
      TOPO_PATH=${TMP}
      LANDUSE_PATH=${TMP}
      HISTORY_PATH[$d]=${TMP}
      RESTART_IN_PATH[$d]=${TMP}
      RESTART_OUT_PATH[$d]=${TMP}
      BOUNDARY_PATH[$d]=${TMP}
      CONSTDB_PATH=$TMPROOT_CONSTDB/dat
      HISTORY_ESFO_PATH=$TMP
    else
      TOPO_PATH="${DATA_TOPO}/const"
      LANDUSE_PATH="${DATA_LANDUSE}/const"
      HISTORY_PATH[$d]=${OUTDIR[$d]}/$time/hist
      HISTORY_ESFO_PATH=${OUTDIR[$d]}/$atime/anal
      if [ $PRESET == 'FUGAKU' ] && (( HIST_LLIO_TMP == 1)) ; then
        HISTORY_PATH[$d]=/local/$time/hist
        if (( atime <= ETIME ));then
          HISTORY_EFSO_PATH=/local/${atime}/gues
        fi
      fi

      loop_prev=$((loop-1))
      loop_next=$((loop+1))
      if ((MAKEINIT == 0)) &&  ((time == STIME)) ; then
        RESTART_IN_PATH[$d]=${INDIR[$d]}/$time/anal
      else
        RESTART_IN_PATH[$d]=${OUTDIR[$d]}/$time/anal
      fi 
      RESTART_OUT_PATH[$d]=${OUTDIR[$d]}/${atime}/anal

      if [ $PRESET == 'FUGAKU' ] && (( ANAL_LLIO_TMP == 1 )) ; then
         if (( loop_prev % ANAL_LLIO_TMP_SKIP != 0 )); then
           RESTART_IN_PATH[$d]=/local/$time/anal
         fi
         if (( loop % ANAL_LLIO_TMP_SKIP != 0 && atime <= ETIME )); then
           RESTART_OUT_PATH[$d]=/local/$atime/anal
         fi
      fi

      BOUNDARY_PATH[$d]=${OUTDIR[$d]}/$time/bdy

      CONSTDB_PATH=$SCALEDIR/scale-rm/test/data
      if [ $PRESET = 'FUGAKU' ] ; then

        CONSTDB_PATH=$TMPROOT_CONSTDB/dat
        if (( BDY_LLIO_TMP == 1 )) ; then
          if ((BDY_ENS ==1));then
            BOUNDARY_PATH[$d]=/local/$time/bdy
          else
            BOUNDARY_PATH[$d]=/share/$time/bdy
          fi
        fi

      fi

    fi
 
    ith=0
    for m in $(seq $mtot); do
      ith=$((ith+1))
###      echo config_file_scale_core $loop $m 
      config_file_scale_core $m &
      if (( ith == SHELL_PROCS )) || ((m == mtot)) ; then 
         wait 
         ith=0
      fi
    done

  #-----------------------------------------------------------------------------
  # letkf
  #-----------------------------------------------------------------------------

  OBS_IN_NAME_LIST=
  OBS_IN_NAME_LIST_EFSO=
  for iobs in $(seq $OBSNUM); do
    if [ "${OBSNAME[$iobs]}" != '' ]; then
      OBSFOOT[$iobs]=${OBSFOOT[$iobs]:-dat}
      if ((DISK_MODE_OBS >= 1)); then
        OBS_IN_NAME_LIST="${OBS_IN_NAME_LIST}'${TMPROOT_OBS}/obs/${OBSNAME[$iobs]}_${atime}.${OBSFOOT[$iobs]}', "
        OBS_IN_NAME_LIST_EFSO="${OBS_IN_NAME_LIST_EFSO}'${TMPROOT_OBS}/obs/${OBSNAME[$iobs]}_${time_efso}.${OBSFOOT[$iobs]}', "
      else
        OBS_IN_NAME_LIST="${OBS_IN_NAME_LIST}'${OBS}/${OBSNAME[$iobs]}_${atime}.${OBSFOOT[$iobs]}', "
        OBS_IN_NAME_LIST_EFSO="${OBS_IN_NAME_LIST_EFSO}'${OBS}/${OBSNAME[$iobs]}_${time_efso}.${OBSFOOT[$iobs]}', "
      fi
    fi
  done

  OBSDA_RUN_LIST=
  for iobs in $(seq $OBSNUM); do
    if [ -n "${OBSOPE_SEPARATE[$iobs]}" ] && ((${OBSOPE_SEPARATE[$iobs]} == 1)); then
      OBSDA_RUN_LIST="${OBSDA_RUN_LIST}.false., "
    else
      OBSDA_RUN_LIST="${OBSDA_RUN_LIST}.true., "
    fi
  done

  DET_RUN_TF='.false.'
  if ((DET_RUN == 1)); then
    DET_RUN_TF='.true.'
  fi

  EFSO_RUN_TF='.false.'
  if (( EFSO_RUN == 1 )); then
    EFSO_RUN_TF='.true.'
  fi

  OBSDA_OUT='.false.'
  if ((OBSOUT_OPT <= 2)); then
    OBSDA_OUT='.true.'
  fi
  SPRD_OUT_TF='.true.'
  if ((SPRD_OUT == 0)); then
    SPRD_OUT_TF='.false.'
  fi
  RTPS_INFL_OUT_TF='.false.'
  if ((RTPS_INFL_OUT == 1)); then
    RTPS_INFL_OUT_TF='.true.'
  fi
  NOBS_OUT_TF='.false.'
  if ((NOBS_OUT == 1)); then
    NOBS_OUT_TF='.true.'
  fi
  INFL_MUL_ADAPTIVE='.false.'
  if ((ADAPTINFL == 1)); then
    INFL_MUL_ADAPTIVE='.true.'
  fi
  INFL_ADD=0.0
  if ((ADDINFL == 1)); then
    INFL_ADD=$INFL_ADD_FACT
  fi
 
  for d in $(seq $DOMNUM); do
    dfmt=$(printf $DOMAIN_FMT $d)

    if ((ISTEP > 3)) ;then 
      for m in $mtot ; do
        if [ "${name_m[$m]}" == 'mgue' ] ; then
          continue
        fi
        cp ${OUTDIR[$d]}/$atime/gues/${name_m[$m]}/* ${OUTDIR[$d]}/$atime/anal/${name_m[$m]}
      done
      cp ${OUTDIR[$d]}/$atime/gues/mean/* ${OUTDIR[$d]}/$atime/gues/sprd
      cp ${OUTDIR[$d]}/$atime/gues/mean/* ${OUTDIR[$d]}/$atime/anal/sprd
    fi

    if ((d == 1)); then
      conf_file_src=$TMP/config.nml.letkf
#      conf_file_src2=$SCRP_DIR/config.nml.scale
      conf_file="$TMP/config/letkf_${atime}.conf"
      conf_file_efso="$TMP/config/efso_${atime}.conf"
    else
      conf_file_src=$TMP/config.nml.letkf.d$d
      #conf_file_src2=$SCRP_DIR/config.nml.scale.d$d
      conf_file="$TMP/config/letkf.d${dfmt}_${atime}.conf"
      conf_file_efso="$TMP/config/efso.d${dfmt}_${atime}.conf"
    fi

    conf_file_src2="$TMP/${name_m[$mmean]}/run.d${dfmt}_${time}.conf"
 
    rm -rf ${OUTDIR[$d]}/$atime/log/letkf
    rm -rf ${OUTDIR[$d]}/$atime/log/efso
    mkdir -p ${OUTDIR[$d]}/$atime/log/letkf
    mkdir -p ${OUTDIR[$d]}/$atime/log/efso

    OBSDEP_OUT_TF=".false."
    OBSDEP_OUT_NC_TF=".false."
    OBSANAL_OUT_TF=".false."
    if (( OBSOUT_OPT < 4 )) ; then
      OBSDEP_OUT_TF=".true."
      OBSDEP_OUT_BASENAME="${OUTDIR[$d]}/obsdep/obsdep_${atime}"
    fi
    OBSDEP_IN_BASENAME="${OUTDIR[$d]}/obsdep/obs/obsdep_${time}" # EFSO
    DEPARTURE_STAT_OUT_BASENAME="${OUTDIR[$d]}/score/score_${atime}"
    OBSNUM_OUT_NC_BASENAME="${OUTDIR[$d]}/score/obsnum_${atime}"
    OBSANAL_IN_BASENAME="${OUTDIR[$d]}/${time}/obs"
    OBSANAL_OUT_BASENAME="${OUTDIR[$d]}/${atime}/obs"
    if ((EFSO_RUN==1));then
      OBSDEP_OUT_NC_TF=".true."
      OBSANAL_OUT_TF=".true."
    fi
    if ((DISK_MODE >= 1)) ;then
#      GUES_IN_BASENAME="${RESTART_OUT_PATH[$d]}/<member>/gues_$(datetime_scale $atime)"
      GUES_IN_BASENAME="${RESTART_OUT_PATH[$d]}/<member>/anal_$(datetime_scale $atime)"
      GUES_MEAN_INOUT_BASENAME="${RESTART_OUT_PATH[$d]}/mean/gues_$(datetime_scale $atime)"
      GUES_SPRD_OUT_BASENAME="${RESTART_OUT_PATH[$d]}/sprd/gues_$(datetime_scale $atime)"
      ANAL_OUT_BASENAME="${RESTART_OUT_PATH[$d]}/<member>/anal_$(datetime_scale $atime)"
      EFSO_ANAL_IN_BASENAME="${RESTART_OUT_PATH[$d]}/mean/anal_$(datetime_scale $atime)"
      RESTART_IN_BASENAME_SCALE="${RESTART_OUT_PATH[$d]}/<member>/gues"
#      EFSO_FCST_FROM_GUES_BASENAME="${HISTORY_EFSO_PATH}/mgue/init_$(datetime_scale $atime)"
#      EFSO_FCST_FROM_ANAL_BASENAME="${HISTORY_EFSO_PATH}/mean/init_$(datetime_scale $atime)"
    else
#      GUES_IN_BASENAME="${RESTART_OUT_PATH[$d]}/../gues/<member>/init_$(datetime_scale $atime)"
      GUES_IN_BASENAME="${RESTART_OUT_PATH[$d]}/../anal/<member>/init_$(datetime_scale $atime)"
      GUES_MEAN_INOUT_BASENAME="${RESTART_OUT_PATH[$d]}/../gues/mean/init_$(datetime_scale $atime)"
      GUES_SPRD_OUT_BASENAME="${RESTART_OUT_PATH[$d]}/../gues/sprd/init_$(datetime_scale $atime)"
      ANAL_OUT_BASENAME="${RESTART_OUT_PATH[$d]}/<member>/init_$(datetime_scale $atime)"
      EFSO_ANAL_IN_BASENAME="${RESTART_OUT_PATH[$d]}/mean/init_$(datetime_scale $atime)"
      EFSO_FCST_FROM_GUES_BASENAME="${HISTORY_EFSO_PATH}/mgue/init_$(datetime_scale $atime)"
      EFSO_FCST_FROM_ANAL_BASENAME="${HISTORY_EFSO_PATH}/mean/init_$(datetime_scale $atime)"
      EFSO_EFCST_FROM_ANAL_BASENAME="${HISTORY_EFSO_PATH}/<member>/init_$(datetime_scale $atime)"
      RESTART_IN_BASENAME_SCALE="${RESTART_OUT_PATH[$d]}/../gues/<member>/init"
    fi

    cat $TMP/config.nml.ensmodel | \
        sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
            -e "/!--CONF_FILES--/a CONF_FILES = \"letkf.d<domain>_${atime}.conf\"," \
            -e "/!--DET_RUN--/a DET_RUN = ${DET_RUN_TF}," \
            -e "/!--EFSO_RUN--/a EFSO_RUN = ${EFSO_RUN_TF}," \
            -e "/!--PPN--/a PPN = $PPN_APPAR," \
            -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
            -e "/!--NUM_DOMAIN--/a NUM_DOMAIN = $DOMNUM," \
            -e "/!--PRC_DOMAINS--/a PRC_DOMAINS = $PRC_DOMAINS_LIST" \
        > ${conf_file}

    cat $conf_file_src | \
        sed -e "/!--OBS_IN_NUM--/a OBS_IN_NUM = $OBSNUM," \
            -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST" \
            -e "/!--OBSDA_RUN--/a OBSDA_RUN = $OBSDA_RUN_LIST" \
            -e "/!--OBSDA_OUT--/a OBSDA_OUT = $OBSDA_OUT" \
            -e "/!--OBSDA_OUT_BASENAME--/a OBSDA_OUT_BASENAME = \"<member>/obsgues.d${dfmt}_${atime}\"," \
            -e "/!--HISTORY_IN_BASENAME--/a HISTORY_IN_BASENAME = \"${HISTORY_PATH[$d]}/<member>/history\"," \
            -e "/!--SLOT_START--/a SLOT_START = $slot_s," \
            -e "/!--SLOT_END--/a SLOT_END = $slot_e," \
            -e "/!--SLOT_BASE--/a SLOT_BASE = $slot_b," \
            -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL = ${LTIMESLOT}.D0," \
            -e "/!--OBSDA_IN--/a OBSDA_IN = .false.," \
            -e "/!--GUES_IN_BASENAME--/a GUES_IN_BASENAME = \"${GUES_IN_BASENAME}\"," \
            -e "/!--GUES_MEAN_INOUT_BASENAME--/a GUES_MEAN_INOUT_BASENAME = \"${GUES_MEAN_INOUT_BASENAME}\"," \
            -e "/!--GUES_SPRD_OUT_BASENAME--/a GUES_SPRD_OUT_BASENAME = \"${GUES_SPRD_OUT_BASENAME}\"," \
            -e "/!--GUES_SPRD_OUT--/a GUES_SPRD_OUT = ${SPRD_OUT_TF}," \
            -e "/!--ANAL_OUT_BASENAME--/a ANAL_OUT_BASENAME = \"${ANAL_OUT_BASENAME}\"," \
            -e "/!--ANAL_SPRD_OUT--/a ANAL_SPRD_OUT = ${SPRD_OUT_TF}," \
            -e "/!--LETKF_TOPOGRAPHY_IN_BASENAME--/a LETKF_TOPOGRAPHY_IN_BASENAME = \"${TOPO_PATH}/topo/topo\"," \
            -e "/!--EFSO_ANAL_IN_BASENAME--/a EFSO_ANAL_IN_BASENAME = \"${EFSO_ANAL_IN_BASENAME}\"," \
            -e "/!--EFSO_FCST_FROM_GUES_BASENAME--/a EFSO_FCST_FROM_GUES_BASENAME = \"${EFSO_FCST_FROM_GUES_BASENAME}\"," \
            -e "/!--EFSO_FCST_FROM_ANAL_BASENAME--/a EFSO_FCST_FROM_ANAL_BASENAME = \"${EFSO_FCST_FROM_ANAL_BASENAME}\"," \
            -e "/!--EFSO_EFCST_FROM_ANAL_BASENAME--/a EFSO_EFCST_FROM_ANAL_BASENAME = \"${EFSO_EFCST_FROM_ANAL_BASENAME}\"," \
            -e "/!--INFL_ADD--/a INFL_ADD = ${INFL_ADD}," \
            -e "/!--INFL_ADD_IN_BASENAME--/a INFL_ADD_IN_BASENAME = \"${DATA_ADDINFL[$d]}/<member>${CONNECTOR}init${sfx}\"," \
            -e "/!--INFL_MUL_ADAPTIVE--/a INFL_MUL_ADAPTIVE = ${INFL_MUL_ADAPTIVE}," \
            -e "/!--INFL_MUL_OUT_BASENAME--/a INFL_MUL_OUT_BASENAME = \"infl.d${dfmt}_$(datetime_scale $atime)\"," \
            -e "/!--RELAX_SPREAD_OUT--/a RELAX_SPREAD_OUT = ${RTPS_INFL_OUT_TF}," \
            -e "/!--RELAX_SPREAD_OUT_BASENAME--/a RELAX_SPREAD_OUT_BASENAME = \"rtpsinfl.d${dfmt}_$(datetime_scale $atime)\"," \
            -e "/!--NOBS_OUT--/a NOBS_OUT = ${NOBS_OUT_TF}," \
            -e "/!--NOBS_OUT_BASENAME--/a NOBS_OUT_BASENAME = \"nobs.d${dfmt}_$(datetime_scale $atime)\"," \
            -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME =  \"${OUTDIR[$d]}/$atime/log/letkf/${name_m[$m]}_LOG\"," \
            -e "/!--DEPARTURE_STAT_OUT_BASENAME--/a DEPARTURE_STAT_OUT_BASENAME = \"${DEPARTURE_STAT_OUT_BASENAME}\"," \
            -e "/!--OBSDEP_OUT--/a OBSDEP_OUT = ${OBSDEP_OUT_TF}," \
            -e "/!--OBSDEP_OUT_NC--/a OBSDEP_OUT_NC = ${OBSDEP_OUT_NC_TF}," \
            -e "/!--OBSANAL_OUT--/a OBSANAL_OUT = ${OBSANAL_OUT_TF}," \
            -e "/!--OBSANAL_IN_BASENAME--/a OBSANAL_IN_BASENAME = \"${OBSANAL_IN_BASENAME}\"," \
            -e "/!--OBSANAL_OUT_BASENAME--/a OBSANAL_OUT_BASENAME = \"${OBSANAL_OUT_BASENAME}\"," \
            -e "/!--OBSDEP_IN_BASENAME--/a OBSDEP_IN_BASENAME = \"${OBSDEP_IN_BASENAME}\"," \
            -e "/!--OBSDEP_OUT_BASENAME--/a OBSDEP_OUT_BASENAME = \"${OBSDEP_OUT_BASENAME}\"," \
            -e "/!--OBSNUM_OUT_NC_BASENAME--/a OBSNUM_OUT_NC_BASENAME = \"${OBSNUM_OUT_NC_BASENAME}\"," \
        >> ${conf_file}

    # Most of these parameters are not important for letkf
    cat $conf_file_src2 | \
        sed -e "s#^RESTART_IN_BASENAME.*#RESTART_IN_BASENAME = \"${RESTART_IN_BASENAME_SCALE}\", #g " \
            -e "s#^TIME_STARTDATE\ =.*#TIME_STARTDATE\ =\ ${atime:0:4},\ ${atime:4:2},\ ${atime:6:2},\ ${atime:8:2},\ ${atime:10:2},\ ${atime:12:2}, #g" \
        >> ${conf_file}


    cat ${conf_file} | \
        sed -e '/^OBS_IN_NAME/d' \
            -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST_EFSO" \
    >> ${conf_file_efso}

#    if ((stage_config == 1)); then
#      echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
#    fi
    if (( PAWR_DECODE == 1 )) ; then
      conf_file_dec_pawr="$TMP/config/dec_pawr_${atime}.conf"
      PAWR_IN_PATH="${PAWR_RAW}/${FNAME_PAWR_RAW}"
      if (( DISK_MODE_OBS >= 1)) ;then
        OUT_PAWR_SUPEROB_PATH="${TMPROOT_OBS}/obs/${OBSNAME[$iobs]}"
        mkdir -p ${TMPROOT_OBS}/obs
      else
        OUT_PAWR_SUPEROB_PATH="${OBS}/${OBSNAME[$iobs]}"
        mkdir -p $OBS
      fi
      cat ${conf_file} | \
      sed -e "/!--PAWR_IN_PATH--/a PAWR_IN_PATH = \"${PAWR_IN_PATH}\"," \
          -e "/!--OUT_PAWR_SUPEROB_PATH--/a OUT_PAWR_SUPEROB_PATH = \"${OUT_PAWR_SUPEROB_PATH}\", " \
          -e "/!--MEMBER_RUN--/a MEMBER_RUN = 1, " \
          -e "s#^MEMBER\ =.*#MEMBER\ =\ 1,#g" \
          -e "s#^CONF_FILES\ =.*#CONF_FILES\ =\ ${conf_file_dec_pawr},#g" \
          -e "s#^TIME_STARTDATE\ =.*#TIME_STARTDATE\ =\ ${atime:0:4},\ ${atime:4:2},\ ${atime:6:2},\ ${atime:8:2},\ ${atime:10:2},\ ${atime:12:2}, #g" \
      >> ${conf_file_dec_pawr}
    fi

   done # [ d in $(seq $DOMNUM) ]

  #-------------------
  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done

echo

#
echo "backup"
cp $TMPS/config/*.conf $OUTDIR/config/
echo

#-------------------------------------------------------------------------------
}

    config_file_init_core (){

    mlocal=$1

      if ((BDY_ENS == 1)); then
        mem_bdy=${name_m[$mlocal]}
      else
        mem_bdy='mean'
      fi

      if ((BDY_FORMAT == 1 || BDY_FORMAT == 2)); then
        FILETYPE_ORG='NetCDF'
      elif ((BDY_FORMAT == 4)); then
        FILETYPE_ORG='GrADS'
      elif ((BDY_FORMAT == 5)); then
        FILETYPE_ORG=
      else
        echo "[Error] $0: Unsupport boundary file types." >&2
        exit 1
      fi
      if ((BDY_FORMAT == 4)); then
        BASENAME_ORG="${TMPROOT_BDYDATA}/${mem_bdy}/gradsbdy_${time}.conf"
      else
        if ((nbdy <= 1)); then
          bdy_no_suffix="_$(printf %05d 0)"
        else
          bdy_no_suffix=
        fi
        BASENAME_ORG="${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}"
      fi

      for d in $(seq $DOMNUM); do
        dfmt=$(printf $DOMAIN_FMT $d)

        if ((d == 1)); then
          conf_file_src=$TMP/config.nml.scale_init
        else
          conf_file_src=$TMP/config.nml.scale_init.d$d
        fi

        if (((loop == 1 && MAKEINIT == 1) || USE_INIT_FROM_BDY == 1)); then
          RESTART_OUTPUT='.true.'
          RESTART_OUT_BASENAME=${BOUNDARY_PATH[$d]}/${mem_bdy}/init_bdy
        else
          RESTART_OUTPUT='.false.'
        fi

        RESTART_OUT_POSTFIX_TIMELABEL_TF=".true."

        if ((mlocal==mmean)); then
          rm -rf  ${OUTDIR[$d]}/$time/log/scale_init
          mkdir -p ${OUTDIR[$d]}/$time/log/scale_init
        fi

      if ((BDY_ENS == 1)) || ((mlocal==mmean)) ; then
        rm -rf ${OUTDIR[$d]}/$time/bdy/$mem_bdy
        mkdir -p ${OUTDIR[$d]}/$time/bdy/$mem_bdy
        mkdir -p ${OUTDIR[$d]}/$time/anal/$mem_bdy
      fi
        conf="$(cat $conf_file_src | \
            sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"${OUTDIR[$d]}/$time/log/scale_init/${name_m[$mlocal]}_LOG\"," \
                -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
                -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
                -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = ${RESTART_OUTPUT}," \
                -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${RESTART_OUT_BASENAME}\"," \
                -e "/!--RESTART_OUT_POSTFIX_TIMELABEL--/a RESTART_OUT_POSTFIX_TIMELABEL = ${RESTART_OUT_POSTFIX_TIMELABEL_TF}," \
                -e "/!--TOPOGRAPHY_IN_BASENAME--/a TOPOGRAPHY_IN_BASENAME = \"${TOPO_PATH}/topo/topo\"," \
                -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${LANDUSE_PATH}/landuse/landuse\"," \
                -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${CONSTDB_PATH}/land/param.bucket.conf\",")"
        conf="$(echo "$conf" | \
          sed -e "/!--BASENAME_ORG--/a BASENAME_ORG = \"${BASENAME_ORG}\"," \
              -e "/!--FILETYPE_ORG--/a FILETYPE_ORG = \"${FILETYPE_ORG}\"," \
              -e "/!--BOUNDARY_UPDATE_DT--/a BOUNDARY_UPDATE_DT = ${BDYINT}.D0,"\
              -e "/!--BASENAME_BOUNDARY--/a BASENAME_BOUNDARY = \"${BOUNDARY_PATH[$d]}/${mem_bdy}/boundary\"," \
              -e "/!--NUMBER_OF_FILES--/a NUMBER_OF_FILES = ${nbdy}," \
              -e "/!--NUMBER_OF_TSTEPS--/a NUMBER_OF_TSTEPS = ${ntsteps}," \
              -e "/!--NUMBER_OF_SKIP_TSTEPS--/a NUMBER_OF_SKIP_TSTEPS = ${ntsteps_skip},")"


        conf_file="$TMP/${name_m[$mlocal]}/init.d${dfmt}_${time}.conf"
        echo "$conf" > ${conf_file}

      done # [ d in $(seq $DOMNUM) ]

      #if ((BDY_FORMAT == 4 && (BDY_ENS == 0 || m == 1))); then
      if ((BDY_FORMAT == 4 )); then
        conf_file="$TMP/${mem_bdy}/gradsbdy_${time}.conf"
        if ((nbdy <= 1)); then
          bdy_no_suffix="_$(printf %05d 0)"
        else
          bdy_no_suffix=
        fi
        cat $TMP/config.nml.grads_boundary | \
            sed -e "s#--DIR--/bdyatm#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_atm_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                -e "s#--DIR--/bdysfc#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_sfc_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                -e "s#--DIR--/bdyland#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_lnd_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
            > ${conf_file}

        #if ((stage_config == 1)); then
        #  #echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
        #fi
      fi # [ BDY_FORMAT == 4 && (BDY_ENS == 0 || m == 1) ]


     }

config_file_scale_core (){

    mlocal=$1

    if ((BDY_ENS == 1)); then
      mem_bdy=${name_m[$mlocal]}
    else
      mem_bdy='mean'
    fi

    for d in $(seq $DOMNUM); do
      dfmt=$(printf $DOMAIN_FMT $d)

      ONLINE_IAM_PARENT=".false."
      if ((d < DOMNUM)); then
        ONLINE_IAM_PARENT=".true."
      fi
      ONLINE_IAM_DAUGHTER=".false."
      if ((d > 1)); then
        ONLINE_IAM_DAUGHTER=".true."
      fi
#      if ((loop == 1 && MAKEINIT == 1)); then
#        RESTART_IN_BASENAME="${name_m[$m]}/init.d${dfmt}"
#      else
#        RESTART_IN_BASENAME="${name_m[$m]}/anal.d${dfmt}"
#      fi
      RESTART_IN_POSTFIX_TIMELABEL_TF=".true."
      RESTART_OUT_POSTFIX_TIMELABEL_TF=".true."

#      if ((loop == 1 )); then
#        RESTART_IN_POSTFIX_TIMELABEL_TF=".false."
#      else
#        RESTART_IN_POSTFIX_TIMELABEL_TF=".true."
#      fi

      RESTART_OUTPUT_TF=".true."

      mkdir -p ${OUTDIR[$d]}/$atime/anal/${name_m[$mlocal]}
      mkdir -p ${OUTDIR[$d]}/$atime/gues/${name_m[$mlocal]}
      mkdir -p ${OUTDIR[$d]}/$time/hist/${name_m[$mlocal]}

      if ((DISK_MODE >= 1)) ;then
        RESTART_IN_BASENAME[$d]="${RESTART_IN_PATH[$d]}/${name_m[$mlocal]}/anal"
        RESTART_OUT_BASENAME[$d]="${RESTART_OUT_PATH[$d]}/${name_m[$mlocal]}/anal"
      else
        if (( MAKEINIT == 1 )) && (( loop == 1 | USE_INIT_FROM_BDY == 1 )) ; then
          RESTART_IN_BASENAME[$d]="${BOUNDARY_PATH[$d]}/${mem_bdy}/init_bdy"
        else
          RESTART_IN_BASENAME[$d]="${RESTART_IN_PATH[$d]}/${name_m[$mlocal]}/init"
        fi
        RESTART_OUT_BASENAME[$d]="${RESTART_OUT_PATH[$d]}/${name_m[$mlocal]}/init"
      fi

      if ((WINDOW_S == LCYCLE && WINDOW_E == LCYCLE));then
        HISTORY_OUT_WAIT=$((LCYCLE+LTIMESLOT)) ### 3D-LETKF : suppress history output
      else
        HISTORY_OUT_WAIT=0
      fi

      if ((d == 1)); then
        conf_file_src=$TMP/config.nml.scale
      else
        conf_file_src=$TMP/config.nml.scale.d$d
      fi
 
      if ((mlocal==mmean)); then
        rm -rf ${OUTDIR[$d]}/$time/log/scale
        mkdir -p ${OUTDIR[$d]}/$time/log/scale
      fi

      conf="$(cat $conf_file_src | \
          sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"${OUTDIR[$d]}/$time/log/scale/${name_m[$mlocal]}_LOG\"," \
              -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
              -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
              -e "/!--TIME_DURATION--/a TIME_DURATION = ${CYCLEFLEN}.D0," \
              -e "/!--TIME_DT_ATMOS_RESTART--/a TIME_DT_ATMOS_RESTART = ${LCYCLE}.D0," \
              -e "/!--TIME_DT_OCEAN_RESTART--/a TIME_DT_OCEAN_RESTART = ${LCYCLE}.D0," \
              -e "/!--TIME_DT_LAND_RESTART--/a TIME_DT_LAND_RESTART = ${LCYCLE}.D0," \
              -e "/!--TIME_DT_URBAN_RESTART--/a TIME_DT_URBAN_RESTART = ${LCYCLE}.D0," \
              -e "/!--ONLINE_DOMAIN_NUM--/a ONLINE_DOMAIN_NUM = ${d}," \
              -e "/!--ONLINE_IAM_PARENT--/a ONLINE_IAM_PARENT = ${ONLINE_IAM_PARENT}," \
              -e "/!--ONLINE_IAM_DAUGHTER--/a ONLINE_IAM_DAUGHTER = ${ONLINE_IAM_DAUGHTER}," \
              -e "/!--RESTART_IN_BASENAME--/a RESTART_IN_BASENAME = \"${RESTART_IN_BASENAME[$d]}\"," \
              -e "/!--RESTART_IN_POSTFIX_TIMELABEL--/a RESTART_IN_POSTFIX_TIMELABEL = ${RESTART_IN_POSTFIX_TIMELABEL_TF}," \
              -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = ${RESTART_OUTPUT_TF}," \
              -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${RESTART_OUT_BASENAME[$d]}\"," \
              -e "/!--RESTART_OUT_POSTFIX_TIMELABEL--/a RESTART_OUT_POSTFIX_TIMELABEL = ${RESTART_OUT_POSTFIX_TIMELABEL_TF}," \
              -e "/!--TOPOGRAPHY_IN_BASENAME--/a TOPOGRAPHY_IN_BASENAME = \"${TOPO_PATH}/topo/topo\"," \
              -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${LANDUSE_PATH}/landuse/landuse\"," \
              -e "/!--FILE_HISTORY_DEFAULT_BASENAME--/a FILE_HISTORY_DEFAULT_BASENAME = \"${HISTORY_PATH[$d]}/${name_m[$mlocal]}/history\"," \
              -e "/!--FILE_HISTORY_DEFAULT_TINTERVAL--/a FILE_HISTORY_DEFAULT_TINTERVAL = ${CYCLEFOUT}.D0," \
              -e "/!--FILE_HISTORY_OUTPUT_WAIT--/a FILE_HISTORY_OUTPUT_WAIT = ${HISTORY_OUT_WAIT}.D0," \
              -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"${OUTDIR[$d]}/$time/log/scale/${name_m[$mlocal]}.d${dfmt}.monitor_${time}\"," \
              -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${CONSTDB_PATH}/land/param.bucket.conf\"," \
              -e "/!--URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME--/a  URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME = \"${CONSTDB_PATH}/urban/param.kusaka01.dat\"," \
              -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${CONSTDB_PATH}/rad/PARAG.29\"," \
              -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${CONSTDB_PATH}/rad/PARAPC.29\"," \
              -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${CONSTDB_PATH}/rad/VARDATA.RM29\"," \
              -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${CONSTDB_PATH}/rad/cira.nc\"," \
              -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${CONSTDB_PATH}/rad/MIPAS\"," \
              -e "/!--ATMOS_PHY_LT_LUT_FILENAME--/a ATMOS_PHY_LT_LUT_FILENAME = \"${CONSTDB_PATH}/lightning/LUT_TK1978_v.txt\"," )" 
      if ((d == 1)); then
        conf="$(echo "$conf" | \
            sed -e "/!--ATMOS_BOUNDARY_IN_BASENAME--/a ATMOS_BOUNDARY_IN_BASENAME = \"${BOUNDARY_PATH[$d]}/${mem_bdy}/boundary\"," )"
      fi

      conf_file="$TMPS/${name_m[$mlocal]}/run.d${dfmt}_${time}.conf"
      echo "$conf" > ${conf_file}

      if ((ENABLE_PARAM_USER == 1)) && [ -e "$TMP/config.nml.scale_user" ]; then
        conf="$(cat $TMP/config.nml.scale_user)"
        if ((OCEAN_INPUT == 1)); then
          if ((OCEAN_FORMAT == 99)); then
#            conf="$(echo "$conf" | \
#                sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"${RESTART_IN_BASENAME[$d]}\",")"
            conf="$(echo "$conf" | \
                 sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"${BOUNDARY_PATH[$d]}/${mem_bdy}/init_bdy_$(datetime_scale $time)\",")"
          fi
        fi
        if ((LAND_INPUT == 1)); then
          if ((LAND_FORMAT == 99)); then
            conf="$(echo "$conf" | \
                sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"${BOUNDARY_PATH[$d]}/${mem_bdy}/init_bdy_$(datetime_scale $time)\",")"
#            conf="$(echo "$conf" | \
#                sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"${RESTART_IN_BASENAME[$d]}\",")"
          fi
        fi
        echo "$conf" >> ${conf_file}
      fi

    done # [ d in $(seq $DOMNUM) ]
}

#===============================================================================

setting () {
#-------------------------------------------------------------------------------
# define steps

nsteps=6
stepname[1]='Run SCALE pp'
stepexecdir[1]="$TMPRUN/scale_pp"
stepexecname[1]="scale-rm_pp_ens"
stepname[2]='Run SCALE init'
stepexecdir[2]="$TMPRUN/scale_init"
stepexecname[2]="scale-rm_init_ens"
stepname[3]='Run ensemble forecasts'
stepexecdir[3]="$TMPRUN/scale"
stepexecname[3]="scale-rm_ens"
if (( OBSOPE_RUN == 0 )) && (( PAWR_DECODE == 1 )) ; then
stepname[4]='Run PAWR decoder'
stepexecdir[4]="$TMPRUN/dec_pawr"
stepexecname[4]="dec_pawr"
elif (( OBSOPE_RUN == 1 )) && (( PAWR_DECODE == 1 )) ; then
  echo "OBSOPE=1 and PAWR_DECODE=1 is not supported. "
  exit 1 
else 
stepname[4]='Run observation operator'
stepexecdir[4]="$TMPRUN/obsope"
stepexecname[4]="obsope"
fi 

stepname[5]='Run LETKF'
stepexecdir[5]="$TMPRUN/letkf"
stepexecname[5]="letkf"

stepname[6]='Run EFSO'
stepexecdir[6]="$TMPRUN/efso"
stepexecname[6]="efso"

if [ PRESET == "FUGAKU" ] && (( USE_LLIO_BIN == 1 )); then
  stepexecbin[1]="$TMP/scale-rm_pp_ens"
  stepexecbin[2]="$TMP/scale-rm_init_ens"
  stepexecbin[3]="$TMP/scale-rm_ens"
  if (( OBSOPE_RUN == 0 )) && (( PAWR_DECODE == 1 )) ; then
    stepexecbin[4]="$DIR/obs/dec_pawr"
  else
    stepexecbin[4]="$DIR/obs/obsope"
  fi
  stepexecbin[5]="$TMP/letkf"
  stepexecbin[6]="$TMP/efso"
else
  for i in `seq $nsteps`; do
     stepexecbin[$i]="./${stepexecname[$i]}"
  done
fi
#-------------------------------------------------------------------------------
# usage help string

USAGE="
[$myname] Run data assimilation cycles.

Configuration files:
  config.main
  config.cycle

Steps:
$(for i in $(seq $nsteps); do echo "  ${i}. ${stepname[$i]}"; done)

Usage: $myname [STIME ETIME ISTEP FSTEP CONF_MODE TIME_LIMIT]

  STIME       Time of the first cycle (format: YYYY[MMDDHHMMSS])
  ETIME       Time of the last  cycle (format: YYYY[MMDDHHMMSS])
               (default: same as STIME)
  ISTEP       The initial step in the first cycle from which this script starts
               (default: the first step)
  FSTEP       The final step in the last cycle by which this script ends
               (default: the last step)
  CONF_MODE   Mode of creating runtime configuration files: 'dynamic' or 'static'
               (default: 'dynamic')
  TIME_LIMIT  Requested time limit (only used when using a job scheduler)
               (default: 30 minutes)
"

#if [ "$1" == '-h' ] || [ "$1" == '--help' ]; then
#  echo "$USAGE" >&2
#  exit 1
#fi

#-------------------------------------------------------------------------------
# set parameters from command line

STIME=${1:-$STIME}; shift
ETIME=${1:-$ETIME}; shift
ISTEP=${1:-$ISTEP}; shift
FSTEP=${1:-$FSTEP}; shift
CONF_MODE=${1:-$CONF_MODE}; shift
TIME_LIMIT="${1:-$TIME_LIMIT}"

#if [ -z "$STIME" ]; then
#  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
#  echo "$USAGE" >&2
#  exit 1
#fi

if ((EFSO_RUN==1)) ;then
  FSTEP=6
fi

#-------------------------------------------------------------------------------
# assign default values to and standardize the parameters

STIME=$(datetime $STIME)
ETIME=$(datetime ${ETIME:-$STIME})
ISTEP=${ISTEP:-1}
FSTEP=${FSTEP:-5}
CONF_MODE=${CONF_MODE:-"static"}
TIME_LIMIT=${TIME_LIMIT:-"0:30:00"}

#-------------------------------------------------------------------------------
# error detection

#if ((MACHINE_TYPE == 10 && ONLINE_STGOUT != 0)); then
#  echo "[Error] $myname: When \$MACHINE_TYPE = 10, \$ONLINE_STGOUT needs to be 0." >&2
#  exit 1
#fi

if [ "$CONF_MODE" != 'static' ] && ((DOMNUM > 1)); then
  echo "[Error] Online nesting with multiple domains is only allowed when the static-config mode is used (\$CONF_MODE = 'static')." 1>&2
  exit 1
fi

if ((RUN_LEVEL == 0)); then
  if ((ENABLE_PARAM_USER == 1)) && [ ! -e "$SCRP_DIR/config.nml.scale_user" ]; then
    echo "[Error] $myname: When \$ENABLE_PARAM_USER = 1, 'config.nml.scale_user' file is required." >&2
    exit 1
  fi
  if ((BDY_FORMAT == 4)) && [ ! -e "$SCRP_DIR/config.nml.grads_boundary" ]; then
    echo "[Error] $myname: When \$BDY_FORMAT = 4, 'config.nml.grads_boundary' file is required." >&2
    exit 1
  fi

  if ((MAKEINIT == 1)); then
    if [ -d "${OUTDIR}/${STIME}/anal" ]; then
      if [ -n "$(ls ${OUTDIR}/${STIME}/anal/*/*.nc 2> /dev/null)" ]; then
        echo "[Error] $myname: Initial ensemble is to be generated (\$MAKEINIT = 1) at \"${OUTDIR}/${STIME}/anal/\", but existing data are found there;" >&2
        echo "        Set \$MAKEINIT = 0 or remove \"${OUTDIR}/${STIME}/anal/*\" before running this job." >&2
        exit 1
      fi
    fi
  fi
fi

#... more detections...

#-------------------------------------------------------------------------------
# common variables

OUT_CYCLE_SKIP=${OUT_CYCLE_SKIP:-1}

CYCLEFLEN=$WINDOW_E     # Model forecast length in a cycle (second)
if (( EFSO_RUN > 0 )) ; then
  EFSO_FCST_LENGTH=${EFSO_FCST_LENGTH:-$WINDOW_E}
  if (( WINDOW_E < EFSO_FCST_LENGTH )) ; then
    CYCLEFLEN=$EFSO_FCST_LENGTH
  fi
fi

if [ -z "$FCSTOUT" ] || ((FCSTOUT >= LTIMESLOT)); then
  CYCLEFOUT=$LTIMESLOT  # Model forecast output interval (second)
elif ((LTIMESLOT % FCSTOUT == 0)); then
  CYCLEFOUT=$FCSTOUT
else
  echo "[Error] If \$FCSTOUT < \$LTIMESLOT, \$LTIMESLOT needs to be an exact multiple of \$FCSTOUT" >&2
  exit 1
fi

if ((BDY_FORMAT >= 1)) && ((BDY_FORMAT <= 4 )) ; then
  if ((BDYCYCLE_INT % BDYINT != 0)); then
    echo "[Error] \$BDYCYCLE_INT needs to be an exact multiple of \$BDYINT" >&2
    exit 1
  fi
  BDY_STARTFRAME_MAX=$((BDYCYCLE_INT / BDYINT))
  if [ -z "$PARENT_REF_TIME" ]; then
    PARENT_REF_TIME=$STIME
    for bdy_startframe in $(seq $BDY_STARTFRAME_MAX); do
      if ((BDY_FORMAT == 1)); then
        BFILE="$DATA_BDY_SCALE/${PARENT_REF_TIME}/${BDY_SCALE_DIR}/${BDY_MEAN}/history${SCALE_SFX_0}"
      elif ((BDY_FORMAT == 2 && BDY_ROTATING == 1)); then
        BFILE="$DATA_BDY_WRF/${PARENT_REF_TIME}/${BDY_MEAN}/wrfout_${PARENT_REF_TIME}" 
      elif ((BDY_FORMAT == 2 && BDY_ROTATING != 1)); then
        BFILE="$DATA_BDY_WRF/${BDY_MEAN}/wrfout_${PARENT_REF_TIME}"
      elif ((BDY_FORMAT == 4 && BDY_ROTATING == 1)); then
        BFILE="$DATA_BDY_GRADS/${PARENT_REF_TIME}/${BDY_MEAN}/atm_${PARENT_REF_TIME}.grd"
      elif ((BDY_FORMAT == 4 && BDY_ROTATING != 1)); then 
        BFILE="$DATA_BDY_GRADS/${BDY_MEAN}/atm_${PARENT_REF_TIME}.grd"
      fi

      if [ -s "$BFILE" ] ; then
        break
      fi 

      if ((bdy_startframe == BDY_STARTFRAME_MAX)); then
        echo "[Error] Cannot find boundary files. "$BFILE >&2
        exit 1
      fi

      PARENT_REF_TIME=$(datetime $PARENT_REF_TIME -${BDYINT} s)
    done
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

print_setting () {
#-------------------------------------------------------------------------------

for vname in DIR DOMAIN @INDIR @OUTDIR @DATA_TOPO DATA_TOPO_BDY_SCALE @DATA_LANDUSE DATA_BDY_SCALE \
             DATA_BDY_SCALE_PREP DATA_BDY_WRF DATA_BDY_NICAM OBS OBSNCEP DET_RUN TOPO_FORMAT \
             LANDUSE_FORMAT LANDUSE_UPDATE BDY_FORMAT BDY_ENS BDYINT BDYCYCLE_INT PARENT_REF_TIME \
             ENABLE_PARAM_USER OCEAN_INPUT OCEAN_FORMAT LAND_INPUT LAND_FORMAT OBSNUM WINDOW_S WINDOW_E \
             LCYCLE LTIMESLOT MEMBER NNODES NNODES_APPAR PPN PPN_APPAR THREADS @SCALE_NP \
             STIME ETIME ISTEP FSTEP CONF_MODE FCSTOUT MAKEINIT OUT_OPT TOPOOUT_OPT \
             LANDUSEOUT_OPT BDYOUT_OPT OBSOUT_OPT LOG_OPT LOG_TYPE; do
  if [ "${vname:0:1}" = '@' ]; then
    for d in $(seq $DOMNUM); do
      vname_d="${vname:1}[$d]"
      printf '  %-20s = %s\n' "$vname_d" "${!vname_d}"
    done
  else
    printf '  %-20s = %s\n' $vname "${!vname}"
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

archive_log () {
#-------------------------------------------------------------------------------

if ((LOG_TYPE >= 3)); then
  time=$STIME
  atime=$(datetime $time $LCYCLE s)
  while ((time <= ETIME)); do
    for d in $(seq $DOMNUM); do
      if ((LOG_OPT <= 2)) && [ -d "${OUTDIR[$d]}/${time}/log/scale_pp" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale_pp.tar scale_pp && rm -fr ${OUTDIR[$d]}/${time}/log/scale_pp ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale_pp.tar.gz scale_pp && rm -fr ${OUTDIR[$d]}/${time}/log/scale_pp ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale_pp.tar scale_pp && rm -fr ${OUTDIR[$d]}/${time}/log/scale_pp
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale_pp.tar.gz scale_pp && rm -fr ${OUTDIR[$d]}/${time}/log/scale_pp
          fi
        fi
      fi

      if ((LOG_OPT <= 2)) && [ -d "${OUTDIR[$d]}/${time}/log/scale_init" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale_init.tar scale_init && rm -fr ${OUTDIR[$d]}/${time}/log/scale_init ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale_init.tar.gz scale_init && rm -fr ${OUTDIR[$d]}/${time}/log/scale_init ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale_init.tar scale_init && rm -fr ${OUTDIR[$d]}/${time}/log/scale_init
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale_init.tar.gz scale_init && rm -fr ${OUTDIR[$d]}/${time}/log/scale_init
          fi
        fi
      fi

      if ((LOG_OPT <= 3)) && [ -d "${OUTDIR[$d]}/${time}/log/scale" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale.tar scale && rm -fr ${OUTDIR[$d]}/${time}/log/scale ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale.tar.gz scale && rm -fr ${OUTDIR[$d]}/${time}/log/scale ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale.tar scale && rm -fr ${OUTDIR[$d]}/${time}/log/scale
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale.tar.gz scale && rm -fr ${OUTDIR[$d]}/${time}/log/scale
          fi
        fi
      fi

      if ((LOG_OPT <= 4)) && [ -d "${OUTDIR[$d]}/${atime}/log/obsope" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${atime}/log -cf ${OUTDIR[$d]}/${atime}/log/obsope.tar obsope && rm -fr ${OUTDIR[$d]}/${atime}/log/obsope ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${atime}/log -czf ${OUTDIR[$d]}/${atime}/log/obsope.tar.gz obsope && rm -fr ${OUTDIR[$d]}/${atime}/log/obsope ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${atime}/log -cf ${OUTDIR[$d]}/${atime}/log/obsope.tar obsope && rm -fr ${OUTDIR[$d]}/${atime}/log/obsope
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${atime}/log -czf ${OUTDIR[$d]}/${atime}/log/obsope.tar.gz obsope && rm -fr ${OUTDIR[$d]}/${atime}/log/obsope
          fi
        fi
      fi

      if ((LOG_OPT <= 4)) && [ -d "${OUTDIR[$d]}/${atime}/log/letkf" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${atime}/log -cf ${OUTDIR[$d]}/${atime}/log/letkf.tar letkf && rm -fr ${OUTDIR[$d]}/${atime}/log/letkf ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${atime}/log -czf ${OUTDIR[$d]}/${atime}/log/letkf.tar.gz letkf && rm -fr ${OUTDIR[$d]}/${atime}/log/letkf ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${atime}/log -cf ${OUTDIR[$d]}/${atime}/log/letkf.tar letkf && rm -fr ${OUTDIR[$d]}/${atime}/log/letkf
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${atime}/log -czf ${OUTDIR[$d]}/${atime}/log/letkf.tar.gz letkf && rm -fr ${OUTDIR[$d]}/${atime}/log/letkf
          fi
        fi
      fi
    done # [ d in $(seq $DOMNUM) ]

    time=$(datetime $time $LCYCLE s)
    atime=$(datetime $time $LCYCLE s)
  done
  if ((TAR_THREAD > 1)); then
    wait
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

obstime () {
#-------------------------------------------------------------------------------
# Determine the observation time slots
#  *Require source 'func_datetime.sh'
#
# Usage: obstime TIME
#
#   TIME  Forecast start time
#
# Other input variables:
#   $LTIMESLOT
#   $WINDOW_S
#   $WINDOW_E
#   $LCYCLE
#
# Return variables:
#   $slot_s
#   $slot_e
#   $slot_b
#   $time_sl[1...$slot_e]
#   $timefmt_sl[1...$slot_e]
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local TIME="$1"

#-------------------------------------------------------------------------------

local otime=$(datetime $TIME)               # HISTORY_OUTPUT_STEP0 = .true.,
#local otime=$(datetime $TIME $LTIMESLOT s)  # HISTORY_OUTPUT_STEP0 = .false.,
local otime_s=$(datetime $TIME $WINDOW_S s)
local otime_e=$(datetime $TIME $WINDOW_E s)
local otime_a=$(datetime $TIME $LCYCLE s)
local is=0
slot_s=0
while ((otime <= otime_e)); do
  is=$((is+1))
  time_sl[$is]=$otime
  timefmt_sl[$is]="$(datetime_fmt ${otime})"
  if ((slot_s == 0 && otime >= otime_s)); then
    slot_s=$is
  fi
  if ((otime == otime_a)); then
    slot_b=$is
  fi
otime=$(datetime $otime $LTIMESLOT s)
done
slot_e=$is

#-------------------------------------------------------------------------------
}

#===============================================================================


