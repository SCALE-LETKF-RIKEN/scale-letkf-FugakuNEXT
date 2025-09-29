#===============================================================================
#
#  Prepare an initial ensemble by perturbing an initial condition.
#  August  2014,              Guo-Yuan Lien
#  October 2014, modified,    Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#
#  Usage:
#    init_perturb [STIME S_PATH]
#
#  STIME   Initial time of the ensemble (format: YYYYMMDDHHMMSS)
#  S_PATH  Source of the initial condition including the basename.
#
#  Use settings:
#    config.main
#
#===============================================================================

cd "$(dirname "$0")"
wkdir=$(pwd)
myname=$(basename "$0")

cd ..

#===============================================================================
# Configuration

. config.main.${PRESET}
res=$? && ((res != 0)) && exit $res

. src/func_datetime.sh
. src/func_util.sh

TMPS=$DIR/tmp/init_perturb

#-------------------------------------------------------------------------------

USAGE="
[$myname] Prepare an initial ensemble by perturbing an initial condition.

Usage: $myname [STIME S_PATH]

  STIME     Initial time of the ensemble (format: YYYYMMDDHHMMSS)
  S_PATH    Source of the initial condition including the basename.
  TIMELABEL Add SCALE timelabel ? (optional : default=0)

"

#-------------------------------------------------------------------------------

if [ "$1" == '-h' ] || [ "$1" == '--help' ]; then
  echo "$USAGE"
  exit 0
fi
if (($# < 2)); then
  echo "$USAGE" >&2
  exit 1
fi

STIME=$(datetime $1)
S_PATH="$2"
TIMELABEL=${3:-0}

if [ $TIMELABEL == 1 ] ; then
  tlabel="_"$(datetime_scale $STIME)
else
  tlabel=""
fi

#-------------------------------------------------------------------------------

cd $wkdir
if [ "${S_PATH:0:1}" != "/" ]; then
  S_PATH=$(cd $(dirname $S_PATH) ; pwd)/$(basename $S_PATH)
fi

if [ ! -s "${S_PATH}${SCALE_SFX_0}" ] ; then
  echo "[Error] $0: Cannot find scale file '$S_PATH${SCALE_SFX_0}'" >&2
  exit 1
fi

#===============================================================================

if [ "$PRESET" == "FUGAKU" ] ; then
  if [ "$(which pip)" == "" ]; then
    echo "loading py-pip..."
    spack load --first py-pip%gcc
  fi 
  if [ "$(pip list | grep netCDF4)" == "" ] || [ "$(pip list | grep numpy)" == ""  ] ; then
    echo "check python packages."
    exit 1
  fi 
  PYTHON="python"
else
  echo "PRESET $PRESET not supported."
  exit 1
fi

#===============================================================================

echo
echo "Prepare output directory..."

#create_outdir

#===============================================================================

echo
echo "Prepare initial members..."

safe_init_tmpdir $TMPS
cd $TMPS

S_basename="$(basename $S_PATH)"
for m in $(seq $MEMBER); do
  mem=$(printf $MEMBER_FMT $m)
  echo "  member $mem"

  cp -f $S_PATH*.nc .
  $PYTHON $SCRP_DIR/init_perturb/init_perturb.py $S_basename
  res=$? && ((res != 0)) && exit $res

  q=0
  while [ -s "$S_basename$(printf $SCALE_SFX $q)" ]; do
    mkdir -p $OUTDIR/$STIME/anal/${mem}
    mv -f "$S_basename$(printf $SCALE_SFX $q)" $OUTDIR/$STIME/anal/${mem}/init${tlabel}$(printf $SCALE_SFX $q)
  q=$((q+1))
  done
done

### mdet 
if [ $DET_RUN == 1 ];then
  mem='mdet'
  echo "  member $mem"
  cp -f $S_PATH*.nc .
  q=0
  while [ -s "$S_basename$(printf $SCALE_SFX $q)" ]; do
    mkdir -p $OUTDIR/$STIME/anal/${mem}
    cp -f $S_basename$(printf $SCALE_SFX $q) $OUTDIR/$STIME/anal/${mem}/init${tlabel}$(printf $SCALE_SFX $q)
    q=$((q+1))
  done
fi

### mean (initial mean = mdet)
mem='mean'
echo "  member $mem"
cp -f $S_PATH*.nc .
q=0
while [ -s "$S_basename$(printf $SCALE_SFX $q)" ]; do
  mkdir -p $OUTDIR/$STIME/anal/${mem}
  cp -f $S_basename$(printf $SCALE_SFX $q) $OUTDIR/$STIME/anal/${mem}/init${tlabel}$(printf $SCALE_SFX $q)
  q=$((q+1))
done



### safe_rm_tmpdir $TMPS

#===============================================================================

echo

exit 0
