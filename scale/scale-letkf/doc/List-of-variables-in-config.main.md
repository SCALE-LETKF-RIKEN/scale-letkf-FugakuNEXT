Last update: 2023.8.6

| Variable | Default value | Explanation |
| --- | --- | --- |
| PRESET |  | Pre-determined combination of settings for particular system environments<ul><li>'FUGAKU': &nbsp;Supercomputer Fugaku</li><li>'Linux_torque': &nbsp;Linux cluster with a Torque job scheduling system (hibuna)</li></ul> |
| RUN_LEVEL |  0 | [SHOULD NOT SET BY USER] Level to control how many parts of the code will be run<ul><li>0: &nbsp;Run everything</li><li>1: &nbsp;Run everything but skipping detecting errors</li><li>2: &nbsp;[NOT IMPLEMENTED] Staging list files are ready; skip generating them</li><li>3: &nbsp;Staging-in has been done; skip staging</li><li>4: &nbsp;Staging-in/out is done outside this script; skiping staging</li></ul> |

***

### Location of the root and input/output directories

| Variable | Default value | Explanation |
| --- | --- | --- |
| DIR |  | Root directory of the SCALE-LETKF source code |
| INDIR |  $OUTDIR | SCALE-LETKF experiment directory for input data<br>(The initial condition data need to exist in this directory unless $MAKEINIT = 1 in config.cycle or config.fcst)<br> |
| OUTDIR |  | SCALE-LETKF experiment directory for output data |

***

### Location of model/data files

| Variable | Default value | Explanation |
| --- | --- | --- |
| SCALEDIR |  | Directory of the SCALE model source code |
| DATADIR |  | Directory of the SCALE database |
| DATA_TOPO |  $INDIR | SCALE-LETKF experiment directory under which prepared topo files exist<br>(effective only when $TOPO_FORMAT = 'prep')  |
| DATA_LANDUSE |  $INDIR | SCALE-LETKF experiment directory under which prepared landuse files exist<br>(effective only when $LANDUSE_FORMAT = 'prep')  |
| DATA_BDY_SCALE_PREP |  $INDIR | SCALE-LETKF experiment directory under which prepared SCALE boundary files exist<br>(effective only when $BDY_FORMAT = 0)  |
| DATA_BDY_SCALE |  | Parent domain's SCALE-LETKF experiment directory for boundary data<br>(effective only when $BDY_FORMAT = 1: **offline-nesting run**) |
| DATA_TOPO_BDY_SCALE |  | Directory of the parent domain's topo files<br>(effective only when $BDY_FORMAT = 1: **offline-nesting run**) |
| DATA_BDY_WRF |  | Directory of boundary data in WRF format<br>(effective only when $BDY_FORMAT = 2) |
| DATA_BDY_GRADS |  | Directory of boundary data in GrADS format<br>(effective only when $BDY_FORMAT = 4) |
| DATA_ADDINFL |  $INDIR | Directory of additive inflation files<br>(effective only when $ADDINFL = 1 in config.cycle)  |
| OBS |  | Directory of observation data files|

***

### Model/data file options

| Variable | Default value | Explanation |
| --- | --- | --- |
| PNETCDF |  0 | Use PnetCDF for single-file I/O?<br>**always set this to 0**<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul>  |
| PNETCDF_BDY_SCALE |  $PNETCDF | Used PnetCDF I/O in the parent domain's SCALE-LETKF experiment?<br>(effective only when $BDY_FORMAT = 1)<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
| DET_RUN |  0 | Enable the deterministic run (Schraff et al. 2016 QJRMS)?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
| TOPO_FORMAT |  | Topography data file format<ul><li>'prep': &nbsp;Use prepared topo files in $DATA_TOPO; this option can be used to save computaitonal time after completing a previous run with the same domain settings and with $TOPOOUT_OPT = 1 in 'config.cycle' or 'config.fcst'</li><li>'GTOPO30': &nbsp;Create topo files using 'GTOPO30' dataset (requires compatible 'config.nml.scale_pp')</li><li>'DEM50M': &nbsp;Create topo files using 'DEM50M' dataset (requires compatible 'config.nml.scale_pp')</li></ul> |
| LANDUSE_FORMAT |  | Land-use data file format<ul><li>'prep': &nbsp;Use prepared landuse files in $DATA_LANDUSE; this option can be used to save computaitonal time after completing a previous run with the same domain settings and with $LANDUSEOUT_OPT = 1 in 'config.cycle' or 'config.fcst'</li><li>'GLCCv2': &nbsp;Create landuse files using 'GLCCv2' dataset (requires compatible 'config.nml.scale_pp')</li><li>'LU100M': &nbsp;Create landuse files using 'LU100M' dataset (requires compatible 'config.nml.scale_pp')</li></ul> |
| LANDUSE_UPDATE |  0 | Use time-variant landuse files<ul><li>0: &nbsp;No (time-invariant landuse files)</li><li>1: &nbsp;Yes</li></ul> |
| BDY_FORMAT |  | Boundary data format<ul><li>0: &nbsp;Use prepared SCALE boundary files with exactly the same domain settings; no additional pre-processing is required (**'scale_init' step will be skipped**); this option can be used to save computaitonal time after completing a previous run with the same domain and cycling/forecast settings and with $BDYOUT_OPT = 1 in 'config.cycle' or 'config.fcst'</li><li>1: &nbsp;SCALE history file: offline-nesting run (requires compatible 'config.nml.scale_init')</li><li>2: &nbsp;WRF (requires compatible 'config.nml.scale_init')</li><li>3: &nbsp;NICAM (requires compatible 'config.nml.scale_init')</li><li>4: &nbsp;GrADS (requires compatible 'config.nml.scale_init')</li><li>5: &nbsp;Idealized experiment without boundary files (**'scale_init' step will be skipped**) </li></ul> |
| BDY_SINGLE_FILE |  0 | Boundary data for each cycle are contained in a single file of a sufficient length?<ul><li>0: &nbsp;No - Assume the length of a single boundary file = $BDYCYCLE_INT which may be shorter than the required forecast length; more than one boundary files may be used for each cycle (e.g., files made by data assimilation cycles)</li><li>1: &nbsp;Yes - Assume the length of a single boundary file >= the required forecast length (i.e., $WINDOW_E for 'cycle' job; $FCSTLEN for 'fcst' job); always only one single boundary file is used for each cycle</li></ul> |
| BDY_SCALE_DIR |  'hist' | Sub-directory name of the SCALE history files in the parent domain's **SCALE-LETKF experiment directory**<br>(effective only when $BDY_FORMAT = 1: **offline-nesting run**)<br> |
| BDY_MEAN |  'mean' | Sub-directory name representing the ensemble mean in the boundary data directory<br> |
| BDY_ENS |  0 | Use ensemble boundary conditions?<ul><li>0: &nbsp;No - Use a fixed boundary condition (under $BDY_MEAN) for all memebers</li><li>1: &nbsp;Yes - Use ensemble boundary conditions</li></ul> |
| BDY_ROTATING |  0 | Use different series of boundary data for different cycles?<ul><li>0: &nbsp;No - Use the same series of boundary data files for all cycles</li><li>1: &nbsp;Yes - Use different series of boundary data files (which can overlap in time) for different cycles</li></ul>See **explanation on "BDY_ROTATING"** for more details.<br>Note that this setting affects the file naming convention in the $DATA_BDY_* boundary data directory; see **File naming convention of data directories**.<br> |
| BDYINT |  $LCYCLE | Time interval (second) of each time frame in the boundary data with multiple time frames in a file.<br>If there is only a single time frame in the boundary data files, set this value to be equal to $BDYCYCLE_INT, the time interval of each boundary data file.<br>In the **offline-nesting run**, this should be equal to the time interval of the history files in the parent domain's experiment.<br> |
| BDYCYCLE_INT |  $BDYINT | Time interval (second) of each boundary data file (regardless of the multiple time frames in a file).<br>In the **offline-nesting run**, this should be equal to $LCYCLE in the parent domain's experiment.<br> |
| PARENT_REF_TIME |  | [OPTIONAL] Any time (format: YYYYMMDDhhmmss) at which a boundary data file exists (for boundary data with multiple time frames in a file, use the start time in a file).<br>This is used as a reference time to search the boundary data files; the program will search available boundary data files by plus/minus increments of $BDYCYCLE_INT based on this reference time.<br>This setting is only important when $LCYCLE < $BDYCYCLE_INT, likely in a **offline-nesing run**. In the case that $LCYCLE >= $BDYCYCLE_INT, this setting can be left blank. |
| ENABLE_PARAM_USER |  0 | Use '&PARAM_USER' namelist section in the 'config.nml.scale_user' configuration file?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes - This setting can only be used with a customized version of SCALE for SCALE-LETKF and it requires 'config.nml.scale_user' configuration file</li></ul> |
| OCEAN_INPUT |  0 | In every cycle, update ocean variables from the ocean input files?<ul><li>0: &nbsp;No - Cycle the ocean variables same as for atmospheric variables</li><li>1: &nbsp;Yes - The list of ocean variables to be updated can be controlled in 'config.nml.scale_user' if using a customized version of SCALE for SCALE-LETKF; otherwise, all ocean variables will be updated</li></ul> |
| OCEAN_FORMAT |  99 | The source of the ocean input files<br>(effective only when $OCEAN_INPUT = 1)<ul><li>0: &nbsp;Use prepared SCALE restart (init) files with exactly the same domain settings; no additional pre-processing is required</li><li>99: &nbsp;Use the same data source of the boundary data (based on $BDY_FORMAT setting)</li></ul> |
| LAND_INPUT |  0 | In every cycle, update land variables from the land input files?<ul><li>0: &nbsp;No - Cycle the land variables same as for atmospheric variables</li><li>1: &nbsp;Yes - The list of land variables to be updated can be controlled in 'config.nml.scale_user' if using a customized version of SCALE for SCALE-LETKF; otherwise, all land variables will be updated</li></ul> |
| LAND_FORMAT |  99 | The source of the land input files<br>(effective only when $LAND_INPUT = 1)<ul><li>0: &nbsp;Use prepared SCALE restart (init) files with exactly the same domain settings; no additional pre-processing is required</li><li>99: &nbsp;Use the same data source of the boundary data (based on $BDY_FORMAT setting)</li></ul> |
| OBSNUM |  | Number of observation data files (in $OBS directory) to be assimilated at the same time |
| OBSNAME[1] |  | Array of the file name prefix for each observation data file [1..$OBSNUM]; see File naming convention of data directories |
| OBSOPE_SEPARATE[1] |  0 | Array of values indicating whether each observation data file should be processed using a separate observation operator program [1..$OBSNUM]<ul><li>0: &nbsp;No - Use built-in observation operator in the LETKF program</li><li>1: &nbsp;Yes - Use a separate observation operator program</li></ul> |
| EFSO_RUN |  0 | Enable the EFSO?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |

***

### Cycling settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| WINDOW_S |  | SCALE forecast time when the assimilation window starts (second)<br>(Set this value to be equal to $LCYCLE for 3D-LETKF) |
| WINDOW_E |  | SCALE forecast time when the assimilation window ends (second)<br>(Set this value to be equal to $LCYCLE for 3D-LETKF) |
| LCYCLE |  | Length of a data assimilation cycle (second) |
| LTIMESLOT |  | Time slot interval for 4D-LETKF (second)<br>(Set this value to be equal to $LCYCLE for 3D-LETKF) |

***

### Parallelization settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| MEMBER |  | Ensemble size |
| NNODES |  | Number of computing nodes used to run the job |
| PPN |  | Number of processes per node (parallelized by MPI with $NNODES*$PPN processes) |
| THREADS |  | Number of threads per MPI process (parallelized by OpenMP or automatic parallelization with this number of threads) |
| SCALE_NP_X |  | Number of subdomains in X direction to run a single member of the SCALE forecast |
| SCALE_NP_Y |  | Number of subdomains in Y direction to run a single member of the SCALE forecast |
| SCALE_NP   |  | SCALE_NP_X \* SCALE_NP_Y = Number of MPI processes to run a single member of the SCALE forecast |

***

### Use of the temporary runtime directories

| Variable | Default value | Explanation |
| --- | --- | --- |
| ONLINE_STAGEOUT |  0 | **Currently not supported** Stage out right after each cycle (do not wait until the end of the job)?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes </li></ul> |
| DISK_ALL_RANK_LOCAL |  0 | [USUALLY SET BY PRESET] Local disks for each MPI process in the same node are independent?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes - Choose this when using the rank-directory on the K computer</li></ul> |
| DISK_MODE |  2 | Type of disk and the way of using the disk for general input/output files (except for the database, boundary, and observation files) in runtime<ul><li>1: &nbsp;Use a shared disk; only create symbolic links of the input/output files to $TMP directory for computation</li><li>2: &nbsp;Use a shared disk; copy the input files to $TMP directory for computation and copy the output files back (staging)</li><li>3: &nbsp;Use local disks; copy the input files to $TMPL directory for computation and copy the output files back (staging)</li></ul> |
| DISK_MODE_CONSTDB |  $DISK_MODE | Type of disk and the way of using the disk for database (topo/landuse) files in runtime<ul><li>1-3: &nbsp;Same meaning as for $DISK_MODE</li></ul> |
| DISK_MODE_BDYDATA |  $DISK_MODE | Type of disk and the way of using the disk for boundary data files in runtime<ul><li>1-3: &nbsp;Same meaning as for $DISK_MODE</li></ul> |
| DISK_MODE_OBS |  $DISK_MODE | Type of disk and the way of using the disk for observation data files in runtime<ul><li>1-3: &nbsp;Same meaning as for $DISK_MODE</li></ul> |
| SYSNAME |  | A string to be used as a part of the temporary directory name |
| TMPSUBDIR |  | A part of the temporary directory name (should be unique in the machine to avoid any potential conflict when running multiple jobs at the same time) |
| TMP |  | Temporary runtime directory path on a shared disk (available on both the head node and all computing nodes)<br>(On the K computer, effective only for the 'micro' job) |
| TMPS |  | Temporary runtime directory path available only on the head node |
| TMPSL |  $TMPS | [OPTIONAL] Temporary runtime directory path available only on the head node and on a local file system<br> |
| TMPL |  | Temporary runtime directory path on local disks (available on all computing nodes; do not need to be available on the head node)<br>(On the K computer, do not need to set this) |
| CLEAR_TMP |  1 | Clear the temporary directories after the completion of the job?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
| USE_LLIO_BIN |  1 | Use LLIO-transfer for binary files? (effective only on Fugaku)<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
| USE_LLIO_DAT |  1 | Use LLIO-transfer for shared data files? (effective only on Fugaku)<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
| USE_SPACK |  1 | Use spack in job scripts to set library paths? (effective only on Fugaku)<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
***

### Environmental settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| MPI_TYPE |  | [USUALLY SET BY PRESET] Type of MPI software<ul><li>'openmpi' &nbsp;Open MPI</li><li>'impi' &nbsp;Intel MPI</li><li>'sgimpt': &nbsp;SGI MPT</li></ul> |
| MPIRUN |  | Path of the 'mpiexec', 'mpirun', or 'impijob' program |
| SCP |  | Command (on computing nodes) used for built-in file staging (to copy runtime files from input/to output directories) |
| SCP_HOSTPREFIX |  | [OPTIONAL] When using scp-like commands for built-in file staging, the prefix for the hostname of the head node |
| STAGE_THREAD |  | Number of threads for built-in parallel file staging |
| TAR_THREAD |  | Number of threads for parallel archiving log files after the completion of the job<br>(effective only when $LOG_TYPE > 3) |
| PYTHON |  | Command for the 'python' program |

***
