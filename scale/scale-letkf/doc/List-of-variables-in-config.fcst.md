Last update: 2023.8.6

| Variable | Default value | Explanation |
| --- | --- | --- |
| STIME |  | Forecast initial time of the first cycle in the job (format: YYYYMMDDhhmmss) |
| ETIME |  | Forecast initial time of the last cycle in the job (format: YYYYMMDDhhmmss) |
| MEMBERS |  | List of members to run the ensemble forecasts, separated by space; e.g., 'mean mdet 1 5', or "$(seq 11 20)" |
| CYCLE |  | [OPTIONAL] Number of cycles to be run in parallel (if possible).<br>This value can be left blank and the program will decide the optimal setting automatically. |
| CYCLE_SKIP |  1 | Interval of cycles at every which the (ensemble) forecasts are conducted<br>($CYCLE_SKIP = 1: run forecasts from every cycle)<br> |
| CONF_MODE |   | Method to generate namelist files for each model/DA program<br>**Currently only the 'static' mode is supported.**<ul><li>'dynamic': &nbsp;Dynamically (and parallelly) generate the namelist files during the job runtime</li><li>'static': &nbsp;Generate all static namelist files on the head node before submitting the computing job. Note that the support to various functions is limited using this mode.</li></ul> |
| TIME_LIMIT |  | The requested "walltime" limit of the computing job |

***

### General settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| MAKEINIT |  | Generate (ensemble) initial conditions as well from the boundary data?<ul><li>0: &nbsp;No - Use the existing initial condition data in $INDIR</li><li>1: &nbsp;Yes - The boundary data in the $DATA_BDY_* directory are not only used for generating boundary conditions, but the same data are also used for generating initial conditions. The results will be then saved (overwritten) to $OUTDIR after the completion of the job.</li></ul> |
| FCSTLEN |  | Forecast length (second) |
| FCSTOUT |  | Output time interval of the forecasts history file (second) |
| RESTARTOUT |  | Output time interval of the forecasts restart file (second) |

***

### Completeness of output files

| Variable | Default value | Explanation |
| --- | --- | --- |
| OUT_OPT |  | Level of completeness of forecast files<br><pre> &nbsp; &nbsp; fcst<br>-- &nbsp; history restart#<br>1: &nbsp; o &nbsp; &nbsp; &nbsp; o<br>2: &nbsp; o<br># = restart file in the end of the forecast</pre> |
| TOPOOUT_OPT |  | Level of completeness of topography files<br><pre>-- &nbsp; topo<br>1: &nbsp; o<br>2: &nbsp; (none)</pre> * The topography file output can be used in a later run with the same domain settings. |
| LANDUSEOUT_OPT |  | Level of completeness of land-use files<br><pre>-- &nbsp; landuse<br>1: &nbsp; o<br>2: &nbsp; (none)</pre> * The land-use file output can be used in a later run with the same domain settings. |
| BDYOUT_OPT |  | Level of completeness of SCALE boundary files<br><pre> &nbsp; &nbsp; bdy<br>-- &nbsp; mean members<br>1: &nbsp; o &nbsp; &nbsp;o<br>2: &nbsp; o<br>3: &nbsp; (none)</pre> * The SCALE boundary file output can be used in a later run with the same domain and forecast settings. |
| LOG_OPT |  | Level of completeness of log files of each program<br><pre> &nbsp; &nbsp; topo landuse bdy scale<br>-- &nbsp; log &nbsp;log &nbsp; &nbsp; log log<br>1-2: o &nbsp; &nbsp;o &nbsp; &nbsp; &nbsp; o &nbsp; o<br>3: &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;o<br>4: &nbsp; (none)</pre> |
| LOG_TYPE |  | Completeness of log file output and the way of saving them<ul><li>1: &nbsp;Only save the log file from the head computing process</li><li>2: &nbsp;Save all log files from all MPI processes</li><li>3: &nbsp;Save all log files from all MPI processes in an archive file (\*.tar)</li><li>4: &nbsp;Save all log files from all MPI processes in a compressed archive file (\*.tar.gz)</li></ul> |

***
