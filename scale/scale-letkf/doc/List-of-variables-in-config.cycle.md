Last update: 2023.8.6

| Variable | Default value | Explanation |
| --- | --- | --- |
| STIME |  | Forecast start time of the first cycle in the job (format: YYYYMMDDhhmmss) |
| ETIME |  | Forecast start time of the last cycle in the job (format: YYYYMMDDhhmmss)<br>(final analysis time: $ETIME + $LCYCLE) |
| ISTEP | 1 | Initial step number of the job within a cycle (1-5)<ul><li>1: &nbsp;scale-rm_pp_ens</li><li>2: &nbsp;scale-rm_init_ens</li><li>3: &nbsp;scale-rm_ens</li><li>4: &nbsp;obsope/obsmake</li><li>5: &nbsp;letkf</li><li>6: &nbsp;EFSO</li></ul>|
| FSTEP | 6 | Final step number of the job within a cycle (1-6)|
| CONF_MODE |   | Method to generate namelist files for each model/DA program<br>**Currently only the 'static' mode is supported.**<ul><li>'dynamic': &nbsp;Dynamically (and parallelly) generate the namelist files during the job runtime</li><li>'static': &nbsp;Generate all static namelist files on the head node before submitting the computing job. Note that the support to various functions is limited using this mode.</li></ul> |
| TIME_LIMIT |  | The requested "walltime" limit of the computing job |

***

### General settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| MAKEINIT |  | Generate (ensemble) initial conditions as well from the boundary data?<ul><li>0: &nbsp;No - Use the existing initial condition data in $INDIR</li><li>1: &nbsp;Yes - The boundary data in the $DATA_BDY_* directory are not only used for generating boundary conditions, but the same data are also used for generating initial conditions. The results will be then saved (overwritten) to $OUTDIR after the completion of the job.</li></ul> |
| FCSTOUT |  | Output time interval of the history files (second).<br>Note that the actual time interval in the history files is min{$LTIMESLOT (in config.main), $FCSTOUT}. |
| ADAPTINFL |  | [NOT WELL TESTED] Enable adaptive multiplicative inflation (Miyoshi 2011 MWR)?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
| ADDINFL |  | Enable additive inflation?<br>(require input additive inflation data in $DATA_ADDINFL)<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |

***

### Completeness of output files

| Variable | Default value | Explanation |
| --- | --- | --- |
| OUT_OPT |  | Choose whether to save the guess/analysis/history files<br><pre> &nbsp; &nbsp; anal &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;gues &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;hist<br>-- &nbsp; mean/ members mean/ members mean members<br>-- &nbsp; sprd &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;sprd<br>1: &nbsp; o &nbsp; &nbsp; o &nbsp; &nbsp; &nbsp; o &nbsp; &nbsp; o &nbsp; &nbsp; &nbsp; o &nbsp; &nbsp;o<br>2: &nbsp; o &nbsp; &nbsp; o &nbsp; &nbsp; &nbsp; o &nbsp; &nbsp; o &nbsp; &nbsp; &nbsp; o<br>3: &nbsp; o &nbsp; &nbsp; o &nbsp; &nbsp; &nbsp; o &nbsp; &nbsp; o<br>4: &nbsp; o &nbsp; &nbsp; o &nbsp; &nbsp; &nbsp; o<br>5: &nbsp; o &nbsp; &nbsp; # &nbsp; &nbsp; &nbsp; o<br>6: &nbsp; o &nbsp; &nbsp; ## &nbsp; &nbsp; &nbsp;o<br>7: &nbsp; o &nbsp; &nbsp; ##<br># = every $OUT_CYCLE_SKIP cycles and the last cycle<br>## = only the last cycle</pre> * The history file output can be used to drive a child domain for **offline nesting**. |
| OUT_CYCLE_SKIP |  | Interval of cycles when the analysis ensemble members are output<br>(effective only when $OUT_OPT = 5) |
| TOPOOUT_OPT |  | Choose whether to save the topography files<br><pre>-- &nbsp; topo<br>1: &nbsp; o<br>2: &nbsp; (none)</pre> * The topography file output can be used in a later run with the same domain settings. |
| LANDUSEOUT_OPT |  | Choose whether to save the land-use files<br><pre>-- &nbsp; landuse<br>1: &nbsp; o<br>2: &nbsp; (none)</pre> * The land-use file output can be used in a later run with the same domain settings. |
| BDYOUT_OPT |  | Choose whether to save the SCALE boundary files<br><pre> &nbsp; &nbsp; bdy<br>-- &nbsp; mean members<br>1: &nbsp; o &nbsp; &nbsp;o<br>2: &nbsp; o<br>3: &nbsp; (none)</pre> * The SCALE boundary file output can be used in a later run with the same domain and cycling settings. |
| OBSOUT_OPT |  | Level of completeness of observation-space diagnostic files<br><pre>-- &nbsp; obsgues# obsanal# obsdep<br>1: &nbsp; o &nbsp; &nbsp; &nbsp; &nbsp;o &nbsp; &nbsp; &nbsp; &nbsp;o<br>2: &nbsp; o &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; o<br>3: &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; o<br>4: &nbsp; (none)<br># = [NOT IMPLEMENTED]</pre> |
| LOG_OPT |  | Level of completeness of log files of each program<br><pre> &nbsp; &nbsp; topo landuse bdy scale obsope letkf<br>-- &nbsp; log &nbsp;log &nbsp; &nbsp; log log &nbsp; log &nbsp; &nbsp;log<br>1-2: o &nbsp; &nbsp;o &nbsp; &nbsp; &nbsp; o &nbsp; o &nbsp; &nbsp; o &nbsp; &nbsp; &nbsp;o<br>3: &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;o &nbsp; &nbsp; o &nbsp; &nbsp; &nbsp;o<br>4: &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;o &nbsp; &nbsp; &nbsp;o<br>5: &nbsp; (none)</pre> |
| LOG_TYPE |  | Completeness of log file output and the way of saving them<ul><li>1: &nbsp;Only save the log file from the head computing process</li><li>2: &nbsp;Save all log files from all MPI processes</li><li>3: &nbsp;Save all log files from all MPI processes in an archive file (\*.tar)</li><li>4: &nbsp;Save all log files from all MPI processes in a compressed archive file (\*.tar.gz)</li></ul> |
| SPRD_OUT |  | Output the ensemble spread file?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
| RTPS_INFL_OUT |  | Output the equivalent multiplicative inflation field when using RTPS inflation?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |
| NOBS_OUT |  | Output the observation number field when using observation number limit?<ul><li>0: &nbsp;No</li><li>1: &nbsp;Yes</li></ul> |

***
