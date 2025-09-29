Last update: 2023.8.6

__Note:__ Variables that are in italic and also have an exclamation mark (__!__) in front in the following tables are those to be automatically determined by the job scripts. Their values __should not be assigned by users__ and they should be written in `config.nml.letkf` as:
```
!--VARIABLE--
```

#### __Special notations__ that can be used in some namelist variables for filename patterns:
- Special notation for member strings: '<member>'

***

#### &PARAM_LOG

Settings related to log files from LETKF

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__LOG_LEVEL__</sub> | 1 | Log message output level: <ul><li>0: &nbsp;Minimum log output</li><li>1: &nbsp;Reduced log output</li><li>2: &nbsp;Normal log output</li><li>3: &nbsp;Verbose log output</li></ul>|
| <sub>__USE_MPI_BARRIER__</sub> | .true. | Enable some MPI_Barrier for better timing measurement? |
| <sub>__LOG_ALL_PRC__</sub> | .false. | Output from all MPI processes? |

***
#### &PARAM_MODEL

Model-related settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__MODEL__</sub> | 'scale-rm' | Model name (always set to 'scale-rm') |
| <sub>__VERIFY_COORD__</sub> | .false. | Verify the vertical coordinate settings with the vertical coordinate values in the input file?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |

***


#### &PARAM_OBSOPE

Observation operator settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__!__*OBS_IN_NUM*</sub> | 1 | Number of input observation data files |
| <sub>__!__*OBS_IN_NAME<br>&nbsp;&nbsp;&nbsp;(nobsfilemax)*</sub> | 'obs.dat' | Array of filenames of each observation file |
| <sub>__OBS_IN_FORMAT<br>&nbsp;&nbsp;&nbsp;(nobsfilemax)__</sub> | 'PREPBUFR' | Array of data format of each observation file<ul><li>'PREPBUFR': &nbsp;Conventional data in LETKF format</li><li>'RADAR': &nbsp;Radar data in LETKF format</li><li>'HIMAWARI8': &nbsp;Himawari-8 data in LETKF format</li></ul> |
| <sub>__!__*OBSDA_RUN<br>&nbsp;&nbsp;&nbsp;(nobsfilemax)*</sub> | .true. | Array setting whether to run observation operator for each observation file<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__OBSDA_OUT__</sub> | .false. | Output observation operator results [i.e., H(x^b)] to files?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes - Can be used as a separate observation operator program</li></ul> |
| <sub>__!__*OBSDA_OUT_BASENAME*</sub> | 'obsda.@@@@' | Base filename pattern of observation operator result outputs ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used) |
| <sub>__!__*HISTORY_IN_BASENAME*</sub> | 'hist.@@@@' | Base filename pattern of input history files for observation operator calculation ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used) |
| <sub>__!__*SLOT_START*</sub> | 1 | Start time-slot number for 4D-LETKF |
| <sub>__!__*SLOT_END*</sub> | 1 | End time-slot number for 4D-LETKF |
| <sub>__!__*SLOT_BASE*</sub> | 1 | Base time-slot number for 4D-LETKF |
| <sub>__!__*SLOT_TINTERVAL*</sub> | 3600.0d0 | Time-slot interval for 4D-LETKF |

***

#### &PARAM_LETKF

General LETKF settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__!__*OBSDA_IN*</sub> | .false. | Input observation operator results [i.e., H(x^b)] from a separate observation operator program |
| <sub>__!__*OBSDA_IN_BASENAME*</sub> | 'obsda.@@@@' | Base filename pattern of input observation operator results from a separate observation operator program ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used) |
| <sub>__!__*GUES_IN_BASENAME*</sub> | 'gues.@@@@' | Base filename pattern of input first-guess files ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used) |
| <sub>__!__*GUES_MEAN_INOUT_BASENAME*</sub> | '' | Base filename of first-guess files for the ensemble mean (may be used as both input and output) |
| <sub>__!__*GUES_SPRD_OUT_BASENAME*</sub> | '' | Base filename of output first-guess ensemble spread |
| <sub>__!__*GUES_SPRD_OUT*</sub> | .true. | Output first-guess ensemble spread?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__!__*ANAL_OUT_BASENAME*</sub> | 'anal.@@@@' | Base filename pattern of output analysis files ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used) |
| <sub>__!__*ANAL_SPRD_OUT*</sub> | .true. | Output analysis ensemble spread?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__!__*LETKF_TOPOGRAPHY_IN_BASENAME*</sub> | 'topo' | Base filename of the input topographic file |
| <sub>__INFL_MUL__</sub> | 1.0d0 | Multiplicative covariance inflation parameter<ul><li>1: &nbsp;Disable the multiplicative inflation</li><li>> 0: &nbsp;Use a global constant inflation parameter</li><li><= 0: &nbsp;Use a 3D inflation field from INFL_MUL_IN_BASENAME file</li></ul> |
| <sub>__INFL_MUL_MIN__</sub> | -1.0d0 | Minimum multiplicative inlfation parameter<ul><li><= 0: &nbsp;No minimum setting</li></ul> |
| <sub>__!__*INFL_MUL_ADAPTIVE*</sub> | .false. | Output adaptively estimated 3D inflation field to INFL_MUL_OUT_BASENAME file?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__!__*INFL_MUL_IN_BASENAME*</sub> | 'infl' | Base filename of input 3D inflation field |
| <sub>__!__*INFL_MUL_OUT_BASENAME*</sub> | 'infl' | Base filename of output (adaptively estimated) 3D inflation field |
| <sub>__INFL_ADD__</sub> | 0.0d0 | Additive covariance inflation parameter; this value will be multiplied to the input additive inflation field when using additive inflation<ul><li>< 0: &nbsp;Disable the additive inflation</li></ul> |
| <sub>__!__*INFL_ADD_IN_BASENAME*</sub> | 'addi.@@@@' | Base filename pattern of the input additive inflation field ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used) |
| <sub>__INFL_ADD_SHUFFLE__</sub> | .false. | Shuffle the ensemble members for additive inflation field?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__INFL_ADD_Q_RATIO__</sub> | .false. | For moisture field, further multiply the additive inflation field by the first-guess ensemble mean values (i.e., use the additive inflation field as ratio)?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__INFL_ADD_REF_ONLY__</sub> | .false. | Apply the additive inflation only around where raining reflectivity (> RADAR_REF_THRES_DBZ) observations exist?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__RELAX_ALPHA__</sub> | 0.0d0 | Relaxation-to-prior-perturbation (RTPP) parameter (Zhang et al. 2004 MWR) |
| <sub>__RELAX_ALPHA_SPREAD__</sub> | 0.0d0 | Relaxation-to-prior-spread (RTPS) parameter (Whitaker and Hamill 2012 MWR) |
| <sub>__!__*RELAX_SPREAD_OUT*</sub> | .false. | Output the equivalent multiplicative inflation field when using RTPS?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__!__*RELAX_SPREAD_OUT_BASENAME*</sub> | 'rtps' | Base filename of equivalent multiplicative inflation field output (when using RTPS) |
| <sub>__RELAX_TO_INFLATED_PRIOR__</sub> | .false. | Choice of using covariance relaxation and multiplicative inflation together<ul><li>.true.: &nbsp;Relaxation to the prior after the multiplicative inflation</li><li>.false.: &nbsp;Relaxation to the original prior before the multiplicative inflation</li></ul> |
| <sub>__GROSS_ERROR__</sub> | 5.0d0 | Threshold of gross error check (times of observation errors) |
| <sub>__GROSS_ERROR_RADAR_REF__</sub> | -1.0d0 | Threshold of gross error check for radar reflectivity data<ul><li>0: &nbsp;Same as GROSS_ERROR</li></ul> |
| <sub>__GROSS_ERROR_RADAR_VR__</sub> | -1.0d0 | Threshold of gross error check for radar radial velocity data<ul><li>0: &nbsp;Same as GROSS_ERROR</li></ul> |
| <sub>__Q_UPDATE_TOP__</sub> | 0.0d0 | Pressure level (Pa) only below which water vapor and hydrometeors are updated |
| <sub>__Q_SPRD_MAX__</sub> | -1.0D0 | Maximum ratio of ensemble spread to ensemble mean for mositure in the analysis; if the analysis ensemble spread is greater than this ratio, scale the ensemble perturbations to reduce the spread to this ratio<ul><li><= 0: &nbsp;Disabled</li></ul> |
| <sub>__BOUNDARY_BUFFER_WIDTH__</sub> | 0.0d0 | Width (m) of the buffer area along the lateral boundary where the analysis increment is gradually reduced to zero |
| <sub>__POSITIVE_DEFINITE_Q__</sub> | .false. | Force setting the negative values in the analysis water vapor field to zero?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__POSITIVE_DEFINITE_QHYD__</sub> | .false. | Force setting the negative values in the analysis hydrometeor fields to zero?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__PS_ADJUST_THRES__</sub> | 100.d0 | Threshold of elevation difference (m) between the station report and the model topography<br>Within the threshold surface pressure observations are assimilated (height adjustment will be performed to compensate this difference); beyond this threshold the surface pressure observations are discarded |
| <sub>__!__*NOBS_OUT*</sub> | .false. | Output the field of actual observation numbers assimilated in each grid when the observation number limit (Hamrud et al. 2015 MWR) is enabled?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__!__*NOBS_OUT_BASENAME*</sub> | 'nobs' | Base filename of the field of actual observation numbers assimilated in each grid when the observation number limit is enabled |

***


#### &PARAM_LETKF_OBS

Observation-specific settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__USE_OBS<br>&nbsp;&nbsp;&nbsp;(nobtype)__</sub> | .true. | Array setting whether each observation report type is used?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__HORI_LOCAL<br>&nbsp;&nbsp;&nbsp;(nobtype)__</sub> | (/500.0d3, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0/) | Array of horizontal covariance localization length scale for each observation report type<ul><li>> 0: &nbsp;Horizontal localization length scale (m)</li><li>0: &nbsp;[NOT IMPLEMENTED] No horizontal localization</li><li>< 0: &nbsp;Same setting as HORI_LOCAL(1)</li></ul> |
| <sub>__VERT_LOCAL<br>&nbsp;&nbsp;&nbsp;(nobtype)__</sub> | (/ 0.4d0,   -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0,   -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, 1000.0d0, -1.0d0, -1.0d0/) | Array of vertical covariance localization length scale for each observation report type<ul><li>> 0: &nbsp;Vertical localization length scale (ln(p) or m depending on the report type - see [note on vertical coordinate](Note-on-vertical-coordinate.md) )</li><li>0: &nbsp;No vertical localization</li><li>< 0: &nbsp;Same setting as VERT_LOCAL(1)</li></ul> |
| <sub>__TIME_LOCAL<br>&nbsp;&nbsp;&nbsp;(nobtype)__</sub> | (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0/) | [NOT IMPLEMENTED] Array of temporal covariabce localization interval for each observation report type<ul><li>> 0: &nbsp;[NOT IMPLEMENTED] Temporal localization interval (sec)</li><li>0: &nbsp;No temporal localization</li><li>< 0: &nbsp;Same setting as TIME_LOCAL(1)</li></ul> |
| <sub>__HORI_LOCAL_RADAR_OBSNOREF__</sub> | -1.0d0 | Horizontal covariance localization length scale (m) for clear-sky radar reflectivity data (<= RADAR_REF_THRES_DBZ)<ul><li>< 0: &nbsp;Same setting as HORI_LOCAL(22) for all radar data</li></ul> |
| <sub>__HORI_LOCAL_RADAR_VR__</sub> | -1.0d0 | Horizontal covariance localization length scale (m) for radar radial velocity data<ul><li>< 0: &nbsp;Same setting as HORI_LOCAL(22) for all radar data</li></ul> |
| <sub>__VERT_LOCAL_RADAR_VR__</sub> | -1.0d0 | Vertical covariance localization length scale (m) for radar radial velocity data<ul><li>< 0: &nbsp;Same setting as HORI_LOCAL(22) for all radar data</li></ul> |
| <sub>__MAX_NOBS_PER_GRID<br>&nbsp;&nbsp;&nbsp;(nobtype)__</sub> | (/ 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/) | Observation number limit: maximum number of observations for each observation report type and each variable assimilated in a grid (Hamrud et al. 2015 MWR)<ul><li>> 0: &nbsp;Enable the observation number limit</li><li>0: &nbsp;No observation number limit</li><li>< 0: &nbsp;Same setting as MAX_NOBS_PER_GRID(1)</li></ul> |
| <sub>__MAX_NOBS_PER_GRID_CRITERION__</sub> | 1 | Criterion to limit the number of observations per grid <ul><li>1: &nbsp;normalized 3D distance (from closest)</li><li>2: &nbsp;localization weight (from largest)</li><li>3: &nbsp;weighted observation error variance (from smallest)</li></ul> |
| <sub>__OBS_MIN_SPACING<br>&nbsp;&nbsp;&nbsp;(nobtype)__</sub> | (/300.0d3, 100.0d3, 100.0d3, 150.0d3, 300.0d3, 150.0d3, 150.0d3, 100.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 300.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3,   1.0d3,  15.0d3,1000.0d3/) | Array of estimates of a typical minimum horizontal observation spacing (m) (in the densest observed area) for each obsetvation report type.<br>* This setting only affects the computational speed but not the analysis results<br>* This setting is used to automatically determine OBS_SORT_GRID_SPACING, effective only when OBS_SORT_GRID_SPACING = 0<ul><li><= 0: &nbsp;Same setting as OBS_MIN_SPACING(1)</li></ul> |
| <sub>__OBS_SORT_GRID_SPACING<br>&nbsp;&nbsp;&nbsp;(nobtype)__</sub> | (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0/) | Array of optimal grid spacing (m) for observation bucket sorting for each observation report type<br>* This setting only affects the computational speed but not the analysis results<ul><li>0: &nbsp;Automatically determined based on HORI_LOCAL, MAX_NOBS_PER_GRID, and OBS_MIN_SPACING</li><li>< 0: &nbsp;Same setting as OBS_SORT_GRID_SPACING(1)</li></ul> |

***

#### &PARAM_LETKF_VAR_LOCAL

Settings of variable localization

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__VAR_LOCAL_UV__</sub> <br> <sub>__VAR_LOCAL_T__</sub> <br> <sub>__VAR_LOCAL_Q__</sub> <br> <sub>__VAR_LOCAL_PS__</sub> <br> <sub>__VAR_LOCAL_RAIN__</sub> <br> <sub>__VAR_LOCAL_TC__</sub> <br> <sub>__VAR_LOCAL_RADAR_REF__</sub> <br> <sub>__VAR_LOCAL_RADAR_VR__</sub> <br> | (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/) | Weighting factors of variable localization from each observed variables. See [variable localization](Variable-localization.md) for details. |

***

#### &PARAM_LETKF_MONITOR

Observation diagnostic settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__DEPARTURE_STAT__</sub> | .true. | Output observation departure statistics (O-B and O-A) for conventional observations?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__DEPARTURE_STAT_RADAR__</sub> | .false. | Output observation departure statistics (O-B and O-A) for radar observations?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__DEPARTURE_STAT_T_RANGE__</sub> | 0.0d0 | Range of time difference to the analysis time, within which observations are considered in the departure statistics<ul><li>0: &nbsp;No time range restriction</li></ul> |
| <sub>__DEPARTURE_STAT_ALL_PROCESSES__</sub> | .true. | Print the departure statistics by all processes?<ul><li>.false.: &nbsp;No - The statistics are only printed by the ensemble mean group, which may save computational time</li><li>.true.: &nbsp;Yes - The same statistics are printed by all processes</li></ul> |
| <sub>__!__*OBSDEP_OUT*</sub> | .true. | Output observation departure (innovation) data for all observations into a binary file?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__!__*OBSDEP_OUT_BASENAME*</sub> | 'obsdep' | Filename of observation departure (innovation) data output |

***

#### &PARAM_LETKF_RADAR

Settings for radar data

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__USE_RADAR_REF__</sub> | .true. | Assimilate radar reflectivity observations?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__USE_RADAR_VR__</sub> | .true. | Assimilate radar radial velocity observations?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__METHOD_REF_CALC__</sub> | 3 | Method to compute the radar reflectivity in the radar observation operator <ul><li>1: &nbsp;WRF method</li><li>2: &nbsp;Tong and Xue 2006, 2008</li><li>3: &nbsp;Observation operator from Xue et al 2007 (coefficients modified by Amemiya et al. 2019)</li></ul> |
| <sub>__USE_METHOD3_REF_MELT__</sub> | .false. | Use radar operator considering melting (Xue et al. 2009QJRMS)?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__USE_T08_RS2014__</sub> | .false. | Use RS2014 in snow obsope (must be consistent with SCALE)<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__USE_TERMINAL_VELOCITY__</sub> | .false. | Consider the terminal velocity of the hydrometeors in the radar observation operator?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__USE_OBSERR_RADAR_REF__</sub> | .false. | Use OBSERR_RADAR_REF for the observation error of radar reflectivity observations instead of that provided in the observation files |
| <sub>__USE_OBSERR_RADAR_VR__</sub> | .false. | Use OBSERR_RADAR_VR for the observation error of radar radial velocity observations instead of that provided in the observation files |
| <sub>__RADAR_OBS_4D__</sub> | .false. | Radar observation data file is in a new format with the "time difference to the analysis time" column, allowing for 4D LETKF?<ul><li>.false.: &nbsp;No - Old-format file (without the "time difference to the analysis time" column)</li><li>.true.: &nbsp;Yes - New-format file (with the "time difference to the analysis time" column)</li></ul> |
| <sub>__RADAR_REF_THRES_DBZ__</sub> | 15.0d0 | Threshold of raining and clear-sky radar reflectivity observations (dBZ) |
| <sub>__MIN_RADAR_REF_MEMBER_OBSNORAIN__</sub> | 1 | Threshold of number of first-guess ensemble members with raining reflectivity to assimilate the *__clear-sky__* radar reflectivity data<br>* The observation data are assimilated only when the number of raining (> RADAR_REF_THRES_DBZ) members is above this threshold. |
| <sub>__MIN_RADAR_REF_MEMBER_OBSRAIN__</sub> | 1 | Same as above, but the threshold to assimilate the *__raining__* reflectivity observations |
| <sub>__MIN_RADAR_REF_DBZ__</sub> | 0.0d0 | Minimum useful radar reflectivity value (dBZ); all reflectivity data below this value are re-assigned to a constant depending on LOW_REF_SHIFT |
| <sub>__LOW_REF_SHIFT__</sub> | 0.0d0 | Shift of the constant refelectivity value for those data smaller than MIN_RADAR_REF_DBZ (dBZ); all reflectivity data below MIN_RADAR_REF_DBZ are set to (MIN_RADAR_REF_DBZ + LOW_REF_SHIFT)<br>* This setting should be zero or negative |
| <sub>__RADAR_ZMAX__</sub> | 99.0d3 | Maximum height level (m) of radar data to be assimilated |
| <sub>__RADAR_ZMIN__</sub> | -99.0d3 | Minimum height level (m) of radar data to be assimilated |

***

#### &PARAM_OBS_ERROR

Observation error settings<br>__* Ususally the observation errors provided in the observation data files, instead of the values here, are used in the assimilation, unless some special options are enabled__

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__OBSERR_RADAR_REF__</sub> | 5.0d0 | Observation error of radar reflectivity observations (dBZ) |
| <sub>__OBSERR_RADAR_VR__</sub> | 3.0d0 | Observation error of radar radial velocity observations (m/s) |

***
