### Analysis and first-guess ensemble mean and spread 

In the LETKF, the estimation of state variables of first guess and analysis are represented by their ensemble mean values. And the uncertainties of them are represented by their ensemble spread (standard deviation). For that, the LETKF output data includes ensemble mean and spread by default.

### File treatment in a DA cycle

The current LETKF code handles those files in a bit tricky way. As it does not have a function to create NetCDF files of SCALE restart file format, it needs to overwrite the files which are prepared beforehand. 

Suppose if you perform a 6-hour DA cycle initialized by the analysis mean at 2024/01/01 00Z. The DA cycle consists of executing 3 binaries; scale-rm_init_ens, scale-rm_ens, and letkf, which respectively correspond to initialization, time integration, and data assimilation steps.   
When the option `OUT_OPT=5` in config.cycle is set, first guess and analysis files are generated in the following manner.  

#### Before scale-rm_ens

| 20240101000000/gues/mean | 20240101000000/gues/\<member\> | 20240101000000/anal/mean | 20240101000000/anal/\<member\> | 
| ------------------------ | ------------------------------ | ------------------------ | ------------------------------ | 
| first guess mean         | None                           | analysis mean            | analysis members               |

#### After scale-rm_ens

| 20240101000000/gues/mean | 20240101000000/gues/\<member\> | 20240101000000/anal/mean | 20240101000000/anal/\<member\> | 
| ------------------------ | ------------------------------ | ------------------------ | ------------------------------ | 
| first guess mean         | None                           | analysis mean            | analysis members               |

| 20240101060000/gues/mean | 20240101060000/gues/\<member\> | 20240101060000/anal/mean | 20240101060000/anal/\<member\> | 
| ------------------------ | ------------------------------ | ------------------------ | ------------------------------ | 
| **6-h fcst from mean (dummy)**   | None                   | **6-h fcst from mean (dummy)** | **first guess** members  |

Note that before letkf step, restart files in mean directories are actually not the ensemble mean, but the forecast from the analysis ensemble mean of the previous analysis time step. This is performed just to create restart files and the values in these files at this step are not used in letkf and just overwritten. 

In the case of 4-D LETKF, you also have an ensemble of history files in separate directories.

If the option `OUT_OPT` is 1, 2, or 3, each member of first guess files for each member are also created in this step. 

#### After letkf

| 20240101000000/gues/mean | 20240101000000/gues/\<member\> | 20240101000000/anal/mean | 20240101000000/anal/\<member\> | 
| ------------------------ | ------------------------------ | ------------------------ | ------------------------------ | 
| first guess mean         | None                           | analysis mean            | analysis members               |

| 20240101060000/gues/mean | 20240101060000/gues/\<member\> | 20240101060000/anal/mean | 20240101060000/anal/\<member\> | 
| ------------------------ | ------------------------------ | ------------------------ | ------------------------------ | 
| **first guess mean**     | None                           | **analysis mean**        | **analysis members**           |

LETKF calculates first guess and analysis mean values and writes them into the files created in the scale-rm_ens step. It also overwrites analysis files of each member.

### Spread file format

Ensemble spread files in gues/sprd, anal/sprd are generated in the similar way with ensemble mean. The spread file has the similar restart file format, but has a different list of variables. Note that the variable set **does not correspond** to the variable name list in the NetCDF file which can be seen by ncdump, as LETKF can't edit NetCDF metadata.

The actual variables in a spread file are as follows. 

| variable name in ncdump | actual variable in sprd file |
| ----------------------- | ---------------------------- |
| DENS                    | ensemble spread of  U        |
| MOMX                    | ensemble spread of  V        | 
| MOMY                    | ensemble spread of  W        | 
| MOMZ                    | ensemble spread of  T        | 
| RHOT                    | ensemble spread of  P        |
| Q\*                     | ensemble spread of  Q\*      |

### Forecast ensemble mean and spread 

Currently this code does not have a function to automatically output mean and spread of ensemble forecast. Ensemble forecast by `scale-rm_ens` just launchs forecasts by `rm_driver` for each member independently. You are supposed to calculate the ensemble mean and spread from history files of the members. Note that the history files in 'mean' directory is **not** an ensemble mean but a single forecast initialized by an analysis ensemble mean.
