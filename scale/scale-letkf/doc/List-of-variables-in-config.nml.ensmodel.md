Last update: 2023.8.6

__Note:__ Variables that are in italic and also have an exclamation mark (__!__) in front in the following tables are those to be automatically determined by the job scripts. Their values __should not be assigned by users__ and they should be written in `config.nml.letkf` as:
```
!--VARIABLE--
```

#### __Special notations__ that can be used in some namelist variables for filename patterns:
- Special notation for member strings: '@@@@'

***

#### &PARAM_ENSEMBLE

Ensemble settings

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__!__*MEMBER*</sub> | 3 | Ensemble size |
| <sub>__!__*MEMBER_RUN*</sub> | 1 | Actual number of ensemble to run in parallel (including mean and mdet) |
| <sub>__!__*CONF_FILES*</sub> | 'run.@@@@.conf' | config files for each member ('@@@@' will be replaced) |
| <sub>__!__*CONF_FILES_SEQNUM*</sub> | .false. |  |
| <sub>__!__*DET_RUN*</sub> | .false. | Enable the deterministic run (Schraff et al. 2016 QJRMS)?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |
| <sub>__DET_RUN_CYCLED__</sub> | .true. | Cycle the deterministic run?<ul><li>.false.: &nbsp;No - Use the forecast from the analysis ensemble mean (in the previous cycle) as the first guess for the deterministic analysis</li><li>.true.: &nbsp;Yes - Use the forecast from the deterministic analysis (in the previous cycle) as the first guess for the deterministic analysis</li></ul> |
| <sub>__!__*EFSO_RUN*</sub> | .false. | Enable the EFSO?<ul><li>.false.: &nbsp;No</li><li>.true.: &nbsp;Yes</li></ul> |

***

#### &PARAM_PROCESSES

Parallelization settings for LETKF

| Variable | Default value | Explanation |
| --- | --- | --- |
| <sub>__!__*PPN*</sub> | 1 | Number of MPI processes used per nodes for the LETKF |
| <sub>__!__*MEM_NODES*</sub> | 1 | Number of nodes used for one ensemble member in the LETKF |
| <sub>__!__*NUM_DOMAIN*</sub> | 1 | Number of domains (for online-nesting) |
| <sub>__!__*PRC_DOMAINS*</sub> | 1 | Number of MPI processes used for one ensemble member in the LETKF |
| <sub>__!__*COLOR_REORDER*</sub> | .false. | coloring reorder for mpi splitting? |

***
