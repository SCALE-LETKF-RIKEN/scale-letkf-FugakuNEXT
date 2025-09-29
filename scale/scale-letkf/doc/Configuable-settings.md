Last update: 2023.8.6

The following files contain the configurable settings for experiments. All of them are located in `letkf/scale/run` directory.

| Configuration file | Note | Type | Used in _cycle_ jobs | Used in _fcst_ jobs |
| --- | --- | --- | --- | --- |
| config.rc | Common configuration variables <br>(no need to change by users) | Bash | ✔ | ✔ |
| [config.main](List-of-variables-in-config.main.md) | Main configuration file | Bash | ✔ | ✔ |
| [config.cycle](List-of-variables-in-config.cycle.md) | Configurations for **cycle** jobs | Bash | ✔ |  |
| [config.fcst](List-of-variables-in-config.fcst.md) | Configurations for **fcst** jobs | Bash |  | ✔ |
| [config.nml.ensmodel](List-of-variables-in-config.nml.ensmodel.md) | Common namelist <br>(no need to change by users) | Fortran namelist | ✔ | ✔ |
| config.nml.scale_pp | Namelist for **scale-rm_pp** | Fortran namelist | (✔) | (✔) |
| config.nml.scale_init | Namelist for **scale-rm_init** | Fortran namelist | (✔) | (✔) |
| config.nml.grads_boundary | Namelist for **scale-rm_init** with GrADS input <br>(only effective when using GrADS-format input files) | Fortran namelist | (✔) | (✔) |
| config.nml.scale | Namelist for **scale-rm** | Fortran namelist | ✔ | ✔ |
| config.nml.scale_user | Namelist for user-customized **scale-rm** <br>(only effective when using a customized SCALE version for LETKF) | Fortran namelist | (✔) | (✔) |
| config.nml.obsope | Namelist for **obsope** <br>(only effective when using a separate observation operator) | Fortran namelist | (✔) |  |
| config.nml.obsmake | Namelist for **obsmake** <br>(only effective when using a separate observation generator) | Fortran namelist | (✔) |  |
| [config.nml.letkf](List-of-variables-in-config.nml.letkf.md) | Namelist for **letkf** | Fortran namelist | ✔ |  |

