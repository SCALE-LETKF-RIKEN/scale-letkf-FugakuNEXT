
Some settings in LETKF about vertical location or length implicitly assume a specific unit. It might be confusing and cause errors which are difficult to trach down. 
This note is made for helping users avoid errors by inconsistency in vertical unit.  

### Vertical coordinate unit of observation data

In [the LETKF observation data format](Observation-file-format.md), vertical location of the observation is indicated in the 4 th element of each record. The coordinate unit depends on the observation data type and element as follows. 

| type # | type name | element # | element name              | unit  |
| ---    | ---       | ---       | ---                       | ---   | 
| 22     | 'PHARAD'  | 4001-4004 | 'REF', 'Vr', 'PRH', 'RE0' | meter |
| 1      | 'ADPUPA'  | 14593     | 'PS'                      | meter |
| else   |      -    | else      | -                         | hPa   | 

### Vertical coordinate unit in localization 

Vertical localization length scales are defined for each observation type in the parameter *VERT_LOCAL* of the namelist [PARAM_LETKF_OBS](List-of-variables-in-config.nml.letkf.md#param_letkf_obs).  
The unit is assumed depending on the observation type as foolows. 

| type # | type name | unit  |
| ---    | ---       | ---   | 
| 22     | 'PHARAD'  | meter |
| else   |      -    | log-p (non-dimensional) |

Therefore, `vert_local(22)` is interpreted as 'm' and the rest of the elements are as 'log-p'. Currently, the LETKF code supports only observation type 1 ('ADPUPA') and 22 ('PHARAD'). In the case if you modify the code to include other type of observation, be sure to set consistent vertical unit in the observation data file and the localization paratemter. 
