# SCALE-LETKF

**SCALE-LETKF** is a data assimilation and numerical weather prediction system utilizing [LETKF](https://github.com/takemasa-miyoshi/letkf). 

## Sources 

### SCALE model

[Scalable Computing for Advanced Library and Environment (SCALE)](https://scale.riken.jp/) - Regional Model (RM)

* [Get source code](https://scale.riken.jp/download/)  

  The current SCALE-LETKF version supports SCALE-5.5.1.

* [Documentation](http://scale.riken.jp/doc/)  
  Model description and user guide.

* [Data files (source of topo/landuse)](https://scale.riken.jp/download/#datasets)

### Database

Input data sets to run the testcases are loacated in the following directories.

**Fugaku:** /share/hp150019/scale_database  
**Hibuna:** /home/amemiya/scale_database  

(**For fugaku users**: Be sure to copy the database to your directory under /vol0003 or /vol0004, as /share is not accesible directly from the compute node.)   

## [Getting started](Getting-started.md) 

## Run testcases

### Ideal experiments

* [Baroclinic waves](Baroclinic-waves.md)

* [Supercell](Supercell.md)

### Real experiments

* [PREPBUFR DA for Japan area](18km_Japan.md)

* [PAWR DA](PAWR-DA.md)

## [Configuable settings](Configuable-settings.md)
