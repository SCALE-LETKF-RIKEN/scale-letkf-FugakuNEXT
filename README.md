# SCALE-LETKF benchmark

Update: Apr 2024

This contains two benchmark experiments as follows. 

1. Idealized supercell experiment
---------------------------------

A simple idealized experiment of a data assimilation cycle.
It consists of two steps corresponding to separate binaries.
- scale-rm_ens : ensemble forecast from 00:25 to 00:30 using SCALE model 
- letkf        : data assimilation of radar observation by LETKF

2. Real-world weather radar data assimilation experiment
--------------------------------------------------------

A real-world case of a data assimilation cycle using a weather radar. (Miyoshi et al. 2023 SC)
It consists of four steps corresponding to separate binaries.
- scale-rm_pp_ens   : preparation of topography and landuse data (necessary only for the first run) 
- scale-rm_init_ens : preparation of initial and boundary condition data using SCALE model 
- scale-rm_ens      : ensemble forecast from 2021-07-30 06:00:00 to 06:00:30 UTC using SCALE model 
- letkf             : data assimilation of radar observation by LETKF

3. Source code
--------------

scale       : version hotfix/5.5.1 (8f7279707115bbfb78d426a7343230f4dd71cd6d) (BSD-2-Clause licence)
scale-letkf : version 5.5.0-v1 (MIT licence)

4. Compilation
--------------

```
$ cd $TOPDIR
$ source setup_spack.sh
$ export SCALE_SYS="FUGAKU"
```

TOPDIR indicates the path of SCALE-LETKF/Full_SCALE-LETKF/CPUGPU directory.

```
$ cd $TOPDIR/scale/scale-rm/src
$ make -j
$ cd $TOPDIR/scale/scale-letkf/scale
$ make
``` 

5. Run SCALE-LETKF
------------------

Get the input files (\*.tar.gz files) from RIKEN Box (see README.md at root directory),
and install to $TOPDIR directory.
If you select the SC23 experiment, you have to expand the file of SCALE-LETKF.dataset-SC23.part1.tar.gz.
For example:

```
$ cd $TOPDIR
$ tar zxvf scale_database.tar.gz
$ tar zxvf SCALE-LETKF.dataset-SC23.part1.tar.gz
$ cd $TOPDIR/test/SC23
$ sh prep.sh 
$ pjsub exec.sh
```

The netcdf files outputs $TOPDIR/result directory,
while the log files outputs the log sub-directory in the test directory.

6. Change settings
------------------

edit prep.sh and execute it to overwrite exec.sh 

GROUP    : your Fugaku user group
MEMBER   : number of ensemble members for LETKF (2-50)
NTHREADS : number of threads per process (1,2,3,4,6,8,12,16,24, or 48 : 12 is recommended for Fugaku)

( for SC23 only )
PRC_NUM_X : number of processes (=subdomains) in X direction (2,4,8,16, or 32)
PRC_NUM_Y : number of processes (=subdomains) in Y direction (2,4,8,16, or 32)

7. Execution time evaluation
----------------------------

The SCALE part records detailed execution times in the standard output, which should be consulted.


The LETKF part does not include a time measurement function;
the execution time of the entire LETKF program is the subject of evaluation.

8. Verification
---------------

Verification of the LETKF program should check the REF statistics in the log file.
The log file (log/NOUT_letkf.4.0) for the first node in the 10-ensemble benchmark test set will output the following:

```
OBSERVATIONAL DEPARTURE STATISTICS [ANALYSIS] (IN THIS SUBDOMAIN):
==========================================================================================
                 U           V           T           Q          PS         REF          Vr
------------------------------------------------------------------------------------------
BIAS           N/A         N/A         N/A         N/A         N/A  -9.547E-01  -1.424E-01
RMSE           N/A         N/A         N/A         N/A         N/A   3.264E+00   1.652E+00
NUMBER           0           0           0           0           0        8390          45
==========================================================================================
OBSERVATIONAL DEPARTURE STATISTICS [ANALYSIS] (GLOBAL):
==========================================================================================
                 U           V           T           Q          PS         REF          Vr
------------------------------------------------------------------------------------------
BIAS           N/A         N/A         N/A         N/A         N/A  -1.331E+00   7.662E-02
RMSE           N/A         N/A         N/A         N/A         N/A   4.719E+00   1.023E+00
NUMBER           0           0           0           0           0      399788      103021
==========================================================================================
 #### TIMER # ...monit_obs_mpi:monit_print:                              0.000252      0.000252
 #### TIMER # ...write_ens_mpi:monit_obs_mpi:                            0.000012      0.000012
 #### TIMER # ...write_ens_mpi:state_trans_inv:                          0.003947      0.003947
```

If these statistics match, there is no problem.

9. Tuning points
----------------

If the number of ensembles is small (<100), the computational cost of LETKF is small.
The computational cost in this case would be accounted for by the weather forecast in SCALE.
If the number of ensembles is sufficiently large (>1000), the computational cost of LETKF also increases relative to the number of ensembles.
In this case, the bottleneck of LETKF is the eigenvalue computation.
