Observation data is processed in the format specific to the SCALE-LETKF. Currently, different formats are used for radar and other observation types. The parameter `OBS_IN_FORMAT` needs to be properly set in order to read the data.  

### Common format 

The common observation format used in SCALE-LETKF consists of 8-record elements of 4-byte floating-point value for each observation. Note that the access type is *sequential* and not *direct*. The following is a part of the subroutine `write_obs` in `scale/common/common_obs_scale.f90` (simplified for brevity). 

```
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obs%nobs
      wk(1) = REAL(obs%elm(n),r_sngl)
      wk(2) = REAL(obs%lon(n),r_sngl)
      wk(3) = REAL(obs%lat(n),r_sngl)
      wk(4) = REAL(obs%lev(n),r_sngl)
      wk(5) = REAL(obs%dat(n),r_sngl)
      wk(6) = REAL(obs%err(n),r_sngl)
      wk(7) = REAL(obs%typ(n),r_sngl)
      wk(8) = REAL(obs%dif(n),r_sngl)
      WRITE(iunit) wk
  END DO
```

The contents are as follows. 

| Record # | Type | Description |
| --- | --- | --- |
| 1 | float | Observation variable (integer code defined in [common_obs_scale.f90](../scale/common/common_obs_scale.f90#L49)) |
| 2 | float | Longitude (degree) |
| 3 | float | Latitude (degree) |
| 4 | float | Vertical level (hPa) , or surface elevation (m) - see [note on vertical coordinate](Note-on-vertical-coordinate.md) |
| 5 | float | Observed value |
| 6 | float | Observation error standard deviation |
| 7 | float | Observation type (integer code defined in [common_obs_scale.f90](../scale/common/common_obs_scale.f90#L89)) |
| 8 | float | Difference of observation time from the target time of data assimilation (second) |


### Radar observation format 

For the radar observation, there is a slight modification in the data format. First, the location of the radar is recorded at the beginning of the file. Second, as it is common to use 3D-LETKF for radar data, the last record describing a time gap can be omitted. See the part of the subroutine `write_obs_radar` (again simplified).

```
IF(RADAR_OBS_4D) THEN
    nrec = 8
  ELSE
    nrec = 7
  END IF

  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  WRITE(iunit) REAL(obs%meta(1),r_sngl)
  WRITE(iunit) REAL(obs%meta(2),r_sngl)
  WRITE(iunit) REAL(obs%meta(3),r_sngl)

  DO n=1,obs%nobs
      wk(1) = REAL(obs%elm(n),r_sngl)
      wk(2) = REAL(obs%lon(n),r_sngl)
      wk(3) = REAL(obs%lat(n),r_sngl)
      wk(4) = REAL(obs%lev(n),r_sngl)
      wk(5) = REAL(obs%dat(n),r_sngl)
      wk(6) = REAL(obs%err(n),r_sngl)
      wk(7) = REAL(obs%typ(n),r_sngl)
      IF(RADAR_OBS_4D) THEN
        wk(8) = REAL(obs%dif(n),r_sngl)
      END IF
      WRITE(iunit) wk(1:nrec)
    end if
  END DO

```

At the beginning of the file, three 4-byte floating point real values are recorded. They are longtitude, latitude, elevation (m) of the radar location, respectively.  
