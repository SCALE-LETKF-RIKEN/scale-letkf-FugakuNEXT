Observation departure statistics output is enabled when the namelist paraemter `OBSDEP_OUT` of PARAM_LETKF_MONITOR is set.   

### Common format 

The common obsdep file format used in SCALE-LETKF consists of 12-record elements of 4-byte floating-point value for each observation. Note that the access type is *sequential* and not *direct*. The following is a part of the subroutine `write_obs_dep` in `scale/common/common_obs_scale.f90` (simplified for brevity). 

```
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obs%nobs
    wk(1) = real(obs(set(n))%elm(idx(n)), r_sngl)
    wk(2) = real(obs(set(n))%lon(idx(n)), r_sngl)
    wk(3) = real(obs(set(n))%lat(idx(n)), r_sngl)
    wk(4) = real(obs(set(n))%lev(idx(n)), r_sngl)
    wk(5) = real(obs(set(n))%dat(idx(n)), r_sngl)
    wk(6) = real(obs(set(n))%err(idx(n)), r_sngl)
    wk(7) = real(obs(set(n))%typ(idx(n)), r_sngl)
    wk(8) = real(obs(set(n))%dif(idx(n)), r_sngl)
    wk(9) = real(qc(n), r_sngl)
    wk(10) = real(omb(n), r_sngl)
    wk(11) = real(oma(n), r_sngl)
    wk(12) = real(spr(n), r_sngl)
    WRITE(iunit) wk
  END DO
```

The contents are as follows. 

| Record # | Type | Description |
| --- | --- | --- |
| 1 | float | Observation variable (integer code defined in [common_obs_scale.f90](../scale/common/common_obs_scale.f90#L49)) |
| 2 | float | Longitude (degree) |
| 3 | float | Latitude (degree) |
| 4 | float | Vertical level (hPa) , or surface elevation (m) in the case of surface pressure observation |
| 5 | float | Observed value |
| 6 | float | Observation error standard deviation |
| 7 | float | Observation type (integer code defined in [common_obs_scale.f90](../scale/common/common_obs_scale.f90#L89)) |
| 8 | float | Difference of observation time from the target time of data assimilation (second) |
| 9 | float | QC flag (integer code defined in [common_obs_scale.f90](../scale/common/common_obs_scale.f90#L144)) |
| 10| float | O-B |
| 11| float | O-A |
| 12| float | First guess ensemble spread in the observation space |

