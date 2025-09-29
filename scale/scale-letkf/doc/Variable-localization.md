### Variable localization

The namelist *PARAM_LETKF_VAR_LOCAL* is dedicated to impose variable localization in the SCALE-LETKF. 
Variable localization may help alleviate the negative impact of sprious correlation due to limited ensemble size, just as spatial and temporal localization. 

Each namelist parameter corresponds to each of observation variable types. Each element of the array corresponds to each model variable to update. 
Default values are all 1, indicating that "all of 11 3-d model variables will be updated equally regardless of an observation variable type."

```
&PARAM_LETKF_VAR_LOCAL
##############            U      V      W      T      P     QV     QC     QR     QI     QS     QG   
VAR_LOCAL_UV        = 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  
VAR_LOCAL_T         = 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  
VAR_LOCAL_Q         = 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  
VAR_LOCAL_PS        = 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  
VAR_LOCAL_RAIN      = 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  
VAR_LOCAL_TC        = 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  
VAR_LOCAL_RADAR_REF = 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  
VAR_LOCAL_RADAR_VR  = 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  
/
```

The parameter is supposed to have values from 0 to 1. If you want to disable analysis update of wind by the assimilation of radar reflectivity, for instance, 
set corresponding elements zero.  
```
VAR_LOCAL_RADAR_REF = 0.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 
```

Or you can partially reduce the influence of the assimilation by using intermediate values.  
In the LETKF, the inverse of the value will be multiplied to the observation error variance along with other localization factors.
If the value is zero, the observation is just skipped and not used for that model variable update. 
```
VAR_LOCAL_RADAR_REF = 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.5d0, 0.5d0, 1.0d0, 0.5d0, 1.0d0, 1.0d0, 
```

Note: Variable localization increase the computational cost of the LETKF, 
as the calculation of transformation matrix at a single grid point needs to be performed more than once
 when model variables use different factors and thus different relative weighting of observations. 

### Reference

“Variable localization” in an ensemble Kalman filter:Application to the carbon cycle data assimilation  
[Kang, Kalnay, Liu, Fung, Miyoshi, and Ide (2011)](https://doi.org/10.1029/2010JD014673)  
