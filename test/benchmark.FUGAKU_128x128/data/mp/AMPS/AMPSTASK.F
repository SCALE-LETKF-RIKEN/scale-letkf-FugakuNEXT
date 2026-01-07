&AMPS_param
! *********************************************************************************
! integer variable of options
! *********************************************************************************
! nbin = 0 means no construction of the group.
!
! 1. level of complexity
! level_com:1
!          :2 a, c axis lengths and dendric length  predicted
!          :3 rosette and irregular length added to 2.
!          :4 aerosol mass components are predicted to 2.
!          :5 aerosol mass components are predicted to 3.
!
! 1. level of complexity
level_comp = 7
! 2. debugging level
debug_level = 1
! 3. what model needs
!    1 -> only tendency is updated.
!    2 -> tendency and variable are updated.
out_type = 2
!
! 4. time period for printing (s)
!T_print_period = 1800
!T_print_period = 1200
T_print_period = 3600000

! 5. output file format
!      defaut : text format, 'binary' : binary format (stream)
output_format = 'binary',
!- Cloud_drop group -
! 5. token (type of hydrometeor)
token_c = 11
! 7. dis_type
dtype_c = 0
! 8. h_break_method
hbreak_c = 0
! 9. flag for prediction (concentration fixed?)
!    1 -> concentration is fixed
!    2 -> concentration is predicted.
flagp_c = 2
!
! - rain group -
! 10. token (type of hydrometeor)
token_r = 1
! 12. dis_type
dtype_r = 1
! 13. type of hydrodynamic breakup
hbreak_r = 2
! 14. flag for prediction (concentration fixed?)
flagp_r = 2
!
! - Solid hydrometeor group -
! 15. token (type of hydrometeor)
token_s = 2
! 17. dis_type
dtype_s = 1
! 18. type of hydrodynamic breakup
hbreak_s = 2
! 19. flag for prediction (concentration fixed?)
flagp_s = 2
!
! - aerosol group -
! 20. token (type of hydrometeor)
token_a = 3
! 21. dis_type for category 1, and 2
dtype_a = 3,4,3,3
! 23. flag for prediction (concentration fixed?)
! 0: no prediction in micro
! 1: predict concentration and diagnose mass content (standard deviation fixed)
! 2: predict concentration and mass content (standard deviation fixed)
! 3: predict concentration and mass content (mean radius fixed)
! -1: initialize aerosol groups, but predict only CCN category and pass it to dynamic model (standard deviation fixed)
! -2: initialize aerosol groups, but predict only IN category and pass it to dynamic model (standard deviation fixed)
! -3: initialize aerosol groups, but do not pass any change of aerosols to dynamic model (standard deviation fixed)
! -4: initialize aerosol groups, but predict only CCN category and pass it to dynamic model (mean radius fixed)
! -5: initialize aerosol groups, but predict only IN category and pass it to dynamic model (mean radius fixed)
! -6: initialize aerosol groups, but do not pass any change of aerosols to dynamic model (mean radius fixed)
!flagp_a = 1
!flagp_a = 2
flagp_a = -3
!
! 24. habit growth determination
!     0: max frequency, 1: random number method
ihabit_gm_random = 1
!
! *********************************************************************************
! real variable of options
! *********************************************************************************
! - Cloud_drop group -
! 1. size_ratio
srat_c = 0.0
! 2. size_addition
sadd_c = 0.0
! 3. minimum_mass
minmass_c = 0.0
! 4. fixed concentration
fcon_c = 25.0
!
! - Rain group -
! 5. size_ratio
!    for n=20
srat_r = 2.540068909
! 6. size_addition
sadd_r = 0.0
! 7. minimum_mass
minmass_r = 4.188790E-09
!    for n=20
sth_r = 1.02
!
! - Solid hydrometeor group -
! 8. size_ratio
! this is for 83 bins from 4.19e-12 to 10.0 g
!!!srat_s = 1.4142135623731
! this is for 20 bins from 4.19e-12 to 10.0 g
srat_s = 4.1581061E+00
! 9. size_addition
sadd_s = 0.0
! 10. minimum_mass
minmass_s = 4.18879020478639E-12
!
! - aerosol group -
! 11. size_ratio
srat_a = 1.0e-6
! 12. size_addition
sadd_a = 0.0
! 13. minimum_mass
minmass_a = 1.0E-21
!
! for physical property of CCN, and IN
!
! 14. molecular weights of soluble part
! NaCl
! M_aps = 58.45,58.45,58.45,58.45
! (NH4)2SO4
! M_aps = 132.14052,132.14052,132.14052,132.14052
! NH4HSO4
M_aps = 115.11,115.11,115.11,115.11,
! 15. molecular weights of insoluble part
! silver iodide (AgI)
! M_api = 234.77,234.77,234.77,234.77
! Kaolinite (Al2Si2O5(OH)4)
M_api = 258.16,258.16,258.16,258.16
! 16. bulk density of soluble material for CCN
! NaCl
! den_aps = 2.16,2.16,2.16,2.16
! (NH4)2SO4
! den_aps = 1.77,1.77,1.77,1.77
! NH4HSO4
den_aps = 1.79,1.79,1.79,1.79,
! 17. bulk density of insoluble material for CCN
! silver iodide (AgI)
!den_api = 5.683,5.683,5.683,5.683
! Kaolinite (Al2Si2O5(OH)4)
den_api = 2.6,2.6,2.6,2.6
! 18. number of ions in dessociated solutes
! NaCl
! nu_aps = 2.0, 2.0, 2.0, 2.0,
! (NH4)2SO4
! nu_aps = 3.0,3.0,3.0,3.0,
! NH4HSO4
nu_aps = 2.0,2.0,2.0,2.0,
! 19. molal coefficient (average)
phi_aps = 0.75, 0.75, 0.75, 0.75,
! 20. mass fraction of soluble material to total mass
eps_ap = 0.7,0.0,0.2,0.7
! 21. log of standard deviation (in natural log of radius in micro m)
! from Jaenicke (1993), marine type, mode II,I,III: log(sig)=0.210, 0.657, 0.396 in diameter
!    log(sig) in radius is exactly the same as those in diameter.
! ap_lnsig = 0.210, 1.0, 0.657, 0.396,
! MPACE (sig=2.04, 2.5)
!ap_lnsig = 0.712949807856125, 1.0, 0.916290731874155, 1.0,
! for W2 in Shipway and Hill (2012)
!ap_lnsig = 0.405465108108164, 1.0, 0.405465108108164, 0.405465108108164,
! for SHEBA
!ap_lnsig = 2.04, 1.0, 2.5, 1.0
ap_lnsig = 0.712949807856125,0.0,0.916290731874155,0.0
! 22. geometrical mean radius (in cm)
!  from Jaenicke marine type
! ap_mean = 0.133e-4,1.0e-4,0.004e-4,0.29e-4
! MPACE
! ap_mean = 0.052e-4, 1.0e-4, 1.3e-4, 0.001e-4,
! SHEBA
ap_mean = 0.052e-4, 1.0e-4, 1.3e-4, 0.0
! for W2 in Shipway and Hill (2012)
!ap_mean = 0.1e-4, 1.0e-4, 0.1e-4, 0.1e-4,
!
! 23. geometrical mean contact parameter (in degree), normal distriubtion
ap_mean_cp = 132.0, 15.5, 132.0, 132.0,
!
! 24. standard deviation of contact parameter (in degree), normal distribution
!
ap_sig_cp = 20.0, 1.4, 20.0, 20.0,
!
! 25. Initial number concentraion of aerosol particles
N_ap_ini = 350.0, 0.00015, 1.8, 0.0,
!
! * character variable of options
! 1. DRCETB
! - directory where the collision coefficiency files are stored.
DRCETB = 'AMPS_DATA/collision_data'
! 2. IFMGTB
! - parameter or database file for mass-geometric relation of ice crystals.
! DRMGTB = 'prmRHI.v1.dat',
! 3. DRAPTB
! - directory where the aerosol activation file is stored.
DRAPTB = 'AMPS_DATA/apact'
! 4. DRSTTB
! - directory where the statistics lookup table is stored.
DRSTTB = 'AMPS_DATA/statpack'
!
! 5. Name of chemicals for aerosols
! APSNAME = 'NACL2','NACL2','NACL2'
! APSNAME = 'NH42SO4','NH42SO4','NH42SO4','NH42SO4'
APSNAME = 'NH4HSO4','NH4HSO4','NH4HSO4','NH4HSO4'
!
! 6. model type of CCN activation
! 1: implicit prediction of vapor (default)
! 2: water saturation adjustment
act_type = 1
!
! sensitivity test
!
coll_level = 1
!
! simulation under fixed Relative Humidity and Temperature
!
! relative humidity over water (%)
!FRHW = 100.0
! relative humidity over ice (ratio)
!FRHI = -999.9
! temperature in Celsius
!FTMP = -5.0
!
! number of time steps relative to one dynamic step
!  collision processes
n_step_cl = 1
! number of time steps relative to collision processes
!  vapor deposition processes
n_step_vp = 10
!
! execution flag for microphysical processes
! on -> 1, off -> 0
!   some processes, other schemes may be chosen with this flag.
!
! index
! 1: printing
! 2: liq-liq collsion-coalescence process
! 3: ice-ice collsion-coalescence process
! 4: liq-ice collsion-coalescence process
! 5: update surface temperature
! 6: vapor deposition on liq
! 7: vapor deposition on ice
! 8: melting shedding process
! 9: hydrodynamic breakup of ice particles
! 10: ice nucleation process
! 11: hydrodynamic breakup of liquid hydrometeors
! 12: autoconversion of cloud droplet bin
! 13: ice nucleation; depositional nucleation
! 14: ice nucleation; contact freezing
! 15: ice nucleation; splinter nulcreation
! 16: ice nucleation; immersion freezing (heterogenious freezing)
! 17: ice nucleation; homogenious freezing
!         0                 1
!         1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
!micexfg = 1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1
micexfg = 1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,2,0,1,0,1

!micexfg = 1,1,1,1,1,1,1,1,0,1,0,0,1,1,0,1,1,1,1,1
!micexfg = 1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,1
!micexfg = 1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,2,0,1,1,1
!micexfg = 1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1
!
/
