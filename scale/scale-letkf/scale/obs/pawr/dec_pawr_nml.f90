module dec_pawr_nml 

  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_obs_scale

  public 
 !--- PARAM_LETKF_RADAR_DECODE
  character(filelenmax) :: PAWR_IN_PATH = "obs.pawr<type>" ! Input raw PAWR file
  real(r_size) :: RADAR_SO_SIZE_HORI = 1000.0d0
  real(r_size) :: RADAR_SO_SIZE_VERT = 1000.0d0
  real(r_size) :: RADAR_MAX_ABS_VR = 100.0d0
  integer :: RADAR_THIN_HORI = 1 ! Thinning horizontal interval (# of grids)
  integer :: RADAR_THIN_VERT = 1 ! Thinning vertical interval (# of grids)
  REAL(r_size) :: MIN_RADAR_REF_DBZ_VR = 5.0d0 ! Minimum reflectivity (dBZ) for Doppler velocity observation
  logical :: RADAR_USE_VR_STD = .true.  !If we are going to use the wind std threshold within each box.
  logical :: RADAR_BIAS_COR_RAIN = .false. ! Simple bias correction for radar obs (rain)
  logical :: RADAR_BIAS_COR_CLR = .false. ! Simple bias correction for radar obs (clear sky)
  real(r_size) :: RADAR_BIAS_RAIN_CONST_DBZ = 0.0_r_size ! Simply bias correction for radar obs (rain)
  real(r_size) :: RADAR_BIAS_CLR_CONST_DBZ = 0.0_r_size ! Simply bias correction for radar obs (clear sky)
  logical :: USE_PAWR_MASK = .false. ! Saitama MP-PAWR shadow mask
  character(filelenmax) :: PAWR_MASK_FILE = '' ! data file for masking
  logical :: OUT_PAWR_GRADS = .false. ! Outut PAWR obs in GrADS format
  character(filelenmax) :: OUT_PAWR_GRADS_PATH = "pawr_ref3d_"   ! Output path
  character(filelenmax) :: OUT_PAWR_SUPEROB_PATH = "radar_" ! Output PAWR obs in LETKF format path
 
  REAL(r_size),save :: MIN_RADAR_REF_VR = 0.0d0    ! Minimum reflectivity (Z) for Doppler velocity observation
  REAL(r_size),save :: RADAR_BIAS_RAIN_CONST = 0.0d0    
  REAL(r_size),save :: RADAR_BIAS_CLR_CONST = 0.0d0     

contains
!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------
subroutine read_nml_letkf_radar_decode
  use scale_io, only: IO_FID_CONF
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_RADAR_DECODE/ &
    PAWR_IN_PATH, &
    RADAR_SO_SIZE_HORI, &
    RADAR_SO_SIZE_VERT, &
    RADAR_THIN_HORI, &
    RADAR_THIN_VERT, &
    MIN_RADAR_REF_DBZ_VR, & 
    RADAR_USE_VR_STD, &
    RADAR_BIAS_COR_RAIN, &
    RADAR_BIAS_COR_CLR, &
    RADAR_BIAS_RAIN_CONST_DBZ, &
    RADAR_BIAS_CLR_CONST_DBZ, &
    RADAR_MAX_ABS_VR, &
    USE_PAWR_MASK, &
    PAWR_MASK_FILE, &
    OUT_PAWR_GRADS, &
    OUT_PAWR_GRADS_PATH, &
    OUT_PAWR_SUPEROB_PATH

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_RADAR_DECODE,iostat=ierr)
  if (ierr < 0) then !--- missing
    if ( LOG_OUT ) write(6,*) '[Warning] /PARAM_LETKF_RADAR_DECODE/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LETKF_RADAR_DECODE. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LETKF_RADAR_DECODE)
  end if

  MIN_RADAR_REF_VR = 10.0d0 ** (MIN_RADAR_REF_DBZ_VR/10.0d0)
  RADAR_BIAS_RAIN_CONST = 10.0d0 ** ( RADAR_BIAS_RAIN_CONST_DBZ / 10.0d0 )
  RADAR_BIAS_CLR_CONST = 10.0d0 ** ( RADAR_BIAS_CLR_CONST_DBZ / 10.0d0 )

return
end subroutine read_nml_letkf_radar_decode


end module
