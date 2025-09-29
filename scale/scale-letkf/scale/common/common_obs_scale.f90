MODULE common_obs_scale
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   07/25/2014 Guo-Yuan Lien     modified for SCALE model
!   .......... See git history for the following revisions
!
!=======================================================================
!
! [LETKF observation format]
!   (In files, all variables are stored in single-precision float)
!
!  column  description
!     (1)  variable type (1..nid_obs; see 'id_*_obs' parameters)
!     (2)  longitude (degree)
!     (3)  latitude (degree)
!     (4)  level/height
!            u,v,t,tv,q,rh: level (hPa)
!            ps: station elevation (m)
!     (5)  observation value
!            wind (m/s)
!            temperature (K)
!            specific humidity (kg/kg)
!            relative humidity (%)
!            surface pressure (hPa)
!     (6)  observation error
!            unit same as observation value
!     (7)  observation platform type (1..nobtype+1; see 'obtypelist' array)
!     (8)  observation time relative to analysis time (sec)
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_nml
  USE common_scale

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: nid_obs_varlocal=8
!  INTEGER,PARAMETER :: nid_obs_varlocal=9 !H08
!
! conventional observations
!
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_tv_obs=3074
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
!
! surface observations codes > 9999
!
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991  ! TC vital
  INTEGER,PARAMETER :: id_tclat_obs=99992  ! TC vital
  INTEGER,PARAMETER :: id_tcmip_obs=99993  ! TC vital
!
! radar observations
!
  INTEGER,PARAMETER :: id_radar_ref_obs=4001
  INTEGER,PARAMETER :: id_radar_ref_zero_obs=4004
  INTEGER,PARAMETER :: id_radar_vr_obs=4002
  INTEGER,PARAMETER :: id_radar_prh_obs=4003
!
! Himawari-8 (H08) observations
!
  INTEGER,PARAMETER :: id_H08IR_obs=8800

  INTEGER,PARAMETER :: elem_uid(nid_obs)= &
     (/id_u_obs, id_v_obs, id_t_obs, id_tv_obs, id_q_obs, id_rh_obs, &
       id_ps_obs, id_rain_obs, id_radar_ref_obs, id_radar_ref_zero_obs, id_radar_vr_obs, id_radar_prh_obs, &
       id_H08IR_obs, id_tclon_obs, id_tclat_obs, id_tcmip_obs/)

  CHARACTER(3),PARAMETER :: obelmlist(nid_obs)= &
     (/'  U', '  V', '  T', ' Tv', '  Q', ' RH', ' PS', 'PRC', 'REF', 'RE0', ' Vr', 'PRH',&
       'H08', 'TCX', 'TCY', 'TCP'/)

  CHARACTER(3),PARAMETER :: obelmlist_varlocal(nid_obs_varlocal)= &
     (/'WND', '  T', 'MOI', ' PS', 'PRC', 'TCV', 'REF', ' Vr'/)
!     (/'WND', '  T', 'MOI', ' PS', 'PRC', 'TCV', 'REF', ' Vr', 'H08'/)

  ! Parameter 'nobtype' is set in common_nml.f90
  CHARACTER(6),PARAMETER :: obtypelist(nobtype)= &
     (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
       'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
       'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
       'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
       'TMPAPR', 'PHARAD', 'H08IRB', 'TCVITL'/) ! H08

  INTEGER,PARAMETER :: max_obs_info_meta = 3 ! maximum array size for type(obs_info)%meta

  TYPE obs_info
    INTEGER :: nobs = 0
    INTEGER,ALLOCATABLE :: elm(:)
    REAL(r_size),ALLOCATABLE :: lon(:)
    REAL(r_size),ALLOCATABLE :: lat(:)
    REAL(r_size),ALLOCATABLE :: lev(:)
    REAL(r_size),ALLOCATABLE :: dat(:)
    REAL(r_size),ALLOCATABLE :: err(:)
    INTEGER,ALLOCATABLE :: typ(:)
    REAL(r_size),ALLOCATABLE :: dif(:)
    REAL(r_size) :: meta(max_obs_info_meta) = undef
    REAL(r_size),ALLOCATABLE :: ri(:)
    REAL(r_size),ALLOCATABLE :: rj(:)
    INTEGER,ALLOCATABLE :: rank(:)
  END TYPE obs_info

  TYPE obs_da_value
    INTEGER :: nobs = 0
    INTEGER :: nobs_in_key = 0
    INTEGER,ALLOCATABLE :: set(:)
    INTEGER,ALLOCATABLE :: idx(:)
    INTEGER,ALLOCATABLE :: key(:)
    REAL(r_size),ALLOCATABLE :: val(:)
    !
    ! obsda%lev array is used only for Himawari-8 assimilation.
    ! This array preserves the most sensitive height derived from transmittance outputs from RTTOV.
    ! For Himawari-8 assimilation, LETKF uses obsda%lev instead of obs%lev.
    ! 
    REAL(r_size),ALLOCATABLE :: ensval(:,:)
    real(r_size), allocatable :: epert(:,:) ! ensemble perturbation (Y18)
    real(r_size), allocatable :: pert(:)    ! ensemble perturbation from each member (Y18)
    REAL(r_size),ALLOCATABLE :: eqv(:,:) ! qv (ensemble)
    REAL(r_size),ALLOCATABLE :: qv(:)    ! qv (mean)
    REAL(r_size),ALLOCATABLE :: tm(:)    ! temp (mean)
    REAL(r_size),ALLOCATABLE :: pm(:)    ! pressure (mean)
    INTEGER,ALLOCATABLE :: qc(:)
  END TYPE obs_da_value

  character(obsformatlenmax), parameter :: obsfmt_prepbufr = 'PREPBUFR'
  character(obsformatlenmax), parameter :: obsfmt_radar    = 'RADAR'
  character(obsformatlenmax), parameter :: obsfmt_h08      = 'HIMAWARI8'
  character(obsformatlenmax), parameter :: obsfmt_radar_nc = 'RADAR-NC'
!  integer, parameter :: nobsformats = 3
!  character(obsformatlenmax), parameter :: obsformat(nobsformats) = &
!    (/obsfmt_prepbufr, obsfmt_radar, obsfmt_h08/)

  INTEGER,PARAMETER :: iqc_good=0
  INTEGER,PARAMETER :: iqc_gross_err=5
  INTEGER,PARAMETER :: iqc_ps_ter=10
  INTEGER,PARAMETER :: iqc_ref_low=11
  INTEGER,PARAMETER :: iqc_ref_mem=12
  INTEGER,PARAMETER :: iqc_radar_vhi=19
  INTEGER,PARAMETER :: iqc_out_vhi=20
  INTEGER,PARAMETER :: iqc_out_vlo=21
  INTEGER,PARAMETER :: iqc_obs_bad=50
  INTEGER,PARAMETER :: iqc_otype=90
  INTEGER,PARAMETER :: iqc_time=97
  INTEGER,PARAMETER :: iqc_out_h=98
  INTEGER,PARAMETER :: iqc_undef=99

  type(obs_info),allocatable,save :: obs(:) ! observation information
  type(obs_da_value),save :: obsda_sort     ! sorted obsda

  integer, save :: obsdep_nobs                     ! obsdep information
  integer, allocatable, save :: obsdep_set(:)      ! 
  integer, allocatable, save :: obsdep_idx(:)      ! 
  integer, allocatable, save :: obsdep_qc(:)       ! 
  real(r_size), allocatable, save :: obsdep_omb(:) ! 
  real(r_size), allocatable, save :: obsdep_oma(:) ! 
  real(r_size), allocatable, save :: obsdep_sprd(:) ! 
  real(r_size), allocatable, save :: obsdep_omb_emean(:) ! 

  REAL(r_size),SAVE :: MIN_RADAR_REF
  REAL(r_size),SAVE :: RADAR_REF_THRES

CONTAINS

!-----------------------------------------------------------------------
! Convert a raw obsID to a sequential obsID (1 - nid_obs)
!-----------------------------------------------------------------------
function uid_obs(id_obs)
  implicit none
  integer :: id_obs
  integer :: uid_obs

  select case(id_obs)
  case(id_u_obs)
    uid_obs = 1
  case(id_v_obs)
    uid_obs = 2
  case(id_t_obs)
    uid_obs = 3
  case(id_tv_obs)
    uid_obs = 4
  case(id_q_obs)
    uid_obs = 5
  case(id_rh_obs)
    uid_obs = 6
  case(id_ps_obs)
    uid_obs = 7
  case(id_rain_obs)
    uid_obs = 8
  case(id_radar_ref_obs)
    uid_obs = 9
  case(id_radar_ref_zero_obs)
    uid_obs = 10
  case(id_radar_vr_obs)
    uid_obs = 11
  case(id_radar_prh_obs)
    uid_obs = 12
  case(id_h08ir_obs) ! H08
    uid_obs = 13     ! H08
  case(id_tclon_obs)
    uid_obs = 14
  case(id_tclat_obs)
    uid_obs = 15
  case(id_tcmip_obs)
    uid_obs = 16
  case default
    uid_obs = -1     ! error
  end select
end function uid_obs
!-----------------------------------------------------------------------
! Convert a raw obsID to a sequential obsID for variable localization (1 - nid_obs_verlocal)
!-----------------------------------------------------------------------
function uid_obs_varlocal(id_obs)
  implicit none
  integer :: id_obs
  integer :: uid_obs_varlocal

  select case(id_obs)
  case(id_u_obs, id_v_obs)
    uid_obs_varlocal = 1
  case(id_t_obs, id_tv_obs)
    uid_obs_varlocal = 2
  case(id_q_obs, id_rh_obs)
    uid_obs_varlocal = 3
  case(id_ps_obs)
    uid_obs_varlocal = 4
  case(id_rain_obs)
    uid_obs_varlocal = 5
  case(id_tclon_obs, id_tclat_obs, id_tcmip_obs)
    uid_obs_varlocal = 6
  case(id_radar_ref_obs, id_radar_ref_zero_obs, id_radar_prh_obs)
    uid_obs_varlocal = 7
  case(id_radar_vr_obs)
    uid_obs_varlocal = 8
  case(id_h08ir_obs)      ! H08
    uid_obs_varlocal = 9  ! H08
  case default
    uid_obs_varlocal = -1 ! error
  end select
end function uid_obs_varlocal

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine set_common_obs_scale
  implicit none

  MIN_RADAR_REF = 10.0d0 ** (MIN_RADAR_REF_DBZ/10.0d0)
  RADAR_REF_THRES = 10.0d0 ** (RADAR_REF_THRES_DBZ/10.0d0)

  return
end subroutine set_common_obs_scale

!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!
! stggrd: grid type of u and v
!  0: non-staggered grid
!  1: staggered grid
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rk,lon,lat,v3d,v2d,yobs,qc,stggrd,typ)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: lon,lat
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  integer, intent(in), optional :: typ
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER,INTENT(OUT) :: qc
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  REAL(r_size) :: u,v,t,q,p,topo

  INTEGER :: stggrd_ = 0
  if (present(stggrd)) stggrd_ = stggrd

  yobs = undef
  qc = iqc_good

  SELECT CASE (elm)
  CASE(id_u_obs,id_v_obs)  ! U,V
    if (stggrd_ == 1) then
      CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5_r_size,rj,u)  !###### should modity itpl_3d to prevent '1.0' problem....??
      CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5_r_size,v)  !######
    else
      CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,u)
      CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,v)
    end if

    if ( present( typ ) ) then
      if ( obtypelist(typ) == 'ASCATW' ) then
        if (stggrd_ == 1) then
          call itpl_2d(v2d(:,:,iv2dd_u10m),ri-0.5_r_size,rj,u)
          call itpl_2d(v2d(:,:,iv2dd_v10m),ri,rj-0.5_r_size,v)
        else
          call itpl_2d(v2d(:,:,iv2dd_u10m),ri,rj,u)
          call itpl_2d(v2d(:,:,iv2dd_v10m),ri,rj,v)
        end if
      endif
    endif

    if (elm == id_u_obs) then
      yobs = u
    else
      yobs = v
    end if
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
  CASE(id_tv_obs)  ! Tv
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,q)
    yobs = yobs * (1.0d0 + fvirt * q)
  CASE(id_q_obs)  ! Q
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,yobs)
  CASE(id_ps_obs) ! PS
    CALL itpl_2d(v2d(:,:,iv2dd_t2m),ri,rj,t)
    CALL itpl_2d(v2d(:,:,iv2dd_q2m),ri,rj,q)
    CALL itpl_2d(v2d(:,:,iv2dd_topo),ri,rj,topo)
    CALL itpl_2d(v2d(:,:,iv2dd_ps),ri,rj,yobs)
    call prsadj(yobs,rk-topo,t,q)
    if (abs(rk-topo) > PS_ADJUST_THRES) then
      if (LOG_LEVEL >= 2) then
        write (6,'(A,F6.1)') '[Warning] PS observation height adjustment exceeds the threshold. dz=', abs(rk-topo)
      end if
      qc = iqc_ps_ter
    end if
!  CASE(id_rain_obs) ! RAIN                        ############# (not finished)
!    CALL itpl_2d(v2d(:,:,iv2dd_rain),ri,rj,yobs) !#############
  CASE(id_rh_obs) ! RH
!    write(6,'(a)') 'Relative humidity obs is not supoorted'
!    stop
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,t)
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,q)
    CALL itpl_3d(v3d(:,:,:,iv3dd_p),rk,ri,rj,p)
    CALL calc_rh(t,q,p,yobs) !!! RH in q/qs, not %
!  CASE(id_tclon_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(1)
!  CASE(id_tclat_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(2)
!  CASE(id_tcmip_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(3)
  CASE DEFAULT
    qc = iqc_otype
  END SELECT

  if ( REJECT_ADPSFC_EXCEPT_PS .and. present( typ ) ) then
    if ( obtypelist(typ) == 'ADPSFC' .and. elm /= id_ps_obs ) then
      qc = iqc_otype
    endif
  endif

  RETURN
END SUBROUTINE Trans_XtoY
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY_radar(elm,radar_lon,radar_lat,radar_z,ri,rj,rk,lon,lat,lev,v3d,v2d,yobs,qc,stggrd,&
                            mv3d,slope3d,ref_add,use_shift)
  use scale_atmos_grid_cartesC_index, only: &
      KHALO
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk,radar_lon,radar_lat,radar_z !!!!! Use only, ri, rj, rk eventually... (radar_lon,lat,z in ri,rj,rk)
  REAL(r_size),INTENT(IN) :: lon,lat,lev
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER,INTENT(OUT) :: qc
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  real(r_size), intent(in), optional :: mv3d(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(in), optional :: slope3d(nlevh,nv3dd)
  real(r_size), intent(out), optional :: ref_add
  logical, intent(in), optional :: use_shift
  INTEGER :: stggrd_ = 0

  REAL(r_size) :: qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr !,rhr
  REAL(r_size) :: dist , dlon , dlat , az , elev , radar_ref,radar_rv

  real(r_size) :: utmp, vtmp

  real(RP) :: lon_RP(1,1), lat_RP(1,1)
  real(RP) :: rotc_RP(2)

  real(r_size) :: pert, mean, slope
  integer :: iv3d, k

  logical :: use_shift_ = .true.

!  integer :: ierr
!  REAL(r_dble) :: rrtimer00,rrtimer
!  rrtimer00 = MPI_WTIME()


  if (present(stggrd)) stggrd_ = stggrd
  if ( present( use_shift ) ) use_shift_ = use_shift

  yobs = undef
  qc = iqc_good

  if (stggrd_ == 1) then
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5_r_size,rj,ur)  !###### should modity itpl_3d to prevent '1.0' problem....??
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5_r_size,vr)  !######
    CALL itpl_3d(v3d(:,:,:,iv3dd_w),rk-0.5_r_size,ri,rj,wr)  !######
  else
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,ur)
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,vr)
    CALL itpl_3d(v3d(:,:,:,iv3dd_w),rk,ri,rj,wr)
  end if
  CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,tr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_p),rk,ri,rj,pr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,qvr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qc),rk,ri,rj,qcr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qr),rk,ri,rj,qrr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qi),rk,ri,rj,qir)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qs),rk,ri,rj,qsr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qg),rk,ri,rj,qgr)
!


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:itpl_3d:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !Compute az and elevation for the current observation.
  !Simple approach (TODO: implement a more robust computation)

  !Azimuth
  dlon=lon-radar_lon
  dlat=lat-radar_lat
  IF ( dlon == 0.0d0 .and. dlat == 0.0d0  )THEN
!      WRITE(6,*)'OOPS',dlon,dlat,lon,lat,radar_lon,radar_lat
    qc = iqc_out_h
    RETURN
  ELSE
    az = rad2deg*atan2(dlon*cos(radar_lat*deg2rad),dlat)
  ENDIF
  !WRITE(6,*)dlon,dlat,lon,lat,radar_lon(ityp),radar_lat(ityp),dlon*cos(radar_lat(ityp))
  !IF( abs(dlon) > maxdlon )maxdlon=abs(dlon)
  !IF( abs(dlat) > maxdlat )maxdlat=abs(dlat)
  !WRITE(6,*)maxdlon,maxdlat
  IF( az < 0) az = 360.0d0 + az
  !elevation
  CALL com_distll_1(lon,lat,radar_lon,radar_lat,dist)
  elev=rad2deg*atan2(lev-radar_z,dist)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:radar_coord:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !Check that the azimuth and elevation angles are within the expected range.
  !Some grid points may be at the radar location.
  !IF( .NOT. ( az .GT. 0.0d0 .AND. az .LT. 360.0d0 ) )RETURN
  !IF( .NOT. ( elev .GT. 0.0d0 .AND. elev .LT. 90.0d0 ))RETURN



  !DEBUG---------------------------------------------------------------
  !WRITE(6,*)'RADAR lat ',radar_lat(ityp),' RADAR lon ',radar_lon(ityp)
  !WRITE(6,*)'OBS   lat ',dlat      ,' OBS   lon ',dlon
  !WRITE(6,*)'AZIMUTH   ',az
  !WRITE(6,*)'RADAR z   ',radar_z ,'  OBS   z   ',lev
  !WRITE(6,*)'ELEVATION ',elev
  !DEGUB---------------------------------------------------------------

  !WRITE(6,*)'BCRV',dlon,dlat,az,elev

  CALL calc_ref_vr(qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr,az,elev,radar_ref,radar_rv)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:calc_ref_vr:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  SELECT CASE (elm)
  CASE(id_radar_ref_obs,id_radar_ref_zero_obs)
!!!!    if (radar_ref < MIN_RADAR_REF) then
!!!!      !In this case we will replace the observation by -RH
!!!!      !This allows us to use pseudo rh observations in some cases.
!!!!      !Later we will take the decision on what to do with these cases...

!!!!      CALL itpl_3d(v3d(:,:,:,iv3dd_rh),rk,ri,rj,rhr)
!!!!!      qc = 
!!!!      yobs = -rhr

!!!!    else                      !!!!!! --------- Pesudo RH: TO BE DONE...
    if (radar_ref < MIN_RADAR_REF) then
      qc = iqc_ref_low
      if ( use_shift_ ) then
        yobs = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT  !!! even if the above qc is bad, still return the value
      else
        if ( radar_ref > 0.0_r_size ) then
          yobs = 10.0d0 * log10(radar_ref)
        else
          yobs = -10.0_r_size
        endif
      endif
    else
      yobs = 10.0d0 * log10(radar_ref)
    end if
!!!!    end if
  CASE(id_radar_vr_obs)
    if (radar_ref < MIN_RADAR_REF) then
      qc = iqc_ref_low
    end if
    yobs = radar_rv  !!! even if the above qc is bad, still return the value
  CASE(id_radar_prh_obs)
    CALL calc_rh(tr,qvr,pr,yobs) !!! RH in q/qs, not %

  CASE DEFAULT
    qc = iqc_otype
  END SELECT

  if ( RADAR_ADDITIVE_Y18 ) then
    if ( present( mv3d ) .and. present( slope3d ) .and. present( ref_add ) ) then
      ref_add = 0.0_r_size
      do iv3d = 1, nv3dd
        call itpl_3d( mv3d(:,:,:,iv3d), rk, ri, rj, mean )
        call itpl_3d(  v3d(:,:,:,iv3d), rk, ri, rj, pert )
        pert = pert - mean
        call itpl_1d( slope3d(:,iv3d), rk, slope )
        ref_add = ref_add + pert * slope
      enddo
    endif
  endif

!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:conversion:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  RETURN
END SUBROUTINE Trans_XtoY_radar
!-----------------------------------------------------------------------
! TC center search
!  [AUTHORS:] T. Miyoshi and M. Kunii
!-----------------------------------------------------------------------
!SUBROUTINE tctrk(ps,t2,ri,rj,trk)
!  IMPLICIT NONE
!  INTEGER,PARAMETER :: isearch = 10 !search radius [grid points]
!  REAL(r_size),INTENT(IN) :: ps(nlon,nlat)
!  REAL(r_size),INTENT(IN) :: t2(nlon,nlat)
!  REAL(r_size),INTENT(IN) :: ri,rj
!  REAL(r_size),INTENT(OUT) :: trk(3) !1:lon, 2:lat, 3:minp
!  REAL(r_size) :: wk(100,3)
!  REAL(r_size) :: slp(nlon,nlat)
!  REAL(r_size) :: p1,p2,p3,p4,p5,a,c,d,xx,yy
!  INTEGER :: i,j,i0,i1,j0,j1,n

!  i0 = MAX(1,FLOOR(ri)-isearch)
!  i1 = MIN(nlon,CEILING(ri)+isearch)
!  j0 = MAX(1,FLOOR(rj)-isearch)
!  j1 = MIN(nlat,CEILING(rj)+isearch)
!  trk = undef

!  DO j=j0,j1
!    DO i=i0,i1
!      slp(i,j) = ps(i,j) * (1.0d0 - 0.0065d0 * phi0(i,j) / &
!        & (t2(i,j) + 0.0065d0 * phi0(i,j))) ** -5.257d0
!    END DO
!  END DO

!  n=0
!  DO j=j0+1,j1-1
!    DO i=i0+1,i1-1
!      IF(slp(i,j) > slp(i  ,j-1)) CYCLE
!      IF(slp(i,j) > slp(i  ,j+1)) CYCLE
!      IF(slp(i,j) > slp(i-1,j  )) CYCLE
!      IF(slp(i,j) > slp(i+1,j  )) CYCLE
!      IF(slp(i,j) > slp(i-1,j-1)) CYCLE
!      IF(slp(i,j) > slp(i+1,j-1)) CYCLE
!      IF(slp(i,j) > slp(i-1,j+1)) CYCLE
!      IF(slp(i,j) > slp(i+1,j+1)) CYCLE
!      p1 = slp(i,j)
!      p2 = slp(i-1,j)
!      p3 = slp(i+1,j)
!      p4 = slp(i,j-1)
!      p5 = slp(i,j+1)
!      c = (p3-p2)*0.5d0
!      d = (p5-p4)*0.5d0
!      a = (p2+p3+p4+p5)*0.25d0 - p1
!      IF(a == 0.0d0) CYCLE
!      xx = -0.5d0 * c / a
!      yy = -0.5d0 * d / a
!      n = n+1
!      wk(n,3) = p1 - a*(xx*xx + yy*yy)
!      wk(n,2) = lat(i,j) * (1.0d0 - yy) + lat(i,j+1) * yy
!      wk(n,1) = lon(i,j) * (1.0d0 - xx) + lon(i+1,j) * xx
!    END DO
!  END DO

!  j=1
!  IF(n > 1) THEN
!    a = wk(1,3)
!    DO i=2,n
!      IF(wk(i,3) < a) THEN
!        a = wk(i,3)
!        j = i
!      END IF
!    END DO
!  END IF
!  trk = wk(j,:)

!END SUBROUTINE tctrk
!-----------------------------------------------------------------------
! Compute relative humidity (RH)
!-----------------------------------------------------------------------
SUBROUTINE calc_rh(t,q,p,rh)
  IMPLICIT NONE
  REAL(r_size),PARAMETER :: t0=273.15d0
  REAL(r_size),PARAMETER :: e0c=6.11d0
  REAL(r_size),PARAMETER :: al=17.3d0
  REAL(r_size),PARAMETER :: bl=237.3d0
  REAL(r_size),PARAMETER :: e0i=6.1121d0
  REAL(r_size),PARAMETER :: ai=22.587d0
  REAL(r_size),PARAMETER :: bi=273.86d0
  REAL(r_size),INTENT(IN) :: t,q,p
  REAL(r_size),INTENT(OUT) :: rh
  REAL(r_size) :: e,es,tc

  e = q * p * 0.01d0 / (0.378d0 * q + 0.622d0)

  tc = t-t0
  IF(tc >= 0.0d0) THEN
    es = e0c * exp(al*tc/(bl+tc))
  ELSE IF(tc <= -15.d0) THEN
    es = e0i * exp(ai*tc/(bi+tc))
  ELSE
    es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
       + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
  END IF

  rh = e/es

  RETURN
END SUBROUTINE calc_rh

!-----------------------------------------------------------------------
! Pressure adjustment for a different height level
!-----------------------------------------------------------------------
SUBROUTINE prsadj(p,dz,t,q)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: p
  REAL(r_size),INTENT(IN) :: dz ! height difference (target - original) [m]
  REAL(r_size),INTENT(IN) :: t  ! temperature [K] at original level
  REAL(r_size),INTENT(IN) :: q  ! humidity [kg/kg] at original level
  REAL(r_size),PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]
  REAL(r_size) :: tv

  IF(dz /= 0) THEN
    tv = t * (1.0d0 + 0.608d0 * q)
    p = p * ((-gamma*dz+tv)/tv)**(gg/(gamma*rd)) !tv is at original level
!    p = p * (tv/(tv+gamma*dz))**(gg/(gamma*rd)) !tv is at target level
  END IF

  RETURN
END SUBROUTINE prsadj

!-----------------------------------------------------------------------
! Compute radar reflectivity and radial wind.
! Radial wind computations for certain methods depend on model reflectivity
! so both functions has been merged into a single one.
! First reflectivity is computed, and the the radial velocity is computed.
!-----------------------------------------------------------------------
SUBROUTINE calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,az,elev,ref,vr)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: qv        !Water vapor
  REAL(r_size), INTENT(IN) :: qc,qr     !Cloud and rain water
  REAL(r_size), INTENT(IN) :: qci,qs,qg !Cloud ice, snow and graupel
  REAL(r_size), INTENT(IN) :: t,p       !Temperature and pressure.
  REAL(r_size), INTENT(IN) :: u,v,w     !velocities with respecto to earth.
  REAL(r_size), INTENT(INOUT) :: ref      !Reflectivity.
  REAL(r_size)              :: ro
  REAL(r_size), INTENT(IN) :: az     !Azimuth respect to the radar.
  REAL(r_size), INTENT(IN) :: elev   !Elevation angle respect to the surface.
  REAL(r_size), INTENT(INOUT) :: vr    !Radial velocity.
  REAL(r_size)  :: qms , qmg !Melting species concentration (METHOD_REF_CALC 3)
  REAL(r_size)  :: qt        !Total condensate mixing ratio (METHOD_REF_CALC 1)
  REAL(r_size)  :: zr , zs , zg !Rain, snow and graupel's reflectivities.
  REAL(r_size)  :: wr , ws , wg !Rain, snow and graupel's mean terminal velocities.
  REAL(r_size)  :: wt           !Total mean terminal velocity.
  REAL(r_size)  :: nor, nos, nog !Rain, snow and graupel's intercepting parameters.
  REAL(r_size)  :: ror, ros, rog , roi !Rain, snow and graupel, ice densities.
  REAL(r_size)  :: a,b,c,d,Cd    !Constant for fall speed computations.
  REAL(r_size)  :: cr_t08,cs_t08,cg_t08,dr_t08,ds_t08,dg_t08    !Constant for fall speed computations (Tomita2008)
  REAL(r_size)  :: cf, pip , roo
  REAL(r_size)  :: ki2 , kr2
  REAL(r_size)  :: lr , ls , lg
  REAL(r_size)  :: tmp_factor , rofactor
  REAL(r_size)  :: p0
  REAL(r_size)  :: Fs, Fg , zms , zmg , fws , fwg !METHOD_REF_CALC 3
  REAL(r_size)  :: qrp , qsp , qgp
  REAL(r_size)  :: maxf                     !Maximum mixture relative concentration. (METHOD_REF_CALC 3)

  REAL(r_size) , parameter :: qeps = 1.0d-20 !Avoid overflow

  REAL(r_size) , parameter :: as_RS14 = 6.9d-2
  REAL(r_size)  :: tc , MOMs_0bs

  !Note: While equivalent reflectivity is assumed to be independent of the radar, in 
  !practice short wavelengths as those associated with K band radars migh frequently
  !experience Mie scattering. In that case, the equivalent reflectivity is not longer
  !radar independent and an appropiate relationship between the forecasted concentrations
  !and the reflectivity should be used.


  !REAL(r_size)  :: trqr,trqs,trqg
  ![P] Pa
  ![T] K
  ![q...] dimensionless
  ![ref] mm^6/m^3
  ![refdb] refdb
  ![U,V,W] m/s
  ![az, elev] degree
  ![vr] m/s

  !This model reflectivity won't be lower than this value.

  !Initialize reflectivities
  zr=0.0d0
  zs=0.0d0
  zg=0.0d0
  zms=0.0d0
  zmg=0.0d0
  ref=0.0d0

  !Compute air density (all methods use this)

  ro =  p / (rd * t)

  !Begin computation of reflectivity and vr

  if (METHOD_REF_CALC == 1) then

    !WRF method: See for example Sugimoto et al. Evaluation of the Performance of Ra
    !dial Velocity Assimilation with WRF-3DVAR System and Simulated Multiple-Doppler
    !Radar Data
    !Only rain is used to estimate the terminal velocity of hidrometeors.
    !Only rain is used to compute equivalent reflectivity.
    !Marshall-Palmer distributions are assumed to find the relationship bestween
    !concentration and equivalent reflectivity.
    ! Sun and Crook 1997 , 1998.
    !Derived for C-band radars.
    !Attenuation is not computed.

    !Reflectivity
    nor=8.0d6      ![m^-4]
    ror=1000.0d0   ![Kg/m3]
    pip=pi ** 1.75 !factor
    cf =10.0d18 * 72 !factor
    p0=1.0d5            !Reference pressure.

    qt=qr + qs + qg  !Assume that the total condensate is rain water
                     !But ignore cloud ice and cloud water

    IF( qt .GT. qeps )THEN
    ref = cf * ( ( ro * qt )**1.75 )
    ref = ref / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
   !ref= 2.04d4 *( ( ro * qt * 1.0d3 ) ** 1.75 ) !Original Sun and Crook expresion.
    ELSE
    ref=0.0d0
    ENDIF

    !Radial wind

    IF ( qt .GT. qeps )THEN
    a=(p0/p)**0.4
    wt = 5.40d0 * a * ( qt ** 0.125 )
    ELSE
    wt=0d0
    ENDIF
   !WRITE(6,*)qr,qs,qg


  else if (METHOD_REF_CALC == 2) then
    !Observation operator from Tong and Xue 2006, 2008 a and b.
    !Based on Smith et al 1975.
    !It includes reflectivity contribution by all the microphisical species.
    !is assumes Marshall and Palmer distributions.
    !Based on C band radars.
    !Attenuation is not computed.
    nor=8.0d6      ![m^-4]
    nos=3.0d6      ![m^-4]
    nog=4.0d4      ![m^-4]
    ror=1000.0d0   ![Kg/m3]
    ros=100.0d0    ![Kg/m3]
    rog=913.0d0    ![Kg/m3] 
    roi=917.0d0    ![Kg/m3]
    roo=1.0d0      ![Kg/m3] Surface air density.
    ki2=0.176d0    !Dielectric factor for ice.
    kr2=0.930d0    !Dielectric factor for water.
    pip=pi ** 1.75 !factor
    cf =1.0d18 * 720 !factor 

    IF( qr .GT. qeps )THEN
    zr= cf * ( ( ro * qr )**1.75 )
    zr= zr / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
    ENDIF
    !The contribution of snow depends on temperature (bright band effect)
    IF( qs .GT. qeps )THEN
    IF ( t <= 273.16 )THEN
     zs = cf * ki2 * ( ros ** 0.25 ) * ( ( ro * qs ) ** 1.75 )
     zs = zs / ( pip * kr2 * ( nos ** 0.75  ) * ( roi ** 2 ) )
    ELSE
    !WARNING: This formulation has to be checked the paper says that 
     !ros instead of roi should be used in thes expresion, but that 
     !leads to unrealistic jumps in the model derived reflectivity.
     zs = cf * ( ( ro * qs ) ** 1.75 )
     zs = zs / ( pip * ( nos ** 0.75 ) * ( roi ** 1.75 ) )
    ENDIF
    ENDIF

    !Only dry graupel contribution is ussed.
    IF( qg .GT. qeps )THEN
    zg= ( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) ** 0.95
    zg= zg * ( ( ro * qg ) ** 1.6625 )

    !zg=(ki2/kr2)*( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) 
    !zg= zg * ( ( ro * qg ) ** 1.75 )

    !if( qg * ro * 1000 .GT. 1.0d0 )THEN
    !  WRITE(6,*)qg, zg, ro, rd , t , p, 10*log10(zg)
    !endif
    ENDIF

    ref = zr + zs + zg

    !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983.
    IF( ref > 0.0d0 )THEN
    !There are hidrometeors, compute their terminal velocity.
    !Change units to be consistent with Lin et al 1983 and
    !to obtain wt in m/s
    nor=nor*1e-3      ![cm^-4]
    nos=nos*1e-3      ![cm^-4]
    nog=nog*1e-3      ![cm^-4]
    ror=ror*1e-3        ![g/cm3]
    ros=ros*1e-3        ![g/cm3]
    rog=rog*1e-3      ![g/cm3] 
    roo=roo*1e-3      ![g/cm3] Surface air density.
    ro= ro*1e-3

    a=2115d0   ![cm**1-b / s]
    b=0.8d0
    c=152.93d0 ![cm**1-b / s]
    d=0.25d0
    Cd=0.6d0

    rofactor= ( roo / ro  ) ** 0.25

    if( qr > qeps )then
      CALL com_gamma( 4.0_r_size + b , tmp_factor )
      lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
      wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
      wr= 1.0d-2*wr * rofactor
    else
      wr = 0.0d0
    endif

    if( qs > qeps )then
      CALL com_gamma( 4.0_r_size + d , tmp_factor )
      ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
      ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
      ws= 1.0d-2*ws * rofactor
    else
      ws = 0.0d0
    endif
 
    if( qg > qeps )then
      CALL com_gamma( 4.5_r_size , tmp_factor )
      lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
      wg= tmp_factor * ( ( ( 4.0d0 * gg * 100.0d0 * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
      wg= 1.0d-2*wg / ( 6.0d0 * ( lg ** 0.5 ) )
    else
      wg = 0.0d0
    endif

    !Reflectivity weighted terminal velocity. 
    wt = ( wr * zr + ws * zs + wg * zg )/ ( zr + zs + zg )

    ELSE

    wt=0.0d0

    ENDIF

  else if (METHOD_REF_CALC == 3) then
    !Observation operator from Xue et al 2007
    !Asumes power law between Z and mass concentration of different 
    !hydrometeor categories. 
    !Derived for X-Band radars tacking into account Mie scattering.
    !Includes a computation of the attenuation. (k)

    MAXF=0.5d0

    !First we need to compute the mixtures between rain, snow and hail
    !Following Jung et al 2007 eq (2) and (3)
    Fg=0.0d0
    Fs=0.0d0
    fwg=0.0d0
    fws=0.0d0
    IF( qr .GT. 0.0d0 .AND. qg .GT. 0.0d0)THEN
      Fg=MAXF * ( min( qr/qg , qg/qr ) )**(1.0d0/3.0d0)
      fwg= qr / ( qr + qg )
    ENDIF
    IF( qr .GT. 0.0d0 .AND. qs .GT. 0.0d0)THEN
      Fs=MAXF * ( min( qr/qs , qs/qr ) )**(1.0d0/3.0d0)
      fws= qr / ( qr + qs )
    ENDIF

    if ( .not. USE_METHOD3_REF_MELT ) then
      Fs = 0.0_r_size
      Fg = 0.0_r_size
    endif

    !Correct the rain, snow and hail mixing ratios assuming
    !that we have a mixture due to melting.

    qrp=(1.0d0-Fs-Fg)*qr

    qsp=(1.0d0-Fs)*qs

    qgp=(1.0d0-Fg)*qg

    !Compute the new species resulting from melting.

    qms=Fs * (qr + qs) !Melting snow concentration.

    qmg=Fg * (qr + qg) !Melting hail concentration.

    !Compute reflectivities for each species including the melting species.

    IF( qrp .GT. qeps )THEN
    zr= 2.53d4 * ( ro * qrp * 1.0d3 )**1.84
    ENDIF
    IF( qsp .GT. qeps )THEN
    zs= 3.48d3 * ( ro * qsp * 1.0d3 )**1.66
    ENDIF
    IF( qgp .GT. qeps )THEN
!!!    zg= 8.18d4 * ( ro * qgp * 1.0d3 )**1.50  !!! hail
    zg= 5.54d3 * ( ro * qgp * 1.0d3 )**1.70   !!! graupel (A. Amemiya 2019.5)
    ENDIF
    IF( qms .GT. qeps )THEN
    zms=( 0.00491 + 5.75*fws - 5.588*(fws**2) )*1.0d5
    zms= zms * ( ro * qms * 1.0d3 )**( 1.67 - 0.202*fws + 0.398*(fws**2) )

    ENDIF
    IF( qmg .GT. qeps )THEN
!!!    zmg=( 0.809 + 10.13*fwg -5.98*(fwg**2) )*1.0d5
!!!    zmg= zmg * ( ro * qmg * 1.0d3 )**( 1.48 + 0.0448*fwg - 0.0313*(fwg**2) ) !!! hail
    zmg=( 0.0358 + 5.27*fwg -9.51*(fwg**2) + 4.68 *(fwg**3) )*1.0d5
    zmg= zmg * ( ro * qmg * 1.0d3 )**( 1.70 + 0.020*fwg + 0.287 * (fwg**2) - 0.186*(fwg**3) ) !!! graupel (A. Amemiya 2020)
    ENDIF

    ref = zr +  zg  + zs + zms + zmg

    !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983. (The distribution parameters are 
    !consistent with the work of Jung et al 2007)

    IF( ref > 0.0d0 )THEN
      !There are hidrometeors, compute their terminal velocity.
      !Units according to Lin et al 1983.
      nor=8.0d-2      ![cm^-4]
      nos=3.0d-2      ![cm^-4]
!      nog=4.0d-4      ![cm^-4]
      nog=4.0d-2      ![cm^-4]   !!! Tomita 2008
      ror=1.0d0        ![g/cm3]
      ros=0.1d0        ![g/cm3]
!      rog=0.917d0      ![g/cm3]
      rog=0.400d0      ![g/cm3]   !!! Tomita 2008
!      roo=0.001d0      ![g/cm3] Surface air density.
      roo=1.28 * 0.001d0      ![g/cm3] Surface air density. !!! Tomita 2008
      ro=1.0d-3 * ro
!      a=2115d0   ![cm**1-b / s]
!      b=0.8d0
!      c=152.93d0 ![cm**1-b / s]
!      d=0.25d0
      Cd=0.6d0  !!! drag_g in SCALE

      cr_t08=0.5d0 ![m**1-b / s]  !!! Tomita 2008 but not used in SCALE
      cr_t08=130.0d0 ![m**1-b / s]     !!! SCALE default
      dr_t08=0.5d0
      cs_t08=4.84d0 ![m**1-b / s]
      ds_t08=0.25d0
      cg_t08=82.5d0 ![m**1-b / s]   !!! Tomita 2008 but not used in SCALE
      dg_t08=0.5d0  !!! SCALE default

      rofactor= ( roo / ro  ) ** 0.5

      IF ( qr .GT. qeps )THEN
      lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25 !!! [cm ^-1]
      CALL com_gamma( 4.0_r_size + dr_t08 , tmp_factor )
      wr= cr_t08 * tmp_factor / ( 6.0d0 * ( ( lr * 1.0e2 ) ** dr_t08 ) ) !!! [m/s]
      wr= wr * rofactor
      ELSE
      wr=0.0d0
      ENDIF

      IF( qs .GT. qeps )THEN
        IF  ( USE_T08_RS2014 ) then
          tc=min(-0.1_r_size, t-273.15_r_size)
          MOMs_0bs = exp(log(10.0_r_size) * loga_(tc,2.0_r_size) + log(ro * qs / as_RS14) * b_(tc,2.0_r_size) )
          ws = cs_t08 * rofactor * exp(log(10.0_r_size) * loga_(tc,2.25_r_size) + log(ro * qs / as_RS14) * b_(tc,2.25_r_size) ) / MOMs_0bs
        ELSE
          ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25 !!! [cm ^-1]
          CALL com_gamma( 4.0_r_size + ds_t08 , tmp_factor )
          ws= cs_t08 * tmp_factor / ( 6.0d0 * ( ( ls * 1.0e2 ) ** ds_t08 ) ) !!! [m/s]
          ws= ws * rofactor
        ENDIF
      ELSE
        ws=0.0d0
      ENDIF


      IF ( qg .GT. qeps )THEN
      lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
      CALL com_gamma( 4.0_r_size + dg_t08 , tmp_factor )
!!!      wg = tmp_factor * ( ( ( 4.0d0 * gg * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
      wg = tmp_factor * ( ( ( 4.0d0 * gg * rog )/( 3.0d0 * Cd * roo ) ) ** 0.5 )  !!! fixed 2021.5.14
      wg = wg / ( 6.0d0 * ( ( lg * 1.0e2 ) ** dg_t08 ) )   !!! [m/s]
!!!      wg= cg_t08 * tmp_factor / ( 6.0d0 * ( ( lg * 1.0e2 ) ** dg_t08 ) ) !!! [m/s]
      wg= wg * rofactor
      ELSE
      wg=0.0d0
      ENDIF

      !Reflectivity weighted terminal velocity. 
      !The melting species are assumed to fail as their non-melting counterparts.
      !however this might not be the case for melting snow.
      wt = ( wr * zr + ws *  zs +  ws * zms + wg *  zg + wg * zmg ) / ( zr + zs + zg + zms + zmg )

    ELSE

      wt=0.0d0

    ENDIF

  else  !IF OVER DIFFERENT OPTIONS

    WRITE(6,*)'[Error] Not recognized method for radar reflectivity and wind computation'
    STOP

  end if ! [METHOD_REF_CALC == ?]


  !Compute radial velocity
  vr = u * cos(elev*deg2rad) * sin(az*deg2rad)
  vr = vr + v * cos(elev*deg2rad) * cos(az*deg2rad)
  IF( USE_TERMINAL_VELOCITY )THEN
    vr = vr + (w - wt)*sin(elev*deg2rad)
  ELSE
    vr = vr + (w)*sin(elev*deg2rad)
  ENDIF

  RETURN

CONTAINS
  function loga_(tems, nm)
    real(r_size) :: loga_
    real(r_size) :: tems
    real(r_size) :: nm
    real(r_size),  parameter   :: coef_a01 =  5.065339_r_size
    real(r_size),  parameter   :: coef_a02 = -0.062659_r_size
    real(r_size),  parameter   :: coef_a03 = -3.032362_r_size
    real(r_size),  parameter   :: coef_a04 =  0.029469_r_size
    real(r_size),  parameter   :: coef_a05 = -0.000285_r_size
    real(r_size),  parameter   :: coef_a06 =  0.31255_r_size
    real(r_size),  parameter   :: coef_a07 =  0.000204_r_size
    real(r_size),  parameter   :: coef_a08 =  0.003199_r_size
    real(r_size),  parameter   :: coef_a09 =  0.0_r_size
    real(r_size),  parameter   :: coef_a10 = -0.015952_r_size
    real(r_size) :: coef_at(4)

    coef_at(1) = coef_a01 + tems * ( coef_a02 + tems * ( coef_a05 + tems * coef_a09 ) )
    coef_at(2) = coef_a03 + tems * ( coef_a04 + tems *   coef_a07 )
    coef_at(3) = coef_a06 + tems *   coef_a08
    coef_at(4) = coef_a10
    loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )

  end function loga_

  function b_(tems, nm)
    real(r_size) :: b_
    real(r_size) :: tems
    real(r_size) :: nm
    real(r_size), parameter   :: coef_b01 =  0.476221_r_size
    real(r_size), parameter   :: coef_b02 = -0.015896_r_size
    real(r_size),  parameter   :: coef_b03 =  0.165977_r_size
    real(r_size),  parameter   :: coef_b04 =  0.007468_r_size
    real(r_size),  parameter   :: coef_b05 = -0.000141_r_size
    real(r_size),  parameter   :: coef_b06 =  0.060366_r_size
    real(r_size),  parameter   :: coef_b07 =  0.000079_r_size
    real(r_size),  parameter   :: coef_b08 =  0.000594_r_size
    real(r_size),  parameter   :: coef_b09 =  0.0_r_size
    real(r_size),  parameter   :: coef_b10 = -0.003577_r_size
    real(r_size) :: coef_bt(4)

    coef_bt(1) = coef_b01 + tems * ( coef_b02 + tems * ( coef_b05 + tems * coef_b09 ) )
    coef_bt(2) = coef_b03 + tems * ( coef_b04 + tems *   coef_b07 )
    coef_bt(3) = coef_b06 + tems *   coef_b08
    coef_bt(4) = coef_b10
    b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )

   end function b_


END SUBROUTINE calc_ref_vr



!-----------------------------------------------------------------------
! Coordinate conversion (find rk in pressure levels)
!
! rk = 0.0d0  : surface observation
!-----------------------------------------------------------------------
!OCL SERIAL
SUBROUTINE phys2ijk(p_full,elem,ri,rj,rlev,rk,qc,typ)
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, KHALO, IA
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: p_full(nlevh,nlonh,nlath)
  INTEGER,INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels (for 3D variable only)
  integer, intent(in), optional :: typ ! observation type
  REAL(r_size),INTENT(OUT) :: rk
  INTEGER,INTENT(OUT) :: qc
  REAL(r_size) :: ak
!  REAL(r_size) :: lnps(nlonh,nlath)
  REAL(r_size) :: lnps(nlevh,nlonh,nlath)
  REAL(r_size) :: plev(nlevh)
  REAL(r_size) :: ptmp
  INTEGER :: i,j,k, ii, jj, ks

  qc = iqc_good
!
! rlev -> rk
!
  if (ri < 1.0d0 .or. ri > nlonh .or. rj < 1.0d0 .or. rj > nlath) then
    if (IA /= 1) then !!! exception: 2D ideal case
      write(6,*) "IA=",IA
      if (LOG_LEVEL >= 1) then
        write (6,'(A)') '[Warning] observation is outside of the horizontal domain'
      end if
      rk = undef
      qc = iqc_out_h
      return
    end if
  end if
  !
  !
  IF(elem > 9999) THEN ! surface observation
    rk = rlev
  ELSE
    if ( present( typ ) ) then
      if ( obtypelist(typ) == 'ASCATW' ) then
        rk = real( KHALO, r_size )
        return
      endif
    endif
    !
    ! horizontal interpolation
    !
    i = CEILING(ri)
    j = CEILING(rj)
    !
    ! Find the lowest valid level
    !
    ks = 1+KHALO
    do jj = j-1, j
      do ii = max(i-1,1), min(i,nlong+IHALO)


!print *, p_full(:,ii,jj)

        DO k=1+KHALO,nlev+KHALO
          if (p_full(k,ii,jj) >= 0.0d0) exit
        END DO
        if (k > ks) ks = k
      end do
    end do
!    DO k=1+KHALO,nlev+KHALO
!!      IF(i <= nlon+IHALO) THEN
!        lnps(i-1:i,j-1:j) = LOG(p_full(k,i-1:i,j-1:j))
!!      ELSE
!!        lnps(i-1,j-1:j) = LOG(p_full(k,i-1,j-1:j))
!!        lnps(1,j-1:j) = LOG(p_full(k,1,j-1:j))
!!      END IF
!      CALL itpl_2d(lnps,ri,rj,plev(k))
!    END DO

    lnps(:,i-1:i,j-1:j) = LOG(p_full(:,i-1:i,j-1:j))
    call itpl_2d_column(lnps,ri,rj,plev)

    !
    ! Log pressure
    !
    rk = LOG(rlev)
    !
    ! determine if rk is within bound.
    !
    IF(rk < plev(nlev+KHALO)) THEN
      call itpl_2d(p_full(nlev+KHALO,:,:),ri,rj,ptmp)
      if (LOG_LEVEL >= 2) then
        write(6,'(A,F12.1,A,F12.1,A,I5)') '[Warning] observation is too high: ptop=', ptmp, ', lev=', rlev, ', elem=', elem
      end if
      rk = undef
      qc = iqc_out_vhi
      RETURN
    END IF
    IF(rk > plev(ks)) THEN
      call itpl_2d(p_full(ks,:,:),ri,rj,ptmp)
!print *, ks, rk, plev(ks)
      if (LOG_LEVEL >= 2) then
        write(6,'(A,F12.1,A,F12.1,A,I5)') '[Warning] observation is too low: pbottom=', ptmp, ', lev=', rlev, ', elem=', elem
      end if
      rk = undef
      qc = iqc_out_vlo

!print *, plev
!print *, elem, ri, rj, rlev, rk, qc

      RETURN
    END IF
    !
    ! find rk
    !
    DO k=ks+1,nlev+KHALO
      IF(plev(k) < rk) EXIT ! assuming descending order of plev
    END DO
    if (k == nlev+KHALO+1) then !tentative
      ak = 0.99
      rk = REAL(k-2,r_size) + ak
    else
      ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
      rk = REAL(k-1,r_size) + ak
    end if
  END IF

  RETURN
END SUBROUTINE phys2ijk
!-----------------------------------------------------------------------
! Coordinate conversion (find rk in height levels)
!
! rk = 0.0d0  : surface observation
!-----------------------------------------------------------------------
!OCL SERIAL
SUBROUTINE phys2ijkz(z_full,ri,rj,rlev,rk,qc)
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, KHALO, IA
!  use common_mpi
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: z_full(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! height levels
  REAL(r_size),INTENT(OUT) :: rk
  INTEGER,INTENT(OUT) :: qc
  REAL(r_size) :: ak
  REAL(r_size) :: zlev(nlevh)
  REAL(r_size) :: ztmp
  INTEGER :: i,j,k, ii, jj, ks

!  integer :: ierr
!  REAL(r_dble) :: rrtimer00,rrtimer
!  rrtimer00 = MPI_WTIME()


  qc = iqc_good
!
! rlev -> rk
!
  if (ri < 1.0d0 .or. ri > nlonh .or. rj < 1.0d0 .or. rj > nlath) then
    if (IA /= 1) then !!! exception: 2D ideal case
      if (LOG_LEVEL >= 1) then
        write (6,'(A)') '[Warning] observation is outside of the horizontal domain'
      end if
      rk = undef
      qc = iqc_out_h
      return
    end if
  end if

!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:check_domain:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! horizontal interpolation
  !
  i = CEILING(ri)
  j = CEILING(rj)
  !
  ! Find the lowest valid level
  !
  ks = 1+KHALO
  do jj = j-1, j
    do ii = max(i-1,1), min(i,nlong+IHALO)
      DO k=1+KHALO,nlev+KHALO
        if (z_full(k,ii,jj) > -300.0d0 .and. z_full(k,ii,jj) < 10000.0d0) exit
      END DO
      if (k > ks) ks = k
    end do
  end do


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:find_lowesr_level:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


!  DO k=1+KHALO,nlev+KHALO
!    CALL itpl_2d(z_full(k,:,:),ri,rj,zlev(k))
!  END DO

  call itpl_2d_column(z_full,ri,rj,zlev)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:itpl_2d:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! determine if rlev is within bound.
  !
  IF(rlev > zlev(nlev+KHALO)) THEN
    call itpl_2d(z_full(nlev+KHALO,:,:),ri,rj,ztmp)
    if (LOG_LEVEL >= 2) then
      write(6,'(A,F8.1,A,F8.1)') '[Warning] observation is too high: ztop=', ztmp, ', lev=', rlev
    end if
    rk = undef
    qc = iqc_out_vhi
    RETURN
  END IF
  IF(rlev < zlev(ks)) THEN
    call itpl_2d(z_full(ks,:,:),ri,rj,ztmp)
    if (LOG_LEVEL >= 2) then
      write(6,'(A,F8.1,A,F8.1)') '[Warning] observation is too low: zbottom=', ztmp, ', lev=', rlev
    end if
    rk = undef
    qc = iqc_out_vlo
    RETURN
  END IF


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:check_height:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! find rk
  !
  DO k=ks+1,nlev+KHALO
    IF(zlev(k) > rlev) EXIT ! assuming ascending order of zlev
  END DO
  if (k == nlev+KHALO+1) then !tentative
    ak = 0.99
    rk = REAL(k-2,r_size) + ak
  else
    ak = (rlev - zlev(k-1)) / (zlev(k) - zlev(k-1))
    rk = REAL(k-1,r_size) + ak
  end if


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:calc_rk:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  RETURN
END SUBROUTINE phys2ijkz
!-----------------------------------------------------------------------
! Coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE phys2ij(rlon,rlat,rig,rjg)
  use scale_atmos_grid_cartesC, only: &
      CXG => ATMOS_GRID_CARTESC_CXG, &
      CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
  use scale_mapprojection, only: &
      MAPPROJECTION_lonlat2xy
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(OUT) :: rig
  REAL(r_size),INTENT(OUT) :: rjg
  real(RP) :: rig_RP
  real(RP) :: rjg_RP
!
! rlon,rlat -> ri,rj
!
  call MAPPROJECTION_lonlat2xy( real(rlon*pi/180.0_r_size, kind=RP), &
                                real(rlat*pi/180.0_r_size, kind=RP), rig_RP, rjg_RP )
  rig = real((rig_RP - CXG(1)) / DX, kind=r_size) + 1.0_r_size
  rjg = real((rjg_RP - CYG(1)) / DY, kind=r_size) + 1.0_r_size

  if (nlonh==1) then !!! adjustment : ideal 2-D case
    write(6,*) "phys2ij : rig adjusted to 1 from",rig
    rig = 1.0_r_size
  end if

  RETURN
END SUBROUTINE phys2ij

SUBROUTINE ij2phys(rig,rjg,rlon,rlat)
  use scale_atmos_grid_cartesC, only: &
      CXG => ATMOS_GRID_CARTESC_CXG, &
      CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rig
  REAL(r_size),INTENT(IN) :: rjg
  REAL(r_size),INTENT(OUT) :: rlon ! (deg)
  REAL(r_size),INTENT(OUT) :: rlat ! (deg)
  real(RP) :: x_RP, y_RP ! (m)
  real(RP) :: rlon_RP, rlat_RP 
!
! ri,rj -> rlon,rlat
!
  x_RP = real((rig - 1.0_r_size), kind=RP)*DX + CXG(1) 
  y_RP = real((rjg - 1.0_r_size), kind=RP)*DY + CYG(1) 

  call MAPPROJECTION_xy2lonlat( x_RP, y_RP, rlon_RP, rlat_RP )

  rlon = real(rlon_RP, kind=r_size)*rad2deg
  rlat = real(rlat_RP, kind=r_size)*rad2deg

  RETURN
END SUBROUTINE ij2phys
!
!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
subroutine itpl_1d(var,rk,var5)
  implicit none
  real(r_size), intent(in) :: var(nlevh)
  real(r_size), intent(in) :: rk
  real(r_size), intent(out) :: var5
  real(r_size) :: ak
  integer :: k

  k = ceiling( rk )
  ak = rk - real(k-1,r_size)

  var5 = var(k-1) * (1 - ak) + var(k) * ak

  return
end subroutine itpl_1d
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  if (nlonh==1) then
    var5 = var(i  ,j-1) * (1-aj) &
       & + var(i  ,j  ) *    aj
  else
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
  end if

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_2d_column(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5(nlevh)
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  if (nlonh==1) then
    var5(:) = var(:,i  ,j-1) * (1-aj) &
          & + var(:,i  ,j  ) *    aj
  else
    var5(:) = var(:,i-1,j-1) * (1-ai) * (1-aj) &
          & + var(:,i  ,j-1) *    ai  * (1-aj) &
          & + var(:,i-1,j  ) * (1-ai) *    aj  &
          & + var(:,i  ,j  ) *    ai  *    aj
  end if

  RETURN
END SUBROUTINE itpl_2d_column

SUBROUTINE itpl_3d(var,rk,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: rk
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj,ak
  INTEGER :: i,j,k

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)
  k = CEILING(rk)
  ak = rk - REAL(k-1,r_size)

  if (nlonh==1) then
    var5 = var(k-1,i  ,j-1) * (1-aj) * (1-ak) &
       & + var(k-1,i  ,j  ) *    aj  * (1-ak) &
       & + var(k,  i  ,j-1) * (1-aj) *    ak  &
       & + var(k,  i  ,j  ) *    aj  *    ak
  else
    var5 = var(k-1,i-1,j-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(k-1,i  ,j-1) *    ai  * (1-aj) * (1-ak) &
       & + var(k-1,i-1,j  ) * (1-ai) *    aj  * (1-ak) &
       & + var(k-1,i  ,j  ) *    ai  *    aj  * (1-ak) &
       & + var(k,  i-1,j-1) * (1-ai) * (1-aj) *    ak  &
       & + var(k,  i  ,j-1) *    ai  * (1-aj) *    ak  &
       & + var(k,  i-1,j  ) * (1-ai) *    aj  *    ak  &
       & + var(k,  i  ,j  ) *    ai  *    aj  *    ak
  end if

  RETURN
END SUBROUTINE itpl_3d
!-----------------------------------------------------------------------
! Monitor observation departure by giving the v3dg,v2dg data
!-----------------------------------------------------------------------
subroutine monit_obs(v3dg,v2dg,topo,nobs,bias,rmse,monit_type,use_key,step)
  use scale_prc, only: &
      PRC_myrank

  implicit none

  REAL(RP),intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),intent(in) :: v2dg(nlon,nlat,nv2d)
  real(r_size),intent(in) :: topo(nlon,nlat)
  INTEGER,INTENT(OUT) :: nobs(nid_obs)
  REAL(r_size),INTENT(OUT) :: bias(nid_obs)
  REAL(r_size),INTENT(OUT) :: rmse(nid_obs)
  LOGICAL,INTENT(OUT) :: monit_type(nid_obs)
  logical,intent(in) :: use_key
  integer,intent(in) :: step

  REAL(r_size) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size) :: v2dgh(nlonh,nlath,nv2dd)
  integer :: nnobs
  integer :: n,nn
  integer :: iset,iidx
  real(r_size) :: ril,rjl,rk,rkz

  real(r_size),allocatable :: oelm(:)
  real(r_size),allocatable :: ohx(:)
  integer,allocatable :: oqc(:)

  integer :: m
  real(r_size) :: obsdep_mean

  call state_to_history(v3dg, v2dg, topo, v3dgh, v2dgh)

  if (use_key) then
    nnobs = obsda_sort%nobs_in_key
  else
    nnobs = obsda_sort%nobs
  end if

  allocate (oelm(nnobs))
  allocate (ohx(nnobs))
  allocate (oqc(nnobs))

#ifdef LETKF_DEBUG
  if (step < 0 .or. step > 2) then
    write (6, *) '[Error] monit_obs: step should be 0, 1, or 2.'
    stop
  end if
#endif
  if (step == 1) then
    obsdep_nobs = nnobs
    if ( .not. allocated( obsdep_set  ) ) allocate( obsdep_set(obsdep_nobs) )
    if ( .not. allocated( obsdep_idx  ) ) allocate( obsdep_idx(obsdep_nobs) )
    if ( .not. allocated( obsdep_qc   ) ) allocate( obsdep_qc (obsdep_nobs) )
    if ( .not. allocated( obsdep_omb  ) ) allocate( obsdep_omb(obsdep_nobs) )
    if ( .not. allocated( obsdep_oma  ) ) allocate( obsdep_oma(obsdep_nobs) )
    if ( .not. allocated( obsdep_sprd ) ) allocate( obsdep_sprd(obsdep_nobs ))
    if ( .not. allocated( obsdep_omb_emean) ) allocate( obsdep_omb_emean(obsdep_nobs) )
  end if

  oqc = -1

!  obs_idx_TCX = -1
!  obs_idx_TCY = -1
!  obs_idx_TCP = -1

!$omp parallel do private(n,nn,iset,iidx,ril,rjl,rk,rkz,m,obsdep_mean)
  do n = 1, nnobs

    if (use_key) then
      nn = obsda_sort%key(n)
    else
      nn = n
    end if

    iset = obsda_sort%set(nn)
    iidx = obsda_sort%idx(nn)

    if (step == 1) then
      obsdep_set(n) = iset
      obsdep_idx(n) = iidx
#ifdef LETKF_DEBUG
    else if (step == 2) then
      if (obsdep_set(n) /= iset) then
        write (6, *) "[Error] 'set' for y_b and y_a are inconsistent!"
        stop
      end if
      if (obsdep_idx(n) /= iidx) then
        write (6, *) "[Error] 'idx' for y_b and y_a are inconsistent!"
        stop
      end if
#endif
    end if

#ifdef LETKF_DEBUG
    if (obsda_sort%qc(nn) /= iqc_good) then
      write (6, *) "[Error] The QC value of this observation provided for monitoring is not good: ", &
                   obsda_sort%qc(nn)
      stop
    end if
#endif

    oelm(n) = obs(iset)%elm(iidx)


    call rij_g2l(PRC_myrank, obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl)
#ifdef LETKF_DEBUG
    if (PRC_myrank /= obs(iset)%rank(iidx) .or. obs(iset)%rank(iidx) == -1) then
      write (6, *) "[Error] This observation provided for monitoring does not reside in my rank: ", &
                   PRC_myrank, obs(iset)%rank(iidx), obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl
      stop
    end if
#endif

    if (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
        abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE) then

      oqc(n) = iqc_otype

      select case (OBS_IN_FORMAT(iset))
      !=========================================================================
      case (obsfmt_prepbufr)
      !-------------------------------------------------------------------------
        call phys2ijk(v3dgh(:,:,:,iv3dd_p),obs(iset)%elm(iidx), &
                      ril,rjl,obs(iset)%lev(iidx),rk,oqc(n),typ=obs(iset)%typ(iidx))
        if (oqc(n) == iqc_good) then
          call Trans_XtoY(obs(iset)%elm(iidx),ril,rjl,rk, &
                          obs(iset)%lon(iidx),obs(iset)%lat(iidx), &
                          v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1,typ=obs(iset)%typ(iidx))
        end if
      !=========================================================================
      case (obsfmt_radar, obsfmt_radar_nc )
      !-------------------------------------------------------------------------
        if (DEPARTURE_STAT_RADAR) then
          call phys2ijkz(v3dgh(:,:,:,iv3dd_hgt),ril,rjl,obs(iset)%lev(iidx),rkz,oqc(n))
          if (oqc(n) == iqc_good) then
            call Trans_XtoY_radar(obs(iset)%elm(iidx),obs(iset)%meta(1), &
                                  obs(iset)%meta(2),obs(iset)%meta(3),ril,rjl,rkz, &
                                  obs(iset)%lon(iidx),obs(iset)%lat(iidx), &
                                  obs(iset)%lev(iidx),v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
            if (oqc(n) == iqc_ref_low) oqc(n) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small
          end if
        end if
      !=========================================================================
      end select

      if (oqc(n) == iqc_good) then
        ohx(n) = obs(iset)%dat(iidx) - ohx(n)
      else
        ohx(n) = undef
      end if

      if (step == 1) then
        obsdep_qc(n) = oqc(n)
        obsdep_omb(n) = ohx(n)
      else if (step == 2) then
        if (obsdep_qc(n) == iqc_good) then ! Use the QC value of y_a only if the QC of y_b is good
          obsdep_qc(n) = oqc(n)            !
        end if                             !
        obsdep_oma(n) = ohx(n)
        obsdep_omb_emean(n) = obsda_sort%val(nn)
      end if

      if ( OBSDEP_OUT_NOQC ) then
        obsdep_qc(n) = iqc_good
      endif

!!! ensemble perturbation output
      obsdep_mean = obsda_sort%ensval(1,nn)
      do m = 2, MEMBER
        obsdep_mean = obsdep_mean + obsda_sort%ensval(m,nn)
      enddo
      obsdep_mean = obsdep_mean / real( MEMBER, kind=r_size )

      obsdep_sprd(n) = ( obsda_sort%ensval(1,nn) - obsdep_mean )**2
      do m = 2, MEMBER
        obsdep_sprd(n) = obsdep_sprd(n) + ( obsda_sort%ensval(m,nn) - obsdep_mean )**2
      enddo
      obsdep_sprd(n) = sqrt( obsdep_sprd(n) / ( MEMBER - 1 ) )

      if (LOG_LEVEL >= 3) then
        write (6, '(2I6,2F8.2,4F12.4,I3)') &
              obs(iset)%elm(iidx), &
              obs(iset)%typ(iidx), &
              obs(iset)%lon(iidx), &
              obs(iset)%lat(iidx), &
              obs(iset)%lev(iidx), &
              obs(iset)%dat(iidx), &
              obs(iset)%err(iidx), &
              ohx(n), &
              oqc(n)
      end if

    end if ! [ DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
           !   abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE ]

  end do ! [ n = 1, nnobs ]
!$omp end parallel do


  call monit_dep(nnobs,oelm,ohx,oqc,nobs,bias,rmse)

  monit_type = .false.
  monit_type(uid_obs(id_u_obs)) = .true.
  monit_type(uid_obs(id_v_obs)) = .true.
  monit_type(uid_obs(id_t_obs)) = .true.
  monit_type(uid_obs(id_tv_obs)) = .true.
  monit_type(uid_obs(id_q_obs)) = .true.
  if (nobs(uid_obs(id_rh_obs)) > nobs(uid_obs(id_q_obs)))then
    monit_type(uid_obs(id_rh_obs)) = .true.
  else
    monit_type(uid_obs(id_q_obs)) = .true.
  end if
  monit_type(uid_obs(id_ps_obs)) = .true.
  if (DEPARTURE_STAT_RADAR) then
    monit_type(uid_obs(id_radar_ref_obs)) = .true.
    monit_type(uid_obs(id_radar_ref_zero_obs)) = .true.
    monit_type(uid_obs(id_radar_vr_obs)) = .true.
!    monit_type(uid_obs(id_radar_prh_obs)) = .true.
  end if

  deallocate (oelm)
  deallocate (ohx)
  deallocate (oqc)

  return
end subroutine monit_obs
!-----------------------------------------------------------------------
! Monitor departure
!  ofmt: output format
!    0: U,V,T(Tv),Q,RH,PS (default)
!    1: U,V,T(Tv),Q,RH,PS,RAIN
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc,nobs,bias,rmse)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  INTEGER,INTENT(OUT) :: nobs(nid_obs)
  REAL(r_size),INTENT(OUT) :: bias(nid_obs)
  REAL(r_size),INTENT(OUT) :: rmse(nid_obs)
  INTEGER :: n,i,ielm

  nobs = 0
  bias = 0.0d0
  rmse = 0.0d0

!!!!!!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i,ielm) REDUCTION(+:nobs,bias,rmse)
  DO n=1,nn
    IF(qc(n) /= iqc_good) CYCLE

    ielm = NINT(elm(n))
    if (ielm == id_tv_obs) then ! compute Tv as T
      ielm = id_t_obs
    end if
    if (ielm == id_radar_ref_zero_obs) then ! compute RE0 as REF
      ielm = id_radar_ref_obs
    end if

    i = uid_obs(ielm)
    nobs(i) = nobs(i) + 1
    bias(i) = bias(i) + dep(n)
    rmse(i) = rmse(i) + dep(n)**2
  END DO

  DO i = 1, nid_obs
    IF(nobs(i) == 0) THEN
      bias(i) = undef
      rmse(i) = undef
    ELSE
      bias(i) = bias(i) / REAL(nobs(i),r_size)
      rmse(i) = SQRT(rmse(i) / REAL(nobs(i),r_size))
    END IF
  END DO

  RETURN
END SUBROUTINE monit_dep
!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_print(nobs,bias,rmse,monit_type)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nobs(nid_obs)
  REAL(r_size),INTENT(IN) :: bias(nid_obs)
  REAL(r_size),INTENT(IN) :: rmse(nid_obs)
  LOGICAL,INTENT(IN),OPTIONAL :: monit_type(nid_obs)

  character(12) :: var_show(nid_obs)
  character(12) :: nobs_show(nid_obs)
  character(12) :: bias_show(nid_obs)
  character(12) :: rmse_show(nid_obs)

  integer :: i, n
  character(4) :: nstr

  logical :: monit_type_(nid_obs)

  monit_type_ = .true.
  if (present(monit_type)) monit_type_ = monit_type

  n = 0
  do i = 1, nid_obs
!!!!!!
!    if (monit_type_(i) .and. i /= uid_obs(id_tv_obs)) then
    if (monit_type_(i) .and. i /= uid_obs(id_tv_obs) .and. i /= uid_obs(id_radar_ref_zero_obs)) then
!!!!!!
      n = n + 1
      write(var_show(n),'(A12)') obelmlist(i)
      write(nobs_show(n),'(I12)') nobs(i)
      if (nobs(i) > 0) then
        write(bias_show(n),'(ES12.3)') bias(i)
        write(rmse_show(n),'(ES12.3)') rmse(i)
      else
        write(bias_show(n),'(A12)') 'N/A'
        write(rmse_show(n),'(A12)') 'N/A'
      end if
    end if
  end do
  write(nstr, '(I4)') n

  if ( LOG_OUT ) then
    WRITE(6,'(A,' // trim(nstr) // "('============'))") '======'
    WRITE(6,'(6x,' // trim(nstr) // 'A)')          var_show(1:n)
    WRITE(6,'(A,' // trim(nstr) // "('------------'))") '------'
    WRITE(6,'(A,' // trim(nstr) // 'A)') 'BIAS  ', bias_show(1:n)
    WRITE(6,'(A,' // trim(nstr) // 'A)') 'RMSE  ', rmse_show(1:n)
    WRITE(6,'(A,' // trim(nstr) // 'A)') 'NUMBER', nobs_show(1:n)
    WRITE(6,'(A,' // trim(nstr) // "('============'))") '======'
  end if

  RETURN
END SUBROUTINE monit_print
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_info_allocate(obs, extended)
  IMPLICIT NONE
  TYPE(obs_info),INTENT(INOUT) :: obs
  logical, optional, intent(in) :: extended

  call obs_info_deallocate(obs)

  ALLOCATE( obs%elm (obs%nobs) )
  ALLOCATE( obs%lon (obs%nobs) )
  ALLOCATE( obs%lat (obs%nobs) )
  ALLOCATE( obs%lev (obs%nobs) )
  ALLOCATE( obs%dat (obs%nobs) )
  ALLOCATE( obs%err (obs%nobs) )
  ALLOCATE( obs%typ (obs%nobs) )
  ALLOCATE( obs%dif (obs%nobs) )

  obs%elm = 0
  obs%lon = 0.0d0
  obs%lat = 0.0d0
  obs%lev = 0.0d0
  obs%dat = 0.0d0
  obs%err = 0.0d0
  obs%typ = 0
  obs%dif = 0.0d0

  if (present(extended)) then
    if (extended) then
      allocate( obs%ri (obs%nobs) )
      allocate( obs%rj (obs%nobs) )
      allocate( obs%rank (obs%nobs) )

      obs%ri = 0.0d0
      obs%rj = 0.0d0
      obs%rank = -1
    end if
  end if

  RETURN
END SUBROUTINE obs_info_allocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_info_deallocate(obs)
  IMPLICIT NONE
  TYPE(obs_info),INTENT(INOUT) :: obs

  IF(ALLOCATED(obs%elm)) DEALLOCATE(obs%elm)
  IF(ALLOCATED(obs%lon)) DEALLOCATE(obs%lon)
  IF(ALLOCATED(obs%lat)) DEALLOCATE(obs%lat)
  IF(ALLOCATED(obs%lev)) DEALLOCATE(obs%lev)
  IF(ALLOCATED(obs%dat)) DEALLOCATE(obs%dat)
  IF(ALLOCATED(obs%err)) DEALLOCATE(obs%err)
  IF(ALLOCATED(obs%typ)) DEALLOCATE(obs%typ)
  IF(ALLOCATED(obs%dif)) DEALLOCATE(obs%dif)

  RETURN
END SUBROUTINE obs_info_deallocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_da_value_allocate(obsda,member)
  IMPLICIT NONE
  TYPE(obs_da_value),INTENT(INOUT) :: obsda
  INTEGER,INTENT(IN) :: member

  call obs_da_value_deallocate(obsda)

  ALLOCATE( obsda%set (obsda%nobs) )
  ALLOCATE( obsda%idx (obsda%nobs) )
  ALLOCATE( obsda%key (obsda%nobs) )
  ALLOCATE( obsda%val (obsda%nobs) )
  ALLOCATE( obsda%qc  (obsda%nobs) )

  obsda%nobs_in_key = 0
  obsda%idx = 0
  obsda%key = 0
  obsda%val = 0.0d0
  obsda%qc = 0

  if (RADAR_PQV) then
    allocate( obsda%tm (obsda%nobs) )
    allocate( obsda%pm (obsda%nobs) )
    allocate( obsda%qv (obsda%nobs) )
    obsda%tm = 0.0d0
    obsda%pm = 0.0d0
    obsda%qv = 0.0d0
  end if

  if (member > 0) then
    ALLOCATE( obsda%ensval (member,obsda%nobs) )
    obsda%ensval = 0.0d0

    if (RADAR_PQV) then
      allocate( obsda%eqv (member,obsda%nobs) )
      obsda%eqv = 0.0d0
    end if
  end if

  if ( RADAR_ADDITIVE_Y18 ) then
    allocate( obsda%pert (obsda%nobs) )
    obsda%pert = 0.0_r_size
    if (member > 0) then
      allocate( obsda%epert (member,obsda%nobs) )
      obsda%epert = 0.0_r_size
    endif
  endif

  RETURN
END SUBROUTINE obs_da_value_allocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_da_value_deallocate(obsda)
  IMPLICIT NONE
  TYPE(obs_da_value),INTENT(INOUT) :: obsda

  obsda%nobs_in_key = 0

  IF(ALLOCATED(obsda%set   )) DEALLOCATE(obsda%set   )
  IF(ALLOCATED(obsda%idx   )) DEALLOCATE(obsda%idx   )
  IF(ALLOCATED(obsda%key   )) DEALLOCATE(obsda%key   )
  IF(ALLOCATED(obsda%val   )) DEALLOCATE(obsda%val   )
  IF(ALLOCATED(obsda%ensval)) DEALLOCATE(obsda%ensval)
  IF(ALLOCATED(obsda%qc    )) DEALLOCATE(obsda%qc    )

  if ( allocated(obsda%eqv ) ) deallocate( obsda%eqv )
  if ( allocated(obsda%pert ) ) deallocate( obsda%pert )
  if ( allocated(obsda%epert ) ) deallocate( obsda%epert )

  RETURN
END SUBROUTINE obs_da_value_deallocate
!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nrec,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nrec
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl),ALLOCATABLE :: wk(:)
  INTEGER :: ios
!  INTEGER :: iu,iv,it,iq,irh,ips,itc
  INTEGER :: iunit
  LOGICAL :: ex
  INTEGER :: sz

  ALLOCATE(wk(nrec))
  nn = 0
!  iu = 0
!  iv = 0
!  it = 0
!  iq = 0
!  irh = 0
!  ips = 0
!  itc = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

! get file size by reading through the entire file... too slow for big files
!-----------------------------
!    DO
!      READ(iunit,IOSTAT=ios) wk
!      IF(ios /= 0) EXIT
!!      SELECT CASE(NINT(wk(1)))
!!      CASE(id_u_obs)
!!        iu = iu + 1
!!      CASE(id_v_obs)
!!        iv = iv + 1
!!      CASE(id_t_obs,id_tv_obs)
!!        it = it + 1
!!      CASE(id_q_obs)
!!        iq = iq + 1
!!      CASE(id_rh_obs)
!!        irh = irh + 1
!!      CASE(id_ps_obs)
!!        ips = ips + 1
!!      CASE(id_tclon_obs)
!!        itc = itc + 1
!!      END SELECT
!      nn = nn + 1
!    END DO

! get file size by INQUIRE statement... may not work for some older fortran compilers
!-----------------------------
    INQUIRE(UNIT=iunit, SIZE=sz)
    IF (MOD(sz, r_sngl * (nrec+2)) /= 0) THEN
!      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
!      RETURN
      WRITE(6,'(3A)') '[Error]',cfile,': unsupported record length'
      STOP
    END IF

    READ(iunit,IOSTAT=ios) wk
    IF(ios /= 0) THEN
      WRITE(6,'(3A,I)') '[Error]',cfile,': Reading error ', ios
      STOP
    ELSE 
      IF (real(int(wk(1)),r_size) /= wk(1)) THEN
        WRITE(6,'(3A,E10.4)') '[Error]',cfile,': Invalid observation value wk(1) = ', wk(1)
        WRITE(6,'(A)')'Check observation data endian'
        STOP
      END IF
    END IF
  
    nn = sz / (r_sngl * (nrec+2))
!-----------------------------

    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
!    WRITE(6,'(A12,I10)') '          U:',iu
!    WRITE(6,'(A12,I10)') '          V:',iv
!    WRITE(6,'(A12,I10)') '      T(Tv):',it
!    WRITE(6,'(A12,I10)') '          Q:',iq
!    WRITE(6,'(A12,I10)') '         RH:',irh
!    WRITE(6,'(A12,I10)') '         Ps:',ips
!    WRITE(6,'(A12,I10)') '   TC TRACK:',itc
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF
  DEALLOCATE(wk)

  RETURN
END SUBROUTINE get_nobs

SUBROUTINE read_obs(cfile,obs)
  use scale_mapprojection, only: &
      MAPPROJECTION_lonlat2xy
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(8)
  real(RP) :: x_RP, y_RP
  INTEGER :: n,iunit

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obs%nobs
    READ(iunit) wk
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_v_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_t_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_tv_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_q_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_ps_obs)
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0 ! hPa -> Pa
    CASE(id_rh_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 0.01 ! percent input
      wk(6) = wk(6) * 0.01 ! percent input
    CASE(id_tcmip_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = real(OBSERR_TCP,kind=r_sngl)
    CASE(id_tclon_obs)
      call MAPPROJECTION_lonlat2xy( real(REAL(wk(2),kind=r_size)*pi/180.0_r_size, kind=RP),&
                                    real(REAL(wk(3),kind=r_size)*pi/180.0_r_size, kind=RP),&
                                    x_RP, y_RP )
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = real(x_RP, kind=r_sngl)
      wk(6) = real(OBSERR_TCX, kind=r_sngl)
    CASE(id_tclat_obs)
      call MAPPROJECTION_lonlat2xy( real(REAL(wk(2),kind=r_size)*pi/180.0_r_size, kind=RP),&
                                    real(REAL(wk(3),kind=r_size)*pi/180.0_r_size, kind=RP),&
                                    x_RP, y_RP )
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = real( y_RP, kind=r_sngl)
      wk(6) = real(OBSERR_TCY, kind=r_sngl)
    END SELECT
    obs%elm(n) = NINT(wk(1))
    obs%lon(n) = REAL(wk(2),r_size)
    obs%lat(n) = REAL(wk(3),r_size)
    obs%lev(n) = REAL(wk(4),r_size)
    obs%dat(n) = REAL(wk(5),r_size)
    obs%err(n) = REAL(wk(6),r_size)
    obs%typ(n) = NINT(wk(7))
    obs%dif(n) = REAL(wk(8),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs

SUBROUTINE write_obs(cfile,obs,append,missing)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL,INTENT(IN),OPTIONAL :: missing
  LOGICAL :: append_
  LOGICAL :: missing_
  REAL(r_sngl) :: wk(8)
  INTEGER :: n,iunit

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  missing_ = .true.
  IF(present(missing)) missing_ = missing

  IF(append_) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',POSITION='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF
  DO n=1,obs%nobs
    if (missing_ .or. abs(obs%dat(n) - undef) > tiny(wk(5))) then
      wk(1) = REAL(obs%elm(n),r_sngl)
      wk(2) = REAL(obs%lon(n),r_sngl)
      wk(3) = REAL(obs%lat(n),r_sngl)
      wk(4) = REAL(obs%lev(n),r_sngl)
      wk(5) = REAL(obs%dat(n),r_sngl)
      wk(6) = REAL(obs%err(n),r_sngl)
      wk(7) = REAL(obs%typ(n),r_sngl)
      wk(8) = REAL(obs%dif(n),r_sngl)
      SELECT CASE(NINT(wk(1)))
      CASE(id_u_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_v_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_t_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_tv_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_q_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_ps_obs)
        wk(5) = wk(5) * 0.01 ! Pa -> hPa
        wk(6) = wk(6) * 0.01 ! Pa -> hPa
      CASE(id_rh_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
        wk(5) = wk(5) * 100.0 ! percent output
        wk(6) = wk(6) * 100.0 ! percent output
      CASE(id_tcmip_obs)
        wk(5) = wk(5) * 0.01 ! Pa -> hPa
        wk(6) = wk(6) * 0.01 ! Pa -> hPa
      END SELECT
      WRITE(iunit) wk
    end if
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs

SUBROUTINE read_obs_da(cfile,obsda,im)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_da_value),INTENT(INOUT) :: obsda
  INTEGER,INTENT(IN) :: im
  REAL(r_sngl) :: wk(4) 
  INTEGER :: n,iunit

!  call obs_da_value_allocate(obsda)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obsda%nobs
    READ(iunit) wk
    obsda%set(n) = NINT(wk(1))
    obsda%idx(n) = NINT(wk(2))  !!!!!! will overflow......
    if (im == 0) then
      obsda%val(n) = REAL(wk(3),r_size)
    else
      obsda%ensval(im,n) = REAL(wk(3),r_size)
    end if
    obsda%qc(n) = NINT(wk(4))
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_da

SUBROUTINE write_obs_da(cfile,obsda,im,append)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_da_value),INTENT(IN) :: obsda
  INTEGER,INTENT(IN) :: im
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL :: append_
  REAL(r_sngl) :: wk(4) 
  INTEGER :: n,iunit

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  IF(append_) THEN
    IF(obsda%nobs <= 0) RETURN
    OPEN(iunit,FILE=cfile,FORM='unformatted',POSITION='append',STATUS='replace')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',STATUS='replace')
  END IF
  DO n=1,obsda%nobs
    wk(1) = REAL(obsda%set(n),r_sngl)
    wk(2) = REAL(obsda%idx(n),r_sngl)  !!!!!! will overflow......
    if (im == 0) then
      wk(3) = REAL(obsda%val(n),r_sngl)
    else
      wk(3) = REAL(obsda%ensval(im,n),r_sngl)
    end if
    wk(4) = REAL(obsda%qc(n),r_sngl)
    WRITE(iunit) wk
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_da

subroutine write_obs_dep(cfile, nobs, set, idx, qc, omb, oma, spr)
  implicit none
  character(*), intent(in) :: cfile
  integer, intent(in) :: nobs
  integer, intent(in) :: set(nobs)
  integer, intent(in) :: idx(nobs)
  integer, intent(in) :: qc(nobs)
  real(r_size), intent(in) :: omb(nobs)
  real(r_size), intent(in) :: oma(nobs)
  real(r_size), intent(in), optional :: spr(nobs)

  real(r_sngl) :: wk(12)
  integer :: n, iunit

  iunit=92

  open (iunit, file=cfile, form='unformatted', access='sequential')
  do n = 1, nobs
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
    select case (nint(wk(1)))
    case (id_u_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_v_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_t_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_tv_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_q_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_ps_obs)
      wk(5) = wk(5) * 0.01 ! Pa -> hPa
      wk(6) = wk(6) * 0.01 ! Pa -> hPa
    case (id_rh_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
      wk(5) = wk(5) * 100.0 ! percent output
      wk(6) = wk(6) * 100.0 ! percent output
    case (id_tcmip_obs)
      wk(5) = wk(5) * 0.01 ! Pa -> hPa
      wk(6) = wk(6) * 0.01 ! Pa -> hPa
    end select
    if (present(spr)) then
      wk(12) = real(spr(n), r_sngl)
      write (iunit) wk
    else
      write (iunit) wk(1:11)
    endif
  end do
  close (iunit)

  return
end subroutine write_obs_dep

SUBROUTINE get_nobs_radar(cfile,nn,radarlon,radarlat,radarz)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
!  REAL(r_sngl) :: wk(8)
  INTEGER :: nrec
  REAL(r_sngl) :: tmp
  INTEGER :: ios
!  INTEGER :: ir,iv
  INTEGER :: iunit
  LOGICAL :: ex
  INTEGER :: sz
  REAL(r_size),INTENT(OUT) :: radarlon,radarlat,radarz

  IF(RADAR_OBS_4D) THEN
    nrec = 8
  ELSE
    nrec = 7
  END IF
  nn = 0
!  iv = 0
!  ir = 0
  iunit=91

  radarlon=0.0d0
  radarlat=0.0d0
  radarz  =0.0d0

  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
!      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
!      RETURN
      WRITE(6,'(3A,I)') '[Error]',cfile,': Reading error', ios
      STOP
     ELSEIF ( tmp < 0.0_r_size .or. tmp > 360.0_r_size ) THEN
        WRITE(6,'(3A,E10.4)') '[Error]',cfile,': Invalid observation value radarlon = ', tmp
        WRITE(6,'(A)')'Check observation data endian'
        STOP
    END IF
    radarlon=REAL(tmp,r_size)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarlat=REAL(tmp,r_size)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarz=REAL(tmp,r_size)

! get file size by reading through the entire file... too slow for big files
!-----------------------------
!    DO
!      READ(iunit,IOSTAT=ios) wk(1:nrec)
!      IF(ios /= 0) EXIT
!!      SELECT CASE(NINT(wk(1)))
!!      CASE(id_radar_ref_obs,id_radar_ref_zero_obs)
!!        ir = ir + 1
!!      CASE(id_radar_vr_obs)
!!        iv = iv + 1
!!      END SELECT
!      nn = nn + 1
!    END DO
!-----------------------------

! get file size by INQUIRE statement... may not work for some older fortran compilers
!-----------------------------
    INQUIRE(UNIT=iunit, SIZE=sz)
    sz = sz - r_sngl * (1+2) * 3 ! substract the radar data header
    IF (MOD(sz, r_sngl * (nrec+2)) /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    nn = sz / (r_sngl * (nrec+2))
!-----------------------------

    WRITE(6,*)' RADAR FILE ', cfile
    WRITE(6,*)' RADAR LON = ',radarlon
    WRITE(6,*)' RADAR LAT = ',radarlat
    WRITE(6,*)' RADAR Z   = ',radarz
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
!    WRITE(6,'(A12,I10)') '   REFLECTIVITY:',ir
!    WRITE(6,'(A12,I10)') 'RADIAL VELOCITY:',iv
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs_radar

SUBROUTINE read_obs_radar(cfile,obs)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(8)
  INTEGER :: nrec
  REAL(r_sngl) :: tmp
  INTEGER :: n,iunit,ios

  IF(RADAR_OBS_4D) THEN
    nrec = 8
  ELSE
    nrec = 7
  END IF
  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  DO n=1,obs%nobs
    READ(iunit) wk(1:nrec)
    obs%elm(n) = NINT(wk(1))
    obs%lon(n) = REAL(wk(2),r_size)
    obs%lat(n) = REAL(wk(3),r_size)
    obs%lev(n) = REAL(wk(4),r_size)
    obs%dat(n) = REAL(wk(5),r_size)
    obs%err(n) = REAL(wk(6),r_size)
!    obs%typ(n) = NINT(wk(7))
    obs%typ(n) = 22
    IF(RADAR_OBS_4D) THEN
      obs%dif(n) = REAL(wk(8),r_size)
    ELSE
      obs%dif(n) = 0.0d0
    END IF
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_radar

SUBROUTINE write_obs_radar(cfile,obs,append,missing)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL,INTENT(IN),OPTIONAL :: missing
  LOGICAL :: append_
  LOGICAL :: missing_
  REAL(r_sngl) :: wk(8)
  INTEGER :: nrec
  INTEGER :: n,iunit

  IF(RADAR_OBS_4D) THEN
    nrec = 8
  ELSE
    nrec = 7
  END IF
  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  missing_ = .true.
  IF(present(missing)) missing_ = missing

  IF(append_) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',POSITION='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF
  WRITE(iunit) REAL(obs%meta(1),r_sngl)
  WRITE(iunit) REAL(obs%meta(2),r_sngl)
  WRITE(iunit) REAL(obs%meta(3),r_sngl)
  DO n=1,obs%nobs
    if (missing_ .or. abs(obs%dat(n) - undef) > tiny(wk(5))) then
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
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_radar

subroutine read_obs_all(obs)
  implicit none

  type(obs_info), intent(out) :: obs(OBS_IN_NUM)
  integer :: iof
  logical :: ex

  do iof = 1, OBS_IN_NUM
    inquire (file=trim(OBS_IN_NAME(iof)), exist=ex)
    if (.not. ex) then
      write(6,*) '[Warning] FILE ',trim(OBS_IN_NAME(iof)),' NOT FOUND'


      obs(iof)%nobs = 0
      call obs_info_allocate(obs(iof), extended=.true.) !!! check why this is necessary !!!


      cycle
    end if

    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call get_nobs(trim(OBS_IN_NAME(iof)),8,obs(iof)%nobs)
    case (obsfmt_radar)
      call get_nobs_radar(trim(OBS_IN_NAME(iof)), obs(iof)%nobs, obs(iof)%meta(1), obs(iof)%meta(2), obs(iof)%meta(3))
    case (obsfmt_radar_nc)
      call get_nobs_radar_nc(trim(OBS_IN_NAME(iof)), obs(iof)%nobs, obs(iof)%meta(1), obs(iof)%meta(2), obs(iof)%meta(3))
    case default
      write(6,*) '[Error] Unsupported observation file format!'
      stop
    end select

    write(6,'(5A,I9,A)') 'OBS FILE [', trim(OBS_IN_NAME(iof)), '] (FORMAT ', &
                         trim(OBS_IN_FORMAT(iof)), '): TOTAL ', &
                         obs(iof)%nobs, ' OBSERVATIONS'

    call obs_info_allocate(obs(iof), extended=.true.)

    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call read_obs(trim(OBS_IN_NAME(iof)),obs(iof))
    case (obsfmt_radar)
      call read_obs_radar(trim(OBS_IN_NAME(iof)),obs(iof))
    case (obsfmt_radar_nc)
      call read_obs_radar_nc( trim( OBS_IN_NAME(iof) ), obs(iof) )
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  return
end subroutine read_obs_all

subroutine write_obs_all(obs, missing, file_suffix)
  implicit none

  type(obs_info), intent(in) :: obs(OBS_IN_NUM)
  logical, intent(in), optional :: missing
  character(len=*), intent(in), optional :: file_suffix
  logical :: missing_
  integer :: iof, strlen1, strlen2
  character(200) :: filestr

  missing_ = .true.
  IF(present(missing)) missing_ = missing

  filestr=""
  do iof = 1, OBS_IN_NUM
    if (present(file_suffix)) then
      strlen1 = len(trim(OBS_IN_NAME(iof)))
      strlen2 = len(trim(file_suffix))
      write (filestr(1:strlen1),'(A)') trim(OBS_IN_NAME(iof))
      write (filestr(strlen1+1:strlen1+strlen2),'(A)') trim(file_suffix)
    else
      filestr = OBS_IN_NAME(iof)
    end if
    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call write_obs(trim(filestr),obs(iof),missing=missing_)
    case (obsfmt_radar)
      call write_obs_radar(trim(filestr),obs(iof),missing=missing_)
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  return
end subroutine write_obs_all

subroutine write_obs_dep_nc( filename, nobs, set, idx, qc, omb, oma, omb_em, sprd, nrank, cnt_rank )
  use netcdf
  use common_ncio
  implicit none

  character(len=*), intent(in) :: filename
  integer, intent(in) :: nobs
  integer, intent(in) :: set(nobs)
  integer, intent(in) :: idx(nobs)
  integer, intent(in) :: qc(nobs)
  real(r_size), intent(in) :: omb(nobs)
  real(r_size), intent(in) :: oma(nobs)
  real(r_size), intent(in) :: omb_em(nobs)
  real(r_size), intent(in) :: sprd(nobs)
  integer, intent(in) :: nrank
  integer, intent(in) :: cnt_rank(nrank)

  integer :: ncid
  integer :: dimid, dimid_rank
  integer :: dim_varid, elm_varid
  integer :: lon_varid, lat_varid
  integer :: lev_varid, dat_varid, qc_varid
  integer :: dif_varid, err_varid
  integer :: omb_varid, oma_varid, omb_em_varid, sprd_varid
  integer :: typ_varid
  integer :: dim_rank_varid, nobs_varid

  character(len=*), parameter :: DIM_NAME = "number"
  character(len=*), parameter :: ELM_NAME = "elm"
  character(len=*), parameter :: LON_NAME = "lon"
  character(len=*), parameter :: LAT_NAME = "lat"
  character(len=*), parameter :: LEV_NAME = "lev"
  character(len=*), parameter :: DAT_NAME = "dat"
  character(len=*), parameter :: QC_NAME  = "qc"
  character(len=*), parameter :: DIF_NAME = "dif"
  character(len=*), parameter :: ERR_NAME = "err"
  character(len=*), parameter :: OMB_NAME = "omb"
  character(len=*), parameter :: OMA_NAME = "oma"
  character(len=*), parameter :: OMB_EM_NAME = "omb_emean"
  character(len=*), parameter :: SPRD_NAME = "sprd"
  character(len=*), parameter :: TYP_NAME = "typ"
  character(len=*), parameter :: DIM_MEM_NAME = "member"
  character(len=*), parameter :: NOBS_NAME = "nobs_rank"

  character(len=*), parameter :: ELM_LONGNAME = "observation id"
  character(len=*), parameter :: LON_LONGNAME = "longitude"
  character(len=*), parameter :: LAT_LONGNAME = "latitude"
  character(len=*), parameter :: LEV_LONGNAME = "level"
  character(len=*), parameter :: DAT_LONGNAME = "observation data"
  character(len=*), parameter :: QC_LONGNAME  = "observation QC flag"
  character(len=*), parameter :: DIF_LONGNAME = "time difference"
  character(len=*), parameter :: ERR_LONGNAME = "observation error"
  character(len=*), parameter :: OMB_LONGNAME = "observation-minus-background"
  character(len=*), parameter :: OMA_LONGNAME = "observation-minus-analysis"
  character(len=*), parameter :: OMB_EM_LONGNAME = "observation-minus-background-ensemble-mean"
  character(len=*), parameter :: SPRD_LONGNAME = "ensemble spread in observation space"
  character(len=*), parameter :: TYP_LONGNAME = "observation platform type"
  character(len=*), parameter :: NOBS_LONGNAME = "number of assimilated observations in each rank"

  integer :: nobs_l(nobs), rank_l(nrank)
  integer :: n

  real(r_sngl) :: elm_l(nobs)
  real(r_sngl) :: lon_l(nobs), lat_l(nobs)
  real(r_sngl) :: lev_l(nobs), dat_l(nobs)
  real(r_sngl) :: dif_l(nobs), err_l(nobs)
  integer :: typ_l(nobs)

  do n = 1, nobs
    nobs_l(n) = n
 
    elm_l(n) = real( obs(set(n))%elm(idx(n)), r_sngl )
    lon_l(n) = real( obs(set(n))%lon(idx(n)), r_sngl )
    lat_l(n) = real( obs(set(n))%lat(idx(n)), r_sngl )
    lev_l(n) = real( obs(set(n))%lev(idx(n)), r_sngl )
    dat_l(n) = real( obs(set(n))%dat(idx(n)), r_sngl )
    dif_l(n) = real( obs(set(n))%dif(idx(n)), r_sngl )
    err_l(n) = real( obs(set(n))%err(idx(n)), r_sngl )

    typ_l(n) = int( obs(set(n))%typ(idx(n)) )
  enddo

  do n = 1, nrank
    rank_l(n) = n - 1
  enddo

  ! Create the file. 
  call ncio_check( nf90_create(trim(filename), nf90_clobber, ncid) )

  ! Define the dimensions. 
  call ncio_check( nf90_def_dim(ncid, DIM_NAME, nobs, dimid) ) 
  call ncio_check( nf90_def_dim(ncid, DIM_MEM_NAME, nrank, dimid_rank) ) 

  ! Define the coordinate variables. 
  call ncio_check( nf90_def_var(ncid, DIM_NAME, NF90_INT, dimid, dim_varid) )
  call ncio_check( nf90_def_var(ncid, DIM_MEM_NAME, NF90_INT, dimid_rank, dim_rank_varid) )

  ! Define the netCDF variables
  call ncio_check( nf90_def_var(ncid, ELM_NAME, NF90_REAL, dimid, elm_varid) )
  call ncio_check( nf90_def_var(ncid, LON_NAME, NF90_REAL, dimid, lon_varid) )
  call ncio_check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, dimid, lat_varid) )
  call ncio_check( nf90_def_var(ncid, LEV_NAME, NF90_REAL, dimid, lev_varid) )
  call ncio_check( nf90_def_var(ncid, DAT_NAME, NF90_REAL, dimid, dat_varid) )
  call ncio_check( nf90_def_var(ncid, QC_NAME,  NF90_INT,  dimid,  qc_varid) )
  call ncio_check( nf90_def_var(ncid, DIF_NAME, NF90_REAL, dimid, dif_varid) )
  call ncio_check( nf90_def_var(ncid, ERR_NAME, NF90_REAL, dimid, err_varid) )
  call ncio_check( nf90_def_var(ncid, TYP_NAME, NF90_INT,  dimid, typ_varid) )

  call ncio_check( nf90_def_var(ncid, OMB_NAME, NF90_REAL, dimid, omb_varid) )
  call ncio_check( nf90_def_var(ncid, OMA_NAME, NF90_REAL, dimid, oma_varid) )
  call ncio_check( nf90_def_var(ncid, OMB_EM_NAME, NF90_REAL, dimid, omb_em_varid) )
  call ncio_check( nf90_def_var(ncid, SPRD_NAME, NF90_REAL, dimid, sprd_varid) )
  call ncio_check( nf90_def_var(ncid, NOBS_NAME, NF90_INT,  dimid_rank, nobs_varid ) )

  ! Add long names for the netCDF variables
  call ncio_check( nf90_put_att(ncid, elm_varid, "long_name", ELM_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lon_varid, "long_name", LON_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lat_varid, "long_name", LAT_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lev_varid, "long_name", LEV_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, dat_varid, "long_name", DAT_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, qc_varid,  "long_name", QC_LONGNAME  ) )
  call ncio_check( nf90_put_att(ncid, dif_varid, "long_name", DIF_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, err_varid, "long_name", ERR_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, typ_varid, "long_name", TYP_LONGNAME ) )

  call ncio_check( nf90_put_att(ncid, omb_varid, "long_name", OMB_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, oma_varid, "long_name", OMA_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, omb_em_varid, "long_name", OMB_EM_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, sprd_varid, "long_name", SPRD_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, nobs_varid, "long_name", NOBS_LONGNAME ) )

  ! Add global attribute
  call ncio_check( nf90_put_att(ncid, NF90_GLOBAL, "total rank for SCALE (nprocs_d)", nrank ) )

  ! End define mode.
  call ncio_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. 
  call ncio_check( nf90_put_var(ncid, dim_varid, nobs_l ) )
  call ncio_check( nf90_put_var(ncid, dim_rank_varid, rank_l ) )

  ! Write the data.
  call ncio_check( nf90_put_var(ncid, elm_varid, elm_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lon_varid, lon_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lat_varid, lat_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lev_varid, lev_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, dat_varid, dat_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, qc_varid,  qc,    start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, dif_varid, dif_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, err_varid, err_l, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, typ_varid, typ_l, start=(/1/), &
                   count=(/nobs/) ) )

  call ncio_check( nf90_put_var(ncid, omb_varid, omb,   start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, oma_varid, oma,   start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, omb_em_varid, omb_em, start=(/1/), &
                   count=(/nobs/) ) )
  call ncio_check( nf90_put_var(ncid, sprd_varid, sprd, start=(/1/), &
                   count=(/nobs/) ) )

  call ncio_check( nf90_put_var(ncid, nobs_varid, cnt_rank, start=(/1/), &
                   count=(/nrank/) ) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine write_obs_dep_nc

subroutine write_monit_nc( nobs, bias, rmse, step, monit_type )
  use netcdf
  use common_ncio
  implicit none

  integer, intent(in) :: nobs(nid_obs)
  real(r_size), intent(in) :: bias(nid_obs)
  real(r_size), intent(in) :: rmse(nid_obs)
  logical, intent(in), optional :: monit_type(nid_obs)
  integer, intent(in) :: step

  real(r_size) :: bias_
  real(r_size) :: rmse_

  integer :: ncid
  integer :: dimid
  integer :: dim_varid

  integer :: nvarids(nid_obs) ! var id for the number of observations
  integer :: rvarids(nid_obs) ! var id for RMSE
  integer :: bvarids(nid_obs) ! var id for bias

  real, parameter :: undef_monit = -9.9999e+30

  character(len=*), parameter :: DIM_NAME = "type" ! 0: guess 
                                                   ! 1: analysis

  logical :: monit_type_(nid_obs)

  integer :: i

  monit_type_ = .true.
  if ( present(monit_type) ) monit_type_ = monit_type

  if ( step == 1 ) then
    ! Create the file. 
    call ncio_check( nf90_create( trim( DEPARTURE_STAT_OUT_BASENAME )//'.nc', nf90_clobber, ncid) )

    ! Define the dimensions. 
    call ncio_check( nf90_def_dim(ncid, DIM_NAME, 2, dimid) ) 
  
    ! Define the coordinate variables. 
    call ncio_check( nf90_def_var(ncid, DIM_NAME, NF90_INT, dimid, dim_varid) )
    call ncio_check( nf90_put_att(ncid, dim_varid, "long_name", "type (1: guess, 2: analysis)" ) )

    ! Define the netCDF variables
    do i = 1, nid_obs
      if ( .not. monit_type_(i) .or. i == uid_obs(id_tv_obs) &
           .or. i == uid_obs(id_radar_ref_zero_obs) ) then
        cycle
      endif

      call ncio_check( nf90_def_var(ncid, 'NOBS_'//trim(adjustl(obelmlist(i))), &
                       NF90_INT,  dimid, nvarids(i) ) )
      call ncio_check( nf90_def_var(ncid, 'RMSE_'//trim(adjustl(obelmlist(i))), &
                       NF90_REAL, dimid, rvarids(i) ) )
      call ncio_check( nf90_def_var(ncid, 'BIAS_'//trim(adjustl(obelmlist(i))), &
                       NF90_REAL, dimid, bvarids(i) ) )

      call ncio_check( nf90_put_att(ncid, rvarids(i) , "missing_value", undef_monit ) )
      call ncio_check( nf90_put_att(ncid, bvarids(i) , "missing_value", undef_monit ) )

      call ncio_check( nf90_put_att(ncid, rvarids(i) , "_FillValue", undef_monit ) )
      call ncio_check( nf90_put_att(ncid, bvarids(i) , "_FillValue", undef_monit ) )

    enddo

    ! Add global attribute
    call ncio_check( nf90_put_att(ncid, NF90_GLOBAL, "title", "Global RMSE and bias scores" ) )

    ! End define mode.
    call ncio_check( nf90_enddef(ncid) )

    ! Write the coordinate variable data. 
    call ncio_check( nf90_put_var(ncid, dim_varid, (/1, 2/) ) )

  else if( step == 2 ) then
    ! Open the file. 
    call ncio_check( nf90_open( trim( DEPARTURE_STAT_OUT_BASENAME )//'.nc', nf90_write, ncid) )

    ! Get variable id.
    do i = 1, nid_obs
      if ( .not. monit_type_(i) .or. i == uid_obs(id_tv_obs) &
           .or. i == uid_obs(id_radar_ref_zero_obs) ) then
        cycle
      endif

      call ncio_check( nf90_inq_varid(ncid, 'NOBS_'//trim(adjustl(obelmlist(i))), &
                       nvarids(i) ) )
      call ncio_check( nf90_inq_varid(ncid, 'RMSE_'//trim(adjustl(obelmlist(i))), &
                       rvarids(i) ) )
      call ncio_check( nf90_inq_varid(ncid, 'BIAS_'//trim(adjustl(obelmlist(i))), &
                       bvarids(i) ) )
    enddo

  endif

  ! Write the data.
  do i = 1, nid_obs
    if ( .not. monit_type_(i) .or. i == uid_obs(id_tv_obs) &
         .or. i == uid_obs(id_radar_ref_zero_obs) ) then
      cycle
    endif
  
    if ( nobs(i) < 0 ) then
      rmse_ = undef_monit
      bias_ = undef_monit
    else
      rmse_ = real( rmse(i), kind=r_sngl )
      bias_ = real( bias(i), kind=r_sngl )
    endif

    call ncio_check( nf90_put_var(ncid, nvarids(i), (/nobs(i)/), &
                     start=(/step/), count=(/1/) ) )
    call ncio_check( nf90_put_var(ncid, rvarids(i), (/rmse_/), &
                     start=(/step/), count=(/1/) ) )
    call ncio_check( nf90_put_var(ncid, bvarids(i), (/bias_/), &
                     start=(/step/), count=(/1/) ) )

  end do

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine write_monit_nc

subroutine write_obsnum_nc( filename, nctype, typ_ctype, elm_u_ctype, num_bqc, num_aqc, tnum_bqc, tnum_aqc )
  use netcdf
  use common_ncio
  implicit none

  character(len=*), intent(in) :: filename
  integer, intent(in) :: nctype
  integer, intent(in) :: typ_ctype(nctype)
  integer, intent(in) :: elm_u_ctype(nctype)
  integer, intent(in) :: num_bqc(nctype)    ! number of each obs
  integer, intent(in) :: num_aqc(nctype)    ! number of each obs
  integer, intent(in) :: tnum_bqc, tnum_aqc ! total numer

  integer :: ncid
  integer :: dimid

  integer :: dim_varid

  character(len=*), parameter :: DIM_NAME = "type" ! 0: before QC
                                                   ! 1: after QC
  integer :: ictype
  integer :: ityp, ielm_u
  integer :: varids(nctype), tvarid

  ! Create the file. 
  call ncio_check( nf90_create(trim(filename), nf90_clobber, ncid) )

  ! Define the dimensions. 
  call ncio_check( nf90_def_dim(ncid, DIM_NAME, 2, dimid) ) 

  ! Define the coordinate variables. 
  call ncio_check( nf90_def_var(ncid, DIM_NAME, NF90_INT, dimid, dim_varid) )
  call ncio_check( nf90_put_att(ncid, dim_varid, "long_name", "type (1: before QC, 2: after QC)" ) )

  ! Define the netCDF variables
  do ictype = 1, nctype
    ityp = typ_ctype(ictype)
    ielm_u = elm_u_ctype(ictype)

    call ncio_check( nf90_def_var(ncid, obtypelist(ityp)//'_'//trim(adjustl(obelmlist(ielm_u))), &
                     NF90_INT,  dimid, varids(ictype) ) )
  end do

  call ncio_check( nf90_def_var(ncid, "TOTAL", &
                   NF90_INT,  dimid, tvarid ) )

  ! Add global attribute
  call ncio_check( nf90_put_att(ncid, NF90_GLOBAL, "title", "Number of observations in the computational domain" ) )

  ! End define mode.
  call ncio_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. 
  call ncio_check( nf90_put_var(ncid, dim_varid, (/1, 2/) ) )

  ! Write the data.
  do ictype = 1, nctype
    call ncio_check( nf90_put_var(ncid, varids(ictype), (/num_bqc(ictype), num_aqc(ictype)/), &
                     start=(/1/), count=(/2/) ) )
  end do
  call ncio_check( nf90_put_var(ncid, tvarid, (/tnum_bqc, tnum_aqc/), &
                   start=(/1/), count=(/2/) ) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine write_obsnum_nc

subroutine get_nobs_radar_nc( cfile, nn, radarlon, radarlat, radarz )
  use netcdf
  use common_ncio
  implicit none

  character(*), intent(in) :: cfile
  integer, intent(out) :: nn
  real(r_size), intent(out) :: radarlon, radarlat, radarz

  integer :: ncid

  ! Open the file. 
  call ncio_check( nf90_open( trim( cfile ), nf90_nowrite, ncid ) ) 

  ! Get obs num
  call ncio_check( nf90_get_att(ncid, nf90_global, "num_valid_obs", nn ))

  ! Get radar location
  call ncio_check( nf90_get_att(ncid, nf90_global, "radar_lon", radarlon ) )
  call ncio_check( nf90_get_att(ncid, nf90_global, "radar_lat", radarlat ) )
  call ncio_check( nf90_get_att(ncid, nf90_global, "radar_hgt", radarz   ) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )
  return
end subroutine get_nobs_radar_nc

subroutine read_obs_radar_nc( cfile, obs )
  use common_nml
  use netcdf
  use common_ncio
  use scale_atmos_grid_cartesC, only: &
      DX, &
      DY
  implicit none

  character(*), intent(in) :: cfile
  type(obs_info), intent(inout) :: obs

  real, allocatable :: obs1d(:)
  real, allocatable :: x1d(:)
  real, allocatable :: y1d(:)
  real, allocatable :: z1d(:)
  real, allocatable :: tdiff1d(:) ! time difference for 4D LETKF
  integer :: nobs ! size 

  integer :: ncid
  integer :: dimid
  integer :: varid
  integer :: varid_x, varid_y, varid_z
  integer :: varid_tdiff

  integer :: start_nc(1)
  integer :: count_nc(1)

  character(len=15), parameter :: VHVAR_NAME = 'RadialVelocity'
  character(len=15), parameter :: ZHVAR_NAME = 'Reflectivity'
  character(len=15) :: VARNAME_SELECT

  integer :: obsid

  integer :: n

  real(r_size) :: radarlon, radarlat, radarz
  real(r_size) :: radar_rig, radar_rjg

  real(r_size) :: rig, rjg
  real(r_size) :: lon, lat
  real(r_size) :: oerr

  ! Open the file. 
  call ncio_check( nf90_open( trim( cfile ), nf90_nowrite, ncid ) ) 

  ! Get dimension size
  call ncio_check( nf90_inq_dimid( ncid, "nobs", dimid ) )
  call ncio_check( nf90_inquire_dimension( ncid, dimid, len=nobs ) )

  ! Get radar location
  call ncio_check( nf90_get_att(ncid, nf90_global, "radar_lon", radarlon ) )
  call ncio_check( nf90_get_att(ncid, nf90_global, "radar_lat", radarlat ) )
  call ncio_check( nf90_get_att(ncid, nf90_global, "radar_hgt", radarz   ) )

  call phys2ij( radarlon, radarlat, radar_rig, radar_rjg )

  ! Get obsid
  call ncio_check( nf90_get_att(ncid, nf90_global, "obsid", obsid ) )

  select case( obsid )
    case ( id_radar_ref_obs )
      VARNAME_SELECT = trim( ZHVAR_NAME )
      oerr = OBSERR_RADAR_REF
    case ( id_radar_vr_obs  )
      VARNAME_SELECT = trim( VHVAR_NAME )
      oerr = OBSERR_RADAR_VR
    case default
      write(6,*) 'Invalid obsid. Stop'
      stop
  end select

  ! Get variable id
  call ncio_check( nf90_inq_varid( ncid, trim( VARNAME_SELECT ), varid ) )
  call ncio_check( nf90_inq_varid( ncid, "x", varid_x ) )
  call ncio_check( nf90_inq_varid( ncid, "y", varid_y ) )
  call ncio_check( nf90_inq_varid( ncid, "z", varid_z ) )
  if ( RADAR_OBS_4D ) then
    call ncio_check( nf90_inq_varid( ncid, "tdiff", varid_tdiff ) )
  endif

  ! Allocate variable
  allocate( obs1d   (nobs) )  
  allocate( x1d     (nobs) )  
  allocate( y1d     (nobs) )  
  allocate( z1d     (nobs) )  
  allocate( tdiff1d (nobs) )  
  
  start_nc = [ 1 ]
  count_nc = [ nobs ]

  call ncio_check( nf90_get_var( ncid, varid, obs1d, &
                      start=start_nc, count=count_nc ) )

  call ncio_check( nf90_get_var( ncid, varid_x, x1d, &
                      start=start_nc, count=count_nc ) )

  call ncio_check( nf90_get_var( ncid, varid_y, y1d, &
                      start=start_nc, count=count_nc ) )

  call ncio_check( nf90_get_var( ncid, varid_z, z1d, &
                      start=start_nc, count=count_nc ) )

  if ( RADAR_OBS_4D ) then
     call ncio_check( nf90_get_var( ncid, varid_tdiff, tdiff1d, &
                         start=start_nc, count=count_nc ) )
  endif

  n = 0
  do n = 1, nobs
    rig = real( x1d(n) / DX, kind=r_size ) + radar_rig
    rjg = real( y1d(n) / DY, kind=r_size ) + radar_rjg

    call ij2phys( rig, rjg, lon, lat )
    obs%elm(n) = obsid
    obs%lon(n) = real( lon, kind=r_size )
    obs%lat(n) = real( lat, kind=r_size )
    obs%lev(n) = real( z1d(n), kind=r_size )
    if ( obsid == id_radar_ref_obs ) then
      if ( RADAR_OBS_IN_DBZ ) then
        obs%dat(n) = real( 10**(obs1d(n)*0.1), kind=r_size ) ! dBZ to dB
      else
        obs%dat(n) = real( obs1d(n), kind=r_size )           ! dB
      endif
    else
      obs%dat(n) = real( obs1d(n), kind=r_size )
    endif
    obs%err(n) = oerr
    obs%typ(n) = 22
    if ( RADAR_OBS_4D ) then
      obs%dif(n) = real( tdiff1d(n), kind=r_size )
    else
      obs%dif(n) = real( 0.0, kind=r_size )
    endif
  enddo

  deallocate( x1d, y1d, z1d )
  deallocate( obs1d )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine read_obs_radar_nc
!---------------------------
subroutine write_obs_anal_rank_nc( filename, ya )
  use netcdf
  use common_ncio
 implicit none

  character(len=*), intent(in) :: filename
  real(r_size), intent(in) :: ya(obsdep_nobs,MEMBER)

  integer :: ncid
  integer :: dimid, dimid_mem
  integer :: dim_varid, dim_varid_mem, elm_varid
  integer :: set_varid, idx_varid
  integer :: lon_varid, lat_varid
  integer :: lev_varid, dat_varid, qc_varid
  integer :: dif_varid, err_varid
  integer :: omb_varid, ombem_varid
  integer :: ya_varid
  integer :: typ_varid

  character(len=*), parameter :: DIM_NAME = "number"
  character(len=*), parameter :: DIM_NAME_MEM = "member"
  character(len=*), parameter :: SET_NAME = "set"
  character(len=*), parameter :: IDX_NAME = "idx"
  character(len=*), parameter :: ELM_NAME = "elm"
  character(len=*), parameter :: LON_NAME = "lon"
  character(len=*), parameter :: LAT_NAME = "lat"
  character(len=*), parameter :: LEV_NAME = "lev"
  character(len=*), parameter :: DAT_NAME = "dat"
  character(len=*), parameter :: QC_NAME  = "qc"
  character(len=*), parameter :: DIF_NAME = "dif"
  character(len=*), parameter :: ERR_NAME = "err"
  character(len=*), parameter :: YA_NAME = "ya"
  character(len=*), parameter :: TYP_NAME = "typ"
  character(len=*), parameter :: OMB_NAME = "omb"
  character(len=*), parameter :: OMBEM_NAME = "omb_emean"

  character(len=*), parameter :: ELM_LONGNAME = "observation id"
  character(len=*), parameter :: SET_LONGNAME = "observation set"
  character(len=*), parameter :: IDX_LONGNAME = "observation index"
  character(len=*), parameter :: LON_LONGNAME = "longitude"
  character(len=*), parameter :: LAT_LONGNAME = "latitude"
  character(len=*), parameter :: LEV_LONGNAME = "level"
  character(len=*), parameter :: DAT_LONGNAME = "observation data"
  character(len=*), parameter :: QC_LONGNAME  = "observation QC flag"
  character(len=*), parameter :: DIF_LONGNAME = "time difference"
  character(len=*), parameter :: ERR_LONGNAME = "observation error"
  character(len=*), parameter :: YA_LONGNAME = "analysis perturbation"
  character(len=*), parameter :: TYP_LONGNAME = "observation platform type"
  character(len=*), parameter :: OMB_LONGNAME = "observation-minus-background"
  character(len=*), parameter :: OMBEM_LONGNAME = "observation-minus-background-ensemble-mean"

  integer :: nobs_l(obsdep_nobs)
  integer :: mem_l(MEMBER)
  integer :: n

  integer :: set_l(obsdep_nobs)
  integer :: idx_l(obsdep_nobs)
  real(r_sngl) :: elm_l(obsdep_nobs)
  real(r_sngl) :: lon_l(obsdep_nobs), lat_l(obsdep_nobs)
  real(r_sngl) :: lev_l(obsdep_nobs), dat_l(obsdep_nobs)
  real(r_sngl) :: dif_l(obsdep_nobs), err_l(obsdep_nobs)
  real(r_sngl) :: ya_l(obsdep_nobs,MEMBER)
  real(r_sngl) :: omb_l(obsdep_nobs)
  real(r_sngl) :: omb_emean_l(obsdep_nobs)
  integer :: typ_l(obsdep_nobs)

  do n = 1, obsdep_nobs
    nobs_l(n) = n
 
    set_l(n) = int( obsdep_set(n) )
    idx_l(n) = int( obsdep_idx(n) )

    elm_l(n) = real( obs(obsdep_set(n))%elm(obsdep_idx(n)), r_sngl )
    lon_l(n) = real( obs(obsdep_set(n))%lon(obsdep_idx(n)), r_sngl )
    lat_l(n) = real( obs(obsdep_set(n))%lat(obsdep_idx(n)), r_sngl )
    lev_l(n) = real( obs(obsdep_set(n))%lev(obsdep_idx(n)), r_sngl )
    dat_l(n) = real( obs(obsdep_set(n))%dat(obsdep_idx(n)), r_sngl )
    dif_l(n) = real( obs(obsdep_set(n))%dif(obsdep_idx(n)), r_sngl )
    err_l(n) = real( obs(obsdep_set(n))%err(obsdep_idx(n)), r_sngl )

    omb_l(n) = real( obsdep_omb(n), r_sngl )
    omb_emean_l(n) = real( obsdep_omb_emean(n), r_sngl )

    ya_l(n,1:MEMBER) = real( ya(n,1:MEMBER), r_sngl )

    typ_l(n) = int( obs(obsdep_set(n))%typ(obsdep_idx(n)) )
  enddo

  do n = 1, MEMBER
    mem_l(n) = n
  enddo

  ! Create the file. 
  call ncio_check( nf90_create(trim(filename), nf90_clobber, ncid) )

  ! Add global attribute
  call ncio_check( nf90_put_att(ncid, NF90_GLOBAL, "title", "Analysis perturbation in the observation space" ) )

  ! Define the dimensions. 
  call ncio_check( nf90_def_dim(ncid, DIM_NAME, obsdep_nobs, dimid) ) 
  call ncio_check( nf90_def_dim(ncid, DIM_NAME_MEM, MEMBER, dimid_mem) ) 

  ! Define the coordinate variables. 
  call ncio_check( nf90_def_var(ncid, DIM_NAME,     NF90_INT, dimid,     dim_varid     ) )
  call ncio_check( nf90_def_var(ncid, DIM_NAME_MEM, NF90_INT, dimid_mem, dim_varid_mem ) )

  ! Define the netCDF variables
  call ncio_check( nf90_def_var(ncid, ELM_NAME, NF90_REAL, dimid, elm_varid) )
  call ncio_check( nf90_def_var(ncid, SET_NAME, NF90_INT, dimid, set_varid) )
  call ncio_check( nf90_def_var(ncid, IDX_NAME, NF90_INT, dimid, idx_varid) )
  call ncio_check( nf90_def_var(ncid, LON_NAME, NF90_REAL, dimid, lon_varid) )
  call ncio_check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, dimid, lat_varid) )
  call ncio_check( nf90_def_var(ncid, LEV_NAME, NF90_REAL, dimid, lev_varid) )
  call ncio_check( nf90_def_var(ncid, DAT_NAME, NF90_REAL, dimid, dat_varid) )
  call ncio_check( nf90_def_var(ncid, QC_NAME,  NF90_INT,  dimid,  qc_varid) )
  call ncio_check( nf90_def_var(ncid, DIF_NAME, NF90_REAL, dimid, dif_varid) )
  call ncio_check( nf90_def_var(ncid, ERR_NAME, NF90_REAL, dimid, err_varid) )
  call ncio_check( nf90_def_var(ncid, TYP_NAME, NF90_INT,  dimid, typ_varid) )
  call ncio_check( nf90_def_var(ncid, OMB_NAME, NF90_REAL, dimid, omb_varid) )
  call ncio_check( nf90_def_var(ncid, OMBEM_NAME, NF90_REAL, dimid, ombem_varid) )

  call ncio_check( nf90_def_var(ncid, YA_NAME, NF90_REAL, (/dimid,dimid_mem/), ya_varid) )

  ! Add long names for the netCDF variables
  call ncio_check( nf90_put_att(ncid, elm_varid, "long_name", ELM_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, set_varid, "long_name", set_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, idx_varid, "long_name", idx_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lon_varid, "long_name", LON_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lat_varid, "long_name", LAT_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, lev_varid, "long_name", LEV_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, dat_varid, "long_name", DAT_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, qc_varid,  "long_name", QC_LONGNAME  ) )
  call ncio_check( nf90_put_att(ncid, dif_varid, "long_name", DIF_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, err_varid, "long_name", ERR_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, typ_varid, "long_name", TYP_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, omb_varid, "long_name", OMB_LONGNAME ) )
  call ncio_check( nf90_put_att(ncid, ombem_varid, "long_name", OMBEM_LONGNAME ) )

  call ncio_check( nf90_put_att(ncid, ya_varid, "long_name", YA_LONGNAME ) )

  ! End define mode.
  call ncio_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. 
  call ncio_check( nf90_put_var(ncid, dim_varid, nobs_l ) )
  call ncio_check( nf90_put_var(ncid, dim_varid_mem, mem_l ) )

  ! Write the data.
  call ncio_check( nf90_put_var(ncid, elm_varid, elm_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, set_varid, set_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, idx_varid, idx_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lon_varid, lon_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lat_varid, lat_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, lev_varid, lev_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, dat_varid, dat_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, qc_varid,  obsdep_qc, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, dif_varid, dif_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, err_varid, err_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, omb_varid, omb_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )
  call ncio_check( nf90_put_var(ncid, ombem_varid, omb_emean_l, start=(/1/), &
                   count=(/obsdep_nobs/) ) )

  call ncio_check( nf90_put_var(ncid, ya_varid, ya_l, start=(/1,1/), &
                   count=(/obsdep_nobs,MEMBER/) ) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine write_obs_anal_rank_nc
!---------------------------
subroutine get_nobs_efso( cfile, nrank, cnt_rank )
  use netcdf
  use common_ncio
  implicit none

  character(*), intent(in) :: cfile
  integer, intent(in) :: nrank
  integer, intent(out) :: cnt_rank(nrank)

  integer :: varid_rank
  integer :: ncid

  ! Open the file. 
  call ncio_check( nf90_open( trim( cfile ), nf90_nowrite, ncid ) )

  ! Get obs num
  call ncio_check( nf90_inq_varid( ncid, "nobs_rank", varid_rank ) )

  call ncio_check( nf90_get_var( ncid, varid_rank, cnt_rank, &
                   start=(/1/), count=(/nrank/) ) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  return
end subroutine get_nobs_efso
!---------------------------
subroutine get_obsdep_efso( cfile, nobs_local, nobs0, set, idx, elm, typ, lon, lat, lev, &
                            dat, err, dif, dep, qc, ya )
  use netcdf
  use common_ncio
  implicit none

  character(*), intent(in) :: cfile
  integer, intent(in) :: nobs_local
  integer, intent(in) :: nobs0

  integer, intent(out) :: set(nobs_local)
  integer, intent(out) :: idx(nobs_local)
  integer, intent(out) :: elm(nobs_local)
  integer, intent(out) :: typ(nobs_local)
  real(r_size), intent(out) :: lon(nobs_local)
  real(r_size), intent(out) :: lat(nobs_local)
  real(r_size), intent(out) :: lev(nobs_local)
  real(r_size), intent(out) :: dat(nobs_local)
  real(r_size), intent(out) :: err(nobs_local)
  real(r_size), intent(out) :: dif(nobs_local)
  real(r_size), intent(out) :: dep(nobs_local)
  integer, intent(out) :: qc(nobs_local)
  real(r_size), intent(out) :: ya(nobs_local,MEMBER)
  integer :: ncid

  integer :: varid_elm, varid_typ
  integer :: varid_set, varid_idx
  integer :: varid_lon, varid_lat, varid_lev
  integer :: varid_dat, varid_err, varid_dif
  integer :: varid_dep, varid_qc
  integer :: varid_ya

  integer :: n

  ! Open the file. 
  call ncio_check( nf90_open( trim( cfile ), nf90_nowrite, ncid ) )

  ! Get variable id
  call ncio_check( nf90_inq_varid( ncid, "elm", varid_elm ) )
  call ncio_check( nf90_inq_varid( ncid, "set", varid_set ) )
  call ncio_check( nf90_inq_varid( ncid, "idx", varid_idx ) )
  call ncio_check( nf90_inq_varid( ncid, "typ", varid_typ ) )
  call ncio_check( nf90_inq_varid( ncid, "lon", varid_lon ) )
  call ncio_check( nf90_inq_varid( ncid, "lat", varid_lat ) )
  call ncio_check( nf90_inq_varid( ncid, "lev", varid_lev ) )
  call ncio_check( nf90_inq_varid( ncid, "dat", varid_dat ) )
  call ncio_check( nf90_inq_varid( ncid, "err", varid_err ) )
  call ncio_check( nf90_inq_varid( ncid, "dif", varid_dif ) )
  call ncio_check( nf90_inq_varid( ncid, "omb", varid_dep ) )
  call ncio_check( nf90_inq_varid( ncid, "qc",  varid_qc  ) )
  call ncio_check( nf90_inq_varid( ncid, "ya", varid_ya ) )

  ! Read variables
  call ncio_check( nf90_get_var( ncid, varid_elm, elm, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_set, set, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_idx, idx, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_typ, typ, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_lon, lon, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_lat, lat, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_lev, lev, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_dat, dat, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_err, err, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_dif, dif, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_dep, dep, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_qc, qc, &
                   start=(/1/), count=(/nobs_local/) ) )

  call ncio_check( nf90_get_var( ncid, varid_ya, ya, &
                   start=(/1,1/), count=(/nobs_local,MEMBER/) ) )

  ! Close the file. 
  call ncio_check( nf90_close(ncid) )

  ! QC should always be good.
  ! qc array stores indices of obs
  do n = 1, nobs_local
    qc(n) = n + nobs0
  enddo

  return
end subroutine get_obsdep_efso
!---------------------------

END MODULE common_obs_scale
