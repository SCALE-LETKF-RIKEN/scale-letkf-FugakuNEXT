MODULE letkf_obs
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009   Takemasa MIYOSHI  created
!   10/04/2012   Guo-Yuan Lien     modified for GFS model
!   08/30/2013   Guo-Yuan Lien     separating obs operator, following changes by Takemasa MIYOSHI
!   01/01/2014   Guo-Yuan Lien     add EFSO, following changes by Daisuke HOTTA
!   October 2014 Guo-Yuan Lien     modified for SCALE model
!   ............ See git history for the following revisions
!
!=======================================================================
!$USE OMP_LIB
  USE common
  use common_nml
  USE common_mpi
  USE common_scale
  USE common_obs_scale
  USE common_mpi_scale
  USE common_letkf

  IMPLICIT NONE
  PUBLIC

  real(r_size),parameter :: dist_zero_fac = 3.651483717        ! SQRT(10.0d0/3.0d0) * 2.0d0
  real(r_size),parameter :: dist_zero_fac_square = 13.33333333 ! dist_zero_fac * dist_zero_fac

  type(obs_da_value),save :: obsda   ! unsortted; sorted data saved to obsda_sort, declared in common_obs_scale

  integer,save :: nobs_extern

  ! combined obs type: {variable type (elm_u), report type (typ)}, allocated only when observations exist
  integer,save :: nctype                        ! number of combined obs type
  integer,save :: ctype_elmtyp(nid_obs,nobtype) ! array of ctype for each combination of (elm_u, typ)
  integer,allocatable,save :: elm_ctype(:)      ! array of elm  for each combined obs type
  integer,allocatable,save :: elm_u_ctype(:)    ! array of elm_u for each combined obs type
  integer,allocatable,save :: typ_ctype(:)      ! array of typ  for each combined obs type
  real(r_size),allocatable,save :: hori_loc_ctype(:) ! array of horizontal localization length for each combined obs type
  real(r_size),allocatable,save :: vert_loc_ctype(:) ! array of vertical localization length for each combined obs type

  integer, parameter :: n_qc_steps = 2
  integer, parameter :: i_before_qc = 1
  integer, parameter :: i_after_qc = 2

  type obs_grid_type
    integer :: ngrd_i
    integer :: ngrd_j
    real(r_size) :: grdspc_i
    real(r_size) :: grdspc_j
    integer :: ngrdsch_i
    integer :: ngrdsch_j
    integer :: ngrdext_i
    integer :: ngrdext_j
    integer, allocatable :: n(:,:,:)
    integer, allocatable :: ac(:,:,:)
    integer, allocatable :: tot(:)
    integer, allocatable :: n_ext(:,:)
    integer, allocatable :: ac_ext(:,:)
    integer :: tot_ext
    integer :: tot_sub(n_qc_steps)    ! only for diagnostic print; 1: before QC; 2: after QC
    integer :: tot_g(n_qc_steps)      ! only for diagnostic print
    integer, allocatable :: next(:,:) ! temporary array
  end type obs_grid_type

  type(obs_grid_type),allocatable,save :: obsgrd(:)

  integer,save :: nobstotalg
  integer,save :: nobstotal
!  integer,save :: maxnobs_sub_per_ctype
  integer,save :: maxnobs_per_ctype
  integer, save :: nobs_local

  integer, allocatable :: obsset(:)
  integer, allocatable :: obsidx(:)
  integer, allocatable :: obselm(:)
  integer, allocatable :: obstyp(:)
  real(r_size), allocatable :: obshdxf(:,:)
  real(r_size), allocatable :: obslon(:)
  real(r_size), allocatable :: obslat(:)
  real(r_size), allocatable :: obslev(:)
  real(r_size), allocatable :: obsdat(:)
  real(r_size), allocatable :: obserr(:)
  real(r_size), allocatable :: obsdif(:)
  real(r_size), allocatable :: obsdep(:)
  integer, allocatable :: obsqc(:)

CONTAINS
!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------
SUBROUTINE set_letkf_obs
  use scale_atmos_grid_cartesC, only: &
    DX, &
    DY

  IMPLICIT NONE
  INTEGER :: n,i,j,ierr,im,iof,iidx

  integer :: n1, n2

  integer :: mem_ref
  real(r_size) :: qvs

  integer :: it,ip
  integer :: ityp,ielm,ielm_u,ictype
  real(r_size) :: target_grdspc

  integer :: myp_i,myp_j
  integer :: ip_i,ip_j

  integer :: nobs_sub(n_qc_steps),nobs_g(n_qc_steps)

  integer :: nobs_elms(nid_obs)
  integer :: nobs_elms_sum(nid_obs)

  integer :: nobs_intern


  character(len=3) :: use_obs_print
  character(4) :: nstr
  


!---
  integer :: cnts
  integer :: cntr(nprocs_d)
  integer :: dspr(nprocs_d)
  integer :: nensobs_div, nensobs_mod
  integer :: im_obs_1, im_obs_2, nensobs_part

  integer :: ns_ext, ne_ext, ns_bufr, ne_bufr
  integer :: ishift, jshift

  type(obs_da_value) :: obsbufs, obsbufr
  integer :: imin1,imax1,jmin1,jmax1,imin2,imax2,jmin2,jmax2
!---




  real(r_size),allocatable :: tmpelm(:)
  INTEGER :: monit_nobs(nid_obs)
  REAL(r_size) :: bias(nid_obs)
  REAL(r_size) :: rmse(nid_obs)

  character(filelenmax) :: obsdafile
  character(11) :: obsda_suffix = '.000000.dat'

  type(obs_da_value) :: obsda_ext

  logical :: ctype_use(nid_obs,nobtype)

!  character(len=timer_name_width) :: timer_str

  call mpi_timer('', 2)

  if ( LOG_OUT ) WRITE(6,'(A)') 'Hello from set_letkf_obs'

  nobs_intern = obsda%nobs - nobs_extern
  if ( LOG_OUT ) then
    WRITE(6,'(A,I10)') 'Internally processed observations: ', nobs_intern
    WRITE(6,'(A,I10)') 'Externally processed observations: ', nobs_extern
    WRITE(6,'(A,I10)') 'Total                observations: ', obsda%nobs
  end if

!-------------------------------------------------------------------------------
! Read externally processed observations
!-------------------------------------------------------------------------------

  if (OBSDA_IN .and. nobs_extern > 0) then

    ! Read externally processed observations
    !---------------------------------------------------------------------------

    n1 = nobs_intern + 1
    n2 = obsda%nobs

    obsda_ext%nobs = nobs_extern
    call obs_da_value_allocate(obsda_ext,0)

    do it = 1, nitmax
      im = myrank_to_mem(it)
      if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
!        write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is reading externally processed observations for member ', &
!              im, ', subdomain id #', myrank_d
        if (im <= MEMBER) then
          obsdafile = OBSDA_IN_BASENAME
          call filename_replace_mem(obsdafile, im)
        else if (im == mmean) then
          obsdafile = OBSDA_MEAN_IN_BASENAME
        else if (im == mmdet) then
          obsdafile = OBSDA_MDET_IN_BASENAME
        end if
        write (obsda_suffix(2:7),'(I6.6)') myrank_d
        if ( LOG_OUT ) write (6,'(A,I6.6,2A)') 'MYRANK ', myrank,' is reading an externally processed obsda file ', trim(obsdafile)//obsda_suffix
        call read_obs_da(trim(obsdafile)//obsda_suffix,obsda_ext,0)

        if (OBSDA_OUT) then
!          write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is appending observations for member ', &
!                im, ', subdomain id #', myrank_d
          if (im <= MEMBER) then
            obsdafile = OBSDA_OUT_BASENAME
            call filename_replace_mem(obsdafile, im)
          else if (im == mmean) then
            obsdafile = OBSDA_MEAN_OUT_BASENAME
          else if (im == mmdet) then
            obsdafile = OBSDA_MDET_OUT_BASENAME
          end if
!          write (obsda_suffix(2:7),'(I6.6)') myrank_d
          if ( LOG_OUT ) write (6,'(A,I6.6,2A)') 'MYRANK ', myrank,' is writing (appending) an obsda file ', trim(obsdafile)//obsda_suffix
          call write_obs_da(trim(obsdafile)//obsda_suffix,obsda_ext,0,append=.true.)
        end if

        ! variables without an ensemble dimension
        if (it == 1) then
          obsda%set(n1:n2) = obsda_ext%set
          obsda%idx(n1:n2) = obsda_ext%idx
#ifdef LETKF_DEBUG
        else
          if (maxval(abs(obsda%set(n1:n2) - obsda_ext%set)) > 0) then
            write (6,'(A)') '[Error] obsda%set are inconsistent among the ensemble'
            stop 99
          end if
          if (maxval(abs(obsda%idx(n1:n2) - obsda_ext%idx)) > 0) then
            write (6,'(A)') '[Error] obsda%idx are inconsistent among the ensemble'
            stop 99
          end if
#endif
        end if

        if ( RADAR_PQV ) then
          call obs_da_value_partial_reduce_iter(obsda, it, n1, n2, obsda_ext%val, obsda_ext%qc, qv=obsda_ext%qv, tm=obsda_ext%tm, pm=obsda_ext%pm )
        elseif ( RADAR_ADDITIVE_Y18 ) then
          call obs_da_value_partial_reduce_iter(obsda, it, n1, n2, obsda_ext%val, obsda_ext%qc, pert=obsda_ext%pert)
        else
          call obs_da_value_partial_reduce_iter(obsda, it, n1, n2, obsda_ext%val, obsda_ext%qc )
        endif

      end if ! [ (im >= 1 .and. im <= MEMBER) .or. im == mmdetin ]
    end do ! [ it = 1, nitmax ]

    call obs_da_value_deallocate(obsda_ext)

    call mpi_timer('set_letkf_obs:read_external_obs:', 2, barrier=MPI_COMM_e)

    ! Broadcast the observation information shared by members (e.g., grid numbers)
    !---------------------------------------------------------------------------

    if (nprocs_e > MEMBER) then
      call MPI_BCAST(obsda%set(n1:n2), nobs_extern, MPI_INTEGER, 0, MPI_COMM_e, ierr)
      call MPI_BCAST(obsda%idx(n1:n2), nobs_extern, MPI_INTEGER, 0, MPI_COMM_e, ierr)
    end if

    call mpi_timer('set_letkf_obs:external_obs_broadcast:', 2, barrier=MPI_COMM_e)
  end if ! [ OBSDA_IN .and. nobs_extern > 0 ]

  ! Allreduce externally processed observations
  !---------------------------------------------------------------------------

  call obs_da_value_allreduce(obsda)

  call mpi_timer('set_letkf_obs:obs_allreduce:', 2)

!-------------------------------------------------------------------------------
! Process observations and quality control (QC)
!-------------------------------------------------------------------------------

  ! Pre-process data
  !-----------------------------------------------------------------------------

  ctype_use(:,:) = .false.
!$OMP PARALLEL PRIVATE(iof,n) REDUCTION(.or.:ctype_use)
  do iof = 1, OBS_IN_NUM
!$OMP DO
    do n = 1, obs(iof)%nobs
      select case (obs(iof)%elm(n))
      case (id_radar_ref_obs)
        if (obs(iof)%dat(n) >= 0.0d0 .and. obs(iof)%dat(n) < 1.0d10) then
          if (obs(iof)%dat(n) < MIN_RADAR_REF) then
            obs(iof)%elm(n) = id_radar_ref_zero_obs
            obs(iof)%dat(n) = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
          else
            obs(iof)%dat(n) = 10.0d0 * log10(obs(iof)%dat(n))
          end if
        else
          obs(iof)%dat(n) = undef
        end if
        if (USE_OBSERR_RADAR_REF) then
          obs(iof)%err(n) = OBSERR_RADAR_REF
        end if
      case (id_radar_ref_zero_obs)
        obs(iof)%dat(n) = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
        if (USE_OBSERR_RADAR_REF) then
          obs(iof)%err(n) = OBSERR_RADAR_REF
        end if
      case (id_radar_vr_obs)
        if (USE_OBSERR_RADAR_VR) then
          obs(iof)%err(n) = OBSERR_RADAR_VR
        end if
      end select

      ! mark (elm, typ) combinations for which observations exist
      ctype_use(uid_obs(obs(iof)%elm(n)), obs(iof)%typ(n)) = .true.
    end do ! [ n = 1, obs(iof)%nobs ]
!$OMP END DO
  end do ! [ iof = 1, OBS_IN_NUM ]
!$OMP END PARALLEL

  ! do this outside of the above obs loop, so these (ctype) arrays can be in ascending order
  nctype = count(ctype_use)
  allocate (elm_ctype     (nctype))
  allocate (elm_u_ctype   (nctype))
  allocate (typ_ctype     (nctype))
  allocate (hori_loc_ctype(nctype))
  allocate (vert_loc_ctype(nctype))
  ictype = 0
  ctype_elmtyp(:,:) = 0
  do ityp = 1, nobtype
    do ielm_u = 1, nid_obs
      if (ctype_use(ielm_u, ityp)) then
        ictype = ictype + 1
        ctype_elmtyp(ielm_u, ityp) = ictype

        elm_ctype(ictype) = elem_uid(ielm_u)
        elm_u_ctype(ictype) = ielm_u
        typ_ctype(ictype) = ityp

        ! horizontal localization
        if (elm_ctype(ictype) == id_radar_ref_zero_obs) then
          hori_loc_ctype(ictype) = HORI_LOCAL_RADAR_OBSNOREF
        else if (elm_ctype(ictype) == id_radar_vr_obs) then
          hori_loc_ctype(ictype) = HORI_LOCAL_RADAR_VR
        else
          hori_loc_ctype(ictype) = HORI_LOCAL(ityp)
        end if
        ! vertical localization
        if (elm_ctype(ictype) == id_radar_vr_obs) then
          vert_loc_ctype(ictype) = VERT_LOCAL_RADAR_VR
        else
          vert_loc_ctype(ictype) = VERT_LOCAL(ityp)
        end if
      end if ! [ ctype_use(ielm_u, ityp) ]
    end do ! [ ielm_u = 1, nid_obs ]
  end do ! [ ityp = 1, nobtype ]

  call mpi_timer('set_letkf_obs:preprocess_data:', 2)

  ! Compute perturbation and departure
  !  -- gross error check
  !  -- QC based on background (radar reflectivity)
  !  -- process Himawari-8 data
  !-----------------------------------------------------------------------------

  allocate(tmpelm(obsda%nobs))

!$OMP PARALLEL PRIVATE(n,i,iof,iidx,mem_ref)
!$OMP DO
  do n = 1, obsda%nobs
    IF(obsda%qc(n) > 0) CYCLE

    iof = obsda%set(n)
    iidx = obsda%idx(n)

    tmpelm(n) = obs(iof)%elm(iidx)



!!!###### RADAR assimilation ######
    if (obs(iof)%elm(iidx) == id_radar_ref_obs .or. obs(iof)%elm(iidx) == id_radar_ref_zero_obs) then
      if (.not. USE_RADAR_REF) then
        obsda%qc(n) = iqc_otype
        cycle
      end if

      if (obs(iof)%dat(iidx) == undef) then
        obsda%qc(n) = iqc_obs_bad
        cycle
      end if

      !!! obsda%ensval: already converted to dBZ
      mem_ref = 0
      do i = 1, MEMBER
        if (obsda%ensval(i,n) > RADAR_REF_THRES_DBZ+1.0d-6) then
          mem_ref = mem_ref + 1
        end if
      end do

      ! Obs: Rain
      if (obs(iof)%dat(iidx) > RADAR_REF_THRES_DBZ+1.0d-6) then
        if (mem_ref < MIN_RADAR_REF_MEMBER_OBSRAIN) then
          obsda%qc(n) = iqc_ref_mem
          if ( LOG_LEVEL >= 3 .and. LOG_OUT ) then
            write (6,'(A)') '* Reflectivity does not fit assimilation criterion'
            write (6,'(A,F6.2,A,F6.2,A,I6,A,F7.3)') &
                  '*  (lon,lat)=(',obs(iof)%lon(iidx),',',obs(iof)%lat(iidx),'), mem_ref=', &
                  mem_ref,', ref_obs=', obs(iof)%dat(iidx)
          end if
         
          !if ( .not. RADAR_PQV ) cycle
          ! When RADAR_PQV=True, pseudo qv obs is assimilated even if mem_ref is
          ! too small
        end if

        if ( RADAR_ADDITIVE_Y18 ) then
          if ( mem_ref < RADAR_ADDITIVE_Y18_MINMEM ) then
            ! Ensemble mean of HX (B in O-B) is assumed to be [MIN_RADAR_REF_DBZ + LOW_REF_SHIFT],
            ! which is a clear sky reflectivity
            ! *Ensemble mean of obsda%epert should be zero
            obsda%ensval(1:MEMBER,n) = obsda%epert(1:MEMBER,n) + MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
            obsda%qc(n) = iqc_good
          endif
        endif

      else
      ! Obs: No rain
        if (mem_ref < MIN_RADAR_REF_MEMBER_OBSNORAIN) then
          obsda%qc(n) = iqc_ref_mem
          if ( LOG_LEVEL >= 3 .and. LOG_OUT ) then
            write (6,'(A)') '* Reflectivity does not fit assimilation criterion'
            write (6,'(A,F6.2,A,F6.2,A,I6,A,F7.3)') &
                  '*  (lon,lat)=(',obs(iof)%lon(iidx),',',obs(iof)%lat(iidx),'), mem_ref=', &
                  mem_ref,', ref_obs=', obs(iof)%dat(iidx)
          end if
          cycle
        end if
      end if
    end if

    if (obs(iof)%elm(iidx) == id_radar_vr_obs) then
      if (.not. USE_RADAR_VR) then
        obsda%qc(n) = iqc_otype
        cycle
      end if
    end if

!    if (obs(iof)%elm(iidx) == id_radar_prh_obs) then
!      if (.not. USE_RADAR_PSEUDO_RH) then
!        obsda%qc(n) = iqc_otype
!        cycle
!      end if
!    end if
!!!###### end RADAR assimilation ######

    obsda%val(n) = obsda%ensval(1,n)
    DO i=2,MEMBER
      obsda%val(n) = obsda%val(n) + obsda%ensval(i,n)
    END DO
    obsda%val(n) = obsda%val(n) / REAL(MEMBER,r_size)

    DO i=1,MEMBER
      obsda%ensval(i,n) = obsda%ensval(i,n) - obsda%val(n) ! Hdx
    END DO
    obsda%val(n) = obs(iof)%dat(iidx) - obsda%val(n) ! y-Hx
    if (DET_RUN) then
      obsda%ensval(mmdetobs,n) = obs(iof)%dat(iidx) - obsda%ensval(mmdetobs,n) ! y-Hx for deterministic run
    end if


    select case (obs(iof)%elm(iidx)) !gross error
    case (id_rain_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RAIN * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_radar_ref_obs,id_radar_ref_zero_obs)

      if( RADAR_PQV .and. obsda%val(n) > RADAR_PQV_OMB ) then
        ! pseudo qv
        obsda%val(n) = obsda%eqv(1,n)
        do i = 2 ,MEMBER
          obsda%val(n) = obsda%val(n) + obsda%eqv(i,n)
        enddo
        obsda%val(n) = obsda%val(n) / real(MEMBER, r_size)

        do i = 1, MEMBER
          obsda%ensval(i,n) = obsda%eqv(i,n) - obsda%val(n) ! Hdx
        enddO

        ! Tetens equation es(Pa)
        qvs = 611.2d0*exp(17.67d0*(obsda%tm(n)-t0c)/(obsda%tm(n) - t0c + 243.5d0))

        ! Saturtion mixing ratio
        qvs = 0.622d0*qvs / ( obsda%pm(n) - qvs )

        obsda%val(n) = qvs - obsda%val(n) ! y-Hx

        if (DET_RUN) then
          obsda%ensval(mmdetobs,n) = qvs - obsda%eqv(mmdetobs,n) ! y-Hx for deterministic run
        end if

        obsda%tm(n) = -1.0d0

      else

        IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_REF * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF

      end if
    case (id_radar_vr_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_VR * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_radar_prh_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_PRH * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case default
      IF(ABS(obsda%val(n)) > GROSS_ERROR * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    end select

    if ( LOG_LEVEL >= 3 .and. LOG_OUT ) then
      write (6, '(2I6,2F8.2,4F12.4,I3)') obs(iof)%elm(iidx), &
                                         obs(iof)%typ(iidx), &
                                         obs(iof)%lon(iidx), &
                                         obs(iof)%lat(iidx), &
                                         obs(iof)%lev(iidx), &
                                         obs(iof)%dat(iidx), &
                                         obs(iof)%err(iidx), &
                                         obsda%val(n), &
                                         obsda%qc(n)
    end if

  END DO ! [ n = 1, obsda%nobs ]
!$OMP END DO
!$OMP END PARALLEL

  call mpi_timer('set_letkf_obs:departure_cal_qc:', 2)

  ! Temporal observation localization !!!!!! not implemented yet !!!!!!
  !-----------------------------------------------------------------------------
!!
!!  DO n=1,nobs
!!    tmperr(n) = tmperr(n) * exp(0.25d0 * (tmpdif(n) / SIGMA_OBST)**2)
!!  END DO
!!
!! SELECT OBS IN THE NODE
!!

  ! Print departure statistics
  !-------------------------------------------------------------------------------

  if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
    write (6, *)
    write (6,'(A,I6,A)') 'OBSERVATIONAL DEPARTURE STATISTICS (IN THIS SUBDOMAIN #', myrank_d, '):'

    call monit_dep(obsda%nobs, tmpelm, obsda%val, obsda%qc, monit_nobs, bias, rmse)
    call monit_print(monit_nobs, bias, rmse)

    call mpi_timer('set_letkf_obs:departure_print:', 2)
  end if

  deallocate(tmpelm)

!-------------------------------------------------------------------------------
! "Bucket sort" of observations of each combined type (with different sorting meshes)
!-------------------------------------------------------------------------------

  allocate (obsgrd(nctype))

  ! Determine mesh size for bucket sort
  !-----------------------------------------------------------------------------

  do ictype = 1, nctype
    ityp = typ_ctype(ictype)

    if (OBS_SORT_GRID_SPACING(ityp) > 0) then
      target_grdspc = OBS_SORT_GRID_SPACING(ityp)
    else if (MAX_NOBS_PER_GRID(ityp) > 0) then
      target_grdspc = 0.1d0 * sqrt(real(MAX_NOBS_PER_GRID(ityp), r_size)) * OBS_MIN_SPACING(ityp) ! need to be tuned
    else
      target_grdspc = hori_loc_ctype(ictype) * dist_zero_fac / 6.0d0                ! need to be tuned
    end if
    obsgrd(ictype)%ngrd_i = min(ceiling(DX * real(nlon,r_size) / target_grdspc), nlon)
    obsgrd(ictype)%ngrd_j = min(ceiling(DY * real(nlat,r_size) / target_grdspc), nlat)
    obsgrd(ictype)%grdspc_i = DX * real(nlon,r_size) / real(obsgrd(ictype)%ngrd_i,r_size)
    obsgrd(ictype)%grdspc_j = DY * real(nlat,r_size) / real(obsgrd(ictype)%ngrd_j,r_size)
    obsgrd(ictype)%ngrdsch_i = ceiling(hori_loc_ctype(ictype) * dist_zero_fac / obsgrd(ictype)%grdspc_i)
    obsgrd(ictype)%ngrdsch_j = ceiling(hori_loc_ctype(ictype) * dist_zero_fac / obsgrd(ictype)%grdspc_j)
    obsgrd(ictype)%ngrdext_i = obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i * 2
    obsgrd(ictype)%ngrdext_j = obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j * 2

    allocate (obsgrd(ictype)%n (  obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j, 0:nprocs_d-1))
    allocate (obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j, 0:nprocs_d-1))
    allocate (obsgrd(ictype)%tot(0:nprocs_d-1))
    allocate (obsgrd(ictype)%n_ext (  obsgrd(ictype)%ngrdext_i, obsgrd(ictype)%ngrdext_j))
    allocate (obsgrd(ictype)%ac_ext(0:obsgrd(ictype)%ngrdext_i, obsgrd(ictype)%ngrdext_j))

    obsgrd(ictype)%n (:,:,:) = 0
    obsgrd(ictype)%ac(:,:,:) = 0
    obsgrd(ictype)%tot(:) = 0
    obsgrd(ictype)%n_ext (:,:) = 0
    obsgrd(ictype)%ac_ext(:,:) = 0
    obsgrd(ictype)%tot_ext = 0
    obsgrd(ictype)%tot_sub(:) = 0
    obsgrd(ictype)%tot_g(:) = 0

    allocate (obsgrd(ictype)%next(obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j))
  end do

  call mpi_timer('set_letkf_obs:bucket_sort_prepare:', 2)

  ! Print observation usage settings
  !-----------------------------------------------------------------------------

  if ( LOG_LEVEL >= 2 .and. LOG_OUT ) then
    write (6, *)
    write (6, '(A)') 'OBSERVATION USAGE SETTINGS (LIST ONLY EXISTING TYPE-VAR):'
    write (6, '(A)') '=================================================================================='
    write (6, '(A)') 'TYPE   VAR  USE HORI_LOC   VERT_LOC TIME_LOC MAX_NOBS MIN_SPAC SORT_MESH_X _MESH_Y'
    write (6, '(A)') '                    (km) (lnP or m)      (s)              (km)        (km)    (km)'
    write (6, '(A)') '----------------------------------------------------------------------------------'
    do ictype = 1, nctype
      ityp = typ_ctype(ictype)
      ielm = elm_ctype(ictype)
      ielm_u = elm_u_ctype(ictype)

      if (USE_OBS(ityp)) then
        if ((ielm == id_radar_ref_obs .or. ielm == id_radar_ref_zero_obs) .and. (.not. USE_RADAR_REF)) then
          use_obs_print = ' No'
        else if (ielm == id_radar_vr_obs .and. (.not. USE_RADAR_VR)) then
          use_obs_print = ' No'
        else if (ielm == id_radar_prh_obs .and. (.not. USE_RADAR_PSEUDO_RH)) then
          use_obs_print = ' No'
        else
          use_obs_print = 'Yes'
        end if
      else
        use_obs_print = ' No'
      end if

      select case (ityp)
      case (22) ! vertical localization in Z
        write (6, '(A6,1x,A3,2x,A3,F9.2,F7.2,A4,F9.2,I9,F9.2,F12.2,F8.2)') obtypelist(ityp), obelmlist(ielm_u), use_obs_print, hori_loc_ctype(ictype)/1000.0d0, &
                  vert_loc_ctype(ictype)/1000.0d0, '[km]', TIME_LOCAL(ityp)/1000.0d0, MAX_NOBS_PER_GRID(ityp), &
                  OBS_MIN_SPACING(ityp)/1000.0d0, obsgrd(ictype)%grdspc_i/1000.0d0, obsgrd(ictype)%grdspc_j/1000.0d0
      case default ! vertical localization in ln(p)
        write (6, '(A6,1x,A3,2x,A3,F9.2,F11.2,F9.2,I9,F9.2,F12.2,F8.2)') obtypelist(ityp), obelmlist(ielm_u), use_obs_print, hori_loc_ctype(ictype)/1000.0d0, &
                  vert_loc_ctype(ictype), TIME_LOCAL(ityp)/1000.0d0, MAX_NOBS_PER_GRID(ityp), &
                  OBS_MIN_SPACING(ityp)/1000.0d0, obsgrd(ictype)%grdspc_i/1000.0d0, obsgrd(ictype)%grdspc_j/1000.0d0
      end select
    end do
    write (6, '(A)') '=================================================================================='

    call mpi_timer('set_letkf_obs:obs_setting_print:', 2)
  end if

  ! First scan: count the observation numbers in each mesh (in each subdomian)
  !-----------------------------------------------------------------------------

  do n = 1, obsda%nobs
    iof = obsda%set(n)
    iidx = obsda%idx(n)
    ictype = ctype_elmtyp(uid_obs(obs(iof)%elm(iidx)), obs(iof)%typ(iidx))

    if (obsda%qc(n) == iqc_good) then
      call ij_obsgrd(ictype, obs(iof)%ri(iidx), obs(iof)%rj(iidx), i, j)
      if (i < 1) i = 1                                          ! Assume the process assignment was correct,
      if (i > obsgrd(ictype)%ngrd_i) i = obsgrd(ictype)%ngrd_i  ! so this correction is only to remedy the round-off problem.
      if (j < 1) j = 1                                          !
      if (j > obsgrd(ictype)%ngrd_j) j = obsgrd(ictype)%ngrd_j  !

      obsgrd(ictype)%n(i,j,myrank_d) = obsgrd(ictype)%n(i,j,myrank_d) + 1
      obsgrd(ictype)%tot_sub(i_after_qc) = obsgrd(ictype)%tot_sub(i_after_qc) + 1 ! only used for diagnostic print (obs number after qc)
    end if

    obsgrd(ictype)%tot_sub(i_before_qc) = obsgrd(ictype)%tot_sub(i_before_qc) + 1 ! only used for diagnostic print (obs number before qc)
  end do

  ! Compute the accumulated numbers in each mesh
  !-----------------------------------------------------------------------------

  do ictype = 1, nctype
    if (ictype > 1) then
      obsgrd(ictype)%ac(0,1,myrank_d) = obsgrd(ictype-1)%ac(obsgrd(ictype-1)%ngrd_i,obsgrd(ictype-1)%ngrd_j,myrank_d)
    end if
    do j = 1, obsgrd(ictype)%ngrd_j
      if (j > 1) then
        obsgrd(ictype)%ac(0,j,myrank_d) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,j-1,myrank_d)
      end if
      do i = 1, obsgrd(ictype)%ngrd_i
        obsgrd(ictype)%ac(i,j,myrank_d) = obsgrd(ictype)%ac(i-1,j,myrank_d) + obsgrd(ictype)%n(i,j,myrank_d)
!        obsgrd(ictype)%next(i,j) = obsgrd(ictype)%ac(i-1,j,myrank_d)
      end do
    end do
    obsgrd(ictype)%next(1:obsgrd(ictype)%ngrd_i,:) = obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i-1,:,myrank_d)
  end do

  call mpi_timer('set_letkf_obs:bucket_sort_first_scan:', 2)

  ! Second scan: save the indices of bucket-sorted observations in obsda%keys(:)
  !-----------------------------------------------------------------------------

  do n = 1, obsda%nobs
    if (obsda%qc(n) == iqc_good) then
      iof = obsda%set(n)
      iidx = obsda%idx(n)
      ictype = ctype_elmtyp(uid_obs(obs(iof)%elm(iidx)), obs(iof)%typ(iidx))

      call ij_obsgrd(ictype, obs(iof)%ri(iidx), obs(iof)%rj(iidx), i, j)
      if (i < 1) i = 1                                         ! Assume the process assignment was correct,
      if (i > obsgrd(ictype)%ngrd_i) i = obsgrd(ictype)%ngrd_i ! so this correction is only to remedy the round-off problem.
      if (j < 1) j = 1                                         !
      if (j > obsgrd(ictype)%ngrd_j) j = obsgrd(ictype)%ngrd_j !

      obsgrd(ictype)%next(i,j) = obsgrd(ictype)%next(i,j) + 1
      obsda%key(obsgrd(ictype)%next(i,j)) = n
    end if
  end do


  call mpi_timer('set_letkf_obs:bucket_sort_second_scan:', 2, barrier=MPI_COMM_d)

  ! ALLREDUCE observation number information from subdomains, and compute total numbers
  !-----------------------------------------------------------------------------

  nobs_sub(:) = 0
  nobs_g(:) = 0
  do ictype = 1, nctype
    if (nprocs_d > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsgrd(ictype)%n, obsgrd(ictype)%ngrd_i*obsgrd(ictype)%ngrd_j*nprocs_d, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i,:,:), (obsgrd(ictype)%ngrd_i+1)*obsgrd(ictype)%ngrd_j*nprocs_d, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
    end if
    call MPI_ALLREDUCE(obsgrd(ictype)%tot_sub, obsgrd(ictype)%tot_g, n_qc_steps, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)

    if (ictype == 1) then
      obsgrd(ictype)%tot(:) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,obsgrd(ictype)%ngrd_j,:)
    else
      obsgrd(ictype)%tot(:) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,obsgrd(ictype)%ngrd_j,:) &
                            - obsgrd(ictype-1)%ac(obsgrd(ictype-1)%ngrd_i,obsgrd(ictype-1)%ngrd_j,:)
    end if
#ifdef LETKF_DEBUG
    if (obsgrd(ictype)%tot(myrank_d) /= obsgrd(ictype)%tot_sub(i_after_qc)) then
      write (6, '(A)') '[Error] Observation counts are inconsistent !!!'
      stop 99
    end if
#endif

    nobs_sub(:) = nobs_sub(:) + obsgrd(ictype)%tot_sub(:)
    nobs_g(:) = nobs_g(:) + obsgrd(ictype)%tot_g(:)

    deallocate (obsgrd(ictype)%next)
  end do

  call mpi_timer('set_letkf_obs:bucket_sort_info_allreduce:', 2)

#ifdef LETKF_DEBUG
  if (nctype > 0) then
    if (obsgrd(nctype)%ac(obsgrd(nctype)%ngrd_i,obsgrd(nctype)%ngrd_j,myrank_d) /= nobs_sub(i_after_qc)) then
      write (6, '(A)') '[Error] Observation counts are inconsistent !!!'
      stop 99
    end if
    if (sum(obsgrd(nctype)%ac(obsgrd(nctype)%ngrd_i,obsgrd(nctype)%ngrd_j,:)) /= nobs_g(i_after_qc)) then
      write (6, '(A)') '[Error] Observation counts are inconsistent !!!'
      stop 99
    end if
  end if
#endif
  nobstotalg = nobs_g(i_after_qc) ! total obs number in the global domain (all types)

  ! Print observation counts for each types
  !-----------------------------------------------------------------------------

  if ( LOG_LEVEL >= 2 .and. LOG_OUT ) then
    write (nstr, '(I4)') nid_obs
    write (6, *)
    write (6, '(A)') 'OBSERVATION COUNTS BEFORE QC (GLOABL):'
    write (6, '(A7,'//nstr//"('========'),A)") '=======', '=========='
    write (6, '(A6,1x,'//nstr//'A8,A10)') 'TYPE  ', obelmlist(:), '     TOTAL'
    write (6, '(A7,'//nstr//"('--------'),A)") '-------', '----------'
    nobs_elms_sum(:) = 0
    do ityp = 1, nobtype
      nobs_elms(:) = 0
      do ielm_u = 1, nid_obs
        if (ctype_elmtyp(ielm_u,ityp) > 0) then
          nobs_elms(ielm_u) = obsgrd(ctype_elmtyp(ielm_u,ityp))%tot_g(i_before_qc)
        end if
      end do
      nobs_elms_sum = nobs_elms_sum + nobs_elms
      write (6, '(A6,1x,'//nstr//'I8,I10)') obtypelist(ityp), nobs_elms(:), sum(nobs_elms(:))
    end do
    write (6, '(A7,'//nstr//"('--------'),A)") '-------', '----------'
    write (6, '(A6,1x,'//nstr//'I8,I10)') 'TOTAL ', nobs_elms_sum(:), nobs_g(i_before_qc)
    write (6, '(A7,'//nstr//"('========'),A)") '=======', '=========='

    write (6, *)
    write (6, '(A)') 'OBSERVATION COUNTS AFTER QC (GLOABL):'
    write (6, '(A7,'//nstr//"('========'),A)") '=======', '=========='
    write (6, '(A6,1x,'//nstr//'A8,A10)') 'TYPE  ', obelmlist(:), '     TOTAL'
    write (6, '(A7,'//nstr//"('--------'),A)") '-------', '----------'
    nobs_elms_sum(:) = 0
    do ityp = 1, nobtype
      nobs_elms(:) = 0
      do ielm_u = 1, nid_obs
        if (ctype_elmtyp(ielm_u,ityp) > 0) then
          nobs_elms(ielm_u) = obsgrd(ctype_elmtyp(ielm_u,ityp))%tot_g(i_after_qc)
        end if
      end do
      nobs_elms_sum = nobs_elms_sum + nobs_elms
      write (6, '(A6,1x,'//nstr//'I8,I10)') obtypelist(ityp), nobs_elms(:), sum(nobs_elms(:))
    end do
    write (6, '(A7,'//nstr//"('--------'),A)") '-------', '----------'
    write (6, '(A6,1x,'//nstr//'I8,I10)') 'TOTAL ', nobs_elms_sum(:), nobs_g(i_after_qc)
    write (6, '(A7,'//nstr//"('========'),A)") '=======', '=========='

    call mpi_timer('set_letkf_obs:obs_count_print_types:', 2)
  end if

  ! Calculate observation numbers in the extended (localization) subdomain,
  ! in preparation for communicating obsetvations in the extended subdomain
  !-----------------------------------------------------------------------------

  call rank_1d_2d(myrank_d, myp_i, myp_j)

  do ictype = 1, nctype
    imin1 = myp_i*obsgrd(ictype)%ngrd_i+1 - obsgrd(ictype)%ngrdsch_i
    imax1 = (myp_i+1)*obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
    jmin1 = myp_j*obsgrd(ictype)%ngrd_j+1 - obsgrd(ictype)%ngrdsch_j
    jmax1 = (myp_j+1)*obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

    do ip = 0, nprocs_d-1
      call rank_1d_2d(ip, ip_i, ip_j)
      imin2 = max(1, imin1 - ip_i*obsgrd(ictype)%ngrd_i)
      imax2 = min(obsgrd(ictype)%ngrd_i, imax1 - ip_i*obsgrd(ictype)%ngrd_i)
      jmin2 = max(1, jmin1 - ip_j*obsgrd(ictype)%ngrd_j)
      jmax2 = min(obsgrd(ictype)%ngrd_j, jmax1 - ip_j*obsgrd(ictype)%ngrd_j)
      if (imin2 > imax2 .or. jmin2 > jmax2) cycle

      ishift = (ip_i - myp_i) * obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
      jshift = (ip_j - myp_j) * obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j
      obsgrd(ictype)%n_ext(imin2+ishift:imax2+ishift, jmin2+jshift:jmax2+jshift) = obsgrd(ictype)%n(imin2:imax2, jmin2:jmax2, ip)
    end do

    if (ictype > 1) then
      obsgrd(ictype)%ac_ext(0,1) = obsgrd(ictype-1)%ac_ext(obsgrd(ictype-1)%ngrdext_i,obsgrd(ictype-1)%ngrdext_j)
    end if
    do j = 1, obsgrd(ictype)%ngrdext_j
      if (j > 1) then
        obsgrd(ictype)%ac_ext(0,j) = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,j-1)
      end if
      do i = 1, obsgrd(ictype)%ngrdext_i
        obsgrd(ictype)%ac_ext(i,j) = obsgrd(ictype)%ac_ext(i-1,j) + obsgrd(ictype)%n_ext(i,j)
      end do
    end do

    if (ictype == 1) then
      obsgrd(ictype)%tot_ext = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,obsgrd(ictype)%ngrdext_j)
    else
      obsgrd(ictype)%tot_ext = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,obsgrd(ictype)%ngrdext_j) &
                             - obsgrd(ictype-1)%ac_ext(obsgrd(ictype-1)%ngrdext_i,obsgrd(ictype-1)%ngrdext_j)
    end if
  end do ! [ ictype = 1, nctype ]

  if (nctype > 0) then
    nobstotal = obsgrd(nctype)%ac_ext(obsgrd(nctype)%ngrdext_i,obsgrd(nctype)%ngrdext_j) ! total obs number in the extended subdomain (all types)

!    maxnobs_sub_per_ctype = maxval(obsgrd(1)%tot(:))
    maxnobs_per_ctype = obsgrd(1)%tot_ext
    do ictype = 2, nctype
!      maxnobs_sub_per_ctype = max(maxnobs_sub_per_ctype, maxval(obsgrd(ictype)%tot(:)))
      maxnobs_per_ctype = max(maxnobs_per_ctype, obsgrd(ictype)%tot_ext)
    end do
  else
    nobstotal = 0
!    maxnobs_sub_per_ctype = 0
    maxnobs_per_ctype = 0
  end if

  call mpi_timer('set_letkf_obs:extdomain_obs_count_cal:', 2)

  ! Construct sorted obsda_sort: 
  !-----------------------------------------------------------------------------

  ! 1) Copy the observation data in own subdomains to send buffer (obsbufs) with
  ! sorted order
  !-----------------------------------------------------------------------------
  if (nctype > 0) then
    cntr(:) = obsgrd(nctype)%ac(obsgrd(nctype)%ngrd_i,obsgrd(nctype)%ngrd_j,:)
    cnts = cntr(myrank_d+1)
    dspr(1) = 0
    do ip = 2, nprocs_d
      dspr(ip) = dspr(ip-1) + cntr(ip-1)
    end do

    nensobs_mod = mod(nensobs, nprocs_e)
    nensobs_div = (nensobs - nensobs_mod) / nprocs_e
    if (myrank_e < nensobs_mod) then
      im_obs_1 = (nensobs_div+1) * myrank_e + 1
      im_obs_2 = (nensobs_div+1) * (myrank_e+1)
      nensobs_part = nensobs_div + 1
    else
      im_obs_1 = (nensobs_div+1) * nensobs_mod + nensobs_div * (myrank_e-nensobs_mod) + 1
      im_obs_2 = (nensobs_div+1) * nensobs_mod + nensobs_div * (myrank_e-nensobs_mod+1)
      nensobs_part = nensobs_div
    end if

    obsbufs%nobs = nobs_sub(i_after_qc)
    call obs_da_value_allocate(obsbufs, nensobs_part)

    do n = 1, nobs_sub(i_after_qc)
      obsbufs%set(n) = obsda%set(obsda%key(n))
      obsbufs%idx(n) = obsda%idx(obsda%key(n))
      obsbufs%val(n) = obsda%val(obsda%key(n))
      if (nensobs_part > 0) then
        obsbufs%ensval(1:nensobs_part,n) = obsda%ensval(im_obs_1:im_obs_2,obsda%key(n))
      end if
      obsbufs%qc(n) = obsda%qc(obsda%key(n))
    end do ! [ n = 1, nobs_sub(i_after_qc) ]

    call mpi_timer('set_letkf_obs:copy_bufs:', 2, barrier=MPI_COMM_d)
  end if ! [ nctype > 0 ]

  call obs_da_value_deallocate(obsda)

  ! 2) Communicate to get global observations;
  !    for variables with an ensemble dimension (ensval),
  !    only obtain data from partial members (nensobs_part) to save memory usage
  !-----------------------------------------------------------------------------
  if (nctype > 0) then
    obsbufr%nobs = nobs_g(i_after_qc)
    call obs_da_value_allocate(obsbufr, nensobs_part)

    call MPI_ALLGATHERV(obsbufs%set, cnts, MPI_INTEGER, obsbufr%set, cntr, dspr, MPI_INTEGER, MPI_COMM_d, ierr)
    call MPI_ALLGATHERV(obsbufs%idx, cnts, MPI_INTEGER, obsbufr%idx, cntr, dspr, MPI_INTEGER, MPI_COMM_d, ierr)
    call MPI_ALLGATHERV(obsbufs%val, cnts, MPI_r_size, obsbufr%val, cntr, dspr, MPI_r_size, MPI_COMM_d, ierr)
    if (nensobs_part > 0) then
      call MPI_ALLGATHERV(obsbufs%ensval, cnts*nensobs_part, MPI_r_size, obsbufr%ensval, cntr*nensobs_part, dspr*nensobs_part, MPI_r_size, MPI_COMM_d, ierr)
    end if
    call MPI_ALLGATHERV(obsbufs%qc, cnts, MPI_INTEGER, obsbufr%qc, cntr, dspr, MPI_INTEGER, MPI_COMM_d, ierr)

    call obs_da_value_deallocate(obsbufs)

    call mpi_timer('set_letkf_obs:mpi_allgatherv:', 2)
  end if ! [ nctype > 0 ]

  ! 3) Copy observation data within the extended (localization) subdomains
  !    from receive buffer (obsbufr) to obsda_sort; rearrange with sorted order
  !-----------------------------------------------------------------------------
  obsda_sort%nobs = nobstotal
  call obs_da_value_allocate(obsda_sort, nensobs)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ip,ip_i,ip_j,ictype,imin1,imax1,jmin1,jmax1,imin2,imax2,jmin2,jmax2,ishift,jshift,j,ns_ext,ne_ext,ns_bufr,ne_bufr)
  do ip = 0, nprocs_d-1
    call rank_1d_2d(ip, ip_i, ip_j)

    do ictype = 1, nctype
      imin1 = myp_i*obsgrd(ictype)%ngrd_i+1 - obsgrd(ictype)%ngrdsch_i
      imax1 = (myp_i+1)*obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
      jmin1 = myp_j*obsgrd(ictype)%ngrd_j+1 - obsgrd(ictype)%ngrdsch_j
      jmax1 = (myp_j+1)*obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

      imin2 = max(1, imin1 - ip_i*obsgrd(ictype)%ngrd_i)
      imax2 = min(obsgrd(ictype)%ngrd_i, imax1 - ip_i*obsgrd(ictype)%ngrd_i)
      jmin2 = max(1, jmin1 - ip_j*obsgrd(ictype)%ngrd_j)
      jmax2 = min(obsgrd(ictype)%ngrd_j, jmax1 - ip_j*obsgrd(ictype)%ngrd_j)
      if (imin2 > imax2 .or. jmin2 > jmax2) cycle

      ishift = (ip_i - myp_i) * obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
      jshift = (ip_j - myp_j) * obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

      do j = jmin2, jmax2
        ns_ext = obsgrd(ictype)%ac_ext(imin2+ishift-1,j+jshift) + 1
        ne_ext = obsgrd(ictype)%ac_ext(imax2+ishift  ,j+jshift)
        if (ns_ext > ne_ext) cycle

        ns_bufr = dspr(ip+1) + obsgrd(ictype)%ac(imin2-1,j,ip) + 1
        ne_bufr = dspr(ip+1) + obsgrd(ictype)%ac(imax2  ,j,ip)

#ifdef LETKF_DEBUG
        if (ne_ext - ns_ext /= ne_bufr - ns_bufr) then
          write (6, '(A)') '[Error] observation grid indices have errors !!!'
          write (6, *) ictype, ip, j, imin1, imax1, jmin1, jmax1, imin2, imax2, jmin2, jmax2, ishift, jshift, &
                       obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j, obsgrd(ictype)%ngrdsch_i, obsgrd(ictype)%ngrdsch_j, ns_ext, ne_ext, ns_bufr, ne_bufr
          stop 99
        end if
#endif

        obsda_sort%set(ns_ext:ne_ext) = obsbufr%set(ns_bufr:ne_bufr)
        obsda_sort%idx(ns_ext:ne_ext) = obsbufr%idx(ns_bufr:ne_bufr)
        obsda_sort%val(ns_ext:ne_ext) = obsbufr%val(ns_bufr:ne_bufr)
        if (nensobs_part > 0) then
          obsda_sort%ensval(im_obs_1:im_obs_2,ns_ext:ne_ext) = obsbufr%ensval(1:nensobs_part,ns_bufr:ne_bufr)
        end if
        obsda_sort%qc(ns_ext:ne_ext) = obsbufr%qc(ns_bufr:ne_bufr)
      end do
    end do ! [ ictype = 1, nctype ]
  end do ! [ ip = 0, nprocs_d-1 ]
!$OMP END PARALLEL DO

  ! Save the keys of observations within the subdomain (excluding the localization buffer area)
  obsda_sort%nobs_in_key = 0
  do ictype = 1, nctype
    imin1 = obsgrd(ictype)%ngrdsch_i + 1
    imax1 = obsgrd(ictype)%ngrdsch_i + obsgrd(ictype)%ngrd_i
    jmin1 = obsgrd(ictype)%ngrdsch_j + 1
    jmax1 = obsgrd(ictype)%ngrdsch_j + obsgrd(ictype)%ngrd_j
    call obs_choose_ext(ictype, imin1, imax1, jmin1, jmax1, obsda_sort%nobs_in_key, obsda_sort%key)

    deallocate (obsgrd(ictype)%n)
    deallocate (obsgrd(ictype)%ac)
    deallocate (obsgrd(ictype)%n_ext)
  end do ! [ ictype = 1, nctype ]

  if (nctype > 0) then
    call obs_da_value_deallocate(obsbufr)

    call mpi_timer('set_letkf_obs:extdomain_copy_bufr:', 2, barrier=MPI_COMM_e)
  end if ! [ nctype > 0 ]

  ! 4) For variables with an ensemble dimension (ensval),
  !    ALLREDUCE among the ensemble dimension to obtain data of all members
  !-----------------------------------------------------------------------------
  if (nprocs_e > 1) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsda_sort%ensval, nensobs*nobstotal, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)

    call mpi_timer('set_letkf_obs:extdomain_allreduce:', 2)
  end if

  if ( LOG_LEVEL >= 3 .and. LOG_OUT ) then
    do n = 1, nobstotal
      write (6, '(I9,2I6,2F8.2,3F12.4,2F10.4,I6,F12.4,I3)') n, &
            obs(obsda_sort%set(n))%elm(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%typ(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%lon(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%lat(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%lev(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%dat(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%err(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%ri(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%rj(obsda_sort%idx(n)), &
            obs(obsda_sort%set(n))%rank(obsda_sort%idx(n)), &
            obsda_sort%val(n), &
            obsda_sort%qc(n)
    end do
  end if

  ! Print observation counts
  !-----------------------------------------------------------------------------

  if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
    write (6, *)
    write (6, '(A,I6,A)') 'OBSERVATION COUNTS (GLOABL AND IN THIS SUBDOMAIN #', myrank_d, '):'
    write (6, '(A)') '====================================================================='
    write (6, '(A)') 'TYPE   VAR      GLOBAL     GLOBAL  SUBDOMAIN  SUBDOMAIN EXT_SUBDOMAIN'
    write (6, '(A)') '             before QC   after QC  before QC   after QC      after QC'
    write (6, '(A)') '---------------------------------------------------------------------'
    do ictype = 1, nctype
      ityp = typ_ctype(ictype)
      ielm_u = elm_u_ctype(ictype)
      write (6, '(A6,1x,A3,1x,4I11,I14)') obtypelist(ityp), obelmlist(ielm_u), obsgrd(ictype)%tot_g(i_before_qc), obsgrd(ictype)%tot_g(i_after_qc), &
                obsgrd(ictype)%tot_sub(i_before_qc), obsgrd(ictype)%tot_sub(i_after_qc), obsgrd(ictype)%tot_ext
    end do
    write (6, '(A)') '---------------------------------------------------------------------'
      write (6, '(A6,5x,4I11,I14)') 'TOTAL ', nobs_g(i_before_qc), nobs_g(i_after_qc), nobs_sub(i_before_qc), nobs_sub(i_after_qc), nobstotal
    write (6, '(A)') '====================================================================='

    call mpi_timer('set_letkf_obs:obs_count_print_qc_extdomain:', 2)
  end if

  if ( OBSNUM_OUT_NC .and. myrank_e == mmean_rank_e .and. myrank_d == 0 ) then
    call write_obsnum_nc( trim( OBSNUM_OUT_NC_BASENAME )//'.nc', nctype, typ_ctype, elm_u_ctype, &
                         obsgrd(1:nctype)%tot_g(i_before_qc), obsgrd(1:nctype)%tot_g(i_after_qc), &
                         nobs_g(i_before_qc), nobs_g(i_after_qc)  )
  endif
  call mpi_timer('set_letkf_obs:write_obs_count_nc:', 2)

  RETURN
END SUBROUTINE set_letkf_obs

!-----------------------------------------------------------------------
! Convert grid (i,j) values to obsgrid (ogi, ogj) sorting mesh
!-----------------------------------------------------------------------
subroutine ij_obsgrd(ctype, ri, rj, ogi, ogj)
  use scale_atmos_grid_cartesC_index, only: &
    IHALO,JHALO
!  use scale_prc, only: &
!    PRC_myrank
  implicit none
  integer, intent(in) :: ctype
  real(r_size), intent(in) :: ri, rj
  integer, intent(out) :: ogi, ogj
  real(r_size) :: ril, rjl

  call rij_g2l(myrank_d, ri, rj, ril, rjl)
  ogi = ceiling((ril - real(IHALO,r_size) - 0.5) * real(obsgrd(ctype)%ngrd_i,r_size) / real(nlon,r_size))
  ogj = ceiling((rjl - real(JHALO,r_size) - 0.5) * real(obsgrd(ctype)%ngrd_j,r_size) / real(nlat,r_size))

  return
end subroutine ij_obsgrd

!-----------------------------------------------------------------------
! Convert grid (i,j) values to obsgrid (ogi, ogj) sorting mesh
! in the extended subdomain
!-----------------------------------------------------------------------
subroutine ij_obsgrd_ext(ctype, ri, rj, ogi, ogj)
  use scale_atmos_grid_cartesC_index, only: &
    IHALO,JHALO
!  use scale_prc, only: &
!    PRC_myrank
  implicit none
  integer, intent(in) :: ctype
  real(r_size), intent(in) :: ri, rj
  integer, intent(out) :: ogi, ogj
  real(r_size) :: ril, rjl

  call rij_g2l(myrank_d, ri, rj, ril, rjl)
  ogi = ceiling((ril - real(IHALO,r_size) - 0.5) * real(obsgrd(ctype)%ngrd_i,r_size) / real(nlon,r_size)) &
      + obsgrd(ctype)%ngrdsch_i
  ogj = ceiling((rjl - real(JHALO,r_size) - 0.5) * real(obsgrd(ctype)%ngrd_j,r_size) / real(nlat,r_size)) &
      + obsgrd(ctype)%ngrdsch_j

  return
end subroutine ij_obsgrd_ext

!-----------------------------------------------------------------------
! Choose observations in a rectangle using the bucket sort results
!-----------------------------------------------------------------------
subroutine obs_choose(ctype, proc, imin, imax, jmin, jmax, nn, nobs_use)
  implicit none
  integer, intent(in) :: ctype
  integer, intent(in) :: proc
  integer, intent(in) :: imin, imax, jmin, jmax
  integer, intent(inout) :: nn
  integer, intent(inout), optional :: nobs_use(:)
  integer :: n, j

  if (imin > imax .or. jmin > jmax) return
  if (obsgrd(ctype)%tot(proc) == 0) return

  do j = jmin, jmax
    if (present(nobs_use)) then
      do n = obsgrd(ctype)%ac(imin-1,j,proc)+1, obsgrd(ctype)%ac(imax,j,proc)
        nn = nn + 1
        nobs_use(nn) = n
      end do
    else
      nn = nn + obsgrd(ctype)%ac(imax,j,proc) - obsgrd(ctype)%ac(imin-1,j,proc)
    end if
  end do

  return
end subroutine obs_choose

!-----------------------------------------------------------------------
! Choose observations in a rectangle using the bucket sort results
! in the extended subdomain
!-----------------------------------------------------------------------
!OCL SERIAL
subroutine obs_choose_ext(ctype, imin, imax, jmin, jmax, nn, nobs_use)
  implicit none
  integer, intent(in) :: ctype
  integer, intent(in) :: imin, imax, jmin, jmax
  integer, intent(inout) :: nn
  integer, intent(inout), optional :: nobs_use(:)
  integer :: n, j

  if (imin > imax .or. jmin > jmax) return
  if (obsgrd(ctype)%tot_ext == 0) return

  do j = jmin, jmax
    if (present(nobs_use)) then
      do n = obsgrd(ctype)%ac_ext(imin-1,j)+1, obsgrd(ctype)%ac_ext(imax,j)
        nn = nn + 1
        nobs_use(nn) = n
      end do
    else
      nn = nn + obsgrd(ctype)%ac_ext(imax,j) - obsgrd(ctype)%ac_ext(imin-1,j)
    end if
  end do

  return
end subroutine obs_choose_ext


!!-----------------------------------------------------------------------
!! Monitor observation diagnostics after the letkf update
!!  Moved from main program code to this subroutine, 2013/12/26 Guo-Yuan Lien
!!-----------------------------------------------------------------------
!SUBROUTINE monit_obs
!  IMPLICIT NONE
!  REAL(r_size),ALLOCATABLE :: ohx(:)
!  REAL(r_size),ALLOCATABLE :: dep(:)
!  INTEGER,ALLOCATABLE :: oqc(:)
!  INTEGER :: l,im
!  CHARACTER(10) :: ombfile='omb.dat'
!  CHARACTER(10) :: omafile='oma.dat'
!  CHARACTER(14) :: obsguesfile='obsguesNNN.dat'
!  CHARACTER(14) :: obsanalfile='obsanalNNN.dat'
!
!  IF(omb_output .AND. myrank == 0) THEN
!    ALLOCATE(ohx(nobs),oqc(nobs),dep(nobs))
!    CALL monit_output('gues',0,ohx,oqc)
!    dep = obsdat - ohx
!    WRITE(6,'(A)') 'OBSERVATIONAL DEPARTURE STATISTICS [gues mean]:'
!    CALL monit_dep(nobs,obselm,dep,oqc)
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',ombfile
!    CALL write_obs2(ombfile,nobs,obselm,obslon,obslat,obslev, &
!                    obsdat,obserr,obstyp,obsdif,dep,obsqc,0)
!  END IF
!  IF(oma_output .AND. myrank == 0) THEN
!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs),dep(nobs))
!    CALL monit_output('anal',0,ohx,oqc)
!    dep = obsdat - ohx
!    WRITE(6,'(A)') 'OBSERVATIONAL DEPARTURE STATISTICS [anal mean]:'
!    CALL monit_dep(nobs,obselm,dep,oqc)
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',omafile
!    CALL write_obs2(omafile,nobs,obselm,obslon,obslat,obslev, &
!                    obsdat,obserr,obstyp,obsdif,dep,obsqc,0)
!  END IF

!  IF(obsgues_output) THEN
!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs))
!    IF(myrank == 0) THEN
!      ohx = obsdat - obsdep
!      oqc = 1
!      WRITE(obsguesfile(8:10),'(A3)') '_me'
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsguesfile
!      CALL write_obs2(obsguesfile,nobs,obselm,obslon,obslat,obslev, &
!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!    END IF
!    l=0
!    DO
!      im = myrank+1 + nprocs * l
!      IF(im > MEMBER) EXIT
!      ohx(:) = obsdat(:) - obsdep(:) + obshdxf(:,im)
!      WRITE(obsguesfile(8:10),'(I3.3)') im
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsguesfile
!      CALL write_obs2(obsguesfile,nobs,obselm,obslon,obslat,obslev, &
!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!      l = l+1
!    END DO
!  END IF

!! This is not an accurate estimate of obsanal.
!! To obtain an accurate estimate of obsanal, use code in [letkf_tools:das_letkf]
!!------
!!  IF(obsanal_output) THEN
!!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs))
!!    CALL monit_output('anal',0,ohx,oqc)
!!    WRITE(obsanalfile(8:10),'(A3)') '_me'
!!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
!!    CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
!!                    obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!!    l=0
!!    DO
!!      im = myrank+1 + nprocs * l
!!      IF(im > MEMBER) EXIT
!!      CALL monit_output('anal',im,ohx,oqc)
!!      WRITE(obsanalfile(8:10),'(I3.3)') im
!!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
!!      CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
!!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!!      l = l+1
!!    END DO
!!  END IF

!  IF(ALLOCATED(ohx)) DEALLOCATE(ohx,oqc)
!  IF(ALLOCATED(dep)) DEALLOCATE(dep)

!  RETURN
!END SUBROUTINE monit_obs
!-----------------------------------------------------------------------
! Monitor h(xb) or h(xa) from a LETKF output file
! Adopted from 'monit_mean' subroutine, 2013/12/24 Guo-Yuan Lien
!-----------------------------------------------------------------------
!  file: 'gues' or 'anal'
!  im:   member # (integer); 0 for ensmean (always called from myrank == 0)
!  ohx:  h(xb) or h(xa)
!  oqc:  quality of h(xb) or h(xa)
!-----------------------------------------------------------------------
!SUBROUTINE monit_output(file,im,ohx,oqc)
!  IMPLICIT NONE
!  CHARACTER(4),INTENT(IN) :: file
!  INTEGER,INTENT(IN) :: im
!  REAL(r_size),INTENT(OUT) :: ohx(nobs)
!  INTEGER,INTENT(OUT) :: oqc(nobs)
!  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3dx)
!  REAL(r_size) :: v2d(nlon,nlat,nv2dx)
!  REAL(r_size) :: v3dtmp(nlon,nlat,nlev,nv3d)
!  REAL(r_size) :: v2dtmp(nlon,nlat,nv2d)
!  REAL(r_size) :: elem
!  REAL(r_size) :: bias_u,bias_v,bias_t,bias_ps,bias_q,bias_rh,bias_rain
!  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_ps,rmse_q,rmse_rh,rmse_rain
!  REAL(r_size) :: hdxf,dep,ri,rj,rk
!  INTEGER :: n,iu,iv,it,iq,ips,irh,irain
!  CHARACTER(11) :: filename='filexxx.grd'
!  INTEGER :: iret
!  REAL(r_size) :: tmpps(nlon*nlat)
!  REAL(r_size) :: tmptv(nlon*nlat,nlev)
!  REAL(r_size) :: tmpp(nlon*nlat,nlev)



!  CHARACTER(4),INTENT(IN) :: file
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,MEMBER,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,MEMBER,nv2d)
!  REAL(RP) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP) :: v2dg(nlon,nlat,nv2d)
!  CHARACTER(9) :: filename='file.0000'
!  integer :: it,im,mstart,mend




!  IF(im == 0) THEN
!    WRITE(filename(1:7),'(A4,A3)') file,'_me'
!  ELSE
!    WRITE(filename(1:7),'(A4,I3.3)') file,im
!  END IF
!  CALL read_grd(filename,v3dtmp,v2dtmp,0) ! read ensemble mean into a temporary array
!  IF(im == 0) THEN
!    WRITE(filename(1:7),'(A4,I3.3)') 'gues',1 ! assume called from myrank == 0
!  ELSE
!    WRITE(filename(1:7),'(A4,I3.3)') 'gues',im
!  END IF
!  CALL read_grdx(filename,v3d,v2d) ! only the orography is used, P will be recalulated






!!-----------------------------------------------------------------------




!      WRITE(filename(1:4),'(A4)') file
!      WRITE(filename(6:9),'(I4.4)') im
!!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',myrank_d,'.nc'
!      call read_restart(filename,v3dg,v2dg)






!!  v3d(:,:,:,iv3d_u) = v3dtmp(:,:,:,iv3d_u)
!!  v3d(:,:,:,iv3d_v) = v3dtmp(:,:,:,iv3d_v)
!!  v3d(:,:,:,iv3d_t) = v3dtmp(:,:,:,iv3d_t)
!!  v3d(:,:,:,iv3d_q) = v3dtmp(:,:,:,iv3d_q)
!!  v3d(:,:,:,iv3d_qc) = v3dtmp(:,:,:,iv3d_qc)
!!  v2d(:,:,iv2d_ps) = v2dtmp(:,:,iv2d_ps)
!!  tmpps = reshape(v2d(:,:,iv2d_ps),(/nlon*nlat/))
!!  tmptv = reshape(v3d(:,:,:,iv3d_t) * (1.0d0 + fvirt * v3d(:,:,:,iv3d_q)),(/nlon*nlat,nlev/))
!!  call sigio_modprd(nlon*nlat,nlon*nlat,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
!!                    gfs_vcoord,iret,tmpps,tmptv,pm=tmpp)
!!  v3d(:,:,:,iv3d_p) = reshape(tmpp,(/nlon,nlat,nlev/))

!!  oqc = 1
!!  DO n=1,nobs
!!    CALL phys2ijk(v3d(:,:,:,iv3d_p),obselm(n),obslon(n),obslat(n),obslev(n),ri,rj,rk)
!!    !
!!    ! For monitoring, don't skip any observation below or above model vertical extent.
!!    ! Just put bad QC but still get estimate.
!!    !
!!    IF(CEILING(rk) > nlev) THEN
!!      rk = REAL(nlev,r_size)
!!      oqc(n) = 0
!!    END IF
!!    IF(CEILING(rk) < 2 .AND. NINT(obselm(n)) /= id_ps_obs) THEN
!!      IF(NINT(obselm(n)) > 9999) THEN
!!        rk = 0.0d0
!!      ELSE IF(NINT(obselm(n)) == id_u_obs .OR. NINT(obselm(n)) == id_v_obs) THEN
!!        rk = 1.00001d0
!!      ELSE
!!        rk = 1.00001d0
!!        oqc(n) = 0
!!      END IF
!!    END IF
!!    IF(NINT(obselm(n)) == id_ps_obs) THEN
!!      CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,rk)
!!      rk = obslev(n) - rk
!!    END IF
!!    IF(NINT(obselm(n)) == id_rain_obs) THEN ! No way to get the accumulated precipitation value
!!      ohx(n) = obsdat(n)
!!      oqc(n) = 0
!!    ELSE
!!      CALL Trans_XtoY(obselm(n),ri,rj,rk,v3d,v2d,ohx(n))
!!    END IF
!!  END DO

!  RETURN
!END SUBROUTINE monit_output
!-----------------------------------------------------------------------
! Read observation diagnostics for EFSO
!  Adapted from Y.Ohta's EFSO code for SPEEDY-LETKF   2013/07/17 D.Hotta
!  Modified, renamed from 'read_monit_obs' to 'set_efso_obs', 2013/12/26 Guo-Yuan Lien
!-----------------------------------------------------------------------
subroutine set_efso_obs
  use scale_atmos_grid_cartesC, only: &
    DX, &
    DY
  implicit none

  real(r_size), allocatable :: tmpdep(:)
  integer, allocatable :: tmpqc0(:,:)
  integer, allocatable :: tmpqc(:)
  integer :: nj(0:nlat-1)
  integer :: njs(1:nlat-1)
  integer :: l,im,i,j,n,nn,ierr
  character(14) :: obsguesfile='obsguesNNN.dat'

  integer :: m
  integer :: ityp, ielm, ielm_u, ictype
  real(r_size) :: target_grdspc
  integer :: iof, iidx

  integer :: nobs_sub(n_qc_steps), nobs_g(n_qc_steps)
  integer :: myp_i, myp_j
  integer :: imin1, imax1, jmin1, jmax1, imin2, imax2, jmin2, jmax2
  integer :: ip_i, ip_j, ip
  integer :: ishift, jshift

  integer :: cnts
  integer :: cntr(nprocs_d)
  integer :: dspr(nprocs_d)
  integer :: nensobs_div, nensobs_mod

  integer :: im_obs_1, im_obs_2, nensobs_part
  type(obs_da_value) :: obsbufs, obsbufr
  integer :: ns_ext, ne_ext, ns_bufr, ne_bufr

  logical :: ctype_use(nid_obs,nobtype)

  integer :: nobs0 ! number of total obs until myrank_d-1

  if ( LOG_OUT ) write(6,'(A)') 'Hello from set_efso_obs'

!  dist_zero = SIGMA_OBS * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zero_rain = SIGMA_OBS_RAIN * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov = SIGMA_OBSV * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov_rain = SIGMA_OBSV_RAIN * SQRT(10.0d0/3.0d0) * 2.0d0

  call get_nobs_efso_mpi( nobstotalg, nobstotal, nobs0 )
  if ( LOG_OUT ) write(6,'(I10,A)') nobstotalg,' TOTAL OBSERVATIONS INPUT'
!
! INITIALIZE GLOBAL VARIABLES
!
  if ( nobstotal > 0 ) then
    allocate( obsset(nobstotal) )
    allocate( obsidx(nobstotal) )
    allocate( obselm(nobstotal) )
    allocate( obstyp(nobstotal) )
    allocate( obslon(nobstotal) )
    allocate( obslat(nobstotal) )
    allocate( obslev(nobstotal) )
    allocate( obsdat(nobstotal) )
    allocate( obserr(nobstotal) )
    allocate( obsdif(nobstotal) )
    allocate( obsdep(nobstotal) )
    allocate( obsqc(nobstotal) )
    allocate( obshdxf(nobstotal,MEMBER) )
  !  allocate( tmpdep(nobstotal) )
  !  allocate( tmpqc0(nobstotal,MEMBER) )
  
!
! Reading background observation data and analysis ensemble observations
!
  
    call get_obsdep_efso_mpi( nobstotal, nobs0, obsset, obsidx, obselm, obstyp, obslon, obslat, obslev, &
                              obsdat, obserr, obsdif, obsdep, obsqc, obshdxf )
  endif
  
  obsda%nobs = nobstotal
  call obs_da_value_allocate( obsda, MEMBER )

  if ( nobstotal > 0 ) then
    obsda%set(:) = obsset
    obsda%idx(:) = obsidx
    obsda%val(:) = obsdep
    obsda%qc(:) = obsqc ! qc should always be good because QC has been applied in a previous cycle
                        ! obsda%qc array stores indices of each obs
  
    do m = 1, MEMBER
      obsda%ensval(m,:) = obshdxf(:,m)
    enddo
  endif

  ctype_use(:,:) = .false.
!#$OMP PARALLEL PRIVATE(iof,n) REDUCTION(.or.:ctype_use)
  do iof = 1, OBS_IN_NUM
!#$OMP DO
    do n = 1, obs(iof)%nobs
      select case (obs(iof)%elm(n))
      case (id_radar_ref_obs)
        if (obs(iof)%dat(n) >= 0.0d0 .and. obs(iof)%dat(n) < 1.0d10) then
          if (obs(iof)%dat(n) < MIN_RADAR_REF) then
            obs(iof)%elm(n) = id_radar_ref_zero_obs
            obs(iof)%dat(n) = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
          else
            obs(iof)%dat(n) = 10.0d0 * log10(obs(iof)%dat(n))
          end if
        else
          obs(iof)%dat(n) = undef
        end if
        if (USE_OBSERR_RADAR_REF) then
          obs(iof)%err(n) = OBSERR_RADAR_REF
        end if
      case (id_radar_ref_zero_obs)
        obs(iof)%dat(n) = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
        if (USE_OBSERR_RADAR_REF) then
          obs(iof)%err(n) = OBSERR_RADAR_REF
        end if
      case (id_radar_vr_obs)
        if (USE_OBSERR_RADAR_VR) then
          obs(iof)%err(n) = OBSERR_RADAR_VR
        end if
      end select

      ! mark (elm, typ) combinations for which observations exist
      ctype_use(uid_obs(obs(iof)%elm(n)), obs(iof)%typ(n)) = .true.
    end do ! [ n = 1, obs(iof)%nobs ]
!#$OMP END DO
  end do ! [ iof = 1, OBS_IN_NUM ]
!#$OMP END PARALLEL


  ! do this outside of the above obs loop, so these (ctype) arrays can be in ascending order
  nctype = count(ctype_use)
  allocate (elm_ctype     (nctype))
  allocate (elm_u_ctype   (nctype))
  allocate (typ_ctype     (nctype))
  allocate (hori_loc_ctype(nctype))
  allocate (vert_loc_ctype(nctype))
  ictype = 0
  ctype_elmtyp(:,:) = 0
  do ityp = 1, nobtype
    do ielm_u = 1, nid_obs
      if (ctype_use(ielm_u, ityp)) then
        ictype = ictype + 1
        ctype_elmtyp(ielm_u, ityp) = ictype

        elm_ctype(ictype) = elem_uid(ielm_u)
        elm_u_ctype(ictype) = ielm_u
        typ_ctype(ictype) = ityp

        ! horizontal localization
        if (elm_ctype(ictype) == id_radar_ref_zero_obs) then
          hori_loc_ctype(ictype) = HORI_LOCAL_RADAR_OBSNOREF
        else if (elm_ctype(ictype) == id_radar_vr_obs) then
          hori_loc_ctype(ictype) = HORI_LOCAL_RADAR_VR
        else
          hori_loc_ctype(ictype) = HORI_LOCAL(ityp)
        end if
        ! vertical localization
        if (elm_ctype(ictype) == id_radar_vr_obs) then
          vert_loc_ctype(ictype) = VERT_LOCAL_RADAR_VR
        else
          vert_loc_ctype(ictype) = VERT_LOCAL(ityp)
        end if
      end if ! [ ctype_use(ielm_u, ityp) ]
    end do ! [ ielm_u = 1, nid_obs ]
  end do ! [ ityp = 1, nobtype ]

!-------------------------------------------------------------------------------
! "Bucket sort" of observations of each combined type (with different sorting meshes)
!-------------------------------------------------------------------------------

  allocate (obsgrd(nctype))

  ! Determine mesh size for bucket sort
  !-----------------------------------------------------------------------------

  do ictype = 1, nctype
    ityp = typ_ctype(ictype)

    if (OBS_SORT_GRID_SPACING(ityp) > 0) then
      target_grdspc = OBS_SORT_GRID_SPACING(ityp)
    else if (MAX_NOBS_PER_GRID(ityp) > 0) then
      target_grdspc = 0.1d0 * sqrt(real(MAX_NOBS_PER_GRID(ityp), r_size)) * OBS_MIN_SPACING(ityp) ! need to be tuned
    else
      target_grdspc = hori_loc_ctype(ictype) * dist_zero_fac / 6.0d0                ! need to be tuned
    end if
    obsgrd(ictype)%ngrd_i = min(ceiling(DX * real(nlon,r_size) / target_grdspc), nlon)
    obsgrd(ictype)%ngrd_j = min(ceiling(DY * real(nlat,r_size) / target_grdspc), nlat)
    obsgrd(ictype)%grdspc_i = DX * real(nlon,r_size) / real(obsgrd(ictype)%ngrd_i,r_size)
    obsgrd(ictype)%grdspc_j = DY * real(nlat,r_size) / real(obsgrd(ictype)%ngrd_j,r_size)
    obsgrd(ictype)%ngrdsch_i = ceiling(hori_loc_ctype(ictype) * dist_zero_fac / obsgrd(ictype)%grdspc_i)
    obsgrd(ictype)%ngrdsch_j = ceiling(hori_loc_ctype(ictype) * dist_zero_fac / obsgrd(ictype)%grdspc_j)
    obsgrd(ictype)%ngrdext_i = obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i * 2
    obsgrd(ictype)%ngrdext_j = obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j * 2

    allocate (obsgrd(ictype)%n (  obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j, 0:nprocs_d-1))
    allocate (obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j, 0:nprocs_d-1))
    allocate (obsgrd(ictype)%tot(0:nprocs_d-1))
    allocate (obsgrd(ictype)%n_ext (  obsgrd(ictype)%ngrdext_i, obsgrd(ictype)%ngrdext_j))
    allocate (obsgrd(ictype)%ac_ext(0:obsgrd(ictype)%ngrdext_i, obsgrd(ictype)%ngrdext_j))

    obsgrd(ictype)%n (:,:,:) = 0
    obsgrd(ictype)%ac(:,:,:) = 0
    obsgrd(ictype)%tot(:) = 0
    obsgrd(ictype)%n_ext (:,:) = 0
    obsgrd(ictype)%ac_ext(:,:) = 0
    obsgrd(ictype)%tot_ext = 0
    obsgrd(ictype)%tot_sub(:) = 0
    obsgrd(ictype)%tot_g(:) = 0

    allocate (obsgrd(ictype)%next(obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j))
  end do

  call mpi_timer('set_efso_obs:bucket_sort_prepare:', 2)

  ! First scan: count the observation numbers in each mesh (in each subdomian)
  !-----------------------------------------------------------------------------

  do n = 1, obsda%nobs
    iof = obsda%set(n)
    iidx = obsda%idx(n)
    ictype = ctype_elmtyp(uid_obs(obs(iof)%elm(iidx)), obs(iof)%typ(iidx))

    !if (obsda%qc(n) == iqc_good) then
      call phys2ij( obs(iof)%lon(iidx), obs(iof)%lat(iidx), obs(iof)%ri(iidx), obs(iof)%rj(iidx) )
      call ij_obsgrd(ictype, obs(iof)%ri(iidx), obs(iof)%rj(iidx), i, j)
      if (i < 1) i = 1                                          ! Assume the process assignment was correct,
      if (i > obsgrd(ictype)%ngrd_i) i = obsgrd(ictype)%ngrd_i  ! so this correction is only to remedy the round-off problem.
      if (j < 1) j = 1                                          !
      if (j > obsgrd(ictype)%ngrd_j) j = obsgrd(ictype)%ngrd_j  !

      obsgrd(ictype)%n(i,j,myrank_d) = obsgrd(ictype)%n(i,j,myrank_d) + 1
      obsgrd(ictype)%tot_sub(i_after_qc) = obsgrd(ictype)%tot_sub(i_after_qc) + 1 ! only used for diagnostic print (obs number after qc)
    !end if

    obsgrd(ictype)%tot_sub(i_before_qc) = obsgrd(ictype)%tot_sub(i_before_qc) + 1 ! only used for diagnostic print (obs number before qc)
  end do

  ! Compute the accumulated numbers in each mesh
  !-----------------------------------------------------------------------------

  do ictype = 1, nctype
    if (ictype > 1) then
      obsgrd(ictype)%ac(0,1,myrank_d) = obsgrd(ictype-1)%ac(obsgrd(ictype-1)%ngrd_i,obsgrd(ictype-1)%ngrd_j,myrank_d)
    end if
    do j = 1, obsgrd(ictype)%ngrd_j
      if (j > 1) then
        obsgrd(ictype)%ac(0,j,myrank_d) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,j-1,myrank_d)
      end if
      do i = 1, obsgrd(ictype)%ngrd_i
        obsgrd(ictype)%ac(i,j,myrank_d) = obsgrd(ictype)%ac(i-1,j,myrank_d) + obsgrd(ictype)%n(i,j,myrank_d)
!        obsgrd(ictype)%next(i,j) = obsgrd(ictype)%ac(i-1,j,myrank_d)
      end do
    end do
    obsgrd(ictype)%next(1:obsgrd(ictype)%ngrd_i,:) = obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i-1,:,myrank_d)
  end do

  call mpi_timer('set_efso_obs:bucket_sort_first_scan:', 2)

  ! Second scan: save the indices of bucket-sorted observations in obsda%keys(:)
  !-----------------------------------------------------------------------------

  do n = 1, obsda%nobs
    !if (obsda%qc(n) == iqc_good) then
      iof = obsda%set(n)
      iidx = obsda%idx(n)
      ictype = ctype_elmtyp(uid_obs(obs(iof)%elm(iidx)), obs(iof)%typ(iidx))

      call ij_obsgrd(ictype, obs(iof)%ri(iidx), obs(iof)%rj(iidx), i, j)
      if (i < 1) i = 1                                         ! Assume the process assignment was correct,
      if (i > obsgrd(ictype)%ngrd_i) i = obsgrd(ictype)%ngrd_i ! so this correction is only to remedy the round-off problem.
      if (j < 1) j = 1                                         !
      if (j > obsgrd(ictype)%ngrd_j) j = obsgrd(ictype)%ngrd_j !

      obsgrd(ictype)%next(i,j) = obsgrd(ictype)%next(i,j) + 1
      obsda%key(obsgrd(ictype)%next(i,j)) = n
    !end if
  end do


  call mpi_timer('set_efso_obs:bucket_sort_second_scan:', 2, barrier=MPI_COMM_d)

  ! ALLREDUCE observation number information from subdomains, and compute total numbers
  !-----------------------------------------------------------------------------

  nobs_sub(:) = 0
  nobs_g(:) = 0
  do ictype = 1, nctype
    if (nprocs_d > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsgrd(ictype)%n, obsgrd(ictype)%ngrd_i*obsgrd(ictype)%ngrd_j*nprocs_d, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i,:,:), (obsgrd(ictype)%ngrd_i+1)*obsgrd(ictype)%ngrd_j*nprocs_d, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
    end if
    call MPI_ALLREDUCE(obsgrd(ictype)%tot_sub, obsgrd(ictype)%tot_g, n_qc_steps, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)

    if (ictype == 1) then
      obsgrd(ictype)%tot(:) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,obsgrd(ictype)%ngrd_j,:)
    else
      obsgrd(ictype)%tot(:) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,obsgrd(ictype)%ngrd_j,:) &
                            - obsgrd(ictype-1)%ac(obsgrd(ictype-1)%ngrd_i,obsgrd(ictype-1)%ngrd_j,:)
    end if
#ifdef LETKF_DEBUG
    if (obsgrd(ictype)%tot(myrank_d) /= obsgrd(ictype)%tot_sub(i_after_qc)) then
      write (6, '(A)') '[Error] Observation counts are inconsistent !!!'
      stop 99
    end if
#endif

    nobs_sub(:) = nobs_sub(:) + obsgrd(ictype)%tot_sub(:)
    nobs_g(:) = nobs_g(:) + obsgrd(ictype)%tot_g(:)

    deallocate (obsgrd(ictype)%next)
  end do

  call mpi_timer('set_efso_obs:bucket_sort_info_allreduce:', 2)

  nobstotalg = nobs_g(i_after_qc) ! total obs number in the global domain (all types)

  ! Calculate observation numbers in the extended (localization) subdomain,
  ! in preparation for communicating obsetvations in the extended subdomain
  !-----------------------------------------------------------------------------

  call rank_1d_2d(myrank_d, myp_i, myp_j)

  do ictype = 1, nctype
    imin1 = myp_i*obsgrd(ictype)%ngrd_i+1 - obsgrd(ictype)%ngrdsch_i
    imax1 = (myp_i+1)*obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
    jmin1 = myp_j*obsgrd(ictype)%ngrd_j+1 - obsgrd(ictype)%ngrdsch_j
    jmax1 = (myp_j+1)*obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

    do ip = 0, nprocs_d-1
      call rank_1d_2d(ip, ip_i, ip_j)
      imin2 = max(1, imin1 - ip_i*obsgrd(ictype)%ngrd_i)
      imax2 = min(obsgrd(ictype)%ngrd_i, imax1 - ip_i*obsgrd(ictype)%ngrd_i)
      jmin2 = max(1, jmin1 - ip_j*obsgrd(ictype)%ngrd_j)
      jmax2 = min(obsgrd(ictype)%ngrd_j, jmax1 - ip_j*obsgrd(ictype)%ngrd_j)
      if (imin2 > imax2 .or. jmin2 > jmax2) cycle

      ishift = (ip_i - myp_i) * obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
      jshift = (ip_j - myp_j) * obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j
      obsgrd(ictype)%n_ext(imin2+ishift:imax2+ishift, jmin2+jshift:jmax2+jshift) = obsgrd(ictype)%n(imin2:imax2, jmin2:jmax2, ip)
    end do

    if (ictype > 1) then
      obsgrd(ictype)%ac_ext(0,1) = obsgrd(ictype-1)%ac_ext(obsgrd(ictype-1)%ngrdext_i,obsgrd(ictype-1)%ngrdext_j)
    end if
    do j = 1, obsgrd(ictype)%ngrdext_j
      if (j > 1) then
        obsgrd(ictype)%ac_ext(0,j) = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,j-1)
      end if
      do i = 1, obsgrd(ictype)%ngrdext_i
        obsgrd(ictype)%ac_ext(i,j) = obsgrd(ictype)%ac_ext(i-1,j) + obsgrd(ictype)%n_ext(i,j)
      end do
    end do
    if (ictype == 1) then
      obsgrd(ictype)%tot_ext = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,obsgrd(ictype)%ngrdext_j)
    else
      obsgrd(ictype)%tot_ext = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,obsgrd(ictype)%ngrdext_j) &
                             - obsgrd(ictype-1)%ac_ext(obsgrd(ictype-1)%ngrdext_i,obsgrd(ictype-1)%ngrdext_j)
    end if
  end do ! [ ictype = 1, nctype ]

  if (nctype > 0) then
    nobstotal = obsgrd(nctype)%ac_ext(obsgrd(nctype)%ngrdext_i,obsgrd(nctype)%ngrdext_j) ! total obs number in the extended subdomain (all types)

!    maxnobs_sub_per_ctype = maxval(obsgrd(1)%tot(:))
    maxnobs_per_ctype = obsgrd(1)%tot_ext
    do ictype = 2, nctype
!      maxnobs_sub_per_ctype = max(maxnobs_sub_per_ctype, maxval(obsgrd(ictype)%tot(:)))
      maxnobs_per_ctype = max(maxnobs_per_ctype, obsgrd(ictype)%tot_ext)
    end do
  else
    nobstotal = 0
!    maxnobs_sub_per_ctype = 0
    maxnobs_per_ctype = 0
  end if

  call mpi_timer('set_efso_obs:extdomain_obs_count_cal:', 2)

  ! Construct sorted obsda_sort: 
  !-----------------------------------------------------------------------------

  ! 1) Copy the observation data in own subdomains to send buffer (obsbufs) with
  ! sorted order
  !-----------------------------------------------------------------------------
  if (nctype > 0) then
    cntr(:) = obsgrd(nctype)%ac(obsgrd(nctype)%ngrd_i,obsgrd(nctype)%ngrd_j,:)
    cnts = cntr(myrank_d+1)
    dspr(1) = 0
    do ip = 2, nprocs_d
      dspr(ip) = dspr(ip-1) + cntr(ip-1)
    end do

    nensobs_mod = mod(nensobs, nprocs_e)
    nensobs_div = (nensobs - nensobs_mod) / nprocs_e
    if (myrank_e < nensobs_mod) then
      im_obs_1 = (nensobs_div+1) * myrank_e + 1
      im_obs_2 = (nensobs_div+1) * (myrank_e+1)
      nensobs_part = nensobs_div + 1
    else
      im_obs_1 = (nensobs_div+1) * nensobs_mod + nensobs_div * (myrank_e-nensobs_mod) + 1
      im_obs_2 = (nensobs_div+1) * nensobs_mod + nensobs_div * (myrank_e-nensobs_mod+1)
      nensobs_part = nensobs_div
    end if

    obsbufs%nobs = nobs_sub(i_after_qc)
    call obs_da_value_allocate(obsbufs, nensobs_part)

    do n = 1, nobs_sub(i_after_qc)
      obsbufs%set(n) = obsda%set(obsda%key(n))
      obsbufs%idx(n) = obsda%idx(obsda%key(n))
      obsbufs%val(n) = obsda%val(obsda%key(n))
      if (nensobs_part > 0) then
        obsbufs%ensval(1:nensobs_part,n) = obsda%ensval(im_obs_1:im_obs_2,obsda%key(n))
      end if
      obsbufs%qc(n) = obsda%qc(obsda%key(n))
    end do ! [ n = 1, nobs_sub(i_after_qc) ]

    call mpi_timer('set_efso_obs:copy_bufs:', 2, barrier=MPI_COMM_d)

  end if ! [ nctype > 0 ]

  call obs_da_value_deallocate(obsda)

  ! 2) Communicate to get global observations;
  !    for variables with an ensemble dimension (ensval),
  !    only obtain data from partial members (nensobs_part) to save memory usage
  !-----------------------------------------------------------------------------
  if (nctype > 0) then
    obsbufr%nobs = nobs_g(i_after_qc)
    call obs_da_value_allocate(obsbufr, nensobs_part)

    call MPI_ALLGATHERV(obsbufs%set, cnts, MPI_INTEGER, obsbufr%set, cntr, dspr, MPI_INTEGER, MPI_COMM_d, ierr)
    call MPI_ALLGATHERV(obsbufs%idx, cnts, MPI_INTEGER, obsbufr%idx, cntr, dspr, MPI_INTEGER, MPI_COMM_d, ierr)
    call MPI_ALLGATHERV(obsbufs%val, cnts, MPI_r_size, obsbufr%val, cntr, dspr, MPI_r_size, MPI_COMM_d, ierr)
    if (nensobs_part > 0) then
      call MPI_ALLGATHERV(obsbufs%ensval, cnts*nensobs_part, MPI_r_size, obsbufr%ensval, cntr*nensobs_part, dspr*nensobs_part, MPI_r_size, MPI_COMM_d, ierr)
    end if
    call MPI_ALLGATHERV(obsbufs%qc, cnts, MPI_INTEGER, obsbufr%qc, cntr, dspr, MPI_INTEGER, MPI_COMM_d, ierr)

    call obs_da_value_deallocate(obsbufs)

    call mpi_timer('set_efso_obs:mpi_allgatherv:', 2)
  end if ! [ nctype > 0 ]

  ! 3) Copy observation data within the extended (localization) subdomains
  !    from receive buffer (obsbufr) to obsda_sort; rearrange with sorted order
  !-----------------------------------------------------------------------------
  obsda_sort%nobs = nobstotal
  call obs_da_value_allocate(obsda_sort, nensobs)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ip,ip_i,ip_j,ictype,imin1,imax1,jmin1,jmax1,imin2,imax2,jmin2,jmax2,ishift,jshift,j,ns_ext,ne_ext,ns_bufr,ne_bufr)
  do ip = 0, nprocs_d-1
    call rank_1d_2d(ip, ip_i, ip_j)

    do ictype = 1, nctype
      imin1 = myp_i*obsgrd(ictype)%ngrd_i+1 - obsgrd(ictype)%ngrdsch_i
      imax1 = (myp_i+1)*obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
      jmin1 = myp_j*obsgrd(ictype)%ngrd_j+1 - obsgrd(ictype)%ngrdsch_j
      jmax1 = (myp_j+1)*obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

      imin2 = max(1, imin1 - ip_i*obsgrd(ictype)%ngrd_i)
      imax2 = min(obsgrd(ictype)%ngrd_i, imax1 - ip_i*obsgrd(ictype)%ngrd_i)
      jmin2 = max(1, jmin1 - ip_j*obsgrd(ictype)%ngrd_j)
      jmax2 = min(obsgrd(ictype)%ngrd_j, jmax1 - ip_j*obsgrd(ictype)%ngrd_j)
      if (imin2 > imax2 .or. jmin2 > jmax2) cycle

      ishift = (ip_i - myp_i) * obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
      jshift = (ip_j - myp_j) * obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

      do j = jmin2, jmax2
        ns_ext = obsgrd(ictype)%ac_ext(imin2+ishift-1,j+jshift) + 1
        ne_ext = obsgrd(ictype)%ac_ext(imax2+ishift  ,j+jshift)
        if (ns_ext > ne_ext) cycle

        ns_bufr = dspr(ip+1) + obsgrd(ictype)%ac(imin2-1,j,ip) + 1
        ne_bufr = dspr(ip+1) + obsgrd(ictype)%ac(imax2  ,j,ip)

        obsda_sort%set(ns_ext:ne_ext) = obsbufr%set(ns_bufr:ne_bufr)
        obsda_sort%idx(ns_ext:ne_ext) = obsbufr%idx(ns_bufr:ne_bufr)
        obsda_sort%val(ns_ext:ne_ext) = obsbufr%val(ns_bufr:ne_bufr)
        if (nensobs_part > 0) then
          obsda_sort%ensval(im_obs_1:im_obs_2,ns_ext:ne_ext) = obsbufr%ensval(1:nensobs_part,ns_bufr:ne_bufr)
        end if
        obsda_sort%qc(ns_ext:ne_ext) = obsbufr%qc(ns_bufr:ne_bufr)
      end do
    end do ! [ ictype = 1, nctype ]
  end do ! [ ip = 0, nprocs_d-1 ]
!$OMP END PARALLEL DO

  ! Save the keys of observations within the subdomain (excluding the localization buffer area)
  obsda_sort%nobs_in_key = 0
  do ictype = 1, nctype
    imin1 = obsgrd(ictype)%ngrdsch_i + 1
    imax1 = obsgrd(ictype)%ngrdsch_i + obsgrd(ictype)%ngrd_i
    jmin1 = obsgrd(ictype)%ngrdsch_j + 1
    jmax1 = obsgrd(ictype)%ngrdsch_j + obsgrd(ictype)%ngrd_j

    call obs_choose_ext(ictype, imin1, imax1, jmin1, jmax1, obsda_sort%nobs_in_key, obsda_sort%key)

    deallocate (obsgrd(ictype)%n)
    deallocate (obsgrd(ictype)%ac)
    deallocate (obsgrd(ictype)%n_ext)
  end do ! [ ictype = 1, nctype ]

  if (nctype > 0) then
    call obs_da_value_deallocate(obsbufr)

    call mpi_timer('set_efso_obs:extdomain_copy_bufr:', 2, barrier=MPI_COMM_e)
  end if ! [ nctype > 0 ]

  ! 4) For variables with an ensemble dimension (ensval),
  !    ALLREDUCE among the ensemble dimension to obtain data of all members
  !-----------------------------------------------------------------------------
  if (nprocs_e > 1) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsda_sort%ensval, nensobs*nobstotal, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)

    call mpi_timer('set_efso_obs:extdomain_allreduce:', 2)
  end if



!!
!! Create observation box
!!
!  nobsgrd = 0
!  nj = 0
!!$OMP PARALLEL PRIVATE(i,j,n,nn)
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    DO n=1,nobs
!      IF(obslat(n) < lat(j) .OR. lat(j+1) <= obslat(n)) CYCLE
!      nj(j) = nj(j) + 1
!    END DO
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    njs(j) = SUM(nj(0:j-1))
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    nn = 0
!    DO n=1,nobs
!      IF(obslat(n) < lat(j) .OR. lat(j+1) <= obslat(n)) CYCLE
!      nn = nn + 1
!    END DO
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    IF(nj(j) == 0) THEN
!      nobsgrd(:,j) = njs(j)
!      CYCLE
!    END IF
!    nn = 0
!    DO i=1,nlon
!      DO n=njs(j)+1,njs(j)+nj(j)
!        IF(i < nlon) THEN
!          IF(obslon(n) < lon(i) .OR. lon(i+1) <= obslon(n)) CYCLE
!        ELSE
!          IF(obslon(n) < lon(nlon) .OR. 360.0d0 <= obslon(n)) CYCLE
!        END IF
!        nn = nn + 1
!      END DO
!      nobsgrd(i,j) = njs(j) + nn
!    END DO
!    IF(nn /= nj(j)) THEN
!!$OMP CRITICAL
!      WRITE(6,'(A,2I10)') 'OBS DATA SORT ERROR: ',nn,nj(j)
!      WRITE(6,'(F6.2,A,F6.2)') lat(j),'< LAT <',lat(j+1)
!      WRITE(6,'(F6.2,A,F6.2)') MINVAL(obslat(njs(j)+1:njs(j)+nj(j))),'< OBSLAT <',MAXVAL(obslat(njs(j)+1:njs(j)+nj(j)))
!!$OMP END CRITICAL
!    END IF
!  END DO
!!$OMP END DO
!!$OMP END PARALLEL

  return
end subroutine set_efso_obs

end module letkf_obs
