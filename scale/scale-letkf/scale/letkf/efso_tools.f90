MODULE efso_tools
!=======================================================================
!
! [PURPOSE:] Module for observation sensitivity calculation
!
! [HISTORY:]
!   07/27/2011 Yoichiro Ohta  created
!   09/29/2011 Yoichiro Ohta  adapted new formulation
!   07/01/2013 Daisuke Hotta  ported to GFS-LETKF system
!   12/19/2013 Guo-Yuan Lien  merged to GFS-LETKF main development
!   01/02/2013 Guo-Yuan Lien  modify output format
!
!=======================================================================
  use common
  use common_scale
  use common_obs_scale
  use letkf_obs
  use common_mpi_scale, only: &
     nij1
!  USE sigio_module

  implicit none

  private
  public lnorm, init_obsense, destroy_obsense !,print_obsense,lnorm,loc_advection
  public obsense, obsense_global, lon2, lat2, nterm

  real(r_size), allocatable :: obsense(:,:)
  real(r_size), allocatable :: obsense_global(:,:)
  real(r_size), ALLOCATABLE :: lon2(:,:)
  real(r_size), ALLOCATABLE :: lat2(:,:)
  integer, parameter :: nterm = 3

contains

subroutine init_obsense( use_global )
  implicit none

  logical, intent(in), optional :: use_global

  logical :: use_global_ = .false.

  if ( present( use_global ) ) use_global_ = use_global

  if ( use_global_ ) then
    allocate( obsense_global(nterm,nobstotalg) )
  else
    allocate( obsense(nterm,nobstotal) )
  endif

  return
end subroutine init_obsense

!-----------------------------------------------------------------------
! Compute norm
! [ref: Eq.(6,7,9), Ota et al. 2013]
!-----------------------------------------------------------------------
! [INPUT]
!  fcst3d,fcst2d: (xmean+X)^f_t  (total field)
!  fcer3d,fcer2d: [1/2(K-1)](e^f_t+e^g_t)
! [OUTPUT]
!  fcst3d,fcst2d: C^(1/2)*X^f_t                    [(J/kg)^(1/2)]
!  fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)  [(J/kg)^(1/2)]
!-----------------------------------------------------------------------
subroutine lnorm(fcst3d,fcst2d,fcer3d,fcer2d)
  use scale_const, only:    &
     GRAV   => CONST_GRAV,  &
     Rdry   => CONST_Rdry,  &
     Rvap   => CONST_Rvap,  &
     CVdry  => CONST_CVdry, &
     LHV    => CONST_LHV0,  &
     PRE00  => CONST_PRE00
  use scale_tracer, only: TRACER_CV
  use scale_atmos_grid_cartesC, only: &
     CDZ => ATMOS_GRID_CARTESC_CDZ 
  use scale_atmos_grid_cartesC_index, only: &
     KHALO
  implicit none

  real(r_size), intent(inout) :: fcst3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(inout) :: fcst2d(nij1,nens,nv2d)
  real(r_size), intent(inout) :: fcer3d(nij1,nlev,nv3d)
  real(r_size), intent(inout) :: fcer2d(nij1,nv2d)

  real(r_size), parameter :: tref = 280.0_r_size
  real(r_size), parameter :: pref = 1000.0e+2_r_size
  real(r_size) :: tmptv(nij1,nlev)
  real(r_size) :: pdelta(nij1,nlev)
  real(r_size) :: weight
  real(r_size) :: rinbv, cptr, qweight, rdtrpr

  real(r_size) :: ps_inv(nij1)

  integer :: iret
  integer :: i, j, k
  integer :: iv3d, iv2d, m
  integer :: ij

  real(r_size) :: qdry, CVtot, Rtot, CPovCV
  real(r_size) :: wmoist

  ! Calculate ensemble mean of forecast
  call ensmean_grd(MEMBER, nens, nij1, fcst3d, fcst2d)

  ! Calculate ensemble forecast perturbations
  do i = 1, MEMBER
    fcst3d(:,:,i,:) = fcst3d(:,:,i,:) - fcst3d(:,:,mmean,:)
    fcst2d(:,i,:) = fcst2d(:,i,:) - fcst2d(:,mmean,:)
  end do

  do ij = 1, nij1
    ps_inv(ij) = 1.0_r_size / fcst3d(ij,1,mmean,iv3d_p)
  enddo

  if ( EFSO_USE_MOIST_ENERGY ) then
    wmoist = 1.0_r_size
  else
    wmoist = 0.0_r_size
  endif

!  rdtrpr = sqrt(rd*tref)/pref
!  ! For surface variables
!  IF(tar_minlev <= 1) THEN
  do iv2d = 1, nv2d
!      IF(i == iv2d_ps) THEN
!        !!! [(Rd*Tr)(dS/4pi)]^(1/2) * (ps'/Pr)
!        fcer2d(:,i) = rdtrpr * wg1(:) * fcer2d(:,i)
!        DO j=1,nbv
!          fcst2d(:,j,i) = rdtrpr * wg1(:) * fcst2d(:,j,i)
!        END DO
!      ELSE
        fcer2d(:,iv2d) = 0.0_r_size
        fcst2d(:,:,iv2d) = 0.0_r_size
!      END IF
  end do
!  ELSE
!    fcer2d(:,:) = 0.0_r_size
!    fcst2d(:,:,:) = 0.0_r_size
!  END IF

!$omp parallel private(k,ij,iv3d,qdry,CVtot,Rtot,CPovCV,weight,m,cptr,qweight)
!$omp do schedule(static) collapse(2)
  do k = 1, nlev
!    IF(k > tar_maxlev .or. k < tar_minlev) THEN
!      fcst3d(:,k,:,:) = 0.0_r_size
!      fcer3d(:,k,:) = 0.0_r_size
!      CYCLE
!    END IF
    do ij = 1, nij1

       qdry  = 1.0_r_size
       CVtot = 0.0_r_size
       do iv3d = iv3d_q, nv3d ! loop over all moisture variables
         qdry  = qdry - fcst3d(ij,k,mmean,iv3d)
         CVtot = CVtot + fcst3d(ij,k,mmean,iv3d) * real( TRACER_CV(iv3d-iv3d_q+1), kind=r_size )
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = real( Rdry, kind=r_size ) * qdry + real( Rvap, kind=r_size ) * fcst3d(ij,k,mmean,iv3d_q)
       CPovCV = ( CVtot + Rtot ) / CVtot

       ! Compute weight
       ! rho * g * dz / p_s
       weight = sqrt( fcst3d(ij,k,mmean,iv3d_p) / ( fcst3d(ij,k,mmean,iv3d_t)  * Rtot ) &
                      * real( GRAV, kind=r_size ) * real( CDZ(k+KHALO), kind=r_size ) * ps_inv(ij) )
       ! Constants
       cptr = sqrt( CPovCV / tref )
       qweight = sqrt( wmoist/( CPovCV*tref ) ) * LHV

      do iv3d = 1, nv3d
        if ( iv3d == iv3d_u .or. iv3d == iv3d_v ) then
          !!! [(dsigma)(dS/4pi)]^(1/2) * u'
          !!! [(dsigma)(dS/4pi)]^(1/2) * v'
          fcer3d(ij,k,iv3d) = weight * fcer3d(ij,k,iv3d)
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d) = weight * fcst3d(ij,k,m,iv3d)
          enddo
        elseif (iv3d == iv3d_t) then
          !!! [(Cp/Tr)(dsigma)(dS/4pi)]^(1/2) * t'
          fcer3d(ij,k,i) = cptr * weight * fcer3d(ij,k,iv3d)
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d) = cptr * weight * fcst3d(ij,k,m,iv3d)
          enddo
        elseif (iv3d == iv3d_q) then
! specific humidity vs vapor concentration?
          !!! [(wg*L^2/Cp/Tr)(dsigma)(dS/4pi)]^(1/2) * q'
          fcer3d(ij,k,iv3d) = qweight * weight * fcer3d(ij,k,iv3d)
          do m = 1, MEMBER
            fcst3d(ij,k,m,iv3d) = qweight * weight * fcst3d(ij,k,m,iv3d)
          enddo
        else
          fcer3d(ij,k,iv3d) = 0.0_r_size
          fcst3d(ij,k,:,iv3d) = 0.0_r_size
        endif
      enddo ! iv3d
    enddo ! ij

  enddo ! k
!$omp end do
!$omp end parallel

  do i = 1, nij1
!    IF(lon1(i) < tar_minlon .or. lon1(i) > tar_maxlon .or. &
!         & lat1(i) < tar_minlat .or. lat1(i) > tar_maxlat) THEN
!      fcer2d(i,:) = 0.0_r_size
!      fcst2d(i,:,:) = 0.0_r_size
!      fcer3d(i,:,:) = 0.0_r_size
!      fcst3d(i,:,:,:) = 0.0_r_size
!    END IF
  end do

  return
end subroutine lnorm

SUBROUTINE loc_advection(ua,va,uf,vf)
  IMPLICIT NONE
  real(r_size),INTENT(IN) :: ua(nij1,nlev)
  real(r_size),INTENT(IN) :: va(nij1,nlev)
  real(r_size),INTENT(IN) :: uf(nij1,nlev)
  real(r_size),INTENT(IN) :: vf(nij1,nlev)
  real(r_size) :: rad2deg, deg2rad
  real(r_size) :: coslat(nij1)
  INTEGER :: i,k
!  ALLOCATE(lon2(nij1,nlev))
!  ALLOCATE(lat2(nij1,nlev))
!  deg2rad = pi/180.0_r_size
!  rad2deg = locadv_rate*eft*3600.0_r_size*180.0_r_size/(pi*re)
!  DO i=1,nij1
!    coslat(i) = 1.0_r_size/cos(lat1(i)*deg2rad)
!  END DO
!  DO k=1,nlev
!    DO i=1,nij1
!      lon2(i,k) = lon1(i) - 0.5_r_size * (ua(i,k) + uf(i,k)) &
!           & * coslat(i) * rad2deg
!      lat2(i,k) = lat1(i) - 0.5_r_size * (va(i,k) + vf(i,k)) &
!           & * rad2deg
!      IF(lat2(i,k) > 90.0_r_size) THEN
!        lat2(i,k) = 180.0_r_size - lat2(i,k)
!        lon2(i,k) = lon2(i,k) + 180.0_r_size
!      ELSE IF(lat2(i,k) < -90.0_r_size) THEN
!        lat2(i,k) = -180.0_r_size - lat2(i,k)
!        lon2(i,k) = lon2(i,k) + 180.0_r_size
!      END IF
!      IF(lon2(i,k) > 360.0_r_size) THEN
!        lon2(i,k) = MOD(lon2(i,k),360.0_r_size)
!      ELSE IF(lon2(i,k) < 0.0_r_size) THEN
!        lon2(i,k) = MOD(lon2(i,k),360.0_r_size) + 360.0_r_size
!      END IF
!    END DO
!  END DO
  RETURN
END SUBROUTINE loc_advection

SUBROUTINE print_obsense
  IMPLICIT NONE
  INTEGER,PARAMETER :: nreg = 3
  INTEGER :: regnh=1, regtr=2, regsh=3          ! Indices for the regions
  real(r_size),PARAMETER :: latbound=20._r_size ! Boundary latitude of TR
  INTEGER :: nobs_sense(nid_obs,nobtype+1,nreg)
  real(r_size) :: sumsense(nid_obs,nobtype+1,nreg)
  real(r_size) :: rate(nid_obs,nobtype+1,nreg)
  INTEGER :: nobs_t
  real(r_size) :: sumsense_t,rate_t
  INTEGER :: nob,oid,otype,ireg,iterm
  CHARACTER(len=2) :: charreg(nreg)
  CHARACTER(len=6) :: charotype
  CHARACTER(len=12) :: ofile(nterm)

!  IF(nobs == 0) RETURN
!  nobs_sense = 0
!  sumsense = 0._r_size
!  rate = 0._r_size
!  charreg(regnh) ='NH'
!  charreg(regtr) ='TR'
!  charreg(regsh) ='SH'
!  ofile(1)='osenseKE.dat'
!  ofile(2)='osensePE.dat'
!  ofile(3)='osenseME.dat'
!
!  ! Binary output (in obs2 format)
!  DO iterm = 1, 3
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',ofile(iterm)
!    CALL write_obs2(ofile(iterm),nobs,obselm,obslon,obslat,obslev, &
!                    obsdat,obserr,obstyp,obsdif,obsense(iterm,:),obsqc,0)
!  END DO
!
!  ! Loop over each observations
!  iterm = 1
!  DO nob=1,nobs
!    ! Select observation elements
!    oid = uid_obs(NINT(obselm(nob)))
!    IF(oid <= 0 .OR. oid > nid_obs) CYCLE
!    otype = NINT(obstyp(nob))
!    ! Select observation types
!    IF(otype <= 0 .OR. otype > nobtype+1) CYCLE
!    ! Select observation regions
!    IF(obslat(nob) > latbound) THEN
!      ireg=regnh
!    ELSE IF(obslat(nob) < -latbound) THEN
!      ireg=regsh
!    ELSE
!      ireg=regtr
!    END IF
!    ! Sum up
!    nobs_sense(oid,otype,ireg) = nobs_sense(oid,otype,ireg) + 1
!    sumsense(oid,otype,ireg) = sumsense(oid,otype,ireg) + obsense(iterm,nob)
!    IF(obsense(iterm,nob) < 0._r_size) THEN
!      rate(oid,otype,ireg) = rate(oid,otype,ireg) + 1._r_size
!    END IF
!  END DO
!
!  WRITE (6, '(A)') '============================================'
!  WRITE (6, '(A,I10)') ' TOTAL NUMBER OF OBSERVATIONS:', nobs
!  WRITE (6, '(A)') '============================================'
!  WRITE (6, '(A)') '              nobs     dJ(KE)       +rate[%]'
!  DO otype = 1,nobtype+1
!    IF(otype <= nobtype) THEN
!      charotype = obtypelist(otype)
!    ELSE
!      charotype = 'OTHERS'
!    END IF
!    nobs_t = SUM(nobs_sense(:,otype,:))
!    IF(nobs_t > 0) THEN
!      sumsense_t = SUM(sumsense(:,otype,:))
!      rate_t = SUM(rate(:,otype,:)) / real(nobs_t,r_size) * 100._r_size
!      WRITE (6, '(A)') '--------------------------------------------'
!      WRITE (6,'(A,1x,A,1x,I8,1x,E12.5,1x,F8.2)') &
!           & charotype,' TOTAL',nobs_t,sumsense_t,rate_t
!    END IF
!    DO ireg = 1,nreg
!      DO oid = 1,nid_obs
!        IF(nobs_sense(oid,otype,ireg) > 0) THEN
!          rate_t = rate(oid,otype,ireg) &
!               & / real(nobs_sense(oid,otype,ireg),r_size) * 100._r_size
!          WRITE (6,'(A,1x,A,1x,A,1x,I8,1x,E12.5,1x,F8.2)') &
!               & charotype,charreg(ireg),obelmlist(oid), &
!               & nobs_sense(oid,otype,ireg), &
!               & sumsense(oid,otype,ireg),   &
!               & rate_t
!        END IF
!      END DO
!    END DO
!  END DO
!  WRITE (6, '(A)') '============================================'

  RETURN
END SUBROUTINE print_obsense

subroutine destroy_obsense
  implicit none

  if( allocated(obsense) ) deallocate(obsense)
  if( allocated(lon2) ) deallocate(lon2)
  if( allocated(lat2) ) deallocate(lat2)

  return
end subroutine destroy_obsense

end module efso_tools
