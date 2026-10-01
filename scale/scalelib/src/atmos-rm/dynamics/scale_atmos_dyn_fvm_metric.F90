!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics FVM metric coefficients
!!
!! @par Description
!!          Time-invariant metric, map-factor and Coriolis coefficients
!!          shared by the FVM short time step schemes
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_fvm_metric
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_atmos_grid_cartesC_index
  use scale_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_FVM_metric_setup
  public :: ATMOS_DYN_FVM_metric_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  ! The following arrays are set in ATMOS_DYN_FVM_metric_setup and are
  ! present on the device (OpenACC) until ATMOS_DYN_FVM_metric_finalize.

  ! map-factor coefficients (always allocated)
  real(RP), public, allocatable :: RMAPF   (:,:,:,:) ! (IA,JA,2,I_XY_MAX) 1 / MAPF

  ! weight of the full level k+1 used to interpolate to the half level k+1/2;
  ! the weight of the level k is 1 - F2H (always allocated)
  real(RP), public, allocatable :: F2H_UYZ (:,:,:)   ! (KA,IA,JA) at (u,y,z)
  real(RP), public, allocatable :: F2H_XVZ (:,:,:)   ! (KA,IA,JA) at (x,v,z)

  ! allocated only if full = .true.
  real(RP), public, allocatable :: MAPF_M12(:,:,:)   ! (IA,JA,I_XY_MAX)   MAPF(1) * MAPF(2)
  real(RP), public, allocatable :: MAPF_R12(:,:,:)   ! (IA,JA,I_XY_MAX)   1 / ( MAPF(1) * MAPF(2) )

  real(RP), public, allocatable :: RGSQRT_XYZ(:,:,:) ! (KA,IA,JA) 1 / GSQRT at (x,y,z)
  real(RP), public, allocatable :: RGSQRT_XYW(:,:,:) ! (KA,IA,JA) 1 / GSQRT at (x,y,w)
  real(RP), public, allocatable :: RGSQRT_UYZ(:,:,:) ! (KA,IA,JA) 1 / GSQRT at (u,y,z)
  real(RP), public, allocatable :: RGSQRT_XVZ(:,:,:) ! (KA,IA,JA) 1 / GSQRT at (x,v,z)

  ! coefficients of the Coriolis and metric terms
  ! COEF_CF_U at (u,y): 0.5*(f(i+1)+f(i)), m1*m2, d(1/m2)/dx, d(1/m1)/dy
  ! COEF_CF_V at (x,v): 0.5*(f(j+1)+f(j)), m1*m2, d(1/m2)/dx, d(1/m1)/dy
  real(RP), public, allocatable :: COEF_CF_U(:,:,:)  ! (4,IA,JA)
  real(RP), public, allocatable :: COEF_CF_V(:,:,:)  ! (4,IA,JA)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
#if 1
#define F2H(k,idx,i,j) (CDZ(k)*GSQRT(k,i,j,idx)/(CDZ(k)*GSQRT(k,i,j,idx)+CDZ(k+1)*GSQRT(k+1,i,j,idx)))
#else
#define F2H(k,idx,i,j) 0.5_RP
#endif
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: allocated_full = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  !! full = .false. sets RMAPF and F2H_* only, which is all the split-explicit
  !! and HIVI schemes need, to avoid allocating unused 3D arrays
  subroutine ATMOS_DYN_FVM_metric_setup( &
       CORIOLI,                &
       MAPF, GSQRT,            &
       CDZ,                    &
       RCDX, RCDY, RFDX, RFDY, &
       full                    )
    implicit none

    real(RP), intent(in) :: CORIOLI(IA,JA)
    real(RP), intent(in) :: MAPF   (IA,JA,2,I_XY_MAX)
    real(RP), intent(in) :: GSQRT  (KA,IA,JA,I_XYZ_MAX)
    real(RP), intent(in) :: CDZ (KA)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in) :: RFDX(IA-1)
    real(RP), intent(in) :: RFDY(JA-1)
    logical,  intent(in) :: full

    integer :: k, i, j, n
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_DYN_FVM_metric_setup",*) 'Setup'

    allocated_full = full

    ! map-factor coefficients
    allocate( RMAPF(IA,JA,2,I_XY_MAX) )
    do n = 1, I_XY_MAX
    do j = 1, JA
    do i = 1, IA
       RMAPF(i,j,1,n) = 1.0_RP / MAPF(i,j,1,n)
       RMAPF(i,j,2,n) = 1.0_RP / MAPF(i,j,2,n)
    enddo
    enddo
    enddo
    !$acc enter data copyin(RMAPF)

    ! F2H weights
    allocate( F2H_UYZ(KA,IA,JA) )
    allocate( F2H_XVZ(KA,IA,JA) )
    do j = 1, JA
    do i = 1, IA
       do k = 1, KA-1
          F2H_UYZ(k,i,j) = F2H(k,I_UYZ,i,j)
          F2H_XVZ(k,i,j) = F2H(k,I_XVZ,i,j)
       enddo
       F2H_UYZ(KA,i,j) = 0.0_RP
       F2H_XVZ(KA,i,j) = 0.0_RP
    enddo
    enddo
    !$acc enter data copyin(F2H_UYZ, F2H_XVZ)

    if ( .not. full ) return

    allocate( MAPF_M12(IA,JA,I_XY_MAX) )
    allocate( MAPF_R12(IA,JA,I_XY_MAX) )
    do n = 1, I_XY_MAX
    do j = 1, JA
    do i = 1, IA
       MAPF_M12(i,j,n) = MAPF(i,j,1,n) * MAPF(i,j,2,n)
       MAPF_R12(i,j,n) = 1.0_RP / MAPF_M12(i,j,n)
    enddo
    enddo
    enddo
    !$acc enter data copyin(MAPF_M12, MAPF_R12)

    allocate( RGSQRT_XYZ(KA,IA,JA) )
    allocate( RGSQRT_XYW(KA,IA,JA) )
    allocate( RGSQRT_UYZ(KA,IA,JA) )
    allocate( RGSQRT_XVZ(KA,IA,JA) )
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RGSQRT_XYZ(k,i,j) = 1.0_RP / GSQRT(k,i,j,I_XYZ)
       RGSQRT_XYW(k,i,j) = 1.0_RP / GSQRT(k,i,j,I_XYW)
       RGSQRT_UYZ(k,i,j) = 1.0_RP / GSQRT(k,i,j,I_UYZ)
       RGSQRT_XVZ(k,i,j) = 1.0_RP / GSQRT(k,i,j,I_XVZ)
    enddo
    enddo
    enddo
    !$acc enter data copyin(RGSQRT_XYZ, RGSQRT_XYW, RGSQRT_UYZ, RGSQRT_XVZ)

    allocate( COEF_CF_U(4,IA,JA) )
    allocate( COEF_CF_V(4,IA,JA) )

    ! at (u, y, z)
    COEF_CF_U(:,:,:) = 0.0_RP
    do j = 2, JA
    do i = 1, IA-1
       COEF_CF_U(1,i,j) = 0.5_RP * ( CORIOLI(i+1,j) + CORIOLI(i,j) )
       COEF_CF_U(2,i,j) = MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY)
       COEF_CF_U(3,i,j) = ( 1.0_RP/MAPF(i+1,j,2,I_XY) - 1.0_RP/MAPF(i,j  ,2,I_XY) ) * RFDX(i)
       COEF_CF_U(4,i,j) = ( 1.0_RP/MAPF(i  ,j,1,I_UV) - 1.0_RP/MAPF(i,j-1,1,I_UV) ) * RCDY(j)
    enddo
    enddo

    ! at (x, v, z)
    COEF_CF_V(:,:,:) = 0.0_RP
    do j = 1, JA-1
    do i = 2, IA
       COEF_CF_V(1,i,j) = 0.5_RP * ( CORIOLI(i,j+1) + CORIOLI(i,j) )
       COEF_CF_V(2,i,j) = MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV)
       COEF_CF_V(3,i,j) = ( 1.0_RP/MAPF(i,j  ,2,I_UV) - 1.0_RP/MAPF(i-1,j,2,I_UV) ) * RCDX(i)
       COEF_CF_V(4,i,j) = ( 1.0_RP/MAPF(i,j+1,1,I_XY) - 1.0_RP/MAPF(i  ,j,1,I_XY) ) * RFDY(j)
    enddo
    enddo
    !$acc enter data copyin(COEF_CF_U, COEF_CF_V)

    return
  end subroutine ATMOS_DYN_FVM_metric_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_DYN_FVM_metric_finalize
    implicit none
    !---------------------------------------------------------------------------

    if ( .not. allocated(RMAPF) ) return

    !$acc exit data delete(RMAPF, F2H_UYZ, F2H_XVZ)
    deallocate( RMAPF   )
    deallocate( F2H_UYZ )
    deallocate( F2H_XVZ )

    if ( allocated_full ) then
       !$acc exit data delete(MAPF_M12, MAPF_R12)
       deallocate( MAPF_M12 )
       deallocate( MAPF_R12 )

       !$acc exit data delete(RGSQRT_XYZ, RGSQRT_XYW, RGSQRT_UYZ, RGSQRT_XVZ)
       deallocate( RGSQRT_XYZ )
       deallocate( RGSQRT_XYW )
       deallocate( RGSQRT_UYZ )
       deallocate( RGSQRT_XVZ )

       !$acc exit data delete(COEF_CF_U, COEF_CF_V)
       deallocate( COEF_CF_U )
       deallocate( COEF_CF_V )
    end if

    allocated_full = .false.

    return
  end subroutine ATMOS_DYN_FVM_metric_finalize

end module scale_atmos_dyn_fvm_metric
