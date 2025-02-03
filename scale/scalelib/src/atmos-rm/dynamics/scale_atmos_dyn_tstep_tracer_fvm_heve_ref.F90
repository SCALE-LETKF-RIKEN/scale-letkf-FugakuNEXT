!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          HEVE FVM scheme for tracer advection in Atmospheric dynamical process with reference
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tstep_tracer_fvm_heve_ref
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer

#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tstep_tracer_fvm_heve_ref

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  ! --> Ong Chia Rui starts ///////////////////////////////////
  ! Rui's code for AMPS non-mass advection
  real(RP), allocatable :: QICONflx_hi(:,:,:,:)
  integer,  allocatable :: QICONupwind_hi(:,:,:,:)
  real(RP), allocatable :: QICON(:,:,:)
  ! --> Ong Chia Rui ends   ///////////////////////////////////
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_Tstep_tracer_fvm_heve_ref( &
       QTRCo, qflx_hi, & ! (out)
       QTRC0, RHOQ_t, & ! (in)
       QTRC_ref, & ! (in)
       DENS0, DENS, & ! (in)
       qflx_ref, num_diff, & ! (in)
       GSQRT, MAPF, & ! (in)
       CDZ, RCDZ, RCDX, RCDY, & ! (in)
       TwoD, & ! (in)
       dtl ) ! (in)
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_fluxZ_XYZ_tracer, &
       ATMOS_DYN_FVM_fluxX_XYZ_tracer, &
       ATMOS_DYN_FVM_fluxY_XYZ_tracer
    implicit none
    real(RP), intent(out) :: QTRCo   (KA,IA,JA)
    real(RP), intent(out) :: qflx_hi (KA,IA,JA,3)
    real(RP), intent(in)  :: QTRC0   (KA,IA,JA)
    real(RP), intent(in)  :: RHOQ_t  (KA,IA,JA)
    real(RP), intent(in)  :: QTRC_ref(KA,IA,JA)
    real(RP), intent(in)  :: DENS0   (KA,IA,JA)
    real(RP), intent(in)  :: DENS    (KA,IA,JA)
    real(RP), intent(in)  :: qflx_ref(KA,IA,JA,3)
    real(RP), intent(in)  :: num_diff(KA,IA,JA,3)
    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,I_XYZ_MAX)
    real(RP), intent(in)  :: MAPF    (IA,JA,2)
    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: RCDZ(KA)
    real(RP), intent(in)  :: RCDX(IA)
    real(RP), intent(in)  :: RCDY(JA)
    logical,  intent(in)  :: TwoD
    real(RP), intent(in)  :: dtl

    real(RP) :: QTRC_N(KA,IA,JA)

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: IIS0, JJS0
    integer  :: i, j, k
    !---------------------------------------------------------------------------

    !$acc data copy(QTRCo) &
    !$acc      copyout(qflx_hi) &
    !$acc      copyin(QTRC0, RHOQ_t, QTRC_ref, DENS0, DENS, qflx_ref, num_diff, &
    !$acc             GSQRT, MAPF, CDZ, RCDZ, RCDX, RCDY) &
    !$acc      create(QTRC_N)

#ifdef DEBUG
    !$acc kernels
    qflx_hi(:,:,:,:) = UNDEF
    !$acc end kernels
#endif

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(IA,JA,KA,QTRC_N,QTRC0,QTRC_ref,EPS)
    !$acc kernels
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA

       if ( QTRC_ref(k,i,j) > EPS ) then
          QTRC_N(k,i,j) = QTRC0(k,i,j) / QTRC_ref(k,i,j)
       else
          QTRC_N(k,i,j) = 0.0_RP
       end if
       ! to detect unphysical values (temptative method)
       ! especially mpace60 faces a problem in daxis in AMPS
       if ( QTRC_N(k,i,j) > 1.e22_RP ) then
          QTRC_N(k,i,j) = 0.0_RP
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! at (x, y, w)
       call ATMOS_DYN_FVM_fluxZ_XYZ_tracer( qflx_hi(:,:,:,ZDIR), & ! (out)
            qflx_ref(:,:,:,ZDIR), QTRC_N(:,:,:), GSQRT(:,:,:,I_XYW), & ! (in)
            num_diff(:,:,:,ZDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, y, z)
       if ( .not. TwoD ) &
       call ATMOS_DYN_FVM_fluxX_XYZ_tracer( qflx_hi(:,:,:,XDIR), & ! (out)
            qflx_ref(:,:,:,XDIR), QTRC_N(:,:,:), GSQRT(:,:,:,I_UYZ), & ! (in)
            num_diff(:,:,:,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, v, z)
       call ATMOS_DYN_FVM_fluxY_XYZ_tracer( qflx_hi(:,:,:,YDIR), & ! (out)
            qflx_ref(:,:,:,YDIR), QTRC_N(:,:,:), GSQRT(:,:,:,I_XVZ), & ! (in)
            num_diff(:,:,:,YDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)


       if ( TwoD ) then
          !$omp parallel do default(none) private(j,k) OMP_SCHEDULE_ &
          !$omp shared(JJS,JJE,IS,KS,KE,QTRCo,QTRC0,DENS0,dtl,qflx_hi,RCDZ,RCDY,MAPF) &
          !$omp shared(GSQRT,RHOQ_t,DENS,I_XYZ)
          !$acc kernels
          do j = JJS, JJE
          do k = KS, KE
             QTRCo(k,IS,j) = ( QTRC0(k,IS,j) * DENS0(k,IS,j) &
                             + dtl * ( - ( ( qflx_hi(k,IS,j,ZDIR) - qflx_hi(k-1,IS,j  ,ZDIR)  ) * RCDZ(k) &
                                         + ( qflx_hi(k,IS,j,YDIR) - qflx_hi(k  ,IS,j-1,YDIR)  ) * RCDY(j) &
                                        ) * MAPF(IS,j,2) / GSQRT(k,IS,j,I_XYZ) &
                            + RHOQ_t(k,IS,j) ) ) / DENS(k,IS,j)
          enddo
          enddo
          !$acc end kernels
       else
          !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JJS,JJE,IIS,IIE,KS,KE,QTRCo,QTRC0,DENS0,dtl,qflx_hi,RCDZ,RCDX,RCDY,MAPF) &
          !$omp shared(GSQRT,RHOQ_t,DENS,I_XYZ)
          !$acc kernels
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             QTRCo(k,i,j) = ( QTRC0(k,i,j) * DENS0(k,i,j) &
                          + dtl * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR)  ) * RCDZ(k) &
                                      + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR)  ) * RCDX(i) &
                                      + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR)  ) * RCDY(j) &
                                  ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j,I_XYZ) &
                          + RHOQ_t(k,i,j) ) ) / DENS(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if

    enddo
    enddo

    !$acc end data
    return
  end subroutine ATMOS_DYN_Tstep_tracer_fvm_heve_ref

end module scale_atmos_dyn_tstep_tracer_fvm_heve_ref
