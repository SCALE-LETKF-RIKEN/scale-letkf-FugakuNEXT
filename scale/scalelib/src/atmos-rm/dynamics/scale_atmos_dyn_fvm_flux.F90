
!-------------------------------------------------------------------------------
!> module scale_atmos_dyn_fvm_flux
!!
!! @par Description
!!          FVM flux scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_fvm_flux
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
  use scale_prc
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
  public :: ATMOS_DYN_FVM_flux_setup

  abstract interface
     subroutine valueW( &
          valW, &
          mflx, val, GSQRT, &
          CDZ )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       implicit none
       real(RP), intent(out) :: valW (KA)
       real(RP), intent(in)  :: mflx (KA)
       real(RP), intent(in)  :: val  (KA)
       real(RP), intent(in)  :: GSQRT(KA)
       real(RP), intent(in)  :: CDZ(KA)
     end subroutine valueW
     subroutine flux_phi( &
          flux, &
          mflx, val, GSQRT, &
          num_diff, &
          CDZ, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       implicit none
       real(RP), intent(inout) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mflx    (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_phi
     subroutine flux_mom( &
          flux, &
          mom, val, DENS, &
          GSQRT, RMAPF, &
          num_diff, &
          CDZ, TwoD, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       implicit none
       real(RP), intent(inout) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mom     (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: RMAPF   (   IA,JA,2) !< 1 / map factor
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       logical,  intent(in)  :: TwoD
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_mom
     subroutine flux_z( &
          flux, &
          mom, val, DENS, &
          GSQRT, J33G, &
          num_diff, &
          CDZ, TwoD, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       implicit none
       real(RP), intent(inout) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mom     (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: J33G
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       logical,  intent(in)  :: TwoD
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_z
     subroutine flux_wz( &
          flux, &
          mom, val, DENS, &
          GSQRT, J33G, &
          num_diff, &
          CDZ, FDZ, &
          dtrk, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       implicit none
       real(RP), intent(inout) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mom     (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: J33G
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       real(RP), intent(in)  :: FDZ(KA-1)
       real(RP), intent(in)  :: dtrk
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_wz
     subroutine flux_j( &
          flux, &
          mom, val, DENS, &
          GSQRT, JG, RMAPF, &
          CDZ, TwoD, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       implicit none
       real(RP), intent(inout) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mom     (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: JG      (KA,IA,JA)
       real(RP), intent(in)  :: RMAPF   (   IA,JA,2) !< 1 / map factor
       real(RP), intent(in)  :: CDZ(KA)
       logical,  intent(in)  :: TwoD
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_j
     subroutine flux_zj_w( &
          flux, fluxJ13, fluxJ23, &
          momz, momx, momy, val, &
          DENS, GSQRT, &
          J33G, J13G, J23G, MAPF, &
          num_diff, &
          CDZ, FDZ, dtrk, TwoD, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       implicit none
       real(RP), intent(inout) :: flux    (KA,IA,JA)
       real(RP), intent(inout) :: fluxJ13 (KA,IA,JA)
       real(RP), intent(inout) :: fluxJ23 (KA,IA,JA)
       real(RP), intent(in)  :: momz    (KA,IA,JA)
       real(RP), intent(in)  :: momx    (KA,IA,JA)
       real(RP), intent(in)  :: momy    (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: J33G
       real(RP), intent(in)  :: J13G    (KA,IA,JA)
       real(RP), intent(in)  :: J23G    (KA,IA,JA)
       real(RP), intent(in)  :: MAPF    (   IA,JA,2) !< map factor
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       real(RP), intent(in)  :: FDZ(KA-1)
       real(RP), intent(in)  :: dtrk
       logical,  intent(in)  :: TwoD
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_zj_w
     subroutine flux_zj( &
          flux, fluxJ13, fluxJ23, &
          momz, momx, momy, val, &
          DENS, GSQRT, F2HW, &
          J33G, J13G, J23G, MAPF, &
          num_diff, &
          CDZ, TwoD, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       implicit none
       real(RP), intent(inout) :: flux    (KA,IA,JA)
       real(RP), intent(inout) :: fluxJ13 (KA,IA,JA)
       real(RP), intent(inout) :: fluxJ23 (KA,IA,JA)
       real(RP), intent(in)  :: momz    (KA,IA,JA)
       real(RP), intent(in)  :: momx    (KA,IA,JA)
       real(RP), intent(in)  :: momy    (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: F2HW    (KA,IA,JA)
       real(RP), intent(in)  :: J33G
       real(RP), intent(in)  :: J13G    (KA,IA,JA)
       real(RP), intent(in)  :: J23G    (KA,IA,JA)
       real(RP), intent(in)  :: MAPF    (   IA,JA,2) !< map factor
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       logical,  intent(in)  :: TwoD
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_zj
  end interface

#ifndef _OPENACC
  procedure(valueW), pointer :: ATMOS_DYN_FVM_flux_valueW_Z => NULL()
#endif
  public :: ATMOS_DYN_FVM_flux_valueW_Z


  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxZ_XYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XYZ
  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxZ_XYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XYZ_tracer

  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxX_XYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XYZ
  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxX_XYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XYZ_tracer

  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxY_XYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XYZ
  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxY_XYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XYZ_tracer


  ! The vertical fluxes of the momenta (the parts of J33G, J13G and J23G at XYW, UYZ and XVZ) are computed
  ! together by ATMOS_DYN_FVM_fluxZJ_XYW, ATMOS_DYN_FVM_fluxZJ_UYZ and ATMOS_DYN_FVM_fluxZJ_XVZ (issue #3754).
  ! With TwoD, fluxJ13 is not computed and momx is not used.
  ! For UYZ and XVZ, GSQRT is used by the part of J33G w, F2HW by the metric parts.
  procedure(flux_zj_w), pointer :: ATMOS_DYN_FVM_fluxZJ_XYW => NULL()
  public :: ATMOS_DYN_FVM_fluxZJ_XYW
  procedure(flux_zj),   pointer :: ATMOS_DYN_FVM_fluxZJ_UYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxZJ_UYZ
  procedure(flux_zj),   pointer :: ATMOS_DYN_FVM_fluxZJ_XVZ => NULL()
  public :: ATMOS_DYN_FVM_fluxZJ_XVZ

  procedure(flux_wz), pointer :: ATMOS_DYN_FVM_fluxZ_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XYW_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_XYW => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XYW
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XYW_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_XYW => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XYW
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XYW_tracer


  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_XYW_tracer

  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_XYW_tracer


  procedure(flux_z), pointer :: ATMOS_DYN_FVM_fluxZ_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_UYZ_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_UYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxX_UYZ
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxX_UYZ_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_UYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxY_UYZ
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxY_UYZ_tracer


  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_UYZ_tracer

  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_UYZ_tracer


  procedure(flux_z), pointer :: ATMOS_DYN_FVM_fluxZ_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XVZ_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_XVZ => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XVZ
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XVZ_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_XVZ => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XVZ
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XVZ_tracer


  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_XVZ_tracer

  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_XVZ_tracer



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
  !
#ifdef _OPENACC

  integer, parameter :: I_UD1 = 1

  integer, parameter :: I_CD2 = 2

  integer, parameter :: I_UD3 = 3

  integer, parameter :: I_UD3KOREN1993 = 4

  integer, parameter :: I_CD4 = 5

  integer, parameter :: I_UD5 = 6

  integer, parameter :: I_CD6 = 7

  integer, parameter :: I_UD7 = 8

  integer, parameter :: I_CD8 = 9

  integer :: i_scheme
  !$acc declare create(i_scheme)
#endif
  !-----------------------------------------------------------------------------
contains
  
  !-----------------------------------------------------------------------------
  !> setup
  subroutine ATMOS_DYN_FVM_flux_setup( &
       scheme, scheme_tracer )
    use scale_prc, only: &
         PRC_abort
    use scale_prc_cartesC, only: &
         PRC_TwoD


   use scale_atmos_dyn_fvm_flux_ud1, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud1, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud1, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud1, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud1, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud1, &
      ATMOS_DYN_FVM_fluxX_XYW_ud1, &
      ATMOS_DYN_FVM_fluxY_XYW_ud1, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud1, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud1, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud1, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud1, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud1, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud1, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud1, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud1, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud1, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud1, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud1, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud1, &
      ATMOS_DYN_FVM_fluxZJ_XYW_ud1, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_ud1, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_ud1

   use scale_atmos_dyn_fvm_flux_cd2, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd2, &
      ATMOS_DYN_FVM_fluxZ_XYZ_cd2, &
      ATMOS_DYN_FVM_fluxX_XYZ_cd2, &
      ATMOS_DYN_FVM_fluxY_XYZ_cd2, &
      ATMOS_DYN_FVM_fluxZ_XYW_cd2, &
      ATMOS_DYN_FVM_fluxX_XYW_cd2, &
      ATMOS_DYN_FVM_fluxY_XYW_cd2, &
      ATMOS_DYN_FVM_fluxJ13_XYW_cd2, &
      ATMOS_DYN_FVM_fluxJ23_XYW_cd2, &
      ATMOS_DYN_FVM_fluxZ_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxX_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxY_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxZ_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxX_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxY_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxZJ_XYW_cd2, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_cd2

   use scale_atmos_dyn_fvm_flux_ud3, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud3, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud3, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud3, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud3, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud3, &
      ATMOS_DYN_FVM_fluxX_XYW_ud3, &
      ATMOS_DYN_FVM_fluxY_XYW_ud3, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud3, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud3, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxZJ_XYW_ud3, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_ud3

   use scale_atmos_dyn_fvm_flux_ud3Koren1993, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxX_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxY_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZJ_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_ud3Koren1993

   use scale_atmos_dyn_fvm_flux_cd4, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd4, &
      ATMOS_DYN_FVM_fluxZ_XYZ_cd4, &
      ATMOS_DYN_FVM_fluxX_XYZ_cd4, &
      ATMOS_DYN_FVM_fluxY_XYZ_cd4, &
      ATMOS_DYN_FVM_fluxZ_XYW_cd4, &
      ATMOS_DYN_FVM_fluxX_XYW_cd4, &
      ATMOS_DYN_FVM_fluxY_XYW_cd4, &
      ATMOS_DYN_FVM_fluxJ13_XYW_cd4, &
      ATMOS_DYN_FVM_fluxJ23_XYW_cd4, &
      ATMOS_DYN_FVM_fluxZ_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxX_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxY_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxZ_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxX_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxY_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxZJ_XYW_cd4, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_cd4

   use scale_atmos_dyn_fvm_flux_ud5, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud5, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud5, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud5, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud5, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud5, &
      ATMOS_DYN_FVM_fluxX_XYW_ud5, &
      ATMOS_DYN_FVM_fluxY_XYW_ud5, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud5, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud5, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxZJ_XYW_ud5, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_ud5

   use scale_atmos_dyn_fvm_flux_cd6, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd6, &
      ATMOS_DYN_FVM_fluxZ_XYZ_cd6, &
      ATMOS_DYN_FVM_fluxX_XYZ_cd6, &
      ATMOS_DYN_FVM_fluxY_XYZ_cd6, &
      ATMOS_DYN_FVM_fluxZ_XYW_cd6, &
      ATMOS_DYN_FVM_fluxX_XYW_cd6, &
      ATMOS_DYN_FVM_fluxY_XYW_cd6, &
      ATMOS_DYN_FVM_fluxJ13_XYW_cd6, &
      ATMOS_DYN_FVM_fluxJ23_XYW_cd6, &
      ATMOS_DYN_FVM_fluxZ_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxX_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxY_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxZ_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxX_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxY_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxZJ_XYW_cd6, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_cd6

   use scale_atmos_dyn_fvm_flux_ud7, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud7, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud7, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud7, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud7, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud7, &
      ATMOS_DYN_FVM_fluxX_XYW_ud7, &
      ATMOS_DYN_FVM_fluxY_XYW_ud7, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud7, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud7, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud7, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud7, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud7, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud7, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud7, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud7, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud7, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud7, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud7, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud7, &
      ATMOS_DYN_FVM_fluxZJ_XYW_ud7, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_ud7, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_ud7

   use scale_atmos_dyn_fvm_flux_cd8, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd8, &
      ATMOS_DYN_FVM_fluxZ_XYZ_cd8, &
      ATMOS_DYN_FVM_fluxX_XYZ_cd8, &
      ATMOS_DYN_FVM_fluxY_XYZ_cd8, &
      ATMOS_DYN_FVM_fluxZ_XYW_cd8, &
      ATMOS_DYN_FVM_fluxX_XYW_cd8, &
      ATMOS_DYN_FVM_fluxY_XYW_cd8, &
      ATMOS_DYN_FVM_fluxJ13_XYW_cd8, &
      ATMOS_DYN_FVM_fluxJ23_XYW_cd8, &
      ATMOS_DYN_FVM_fluxZ_UYZ_cd8, &
      ATMOS_DYN_FVM_fluxX_UYZ_cd8, &
      ATMOS_DYN_FVM_fluxY_UYZ_cd8, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_cd8, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_cd8, &
      ATMOS_DYN_FVM_fluxZ_XVZ_cd8, &
      ATMOS_DYN_FVM_fluxX_XVZ_cd8, &
      ATMOS_DYN_FVM_fluxY_XVZ_cd8, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_cd8, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_cd8, &
      ATMOS_DYN_FVM_fluxZJ_XYW_cd8, &
      ATMOS_DYN_FVM_fluxZJ_UYZ_cd8, &
      ATMOS_DYN_FVM_fluxZJ_XVZ_cd8

    implicit none
    character(len=*), intent(in) :: scheme
    character(len=*), intent(in) :: scheme_tracer
    !---------------------------------------------------------------------------

    nullify( ATMOS_DYN_FVM_fluxZJ_XYW, ATMOS_DYN_FVM_fluxZJ_UYZ, ATMOS_DYN_FVM_fluxZJ_XVZ )

    select case( scheme )

    case( "UD1" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud1 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 1
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud1

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud1

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud1

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud1

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud1

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud1

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud1

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud1

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud1

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud1

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_ud1
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_ud1
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_ud1

      if ( ( .not. PRC_TwoD ) .and. IHALO < 1 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 1
         call PRC_abort
      end if
      if ( JHALO < 1 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 1
         call PRC_abort
      end if


    case( "CD2" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the cd2 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 2
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_cd2

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_cd2

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_cd2

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_cd2

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_cd2

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_cd2

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_cd2

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_cd2

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_cd2

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_cd2

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_cd2
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_cd2
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_cd2

      if ( ( .not. PRC_TwoD ) .and. IHALO < 1 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 1
         call PRC_abort
      end if
      if ( JHALO < 1 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 1
         call PRC_abort
      end if


    case( "UD3" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud3 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 3
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud3

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud3

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud3

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud3

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud3

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud3

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud3

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud3

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud3

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud3

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_ud3
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_ud3
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_ud3

      if ( ( .not. PRC_TwoD ) .and. IHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 2
         call PRC_abort
      end if
      if ( JHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 2
         call PRC_abort
      end if


    case( "UD3KOREN1993" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud3Koren1993 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 4
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud3Koren1993

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_ud3Koren1993
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_ud3Koren1993
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_ud3Koren1993

      if ( ( .not. PRC_TwoD ) .and. IHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 2
         call PRC_abort
      end if
      if ( JHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 2
         call PRC_abort
      end if


    case( "CD4" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the cd4 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 5
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_cd4

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_cd4

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_cd4

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_cd4

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_cd4

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_cd4

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_cd4

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_cd4

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_cd4

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_cd4

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_cd4
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_cd4
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_cd4

      if ( ( .not. PRC_TwoD ) .and. IHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 2
         call PRC_abort
      end if
      if ( JHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 2
         call PRC_abort
      end if


    case( "UD5" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud5 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 6
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud5

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud5

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud5

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud5

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud5

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud5

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud5

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud5

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud5

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud5

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_ud5
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_ud5
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_ud5

      if ( ( .not. PRC_TwoD ) .and. IHALO < 3 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 3
         call PRC_abort
      end if
      if ( JHALO < 3 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 3
         call PRC_abort
      end if


    case( "CD6" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the cd6 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 7
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_cd6

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_cd6

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_cd6

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_cd6

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_cd6

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_cd6

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_cd6

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_cd6

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_cd6

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_cd6

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_cd6
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_cd6
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_cd6

      if ( ( .not. PRC_TwoD ) .and. IHALO < 3 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 3
         call PRC_abort
      end if
      if ( JHALO < 3 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 3
         call PRC_abort
      end if


    case( "UD7" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud7 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 8
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud7

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud7

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud7

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud7

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud7

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud7

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud7

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud7

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud7

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud7

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_ud7
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_ud7
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_ud7

      if ( ( .not. PRC_TwoD ) .and. IHALO < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 4
         call PRC_abort
      end if
      if ( JHALO < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 4
         call PRC_abort
      end if

      if ( KMAX >= 3 .and. KMAX < 5 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'KMAX must be 1, 2, or >= ', 5, ' for UD7: ', KMAX
         call PRC_abort
      end if


    case( "CD8" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the cd8 scheme is used for flux calculation'

#ifdef _OPENACC
      i_scheme = 9
#else

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_cd8

#endif

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_cd8

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_cd8

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_cd8

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_cd8

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_cd8

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_cd8

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_cd8

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_cd8

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_cd8

      ATMOS_DYN_FVM_fluxZJ_XYW => ATMOS_DYN_FVM_fluxZJ_XYW_cd8
      ATMOS_DYN_FVM_fluxZJ_UYZ => ATMOS_DYN_FVM_fluxZJ_UYZ_cd8
      ATMOS_DYN_FVM_fluxZJ_XVZ => ATMOS_DYN_FVM_fluxZJ_XVZ_cd8

      if ( ( .not. PRC_TwoD ) .and. IHALO < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 4
         call PRC_abort
      end if
      if ( JHALO < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 4
         call PRC_abort
      end if

      if ( KMAX >= 3 .and. KMAX < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'KMAX must be 1, 2, or >= ', 4, ' for CD8: ', KMAX
         call PRC_abort
      end if


    case default
       LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'scheme is invalid: ', scheme
       call PRC_abort
    end select

    !$acc update device(i_scheme)

    if ( .not. ( associated( ATMOS_DYN_FVM_fluxZJ_XYW ) .and. &
                 associated( ATMOS_DYN_FVM_fluxZJ_UYZ ) .and. &
                 associated( ATMOS_DYN_FVM_fluxZJ_XVZ ) ) ) then
       LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'the vertical momentum fluxes are not set for ', trim(scheme)
       call PRC_abort
    end if

    select case( scheme_tracer )

    case( "UD1" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud1 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud1

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud1

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud1

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud1

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud1

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud1

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud1

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud1

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud1

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud1

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud1

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud1

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud1

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud1

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud1

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud1

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud1

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud1


      if ( ( .not. PRC_TwoD ) .and. IHALO < 1 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 1
         call PRC_abort
      end if
      if ( JHALO < 1 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 1
         call PRC_abort
      end if


    case( "CD2" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the cd2 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_cd2

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_cd2

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_cd2

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_cd2

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_cd2

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_cd2

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_cd2

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_cd2

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_cd2

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_cd2

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_cd2

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_cd2

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_cd2

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_cd2

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_cd2

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_cd2

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_cd2

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_cd2


      if ( ( .not. PRC_TwoD ) .and. IHALO < 1 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 1
         call PRC_abort
      end if
      if ( JHALO < 1 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 1
         call PRC_abort
      end if


    case( "UD3" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud3 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud3

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud3

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud3

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud3

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud3

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud3

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud3

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud3

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud3

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud3

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud3

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud3

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud3

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud3

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud3

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud3

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud3

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud3


      if ( ( .not. PRC_TwoD ) .and. IHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 2
         call PRC_abort
      end if
      if ( JHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 2
         call PRC_abort
      end if


    case( "UD3KOREN1993" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud3Koren1993 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud3Koren1993


      if ( ( .not. PRC_TwoD ) .and. IHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 2
         call PRC_abort
      end if
      if ( JHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 2
         call PRC_abort
      end if


    case( "CD4" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the cd4 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_cd4

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_cd4

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_cd4

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_cd4

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_cd4

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_cd4

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_cd4

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_cd4

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_cd4

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_cd4

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_cd4

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_cd4

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_cd4

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_cd4

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_cd4

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_cd4

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_cd4

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_cd4


      if ( ( .not. PRC_TwoD ) .and. IHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 2
         call PRC_abort
      end if
      if ( JHALO < 2 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 2
         call PRC_abort
      end if


    case( "UD5" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud5 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud5

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud5

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud5

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud5

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud5

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud5

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud5

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud5

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud5

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud5

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud5

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud5

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud5

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud5

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud5

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud5

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud5

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud5


      if ( ( .not. PRC_TwoD ) .and. IHALO < 3 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 3
         call PRC_abort
      end if
      if ( JHALO < 3 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 3
         call PRC_abort
      end if


    case( "CD6" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the cd6 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_cd6

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_cd6

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_cd6

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_cd6

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_cd6

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_cd6

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_cd6

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_cd6

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_cd6

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_cd6

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_cd6

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_cd6

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_cd6

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_cd6

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_cd6

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_cd6

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_cd6

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_cd6


      if ( ( .not. PRC_TwoD ) .and. IHALO < 3 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 3
         call PRC_abort
      end if
      if ( JHALO < 3 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 3
         call PRC_abort
      end if


    case( "UD7" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the ud7 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud7

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud7

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud7

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud7

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud7

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud7

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud7

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud7

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud7

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud7

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud7

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud7

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud7

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud7

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud7

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud7

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud7

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud7


      if ( ( .not. PRC_TwoD ) .and. IHALO < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 4
         call PRC_abort
      end if
      if ( JHALO < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 4
         call PRC_abort
      end if


    case( "CD8" )
      LOG_INFO("ATMOS_DYN_FVM_flux_setup",*) 'the cd8 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_cd8

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_cd8

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_cd8

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_cd8

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_cd8

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_cd8

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_cd8

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_cd8

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_cd8

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_cd8

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_cd8

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_cd8

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_cd8

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_cd8

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_cd8

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_cd8

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_cd8

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_cd8


      if ( ( .not. PRC_TwoD ) .and. IHALO < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'IHALO must be >= ', 4
         call PRC_abort
      end if
      if ( JHALO < 4 ) then
         LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'JHALO must be >= ', 4
         call PRC_abort
      end if


    case default
       LOG_ERROR("ATMOS_DYN_FVM_flux_setup",*) 'scheme is invalid: ', scheme_tracer
       call PRC_abort
    end select

  end subroutine ATMOS_DYN_FVM_flux_setup

#ifdef _OPENACC

  subroutine ATMOS_DYN_FVM_flux_valueW_Z( valW, mflx, val, GSQRT, CDZ )
    !$acc routine vector

   use scale_atmos_dyn_fvm_flux_ud1, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud1

   use scale_atmos_dyn_fvm_flux_cd2, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd2

   use scale_atmos_dyn_fvm_flux_ud3, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud3

   use scale_atmos_dyn_fvm_flux_ud3Koren1993, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud3Koren1993

   use scale_atmos_dyn_fvm_flux_cd4, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd4

   use scale_atmos_dyn_fvm_flux_ud5, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud5

   use scale_atmos_dyn_fvm_flux_cd6, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd6

   use scale_atmos_dyn_fvm_flux_ud7, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud7

   use scale_atmos_dyn_fvm_flux_cd8, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd8

    implicit none

    real(RP), intent(out) :: valW (KA)

    real(RP), intent(in)  :: mflx (KA)

    real(RP), intent(in)  :: val  (KA)

    real(RP), intent(in)  :: GSQRT(KA)

    real(RP), intent(in)  :: CDZ  (KA)


    select case ( i_scheme )

    case( 1 )
      call ATMOS_DYN_FVM_flux_valueW_Z_ud1( valW, mflx, val, GSQRT, CDZ )

    case( 2 )
      call ATMOS_DYN_FVM_flux_valueW_Z_cd2( valW, mflx, val, GSQRT, CDZ )

    case( 3 )
      call ATMOS_DYN_FVM_flux_valueW_Z_ud3( valW, mflx, val, GSQRT, CDZ )

    case( 4 )
      call ATMOS_DYN_FVM_flux_valueW_Z_ud3Koren1993( valW, mflx, val, GSQRT, CDZ )

    case( 5 )
      call ATMOS_DYN_FVM_flux_valueW_Z_cd4( valW, mflx, val, GSQRT, CDZ )

    case( 6 )
      call ATMOS_DYN_FVM_flux_valueW_Z_ud5( valW, mflx, val, GSQRT, CDZ )

    case( 7 )
      call ATMOS_DYN_FVM_flux_valueW_Z_cd6( valW, mflx, val, GSQRT, CDZ )

    case( 8 )
      call ATMOS_DYN_FVM_flux_valueW_Z_ud7( valW, mflx, val, GSQRT, CDZ )

    case( 9 )
      call ATMOS_DYN_FVM_flux_valueW_Z_cd8( valW, mflx, val, GSQRT, CDZ )

    end select

    return
  end subroutine ATMOS_DYN_FVM_flux_valueW_Z

#endif

end module scale_atmos_dyn_fvm_flux

!--
! vi:set readonly sw=4 ts=8
!
!Local Variables:
!mode: f90
!buffer-read-only: t
!End:
!
!++
