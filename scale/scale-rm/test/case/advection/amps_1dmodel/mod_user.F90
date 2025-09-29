!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_finalize
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update

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
  logical, private :: USER_do = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config before setup of tracers
  subroutine USER_tracer_setup
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    ! if you want to add tracers, call the TRACER_regist subroutine.
    ! e.g.,
!    integer, parameter     :: NQ = 1
!    integer                :: QS
!    character(len=H_SHORT) :: NAME(NQ)
!    character(len=H_MID)   :: DESC(NQ)
!    character(len=H_SHORT) :: UNIT(NQ)
!
!    data NAME (/ 'name' /)
!    data DESC (/ 'tracer name' /)
!    data UNIT (/ 'kg/kg' /)
    !---------------------------------------------------------------------------

!    call TRACER_regist( QS,   & ! [OUT]
!                        NQ,   & ! [IN]
!                        NAME, & ! [IN]
!                        DESC, & ! [IN]
!                        UNIT  ) ! [IN]

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'This module is dummy.'

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Finalization
  subroutine USER_finalize
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_finalize

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    use scale_atmos_phy_mp_amps, only: &
       ATMOS_PHY_MP_amps_init_qtrc_1DMODEL
    use scale_const, only: &
       CONST_GRAV
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    implicit none
    !---------------------------------------------------------------------------

    real(RP) :: RHO
    real(RP) :: VELX
    real(RP) :: VELY
    real(RP) :: POTT

    !real(RP) :: QHYD(KA,IA,JA,N_HYD)
    !real(RP) :: QNUM(KA,IA,JA,N_HYD)

    real(RP) :: MK1DMODEL_rho = 0.0D0
    real(RP) :: MK1DMODEL_v = 0.0D0
    real(RP) :: MK1DMODEL_pot = 0.0D0
    real(RP) :: MK1DMODEL_qv = 0.0D0
    real(RP) :: MK1DMODEL_qc = 0.0D0
    real(RP) :: MK1DMODEL_ni_ylimit = 0.0D0
    integer  :: MK1DMODEL_case = 1
    namelist / PARAM_USER_MKINIT / &
       MK1DMODEL_rho, &
       MK1DMODEL_v, &
       MK1DMODEL_pot, &
       MK1DMODEL_qv, &
       MK1DMODEL_qc, &
       MK1DMODEL_ni_ylimit, &
       MK1DMODEL_case

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! REMEMBER TO SWITCH OFF AMPS MICROPHYSICS CALCULATION

    LOG_NEWLINE
    LOG_INFO("MKINIT_1DMODEL",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_1DMODEL",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER_MKINIT,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_1DMODEL",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_1DMODEL",*) 'Not appropriate names in namelist PARAM_MKINIT_1DMODEL. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER_MKINIT)

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = MK1DMODEL_rho
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = MK1DMODEL_rho * MK1DMODEL_v

       RHOT(k,i,j) = MK1DMODEL_rho * ( MK1DMODEL_pot )

       QTRC(k,i,j,1) = MK1DMODEL_qv

    enddo
    enddo
    enddo

    CONST_GRAV = 0.0_RP

    call ATMOS_PHY_MP_amps_init_qtrc_1DMODEL( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                              DENS(:,:,:), &
                                              QTRC(:,:,:,QS_MP+1:QE_MP), &
                                              MK1DMODEL_ni_ylimit, &
                                              MK1DMODEL_case )

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calculation tendency
  subroutine USER_calc_tendency
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_update
    use scale_prc, only: &
       PRC_abort
    implicit none
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       call PRC_abort
    endif

    return
  end subroutine USER_update

end module mod_user
