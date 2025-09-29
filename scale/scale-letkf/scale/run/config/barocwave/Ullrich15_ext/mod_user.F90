!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Set boundary conditions for baroclinic wave in a channel based on Ullrich et al. (2015).
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

  use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       PRES

  use scale_prc

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
  public :: USER_resume

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,  public, save :: USER_resume_do = .false. !< do user step?
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private, save :: USER_do = .false. !< do user step?

  character(len=H_LONG),private, save :: REFSTATE_IN_BASENAME=''
  real(RP), private, save :: NUDGE_UV_DAY_TOP = -999.9D0
  real(RP), private, save :: NUDGE_UV_DAY_BOT = -999.9D0
  real(RP), private, save :: NUDGE_UV_TRANS_HEIGHT = 3.0D3
  real(RP), private, save :: NUDGE_UV_TRANS_WIDTH  = 1.0D3
  real(RP), private, save :: NUDGE_T_DAY_TOP = -999.9D0
  real(RP), private, save :: NUDGE_T_DAY_BOT = -999.9D0
  real(RP), private, save :: NUDGE_T_TRANS_HEIGHT = 3.0D3
  real(RP), private, save :: NUDGE_T_TRANS_WIDTH  = 1.0D3

  real(RP), private, allocatable, save :: ALPHA_NUDGE_UV(:,:,:)
  real(RP), private, allocatable, save :: ALPHA_NUDGE_T(:,:,:)

  real(RP), private, allocatable, save :: MOMX_init(:,:,:)   ! momentum x  [kg/m2/s]
  real(RP), private, allocatable, save :: MOMY_init(:,:,:)   ! momentum y  [kg/m2/s]
  real(RP), private, allocatable, save :: RHOT_init(:,:,:)   ! DENS * POTT [K*kg/m3]
 
  logical, private, save :: flag_init = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine USER_tracer_setup

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ
    use scale_atmos_grid_cartesC_index, only: &
       KA, IA, JA 
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_close, &
       FILE_CARTESC_read 
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       REFSTATE_IN_BASENAME, &
       NUDGE_UV_DAY_TOP, & 
       NUDGE_UV_DAY_BOT, & 
       NUDGE_UV_TRANS_HEIGHT, & 
       NUDGE_UV_TRANS_WIDTH, & 
       NUDGE_T_DAY_TOP, & 
       NUDGE_T_DAY_BOT, & 
       NUDGE_T_TRANS_HEIGHT, & 
       NUDGE_T_TRANS_WIDTH 

    integer :: ierr
    integer :: i, j, k
    integer :: refstate_fid
    real(RP) :: alpha_uv_top=0.0D0, alpha_uv_bot=0.0D0, alpha_t_top=0.0D0, alpha_t_bot=0.0D0
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'
    LOG_INFO("USER_setup",*) 'User procedure in test/case/barocwave/Ullrich15'

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

    allocate( ALPHA_NUDGE_UV(KA,IA,JA) )
    allocate( ALPHA_NUDGE_T(KA,IA,JA)  )
    allocate( MOMX_init(KA,IA,JA)      )
    allocate( MOMY_init(KA,IA,JA)      )
    allocate( RHOT_init(KA,IA,JA)      )

    if (NUDGE_UV_DAY_TOP.gt.0.0) alpha_uv_top = 1.0 / ( NUDGE_UV_DAY_TOP * 86400.0 )
    if (NUDGE_UV_DAY_BOT.gt.0.0) alpha_uv_bot = 1.0 / ( NUDGE_UV_DAY_BOT * 86400.0 )
    if (NUDGE_T_DAY_TOP .gt.0.0) alpha_t_top  = 1.0 / ( NUDGE_T_DAY_TOP  * 86400.0 )
    if (NUDGE_T_DAY_BOT .gt.0.0) alpha_t_bot  = 1.0 / ( NUDGE_T_DAY_BOT  * 86400.0 )
  
    do k = 1, KA; do i = 1, IA ; do j = 1, JA
      ALPHA_NUDGE_UV(k, i, j) =  0.5 * &
         (   ( alpha_uv_top - alpha_uv_bot ) * tanh( ( REAL_CZ(k, i, j) - NUDGE_UV_TRANS_HEIGHT ) / NUDGE_UV_TRANS_WIDTH ) & 
           + ( alpha_uv_top + alpha_uv_bot ) )
      ALPHA_NUDGE_T(k, i, j)  =  0.5 * &
         (   ( alpha_t_top  - alpha_t_bot  ) * tanh( ( REAL_CZ(k, i, j) - NUDGE_T_TRANS_HEIGHT  ) / NUDGE_T_TRANS_WIDTH  ) & 
           + ( alpha_t_top  + alpha_t_bot  ) )
    end do; end do; end do

    if (trim(REFSTATE_IN_BASENAME) /= '') then
      flag_init = .true.
      call FILE_CARTESC_open( REFSTATE_IN_BASENAME, refstate_fid )
      call FILE_CARTESC_read( refstate_fid, 'MOMX', 'ZXHY', & ! [IN]
                               MOMX_init(:,:,:)                                ) ! [OUT]
      call FILE_CARTESC_read( refstate_fid, 'MOMY', 'ZXYH', & ! [IN]
                               MOMY_init(:,:,:)                                ) ! [OUT]
      call FILE_CARTESC_read( refstate_fid, 'RHOT', 'ZXY',  & ! [IN]
                               RHOT_init(:,:,:)                                ) ! [OUT]
      call FILE_CARTESC_close( refstate_fid )
    end if    

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
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Dummy
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine USER_calc_tendency
    use scale_atmos_grid_cartesC_index, only: &
       KA, IA, JA
    use mod_atmos_vars, only: &
       RHOU_tp, & 
       RHOV_tp, & 
       RHOT_tp  
     implicit none

     integer :: k, i, j
    !---------------------------------------------------------------------------

    if (flag_init) then
      do k = 1, KA; do i = 1, IA ; do j = 1, JA
        RHOU_tp(k, i, j) = RHOU_tp(k, i, j) - ALPHA_NUDGE_UV(k, i, j) * ( MOMX(k, i, j) - MOMX_init(k, i, j) ) 
        RHOV_tp(k, i, j) = RHOV_tp(k, i, j) - ALPHA_NUDGE_UV(k, i, j) * ( MOMY(k, i, j) - MOMY_init(k, i, j) ) 
        RHOT_tp(k, i, j) = RHOT_tp(k, i, j) - ALPHA_NUDGE_T(k, i, j)  * ( RHOT(k, i, j) - RHOT_init(k, i, j) ) 
      end do; end do; end do
    else
      flag_init=.true.
      MOMX_init = MOMX
      MOMY_init = MOMY
      RHOT_init = RHOT   
    end if

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_update
    use scale_prc_cartesC, only: &
       PRC_HAS_N, &
       PRC_HAS_S
    implicit none

    integer :: k, i, j

    !---------------------------------------------------------------------------

    if ( .not. USER_do ) return

    ! Apply the boundary condition at y=+Ly and y=-Ly

    if ( .NOT. PRC_HAS_N ) then
       MOMY(:,:,JE)   = 0.0_RP
       do j = 1, JHALO
          MOMY(:,:,JE+j) = - MOMY(:,:,JE-j  )
          DENS(:,:,JE+j) = + DENS(:,:,JE-j+1)
          MOMX(:,:,JE+j) = + MOMX(:,:,JE-j+1)
          MOMZ(:,:,JE+j) = + MOMZ(:,:,JE-j+1)
          RHOT(:,:,JE+j) = + RHOT(:,:,JE-j+1)
       enddo
    end if

    if ( .NOT. PRC_HAS_S ) then
       MOMY(:,:,JS-1) = 0.0_RP
       do j = 1, JHALO
          if ( j < JHALO ) MOMY(:,:,JS-j-1) = - MOMY(:,:,JS+j-1)
          DENS(:,:,JS-j) = + DENS(:,:,JS+j-1)
          MOMX(:,:,JS-j) = + MOMX(:,:,JS+j-1)
          MOMZ(:,:,JS-j) = + MOMZ(:,:,JS+j-1)
          RHOT(:,:,JS-j) = + RHOT(:,:,JS+j-1)
       enddo
    end if

    return
  end subroutine USER_update

end module mod_user
