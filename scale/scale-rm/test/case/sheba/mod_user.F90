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

  character(len=H_SHORT) :: USER_experiment

  ! *********************************************************************
  ! -- these are defined for SHEBA and MPACE run
  ! *********************************************************************
  logical :: SHEBA_SWITCH_ACCE = .false.
  logical :: SHEBA_SWITCH_MOMZ = .false.
  logical :: SHEBA_SWITCH_RHOU = .false.
  logical :: SHEBA_SWITCH_RHOV = .false.
  logical :: SHEBA_SWITCH_DENS = .false.
  logical :: SHEBA_SWITCH_QVAP = .false.
  logical :: SHEBA_SWITCH_RHOT = .false.
  logical :: SHEBA_SWITCH_TEMP = .false.
  logical :: SHEBA_SWITCH_RHOQ = .false.

  real(RP), allocatable :: largeScaleTTendency(:) ! large-scale temperature forcing
  real(RP), allocatable :: largeScaleQTendency(:) ! large-scale vapor forcing
  real(RP), allocatable :: WLS(:) ! large-scale sinking
  real(RP) :: sfc_largeScaleTTendency, sfc_largeScaleQTendency, sfc_wls

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
    use scale_const, only: &
       PI => CONST_PI
    use scale_atmos_grid_cartesC, only: &
       CZ  => ATMOS_GRID_CARTESC_CZ, &
       FZ  => ATMOS_GRID_CARTESC_FZ
    implicit none

    integer, parameter :: EXP_klim = 501 ! there was a bug: EXP_klim = 100 (3/4/2020), after
                                         ! SHEBA has been run
    integer            :: EXP_kmax

    logical  :: USER_const = .true.

    real(RP) :: EXP_z   (EXP_klim+1) ! height      [m]

    real(RP) :: SFC_TTND             ! surface large-scale temperature tendency [K/s]
    real(RP) :: EXP_ttnd(EXP_klim+1) ! large-scale temperature tendency [K/s]
    real(RP) :: SFC_QTND             ! surface large-scale vapor tendency [kg/kg/s]
    real(RP) :: EXP_qtnd(EXP_klim+1) ! large-scale vapor tendency [kg/kg/s]
    real(RP) :: SFC_WSIK             ! surface large-scale sinking [m/s]
    real(RP) :: EXP_Wsik(EXP_klim+1) ! large-scale sinking [m/s]

    character(len=H_LONG) :: USER_file = ''

    real(RP) :: T_tend = 0.0_RP
    real(RP) :: Q_tend = 0.0_RP
    real(RP) :: W_sink = 0.0_RP

    real(RP) :: fact1, fact2, ttnd(KA), qtnd(KA), wlse(KA)

    integer :: k, kref
    integer :: fid
    integer :: ierr

    namelist / PARAM_USER / &
       USER_do, &
       USER_file, &
       USER_const, &
       USER_experiment, &
       SHEBA_SWITCH_ACCE, &
       SHEBA_SWITCH_MOMZ, &
       SHEBA_SWITCH_RHOU, &
       SHEBA_SWITCH_RHOV, &
       SHEBA_SWITCH_DENS, &
       SHEBA_SWITCH_QVAP, &
       SHEBA_SWITCH_RHOT, &
       SHEBA_SWITCH_TEMP, &
       SHEBA_SWITCH_RHOQ

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

    ! initialization of local arrays
    allocate(largeScaleTTendency(KA),largeScaleQTendency(KA),WLS(KA))
    largeScaleTTendency(:) = 0.0_RP
    largeScaleQTendency(:) = 0.0_RP
    WLS(:) = 0.0_RP

    ! temporarily set to zero
    if ( USER_const .eqv. .true. ) then

       largeScaleTTendency(KS:KE) = T_tend
       largeScaleQTendency(KS:KE) = Q_tend
       WLS(KS:KE) = W_sink

       sfc_largeScaleTTendency = largeScaleTTendency(KS)
       sfc_largeScaleQTendency = largeScaleQTendency(KS)
       sfc_wls = WLS(KS)

       WLS(KE+1) = WLS(KE)
       WLS(KS-1) = WLS(KS)

    else

       fid = IO_get_available_fid()
       open( fid,                                 &
             file   = trim(USER_file), &
             form   = 'formatted',                &
             status = 'old',                      &
             iostat = ierr                        )

       if ( ierr /= 0 ) then
          LOG_ERROR("read_largescale_sheba",*) '[user_setup/read_largescale] Input file not found!'
          call PRC_abort
       endif

       !--- read sounding file till end
       read(fid,*) SFC_TTND, SFC_QTND, SFC_WSIK

       LOG_INFO("USER_setup",*) '+ Surface large-scale temperature tendency [K/s]', SFC_TTND

       sfc_largeScaleTTendency = SFC_TTND
       sfc_largeScaleQTendency = SFC_QTND
       sfc_wls = SFC_WSIK

       do k = 2, EXP_klim
          read(fid,*,iostat=ierr) EXP_z(k), EXP_ttnd(k), EXP_qtnd(k), EXP_wsik(k)
          if ( ierr /= 0 ) exit
       enddo

       EXP_kmax = k - 1
       close(fid)

       ! Boundary
       EXP_z   (1)          = 0.0_RP
       EXP_ttnd(1)          = SFC_TTND
       EXP_qtnd(1)          = SFC_QTND
       EXP_wsik(1)          = SFC_WSIK
       EXP_z   (EXP_kmax+1) = 100.E3_RP
       EXP_ttnd(EXP_kmax+1) = EXP_ttnd(EXP_kmax)
       EXP_qtnd(EXP_kmax+1) = EXP_qtnd(EXP_kmax)
       EXP_wsik(EXP_kmax+1) = EXP_wsik(EXP_kmax)

       !--- linear interpolate to model grid
       do k = KS, KE
          do kref = 2, EXP_kmax+1
             if (       CZ(k) >  EXP_z(kref-1) &
                  .AND. CZ(k) <= EXP_z(kref  ) ) then

                fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
                fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

                ttnd(k) = EXP_ttnd(kref-1) * fact1 &
                        + EXP_ttnd(kref  ) * fact2
                qtnd(k) = EXP_qtnd(kref-1) * fact1 &
                        + EXP_qtnd(kref  ) * fact2
                wlse(k) = EXP_wsik(kref-1) * fact1 &
                        + EXP_wsik(kref  ) * fact2

             endif
          enddo
       enddo

       do k = KS, KE
          largeScaleTTendency(k) = ttnd(k)
          largeScaleQTendency(k) = qtnd(k)
          WLS(k) = wlse(k)
       enddo

       WLS(KE+1) = WLS(KE)
       WLS(KS-1) = sfc_wls
    endif

    LOG_NEWLINE
    LOG_INFO("USER_setup",'(1x,A)') 'Large-scale advective tendency of temperature, water vapor, and subsidence'
    LOG_INFO_CONT('(1x,A)') '====================================================='
    LOG_INFO_CONT('(1x,A)') '      GRID CENTER         Tadv            Qadv             W'
    do k = KS-1, KE
       LOG_INFO_CONT('(1x,A,ES15.5,A,ES15.5,A,ES15.5,A,ES15.5)') '    ', CZ(k), '   ', largeScaleTTendency(k), '   ', largeScaleQTendency(k), '   ', WLS(k)
    enddo
    LOG_INFO_CONT('(1x,A)') '====================================================='

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Finalization
  subroutine USER_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_finalize",*) 'Finalize'

    if (allocated(largeScaleTTendency)) deallocate(largeScaleTTendency)
    if (allocated(largeScaleQTendency)) deallocate(largeScaleQTendency)
    if (allocated(WLS)) deallocate(WLS)
    if (allocated(largeScaleTTendency)) deallocate(largeScaleTTendency)

    return
  end subroutine USER_finalize

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       N_HYD, &
       I_HC
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
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    implicit none
    !---------------------------------------------------------------------------

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV1D(KA)
    real(RP) :: QCI1D(KA)
    real(RP) :: QNUM1D(KA)

    real(RP) :: qv(KA,IA,JA)
    real(RP) :: QHYD(KA,IA,JA,N_HYD)
    real(RP) :: QNUM(KA,IA,JA,N_HYD)

    real(RP) :: bubbles(KA,JA), temp

    integer  :: SHEBA_fluctuationNumberLayers = 10
    real(RP) :: SHEBA_fluctuation = 0.5D0
    namelist / PARAM_USER_MKINIT / &
       SHEBA_fluctuation, &
       SHEBA_fluctuationNumberLayers

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_SHEBA",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_SHEBA",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER_MKINIT,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_SHEBA",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_SHEBA",*) 'Not appropriate names in namelist PARAM_MKINIT_SHEBA. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER_MKINIT)

    call read_sounding_sheba( RHO, VELX, VELY, POTT, QV1D, QCI1D, QNUM1D ) ! (out)

    !$acc data copyin(RHO,VELX,VELY,POTT,QV1D,QCI1D,QNUM1D) &
    !$acc      copyout(DENS,MOMZ,MOMX,MOMY,RHOT,QTRC(:,:,:,QS_MP:QE_MP)) &
    !$acc      create(bubbles,QHYD,QNUM)

    ! initiate small fluctuation in potential temperature
    !$omp parallel do private(temp)
    !$acc kernels
    do j = JSB, JEB
    do k = KS, KS+SHEBA_fluctuationNumberLayers-1
       call random_number(temp)
       bubbles(k,j) = SHEBA_fluctuation * temp
    enddo
    enddo
    !$acc end kernels
    !$omp parallel do
    !$acc kernels
    do j = JSB, JEB
    do k = KS+SHEBA_fluctuationNumberLayers, KE
       bubbles(k,j) = 0.0_RP
    enddo
    enddo
    !$acc end kernels

    !$omp workshare
    !$acc kernels
    QHYD(:,:,:,:) = 0.0_RP
    QNUM(:,:,:,:) = 0.0_RP
    !$acc end kernels
    !$omp end workshare


    !$omp parallel do
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = RHO(k) * bubbles(k,j)
       ! rotated 90 deg
       MOMX(k,i,j) = RHO(k) * (-VELY(k))
       MOMY(k,i,j) = RHO(k) *   VELX(k)

       RHOT(k,i,j) = RHO(k) * ( POTT(k) )

       qv  (k,i,j) = QV1D(k)

       QHYD(k,i,j,I_HC) = QCI1D(k) ! only the smallest (bin) cloud ice are initialized
       if (QNUM1D(k) /= 0.0_RP) then
          QNUM(k,i,j,I_HC) = QNUM1D(k)*1000000.0_RP
       else
          ! number concentration (cm-3->m-3) of smallest droplets
          ! assuming that each droplet mass is 1 x 10^-14 g (smallest bin)
          QNUM(k,i,j,I_HC) = QCI1D(k)*RHO(k)/1.E-17_RP*1.E-6_RP *1000000.0_RP
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                        qv(:,:,:), QHYD(:,:,:,:), & ! [IN]
                                        DENS(:,:,:),              & ! [IN]
                                        QTRC(:,:,:,QS_MP:QE_MP),  & ! [OUT]
                                        QNUM=QNUM(:,:,:,:)        ) ! [IN]

    !$acc end data

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calculation tendency
  subroutine USER_calc_tendency
    use scale_const, only: &
       CONST_GRAV
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_atmos_vars, only: &
       DENS,  &
       MOMZ   => MOMZ_av, &
       U, &
       V, &
       W, &
       QTRC,  &
       PRES,  &
       TEMP,  &
       POTT,  &
       RHOT,  &
       DENS_t => DENS_tp, & ! cell center
       RHOU_t => RHOU_tp, & ! cell center
       RHOV_t => RHOV_tp, & ! cell center
       MOMZ_t => MOMZ_tp, & ! cell face
       RHOT_t => RHOT_tp, & ! cell center
       RHOQ_t => RHOQ_tp    ! cell center
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none
    !---------------------------------------------------------------------------

    real(RP) :: MOMZ_t_USER(KA,IA,JA)
    real(RP) :: RHOU_t_USER(KA,IA,JA)
    real(RP) :: RHOV_t_USER(KA,IA,JA)
    real(RP) :: DENS_t_USER(KA,IA,JA)
    real(RP) :: RHOT_t_USER(KA,IA,JA)
    real(RP) :: RHOH_t_USER(KA,IA,JA)
    real(RP) :: TEMP_t_USER(KA,IA,JA)
    real(RP) :: RHOQ_t_USER(KA,IA,JA,QS_MP:QE_MP)

    real(RP) :: subsidence_sink(KA,IA,JA)

    real(RP) :: DENS_column(KA), RHOT_column(KA), U_column(KA), V_column(KA), W_column(KA), POTT_column(KA)
    real(RP) :: QTRC_column(KA,QA), MOMZ_column(KA), TEMP_column(KA)
    real(RP) :: FZ(KA), FDZ(KA), RFDZ(KA), RCDZ(KA)

    real(RP) :: SINK_DUP, SINK_UP, SINK_CEN


    integer  :: k, i, j, iq

    if ( .not. USER_do ) then
       return
    endif

    MOMZ_t_USER(:,:,:) = 0.0_RP
    RHOU_t_USER(:,:,:) = 0.0_RP
    RHOV_t_USER(:,:,:) = 0.0_RP
    DENS_t_USER(:,:,:) = 0.0_RP
    RHOT_t_USER(:,:,:) = 0.0_RP
    RHOH_t_USER(:,:,:) = 0.0_RP
    TEMP_t_USER(:,:,:) = 0.0_RP

    RHOQ_t_USER(:,:,:,:) = 0.0_RP

    subsidence_sink = 0.0_RP

    !$omp parallel do &
    !$omp private(FZ,FDZ,RFDZ,RCDZ,DENS_column,TEMP_column,POTT_column,U_column,V_column,W_column,RHOT_column,MOMZ_column,QTRC_column, &
    !$omp         SINK_CEN,SINK_UP,SINK_DUP)
    do j = JS, JE
    do i = IS, IE

       ! define grid interval
       FZ(1:KA) = REAL_FZ(1:KA,i,j)
       FDZ(KS-1) = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j)
       RFDZ(KS-1) = 1.0_RP / FDZ(KS-1)
       do k = KS, KE
          FDZ(k) = REAL_CZ(k+1,i,j) - REAL_CZ(k  ,i,j)
          RFDZ(k) = 1.0_RP / FDZ(k)
          RCDZ(k) = 1.0_RP / ( REAL_FZ(k  ,i,j) - REAL_FZ(k-1,i,j) )
       enddo

       ! define columnar prognotic variables with boundary condition (upper and lower)
       do k = KS, KE
          DENS_column(k) = DENS(k,i,j)
          TEMP_column(k) = TEMP(k,i,j)
          POTT_column(k) = POTT(k,i,j)
          U_column(k) = U(k,i,j)
          V_column(k) = V(k,i,j)
          W_column(k) = W(k,i,j)
          RHOT_column(k) = RHOT(k,i,j)
          MOMZ_column(k) = MOMZ(k,i,j)
          do iq = 1, QA
             QTRC_column(k,iq) = QTRC(k,i,j,iq)
          enddo
       enddo
       ! the upper boundary condition by linear extrapolation
       DENS_column(KE+1) = 2.0_RP*DENS_column(KE) - DENS_column(KE-1)
       TEMP_column(KE+1) = 2.0_RP*TEMP_column(KE) - TEMP_column(KE-1)
       POTT_column(KE+1) = 2.0_RP*POTT_column(KE) - POTT_column(KE-1)
       RHOT_column(KE+1) = 2.0_RP*RHOT_column(KE) - RHOT_column(KE-1)
       U_column(KE+1)    = 2.0_RP*U_column(KE)    - U_column(KE-1)
       V_column(KE+1)    = 2.0_RP*V_column(KE)    - V_column(KE-1)
       W_column(KE+1)    = -W_column(KE)
       MOMZ_column(KE:KE+1) = 0.0_RP
       do iq = 1, QA
          QTRC_column(KE+1,iq) = 2.0_RP*QTRC_column(KE,iq) - QTRC_column(KE-1,iq)
       enddo
       ! the lower boundary condition assumes constant extrapolation
       DENS_column(KS-1) = DENS_column(KS)
       TEMP_column(KS-1) = TEMP_column(KS)
       POTT_column(KS-1) = POTT_column(KS)
       U_column(KS-1)    = U_column(KS)
       V_column(KS-1)    = V_column(KS)
       W_column(KS-1)    = W_column(KS)
       RHOT_column(KS-1) = RHOT_column(KS)
       MOMZ_column(KS-1) = MOMZ_column(KS)
       MOMZ_column(KS-2) = MOMZ_column(KS) ! z momentum needs two halo values
       do iq = 1, QA
          QTRC_column(KS-1,iq) = QTRC_column(KS,iq)
       enddo

       if ( trim(USER_experiment) == 'SHEBA' ) then

          ! compute large-scale forcing
          do k = KS, KE

             ! large-scale cooling forcing
             TEMP_t_USER(k,i,j) = TEMP_t_USER(k,i,j) + largeScaleTTendency(k) * RHOT(k,i,j) / TEMP(k,i,j)

             ! large-scale vapor forcing
             RHOQ_t_USER(k,i,j,QS_MP) = RHOQ_t_USER(k,i,j,QS_MP) + largeScaleQTendency(k) * DENS(k,i,j)

             ! large-scale vapor forcing
             DENS_t_USER(k,i,j) = DENS_t_USER(k,i,j) + largeScaleQTendency(k) * DENS(k,i,j)

             ! large-scale sinking forcing, the simple first-order upwind advection scheme is used
             if ( PRES(k,i,j) >= 95700.0_RP ) then
                SINK_CEN = -(-0.000008233_RP*PRES(k,i,j) + 0.8379_RP) / &
                                 (CONST_GRAV*DENS(k,i,j))
             else if ( PRES(k,i,j) >= 60000.0_RP ) then
                SINK_CEN = (-0.05_RP) / (CONST_GRAV*DENS(k,i,j))
             else
                SINK_CEN = 0.0_RP
             endif
             if ( k == KE ) then
                if ( PRES(k-1,i,j) >= 95700.0_RP ) then
                   SINK_UP = -(-0.000008233_RP*PRES(k-1,i,j) + 0.8379_RP) / &
                                   (CONST_GRAV*DENS(k-1,i,j))
                else if ( PRES(k-1,i,j) >= 60000.0_RP ) then
                   SINK_UP = (-0.05_RP) / (CONST_GRAV*DENS(k-1,i,j))
                else
                   SINK_UP = 0.0_RP
                endif
                SINK_UP = 2.0_RP*SINK_CEN - SINK_UP
             else
                if ( PRES(k+1,i,j) >= 95700.0_RP ) then
                   SINK_UP = -(-0.000008233_RP*PRES(k+1,i,j) + 0.8379_RP) / &
                                   (CONST_GRAV*DENS(k+1,i,j))
                else if ( PRES(k+1,i,j) >= 60000.0_RP ) then
                   SINK_UP = (-0.05_RP) / (CONST_GRAV*DENS(k+1,i,j))
                else
                   SINK_UP = 0.0_RP
                endif
             endif
             if ( k == KS ) then
                SINK_DUP = -SINK_CEN
             else
                if ( PRES(k-1,i,j) >= 95700.0_RP ) then
                   SINK_DUP = -(-0.000008233_RP*PRES(k-1,i,j) + 0.8379_RP) / &
                                   (CONST_GRAV*DENS(k-1,i,j))
                else if ( PRES(k-1,i,j) >= 60000.0_RP ) then
                   SINK_DUP = (-0.05_RP) / (CONST_GRAV*DENS(k-1,i,j))
                else
                   SINK_DUP = 0.0_RP
                endif
             endif

             subsidence_sink(k,i,j) = SINK_CEN

             ! -- x momentum --
             RHOU_t_USER(k,i,j) = RHOU_t_USER(k,i,j) - &
                                  0.5_RP*(SINK_UP + SINK_CEN)* &
                                  !0.5_RP*(DENS_column(k+1)*SINK_UP + &
                                  !        DENS_column(k  )*SINK_CEN)* &
                                  (U_column(k+1) - U_column(k))*RFDZ(k)

             ! -- y momentum --
             RHOV_t_USER(k,i,j) = RHOV_t_USER(k,i,j) - &
                                  0.5_RP*(SINK_UP + SINK_CEN)* &
                                  !0.5_RP*(DENS_column(k+1)*SINK_UP + &
                                  !        DENS_column(k  )*SINK_CEN)* &
                                  (V_column(k+1) - V_column(k))*RFDZ(k)

             ! -- mixing ratio --
             if ( SHEBA_SWITCH_QVAP ) then
                iq = QS_MP
                RHOQ_t_USER(k,i,j,iq) = RHOQ_t_USER(k,i,j,iq) &
                                      - 0.5_RP*(SINK_UP + SINK_CEN) &
                                      * (QTRC_column(k+1,iq) - QTRC_column(k,iq))*RFDZ(k)
             else
                do iq = QS_MP, QE_MP
                   RHOQ_t_USER(k,i,j,iq) = RHOQ_t_USER(k,i,j,iq) &
                                         - 0.5_RP*(SINK_UP + SINK_CEN) &
                                         * (QTRC_column(k+1,iq) - QTRC_column(k,iq))*RFDZ(k)
                enddo
             endif

             ! -- energy --
             RHOT_t_USER(k,i,j) = RHOT_t_USER(k,i,j) - & ! DENS_column(k)* &
                                  0.5_RP*(SINK_UP + SINK_CEN)* &
                                  (POTT_column(k+1) - POTT_column(k))*RFDZ(k)

          enddo

          ! large-scale sinking forcing, the simple first-order upwind advection scheme is used
          do k = KS, KE-1
             if ( PRES(k,i,j) >= 95700.0_RP ) then
                SINK_CEN = -(-0.000008233_RP*PRES(k,i,j) + 0.8379_RP) / &
                                 (CONST_GRAV*DENS(k,i,j))
             else if ( PRES(k,i,j) >= 60000.0_RP ) then
                SINK_CEN = (-0.05_RP) / (CONST_GRAV*DENS(k,i,j))
             else
                SINK_CEN = 0.0_RP
             endif
             if ( PRES(k+1,i,j) >= 95700.0_RP ) then
                SINK_UP = -(-0.000008233_RP*PRES(k+1,i,j) + 0.8379_RP) / &
                                (CONST_GRAV*DENS(k+1,i,j))
             else if ( PRES(k+1,i,j) >= 60000.0_RP ) then
                SINK_UP = (-0.05_RP) / (CONST_GRAV*DENS(k+1,i,j))
             else
                SINK_UP = 0.0_RP
             endif
             if ( k == KE-1 ) then
                SINK_DUP = 2.0_RP*SINK_UP - SINK_CEN
             else
                if ( PRES(k+2,i,j) >= 95700.0_RP ) then
                   SINK_DUP = -(-0.000008233_RP*PRES(k+2,i,j) + 0.8379_RP) / &
                                    (CONST_GRAV*DENS(k+2,i,j))
                else if ( PRES(k+2,i,j) >= 60000.0_RP ) then
                   SINK_DUP = (-0.05_RP) / (CONST_GRAV*DENS(k+2,i,j))
                else
                   SINK_DUP = 0.0_RP
                endif
             endif

             ! -- z momentum --
             if (SHEBA_SWITCH_ACCE) then
                MOMZ_t_USER(k,i,j) = MOMZ_t_USER(k,i,j) &
                                   + 0.5_RP*(SINK_CEN + SINK_UP) &
                                   * 0.5_RP*(DENS_column(k) + DENS_column(k+1))
             else
                MOMZ_t_USER(k,i,j) = MOMZ_t_USER(k,i,j) &
                                   - SINK_UP &
                                   !* DENS_column(k+1)*SINK_UP &
                                   * ( MOMZ_column(k+1)*2.0_RP/( DENS_column(k+2) &
                                                             + DENS_column(k+1)) &
                                     - MOMZ_column(k  )*2.0_RP/( DENS_column(k+1) &
                                                               + DENS_column(k)))*RCDZ(k+1)
             endif
          enddo

       else if ( trim(USER_experiment) == 'MPACE' .or. trim(USER_experiment) == 'ISDAC' ) then

          ! compute large-scale forcing
          do k = KS, KE

             ! large-scale cooling forcing
             TEMP_t_USER(k,i,j) = TEMP_t_USER(k,i,j) + largeScaleTTendency(k) &
                                                     * DENS(k,i,j) * POTT(k,i,j) / TEMP(k,i,j)

             ! large-scale vapor forcing
             RHOQ_t_USER(k,i,j,QS_MP) = RHOQ_t_USER(k,i,j,QS_MP) + largeScaleQTendency(k) * DENS(k,i,j)

             ! large-scale vapor forcing
             DENS_t_USER(k,i,j) = DENS_t_USER(k,i,j) + largeScaleQTendency(k) * DENS(k,i,j)

             SINK_CEN = WLS(k)
             SINK_UP = WLS(k+1)
             SINK_DUP = WLS(k-1)

             subsidence_sink(k,i,j) = SINK_CEN

             ! -- x momentum --
             RHOU_t_USER(k,i,j) = RHOU_t_USER(k,i,j) &
                                - 0.5_RP*(SINK_UP + SINK_CEN) &
                                * (U_column(k+1) - U_column(k))*RFDZ(k)

             ! -- y momentum --
             RHOV_t_USER(k,i,j) = RHOV_t_USER(k,i,j) &
                                - 0.5_RP*(SINK_UP + SINK_CEN) &
                                * (V_column(k+1) - V_column(k))*RFDZ(k)

             ! -- mixing ratio --
             if ( SHEBA_SWITCH_QVAP ) then
                iq = QS_MP
                RHOQ_t_USER(k,i,j,iq) = RHOQ_t_USER(k,i,j,iq) &
                                      - 0.5_RP*(SINK_UP + SINK_CEN) &
                                      * (QTRC_column(k+1,iq) - QTRC_column(k,iq))*RFDZ(k)
             else
                do iq = QS_MP, QE_MP
                   RHOQ_t_USER(k,i,j,iq) = RHOQ_t_USER(k,i,j,iq) &
                                         - 0.5_RP*(SINK_UP + SINK_CEN) &
                                         * (QTRC_column(k+1,iq) - QTRC_column(k,iq))*RFDZ(k)
                enddo
             endif

             ! -- energy --
             RHOT_t_USER(k,i,j) = RHOT_t_USER(k,i,j) &
                                - 0.5_RP*(SINK_UP + SINK_CEN) &
                                * (POTT_column(k+1) - POTT_column(k))*RFDZ(k)

          enddo

          ! large-scale sinking forcing, the simple first-order upwind advection scheme is used
          do k = KS, KE-1
             SINK_CEN = WLS(k)
             SINK_UP = WLS(k+1)
             SINK_DUP = WLS(k+2)

             ! -- z momentum --
             if (SHEBA_SWITCH_ACCE) then
                MOMZ_t_USER(k,i,j) = MOMZ_t_USER(k,i,j) &
                                   + 0.5_RP*(SINK_CEN + SINK_UP) &
                                   * 0.5_RP*(DENS_column(k) + DENS_column(k+1))
             else
                MOMZ_t_USER(k,i,j) = MOMZ_t_USER(k,i,j) &
                                   - SINK_UP &
                                   * ( MOMZ_column(k+1)*2.0_RP/( DENS_column(k+2) &
                                                               + DENS_column(k+1)) &
                                     -  MOMZ_column(k  )*2.0_RP/( DENS_column(k+1) &
                                                                + DENS_column(k  )))*RCDZ(k+1)
             endif
          enddo

       endif ! end of experiment

    enddo
    enddo

    call FILE_HISTORY_in( MOMZ_t_USER(:,:,:), 'MOMZ_t_USER', 'large-scale momz sinking',            'kg/m2/s2',  fill_halo=.true. )

    call FILE_HISTORY_in( RHOU_t_USER(:,:,:), 'RHOU_t_USER', 'large-scale rhou sinking',            'kg/m2/s2',  fill_halo=.true. )

    call FILE_HISTORY_in( RHOV_t_USER(:,:,:), 'RHOV_t_USER', 'large-scale rhov sinking',            'kg/m2/s2',  fill_halo=.true. )

    call FILE_HISTORY_in( DENS_t_USER(:,:,:), 'DENS_t_USER', 'large-scale rho sinking',               'kg/m3/s',   fill_halo=.true. )

    call FILE_HISTORY_in( RHOT_t_USER(:,:,:), 'RHOT_t_USER', 'large-scale pt sinking', 'K*kg/m3/s', fill_halo=.true. )

    call FILE_HISTORY_in( TEMP_t_USER(:,:,:), 'TEMP_t_USER', 'large-scale pt cooling', 'K*kg/m3/s', fill_halo=.true. )

    iq = QS_MP
    call FILE_HISTORY_in( RHOQ_t_USER(:,:,:,iq), 'RHOQ_t_USER', 'large-scale qv sinking',          'kg/m3/s',   fill_halo=.true. )

    call FILE_HISTORY_in( subsidence_sink(:,:,:), 'subsidence_sink', 'large-scale sinking', 'm/s', fill_halo=.true. )

    if ( SHEBA_SWITCH_MOMZ ) then
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE
          MOMZ_t(k,i,j) = MOMZ_t(k,i,j) + MOMZ_t_USER(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( SHEBA_SWITCH_RHOU ) then
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE
          RHOU_t(k,i,j) = RHOU_t(k,i,j) + RHOU_t_USER(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( SHEBA_SWITCH_RHOV ) then
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE
          RHOV_t(k,i,j) = RHOV_t(k,i,j) + RHOV_t_USER(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( SHEBA_SWITCH_DENS ) then
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE
          DENS_t(k,i,j) = DENS_t(k,i,j) + DENS_t_USER(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( SHEBA_SWITCH_RHOT ) then
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE
          RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_USER(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( SHEBA_SWITCH_TEMP ) then
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE
          RHOT_t(k,i,j) = RHOT_t(k,i,j) + TEMP_t_USER(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( SHEBA_SWITCH_RHOQ ) then
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE
          do iq = 1, QA
             RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_USER(k,i,j,iq)
          enddo
       enddo
       enddo
       enddo
    endif

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
       !call PRC_abort
    endif

    return
  end subroutine USER_update

  !-----------------------------------------------------------------------------
  !> Read sounding data from file specialized for SHEBA project
  subroutine read_sounding_sheba( &
       DENS, VELX, VELY, POTT, QV, QCI, QNCI )
    use scale_const, only: &
       P00   => CONST_PRE00, &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap, &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ, &
       FZ => ATMOS_GRID_CARTESC_FZ
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    implicit none

    real(RP), intent(out) :: DENS(KA)
    real(RP), intent(out) :: VELX(KA)
    real(RP), intent(out) :: VELY(KA)
    real(RP), intent(out) :: POTT(KA)
    real(RP), intent(out) :: QV  (KA)
    real(RP), intent(out) :: QCI (KA)
    real(RP), intent(out) :: QNCI(KA)

    real(RP) :: TEMP(KA)
    real(RP) :: PRES(KA)
    real(RP) :: QC  (KA)

    real(RP) :: QV_ (KA)

    character(len=H_LONG) :: ENV_IN_SOUNDING_file = ''
    logical :: USE_HYDROSTATIC

    integer, parameter :: EXP_klim = 501 ! there was a bug: EXP_klim = 100 (30/3/2020)
    integer            :: EXP_kmax

    real(RP) :: SFC_THETA            ! surface potential temperature [K]
    real(RP) :: SFC_TEMP             ! temperature [K]
    real(RP) :: SFC_PRES             ! surface pressure [hPa]
    real(RP) :: SFC_RH               ! surface relative humidity
    real(RP) :: SFC_QV               ! surface watervapor [g/kg]
    real(RP) :: SFC_QCI              ! surface liquid and ice [g/kg]
    real(RP) :: SFC_QNCI             ! surface liquid and ice conc. [/cm3]

    real(RP) :: pres_sfc
    real(RP) :: pott_sfc
    real(RP) :: temp_sfc
    real(RP) :: qv_sfc
    real(RP) :: qc_sfc

    real(RP) :: EXP_z   (EXP_klim+1) ! height      [m]
    real(RP) :: EXP_pres(EXP_klim+1) ! pressure    [Pa]
    real(RP) :: EXP_temp(EXP_klim+1) ! temperature [K]
    real(RP) :: EXP_pott(EXP_klim+1) ! potential temperature [K]
    real(RP) :: EXP_rh  (EXP_klim+1) ! relative humidity
    real(RP) :: EXP_qv  (EXP_klim+1) ! water vapor [g/kg]
    real(RP) :: EXP_qci (EXP_klim+1) ! liquid and ice [g/kg]
    real(RP) :: EXP_qnci(EXP_klim+1) ! liquid and ice conc. [/cm3]
    real(RP) :: EXP_u   (EXP_klim+1) ! velocity u  [m/s]
    real(RP) :: EXP_v   (EXP_klim+1) ! velocity v  [m/s]

    real(RP) :: fact1, fact2, R_gas, CP_gas, qd_assumption

    logical :: converged

    integer :: k, kref
    integer :: fid
    integer :: ierr

    namelist / PARAM_MKINIT_SOUNDING / &
       ENV_IN_SOUNDING_file, &
       USE_HYDROSTATIC

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SOUNDING,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("read_sounding_sheba",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("read_sounding_sheba",*) 'Not appropriate names in namelist PARAM_MKINIT_SOUNDING. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SOUNDING)

    !--- prepare sounding profile
    LOG_INFO("read_sounding_sheba",*) 'Input sounding file:', trim(ENV_IN_SOUNDING_file)
    fid = IO_get_available_fid()
    open( fid,                                 &
          file   = trim(ENV_IN_SOUNDING_file), &
          form   = 'formatted',                &
          status = 'old',                      &
          iostat = ierr                        )

    if ( ierr /= 0 ) then
       LOG_ERROR("read_sounding_sheba",*) '[mod_mkinit/read_sounding] Input file not found!'
    endif


    !--- read sounding file till end
    read(fid,*) SFC_PRES, SFC_TEMP, SFC_THETA, SFC_QV, SFC_QCI, SFC_QNCI

    LOG_INFO("read_sounding_sheba",*) '+ Surface pressure [hPa]',                         SFC_PRES
    LOG_INFO("read_sounding_sheba",*) '+ Surface temperature [K]',                        SFC_TEMP
    LOG_INFO("read_sounding_sheba",*) '+ Surface pot. temp  [K]',                         SFC_THETA
    LOG_INFO("read_sounding_sheba",*) '+ Surface water vapor [g/kg]',                     SFC_QV
    LOG_INFO("read_sounding_sheba",*) '+ Surface liquid and ice [g/kg]',                  SFC_QCI
    LOG_INFO("read_sounding_sheba",*) '+ Surface liquid and ice conc. [/cm3]',            SFC_QNCI

    do k = 2, EXP_klim
       read(fid,*,iostat=ierr) EXP_z(k), EXP_pres(k), EXP_temp(k), EXP_pott(k), EXP_qv(k), EXP_u(k), EXP_v(k), EXP_qci(k), EXP_qnci(k)
       if ( ierr /= 0 ) exit
    enddo

    EXP_kmax = k - 1
    close(fid)

    ! Boundary
    EXP_z   (1)          = 0.0_RP
    EXP_pres(1)          = SFC_PRES
    EXP_temp(1)          = SFC_TEMP
    EXP_pott(1)          = SFC_THETA
    EXP_qv  (1)          = SFC_QV
    EXP_qci (1)          = SFC_QCI
    EXP_qnci(1)          = SFC_QNCI
    EXP_u   (1)          = EXP_u   (2)
    EXP_v   (1)          = EXP_v   (2)
    EXP_z   (EXP_kmax+1) = 100.E3_RP
    EXP_pres(EXP_kmax+1) = EXP_pres(EXP_kmax)
    EXP_temp(EXP_kmax+1) = EXP_temp(EXP_kmax)
    EXP_pott(EXP_kmax+1) = EXP_pott(EXP_kmax)
    EXP_qv  (EXP_kmax+1) = EXP_qv  (EXP_kmax)
    EXP_qci (EXP_kmax+1) = EXP_qci (EXP_kmax)
    EXP_qnci(EXP_kmax+1) = EXP_qnci(EXP_kmax)
    EXP_u   (EXP_kmax+1) = EXP_u   (EXP_kmax)
    EXP_v   (EXP_kmax+1) = EXP_v   (EXP_kmax)

    do k = 1, EXP_kmax+1
       EXP_qv(k)   = EXP_qv(k)   * 1.E-3_RP ! [g/kg]->[kg/kg]
       EXP_qci(k)  = EXP_qci(k)  * 1.E-3_RP ! [g/kg]->[kg/kg]
       EXP_pres(k) = EXP_pres(k) * 1.E2_RP  ! [hPa]->[Pa]
    enddo

    ! calc in dry condition
    pres_sfc = SFC_PRES * 1.E2_RP ! [hPa]->[Pa]
    pott_sfc = SFC_THETA
    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       qv_sfc   = SFC_QV  * 1.E-3_RP ! [g/kg]->[kg/kg]
       qc_sfc   = SFC_QCI * 1.E-3_RP ! [g/kg]->[kg/kg]
    end if

    !--- linear interpolate to model grid
    do k = KS, KE
       do kref = 2, EXP_kmax+1
          if (       CZ(k) >  EXP_z(kref-1) &
               .AND. CZ(k) <= EXP_z(kref  ) ) then

             fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
             fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

             PRES(k) = EXP_pres(kref-1) * fact1 &
                     + EXP_pres(kref  ) * fact2
             TEMP(k) = EXP_temp(kref-1) * fact1 &
                     + EXP_temp(kref  ) * fact2
             POTT(k) = EXP_pott(kref-1) * fact1 &
                     + EXP_pott(kref  ) * fact2
             QV  (k) = EXP_qv  (kref-1) * fact1 &
                     + EXP_qv  (kref  ) * fact2
             VELX(k) = EXP_u   (kref-1) * fact1 &
                     + EXP_u   (kref  ) * fact2
             VELY(k) = EXP_v   (kref-1) * fact1 &
                     + EXP_v   (kref  ) * fact2
             QCI (k) = EXP_qci (kref-1) * fact1 &
                     + EXP_qci (kref  ) * fact2
             QNCI(k) = EXP_qnci(kref-1) * fact1 &
                     + EXP_qnci(kref  ) * fact2
          endif
       enddo
    enddo

    if ( ATMOS_HYDROMETEOR_dry ) QV(:) = 0.0_RP

    if ( USE_HYDROSTATIC ) then
       ! make density & pressure profile in moist condition
       call HYDROSTATIC_buildrho( KA, KS, KE, &
                                  POTT(:), QV(:), QCI(:),              & ! [IN]
                                  pres_sfc, pott_sfc, qv_sfc, qc_sfc,  & ! [IN]
                                  CZ(:), FZ(:),                        & ! [IN]
                                  DENS(:), TEMP(:), PRES(:), temp_sfc, & ! [OUT]
                                  converged                            ) ! [OUT]
    else
       do k = KS, KE
          ! use gas law to deduce density
          ! ASSUMPTION: rho_d = rho
          qd_assumption = 1.0_RP - QV(k) - QCI(k)
          R_gas = qd_assumption * Rdry  + QV(k) * Rvap
          CP_gas = qd_assumption * CPdry + QV(k) * CPvap + QCI(k) * CL
          !R_gas = (1.0_RP - QV(k)) * Rdry  + QV(k) * Rvap
          DENS(k) = PRES(k) / R_gas / TEMP(k)
          POTT(k) = TEMP(k) * (P00 / PRES(k))**(R_gas/CP_gas)
       enddo
    endif

    LOG_INFO("read_sounding_sheba",*) 'Checking initial condition'
    LOG_INFO("read_sounding_sheba",*) 'Z          D          T          P          O          V          C'
    do k = KS, KE
       LOG_INFO("read_sounding_sheba",'(7ES15.6)') CZ(k), DENS(k), TEMP(k), PRES(k), POTT(k),  QV(k), QCI(k)
    enddo

    return
  end subroutine read_sounding_sheba

end module mod_user
