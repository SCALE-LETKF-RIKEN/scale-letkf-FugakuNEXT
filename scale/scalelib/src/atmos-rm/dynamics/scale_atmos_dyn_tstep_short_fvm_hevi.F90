!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          HEVI FVM scheme for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#ifdef _OPENACC
#define HEVI_FISSION 1
#endif

#ifdef PROFILE_FAPP
#define PROFILE_START(name) call fapp_start(name, 1, 1)
#define PROFILE_STOP(name)  call fapp_stop (name, 1, 1)
#else
#define PROFILE_START(name)
#define PROFILE_STOP(name)
#endif

! NVTX is provided by the NVHPC compiler (-cudalib=nvtx)
#if defined(_OPENACC) && defined(NVIDIA)
#define NVTX_PUSH(name) nvtx_level = nvtxRangePush(name)
#define NVTX_POP()      nvtx_level = nvtxRangePop()
#else
#define NVTX_PUSH(name)
#define NVTX_POP()
#endif

#include "scalelib.h"

! the vector length macro is defined in scalelib.h,
! so this check must come after the include
#if defined(HEVI_FISSION) && LSIZE > 1
#error "HEVI_FISSION is not supported with LSIZE > 1"
#endif

module scale_atmos_dyn_tstep_short_fvm_hevi
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
#if defined DEBUG || defined QUICKDEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
#if defined(_OPENACC) && defined(NVIDIA)
  use nvtx
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tstep_short_fvm_hevi_regist
  public :: ATMOS_DYN_Tstep_short_fvm_hevi_setup
  public :: ATMOS_DYN_Tstep_short_fvm_hevi

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
#if 1
! weight for full-level k+1; weight for k is 1 - F2H(k,idx)
#define F2H(k,idx) (CDZ(k)*GSQRT(k,i,j,idx)/(CDZ(k)*GSQRT(k,i,j,idx)+CDZ(k+1)*GSQRT(k+1,i,j,idx)))
#else
#define F2H(k,idx) 0.5_RP
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: NB = 1
  integer,  private, parameter :: VA_FVM_HEVI = 0
  integer                      :: IFS_OFF
  integer                      :: JFS_OFF

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_Tstep_short_fvm_hevi_regist( &
       ATMOS_DYN_TYPE, &
       VA_out,         &
       VAR_NAME,       &
       VAR_DESC,       &
       VAR_UNIT        )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*),       intent(in)  :: ATMOS_DYN_TYPE
    integer,                intent(out) :: VA_out         !< number of prognostic variables
    character(len=H_SHORT), intent(out) :: VAR_NAME(:)    !< name   of the variables
    character(len=H_MID),   intent(out) :: VAR_DESC(:)    !< desc.  of the variables
    character(len=H_SHORT), intent(out) :: VAR_UNIT(:)    !< unit   of the variables
    !---------------------------------------------------------------------------

    if ( ATMOS_DYN_TYPE /= 'FVM-HEVI' .AND. ATMOS_DYN_TYPE /= 'HEVI' ) then
       LOG_ERROR("ATMOS_DYN_Tstep_short_fvm_hevi_regist",*) 'ATMOS_DYN_TYPE is not FVM-HEVI. Check!'
       call PRC_abort
    endif

    VA_out      = VA_FVM_HEVI
    VAR_NAME(:) = ""
    VAR_DESC(:) = ""
    VAR_UNIT(:) = ""

    LOG_NEWLINE
    LOG_INFO("ATMOS_DYN_Tstep_short_fvm_hevi_regist",*) 'Register additional prognostic variables (HEVI)'
    if ( VA_out < 1 ) then
       LOG_INFO_CONT(*) '=> nothing.'
    endif

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_hevi_regist

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tstep_short_fvm_hevi_setup
    implicit none
    !---------------------------------------------------------------------------

    LOG_INFO("ATMOS_DYN_Tstep_short_fvm_hevi_setup",*) 'HEVI Setup'

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_hevi_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_Tstep_short_fvm_hevi( &
       DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
       PROG_RK,                                     &
       mflx_hi, tflx_hi,                            &
       DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
       DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
       PROG0, PROG,                                 &
       DPRES0, RT2P, CORIOLI,                       &
       num_diff, wdamp_coef, divdmp_coef, DDIV,     &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T,               &
       FLAG_FCT_ALONG_STREAM,                       &
       CDZ, FDZ, FDX, FDY,                          &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
       PHI, GSQRT, J13G, J23G, J33G, MAPF,          &
       REF_dens, REF_rhot,                          &
       BND_W, BND_E, BND_S, BND_N, TwoD,            &
       dtrk, last                                   )
    use scale_atmos_grid_cartesC_index
    use scale_const, only: &
#ifdef DRY
       Rdry   => CONST_Rdry,  &
       CVdry  => CONST_CVdry, &
       CPdry  => CONST_CPdry, &
#endif
       EPS    => CONST_EPS,   &
       GRAV   => CONST_GRAV,  &
       P00    => CONST_PRE00
    use scale_atmos_dyn_fvm_fct, only: &
       ATMOS_DYN_FVM_fct
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_flux_valueW_Z, &
       ATMOS_DYN_FVM_fluxZ_XYZ,     &
       ATMOS_DYN_FVM_fluxX_XYZ,     &
       ATMOS_DYN_FVM_fluxY_XYZ,     &
       ATMOS_DYN_FVM_fluxZ_XYW,     &
       ATMOS_DYN_FVM_fluxJ13_XYW,   &
       ATMOS_DYN_FVM_fluxJ23_XYW,   &
       ATMOS_DYN_FVM_fluxX_XYW,     &
       ATMOS_DYN_FVM_fluxY_XYW,     &
       ATMOS_DYN_FVM_fluxZ_UYZ,     &
       ATMOS_DYN_FVM_fluxJ13_UYZ,   &
       ATMOS_DYN_FVM_fluxJ23_UYZ,   &
       ATMOS_DYN_FVM_fluxX_UYZ,     &
       ATMOS_DYN_FVM_fluxY_UYZ,     &
       ATMOS_DYN_FVM_fluxZ_XVZ,     &
       ATMOS_DYN_FVM_fluxJ13_XVZ,   &
       ATMOS_DYN_FVM_fluxJ23_XVZ,   &
       ATMOS_DYN_FVM_fluxX_XVZ,     &
       ATMOS_DYN_FVM_fluxY_XVZ
#ifdef HIST_TEND
    use scale_file_history, only: &
       FILE_HISTORY_in
#endif
    use scale_matrix, only: &
       MATRIX_SOLVER_TRIDIAGONAL, &
       MATRIX_SOLVER_TRIDIAGONAL_1D_CR
#if defined(HEVI_FISSION) && defined(USE_CUDALIB)
    use scale_prc, only: &
       PRC_abort
#endif
    implicit none

    real(RP), intent(out) :: DENS_RK(KA,IA,JA)   ! prognostic variables
    real(RP), intent(out) :: MOMZ_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMX_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMY_RK(KA,IA,JA)   !
    real(RP), intent(out) :: RHOT_RK(KA,IA,JA)   !

    real(RP), intent(out) :: PROG_RK(KA,IA,JA,VA)  !

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3) ! rho * vel(x,y,z)
    real(RP), intent(out)   :: tflx_hi(KA,IA,JA,3) ! rho * theta * vel(x,y,z)

    real(RP), intent(in),target :: DENS0(KA,IA,JA) ! prognostic variables
    real(RP), intent(in),target :: MOMZ0(KA,IA,JA) ! at previous dynamical time step
    real(RP), intent(in),target :: MOMX0(KA,IA,JA) !
    real(RP), intent(in),target :: MOMY0(KA,IA,JA) !
    real(RP), intent(in),target :: RHOT0(KA,IA,JA) !

    real(RP), intent(in) :: DENS(KA,IA,JA)   ! prognostic variables
    real(RP), intent(in) :: MOMZ(KA,IA,JA)   ! at previous RK step
    real(RP), intent(in) :: MOMX(KA,IA,JA)   !
    real(RP), intent(in) :: MOMY(KA,IA,JA)   !
    real(RP), intent(in) :: RHOT(KA,IA,JA)   !

    real(RP), intent(in) :: DENS_t(KA,IA,JA)
    real(RP), intent(in) :: MOMZ_t(KA,IA,JA)
    real(RP), intent(in) :: MOMX_t(KA,IA,JA)
    real(RP), intent(in) :: MOMY_t(KA,IA,JA)
    real(RP), intent(in) :: RHOT_t(KA,IA,JA)

    real(RP), intent(in)  :: PROG0(KA,IA,JA,VA)
    real(RP), intent(in)  :: PROG (KA,IA,JA,VA)

    real(RP), intent(in) :: DPRES0(KA,IA,JA)
    real(RP), intent(in) :: RT2P(KA,IA,JA)
    real(RP), intent(in) :: CORIOLI(IA,JA)
    real(RP), intent(in) :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in) :: wdamp_coef(KA)
    real(RP), intent(in) :: divdmp_coef
    real(RP), intent(in) :: DDIV(KA,IA,JA)

    logical,  intent(in) :: FLAG_FCT_MOMENTUM
    logical,  intent(in) :: FLAG_FCT_T
    logical,  intent(in) :: FLAG_FCT_ALONG_STREAM

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: FDZ(KA-1)
    real(RP), intent(in) :: FDX(IA-1)
    real(RP), intent(in) :: FDY(JA-1)
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in) :: RFDZ(KA-1)
    real(RP), intent(in) :: RFDX(IA-1)
    real(RP), intent(in) :: RFDY(JA-1)

    real(RP), intent(in)  :: PHI     (KA,IA,JA)           !< geopotential
    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,I_XYZ_MAX) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,I_XYZ_MAX) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,I_XYZ_MAX) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                         !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF    (IA,JA,2,I_XY_MAX)   !< map factor
    real(RP), intent(in)  :: REF_dens(KA,IA,JA)           !< reference density
    real(RP), intent(in)  :: REF_rhot(KA,IA,JA)

    logical,  intent(in)  :: BND_W
    logical,  intent(in)  :: BND_E
    logical,  intent(in)  :: BND_S
    logical,  intent(in)  :: BND_N
    logical,  intent(in)  :: TwoD

    real(RP), intent(in)  :: dtrk
    logical,  intent(in)  :: last


    ! diagnostic variables (work space)
    real(RP) :: POTT(KA,IA,JA) ! potential temperature [K]
    real(RP) :: DPRES(KA,IA,JA) ! pressure deviation from reference pressure

    real(RP) :: qflx_hi (KA,IA,JA,3)
    real(RP) :: qflx_J13(KA,IA,JA)
    real(RP) :: qflx_J23(KA,IA,JA)

    real(RP) :: advch ! horizontal advection
    real(RP) :: advcv ! vertical advection
    real(RP) :: wdmp  ! Raileight damping
    real(RP) :: div   ! divergence damping
    real(RP) :: pg    ! pressure gradient force
    real(RP) :: cf    ! colioris force
    real(RP) :: advc  ! advection
    real(RP) :: momy_u ! momentum y at u point
    real(RP) :: momx_v ! momentum x at v point
    real(RP) :: dpresm ! pressure deviation from reference pressure at k-1
    real(RP) :: dpresk ! pressure deviation from reference pressure at k
    real(RP) :: dpresp ! pressure deviation from reference pressure at k+1
    real(RP) :: f2h1, f2h2   ! F2H weights at k (f2h2 = 1 - f2h1)
    real(RP) :: f2h1m, f2h2m ! F2H weights at k-1
#ifdef HIST_TEND
    real(RP) :: advch_t(KA,IA,JA,5)
    real(RP) :: advcv_t(KA,IA,JA,5)
    real(RP) :: wdmp_t(KA,IA,JA)
    real(RP) :: ddiv_t(KA,IA,JA,3)
    real(RP) :: pg_t(KA,IA,JA,3)
    real(RP) :: cf_t(KA,IA,JA,2)
    logical  :: lhist
#endif

    ! for implicit solver
#ifdef _OPENACC
    real(RP) :: A0, A1
#else
    real(RP) :: A(KA)
#endif
    real(RP) :: B
    real(RP) :: fact
    real(RP) :: Sr(KA,IA,JA)
    real(RP) :: Sw(KA,IA,JA)
    real(RP) :: St(KA,IA,JA)

#ifdef HEVI_FISSION
    ! work arrays to pass intermediate values between the split kernels
    real(RP) :: pg_work(KA,IA,JA)
    real(RP) :: cf_work(KA,IA,JA)
    real(RP) :: advc_work(KA,IA,JA)
    real(RP) :: div_work(KA,IA,JA)
#endif

#ifdef HEVI_FISSION
    real(RP), allocatable :: PT_work(:,:,:)
    real(RP), allocatable :: Ci_work(:,:,:)
    real(RP), allocatable :: Co_work(:,:,:)
    real(RP), allocatable :: F1_work(:,:,:)
    real(RP), allocatable :: F2_work(:,:,:)
    real(RP), allocatable :: F3_work(:,:,:)
    ! The i and j indices are taken implicitly from the enclosing loop, so these
    ! macros must only be used inside the body of an (i,j) loop nest.
    ! They are #undef'ed just after the end of this subroutine.
#define PT(k,l) PT_work(k,i,j)
#define Ci(k,l) Ci_work(k,i,j)
#define Co(k,l) Co_work(k,i,j)
#define F1(k,l) F1_work(k,i,j)
#define F2(k,l) F2_work(k,i,j)
#define F3(k,l) F3_work(k,i,j)
#else
    real(RP) :: PT(KA,LSIZE)
    real(RP) :: Ci(KS:KE-1,LSIZE)
    real(RP) :: Co(KS:KE-1,LSIZE)
    real(RP) :: F1(KS:KE-1,LSIZE)
    real(RP) :: F2(KS:KE-1,LSIZE)
    real(RP) :: F3(KS:KE-1,LSIZE)
#endif

#ifdef _OPENACC
    real(RP) :: work(KMAX-1,4) ! for CR
#endif
#if defined(_OPENACC) && defined(NVIDIA)
    integer :: nvtx_level
#endif

    ! for temporary variables
    real(RP) :: tmp

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j, ii
    integer :: iss, iee
#if LSIZE == 1
    integer, parameter :: l = 1
#else
    integer  :: l
#endif


#ifdef DEBUG
    POTT(:,:,:) = UNDEF
    DPRES(:,:,:) = UNDEF

    qflx_hi (:,:,:,:) = UNDEF
    qflx_J13(:,:,:)   = UNDEF
    qflx_J23(:,:,:)   = UNDEF
#endif

#if defined DEBUG || defined QUICKDEBUG
    DENS_RK(   1:KS-1,:,:)   = UNDEF
    DENS_RK(KE+1:KA  ,:,:)   = UNDEF
    MOMZ_RK(   1:KS-1,:,:)   = UNDEF
    MOMZ_RK(KE+1:KA  ,:,:)   = UNDEF
    MOMX_RK(   1:KS-1,:,:)   = UNDEF
    MOMX_RK(KE+1:KA  ,:,:)   = UNDEF
    MOMY_RK(   1:KS-1,:,:)   = UNDEF
    MOMY_RK(KE+1:KA  ,:,:)   = UNDEF
    RHOT_RK(   1:KS-1,:,:)   = UNDEF
    RHOT_RK(KE+1:KA  ,:,:)   = UNDEF
    PROG_RK(   1:KS-1,:,:,:) = UNDEF
    PROG_RK(KE+1:KA  ,:,:,:) = UNDEF
#endif

#ifdef HIST_TEND
    lhist = last
#endif

#ifdef PROFILE_FIPP
    call fipp_start()
#endif

#if defined(HEVI_FISSION) && defined(USE_CUDALIB)
    ! The cuSPARSE branch of MATRIX_SOLVER_tridiagonal_3D solves the whole
    ! (IA,JA) plane in one batched call and ignores the IS:IE / JS:JE loop
    ! range, so with cache blocking it would repeat the full-plane solve for
    ! every block and read the work arrays outside the current block, which
    ! are not filled yet.
    if ( IBLOCK /= IMAX .or. JBLOCK /= JMAX ) then
       LOG_ERROR("ATMOS_DYN_Tstep_short_fvm_hevi",*) 'cache blocking is not supported with cuSPARSE. Set IBLOCK = IMAX and JBLOCK = JMAX. Check!'
       call PRC_abort
    end if
#endif

    IFS_OFF = 1
    JFS_OFF = 1
    if ( BND_W ) IFS_OFF = 0
    if ( BND_S ) JFS_OFF = 0

#ifdef HEVI_FISSION
    ! allocated here (not inside the block loop) to avoid repeated device
    ! allocation, and freed on return so that no memory is held outside dynamics
    allocate( PT_work(KA,      IS:IE, JS:JE) )
    allocate( Ci_work(KS:KE-1, IS:IE, JS:JE) )
    allocate( Co_work(KS:KE-1, IS:IE, JS:JE) )
    allocate( F1_work(KS:KE-1, IS:IE, JS:JE) )
    allocate( F2_work(KS:KE-1, IS:IE, JS:JE) )
    allocate( F3_work(KS:KE-1, IS:IE, JS:JE) )
#endif

    !$acc data &
    !$acc copy(mflx_hi) &
    !$acc copyout(DENS_RK,MOMZ_RK,MOMX_RK,MOMY_RK,RHOT_RK,PROG_RK, &
    !$acc         tflx_hi) &
#ifdef HIST_TEND
    !$acc copyout(advch_t,advcv_t,wdmp_t,ddiv_t,pg_t,cf_t) &
#endif
    !$acc copyin(DENS0,MOMZ0,MOMX0,MOMY0,RHOT0, &
    !$acc        DENS,MOMZ,MOMX,MOMY,RHOT,DENS_t,MOMZ_t,MOMX_t,MOMY_t,RHOT_t, &
    !$acc        PROG0,PROG, &
    !$acc        DPRES0,RT2P,CORIOLI,num_diff,wdamp_coef,divdmp_coef,DDIV, &
    !$acc        FLAG_FCT_MOMENTUM,FLAG_FCT_T,FLAG_FCT_ALONG_STREAM, &
    !$acc        CDZ,FDZ,FDX,FDY,RCDZ,RCDX,RCDY,RFDZ,RFDX,RFDY, &
    !$acc        PHI,GSQRT,J13G,J23G,J33G,MAPF,REF_dens,REF_rhot, &
    !$acc        BND_W,BND_E,BND_S,BND_N,TwoD,dtrk,last) &
    !$acc create(POTT,DPRES, &
    !$acc        qflx_hi,qflx_J13,qflx_J23, &
#ifdef HEVI_FISSION
    !$acc        pg_work,cf_work,advc_work,div_work, &
    !$acc        PT_work,Ci_work,Co_work,F1_work,F2_work,F3_work, &
#endif
    !$acc        Sr,Sw,St)

#ifdef HEVI_FISSION
    if ( divdmp_coef .le. 0.0_RP ) then
       !$acc kernels
       div_work(:,:,:) = 0.0_RP
       !$acc end kernels
    end if
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       PROFILE_START("hevi_pres")
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JJS,JJE,IIS,IIE,IA,KS,KE,DPRES0,RT2P,RHOT,REF_rhot,DPRES,DENS,PHI)
#ifdef HEVI_FISSION
       !$acc kernels async(0)
#else
       !$acc kernels
#endif
       do j = JJS, JJE+1
       do i = IIS, min(IIE+1,IA)
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, DPRES0(k,i,j) )
             call CHECK( __LINE__, RT2P(k,i,j) )
             call CHECK( __LINE__, RHOT(k,i,j) )
             call CHECK( __LINE__, REF_rhot(k,i,j) )
#endif
             DPRES(k,i,j) = DPRES0(k,i,j) + RT2P(k,i,j) * ( RHOT(k,i,j) - REF_rhot(k,i,j) )
          enddo
#ifdef HEVI_FISSION
       enddo
       enddo
       !$acc end kernels

       !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
       !$omp private(i,j)
       !$acc kernels async(1)
       do j = JJS, JJE+1
       do i = IIS, min(IIE+1,IA)
#endif
          DPRES(KS-1,i,j) = DPRES0(KS-1,i,j) - DENS(KS,i,j) * ( PHI(KS-1,i,j) - PHI(KS+1,i,j) )
          DPRES(KE+1,i,j) = DPRES0(KE+1,i,j) - DENS(KE,i,j) * ( PHI(KE+1,i,j) - PHI(KE-1,i,j) )
       enddo
       enddo
       !$acc end kernels
       PROFILE_STOP("hevi_pres")

       !##### continuity equation #####

       PROFILE_START("hevi_mflx_z")
       ! at (x, y, w)
       if ( TwoD ) then
          !$omp parallel do default(none) private(j,k) OMP_SCHEDULE_ &
          !$omp shared(JJS,JJE,IS,KS,KE,GSQRT,MOMY,J23G,mflx_hi,MAPF,num_diff)
          !$acc kernels
          do j = JJS, JJE
             mflx_hi(KS-1,IS,j,ZDIR) = 0.0_RP
             do k = KS, KE-1
#ifdef DEBUG
                call CHECK( __LINE__, MOMY(k+1,IS,j) )
                call CHECK( __LINE__, MOMY(k+1,IS,j-1) )
                call CHECK( __LINE__, MOMY(k  ,IS,j) )
                call CHECK( __LINE__, MOMY(k  ,IS,j-1) )
#endif
                mflx_hi(k,IS,j,ZDIR) = J23G(k,IS,j,I_XYW) * 0.25_RP &
                                     * ( MOMY(k+1,IS,j)+MOMY(k+1,IS,j-1) &
                                       + MOMY(k  ,IS,j)+MOMY(k  ,IS,j-1) ) &
                                      + GSQRT(k,IS,j,I_XYW) / MAPF(IS,j,2,I_XY) * num_diff(k,IS,j,I_DENS,ZDIR)
             enddo
             mflx_hi(KE,IS,j,ZDIR) = 0.0_RP
          enddo
          !$acc end kernels
       else
          !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JJS,JJE,IIS,IIE,KS,KE,GSQRT,MOMX,MOMY,J13G,J23G,mflx_hi,MAPF,num_diff)
#ifdef HEVI_FISSION
          !$acc kernels async(0)
#else
          !$acc kernels
#endif
          do j = JJS, JJE
          do i = IIS-1, IIE
             do k = KS, KE-1
#ifdef DEBUG
                call CHECK( __LINE__, MOMX(k+1,i  ,j) )
                call CHECK( __LINE__, MOMX(k+1,i-1,j) )
                call CHECK( __LINE__, MOMX(k  ,i  ,j) )
                call CHECK( __LINE__, MOMX(k  ,i+1,j) )
                call CHECK( __LINE__, MOMY(k+1,i,j) )
                call CHECK( __LINE__, MOMY(k+1,i,j-1) )
                call CHECK( __LINE__, MOMY(k  ,i,j) )
                call CHECK( __LINE__, MOMY(k  ,i,j-1) )
#endif
                mflx_hi(k,i,j,ZDIR) = J13G(k,i,j,I_XYW) * 0.25_RP / MAPF(i,j,2,I_XY) &
                                    * ( MOMX(k+1,i,j)+MOMX(k+1,i-1,j) &
                                      + MOMX(k  ,i,j)+MOMX(k  ,i-1,j) ) &
                                    + J23G(k,i,j,I_XYW) * 0.25_RP / MAPF(i,j,1,I_XY) &
                                    * ( MOMY(k+1,i,j)+MOMY(k+1,i,j-1) &
                                      + MOMY(k  ,i,j)+MOMY(k  ,i,j-1) ) &
                                     + GSQRT(k,i,j,I_XYW) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) ) * num_diff(k,i,j,I_DENS,ZDIR)
             enddo
#ifdef HEVI_FISSION
          enddo
          enddo
          !$acc end kernels

          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j)
          !$acc kernels async(1)
          do j = JJS, JJE
          do i = IIS-1, IIE
#endif
             mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
             mflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
          !$acc end kernels
       end if
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_mflx_z")

       if ( .not. TwoD ) then
          PROFILE_START("hevi_mflx_x")
          iss = max(IIS-1,IS-IFS_OFF)
          iee = min(IIE,IEH)
          ! at (u, y, z)
          !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JJS,JJE,iss,iee,KS,KE,GSQRT,I_UYZ,MOMX,num_diff,mflx_hi,MAPF,I_UY)
          !$acc kernels
          do j = JJS, JJE
          do i = iss, iee
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
             call CHECK( __LINE__, MOMX(k,i,j) )
             call CHECK( __LINE__, num_diff(k,i,j,I_DENS,XDIR) )
#endif
             mflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY) &
                  * ( MOMX(k,i,j) + num_diff(k,i,j,I_DENS,XDIR) )
          enddo
          enddo
          enddo
          !$acc end kernels
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          PROFILE_STOP("hevi_mflx_x")
       end if

       ! at (x, v, z)
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JJS,JS,JFS_OFF,JJE,JEH,IIS,IIE,KS,KE,GSQRT,MOMY,num_diff,mflx_hi,MAPF)
       !$acc kernels
       do j = max(JJS-1,JS-JFS_OFF), min(JJE,JEH)
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
          call CHECK( __LINE__, MOMY(k,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,YDIR) )
#endif
          mflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV) &
               * ( MOMY(k,i,j) + num_diff(k,i,j,I_DENS,YDIR) )
       enddo
       enddo
       enddo
       !$acc end kernels
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update density
       PROFILE_START("hevi_sr")
#ifdef HEVI_FISSION
       ! mflx_hi(:,:,:,ZDIR) is computed asynchronously
       !$acc wait
#endif
       if ( TwoD ) then
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(j,k,advcv,advch,tmp) &
#ifdef HIST_TEND
          !$omp shared(advch_t,lhist) &
#endif
          !$omp shared(JJS,JJE,IS,KS,KE) &
          !$omp shared(DENS0,Sr,mflx_hi,DENS_t) &
          !$omp shared(RCDZ,RCDY,MAPF,GSQRT)
          !$acc kernels
          do j = JJS, JJE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, DENS0(k,IS,j) )
             call CHECK( __LINE__, mflx_hi(k  ,IS,j  ,YDIR) )
             call CHECK( __LINE__, mflx_hi(k  ,IS,j-1,YDIR) )
             call CHECK( __LINE__, DENS_t(k,IS,j) )
#endif
             advcv = - ( mflx_hi(k,IS,j,ZDIR)-mflx_hi(k-1,IS,j,  ZDIR) ) * RCDZ(k)
             advch = - ( mflx_hi(k,IS,j,YDIR)-mflx_hi(k  ,IS,j-1,YDIR) ) * RCDY(j)
             tmp = ( advcv + advch ) * MAPF(IS,j,2,I_XY) / GSQRT(k,IS,j,I_XYZ)
             Sr(k,IS,j) = tmp + DENS_t(k,IS,j)
#ifdef HIST_TEND
             if ( lhist ) then
                advch_t(k,IS,j,I_DENS) = tmp
             end if
#endif
          enddo
          enddo
          !$acc end kernels
       else
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,advcv,advch,tmp) &
#ifdef HIST_TEND
          !$omp shared(advch_t,lhist) &
#endif
          !$omp shared(JJS,JJE,IIS,IIE,KS,KE) &
          !$omp shared(DENS0,Sr,mflx_hi,DENS_t) &
          !$omp shared(RCDZ,RCDX,RCDY,MAPF,GSQRT)
          !$acc kernels
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, DENS0(k,i,j) )
             call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, mflx_hi(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, mflx_hi(k  ,i  ,j-1,YDIR) )
             call CHECK( __LINE__, DENS_t(k,i,j) )
#endif
             advcv = -   ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i  ,j,  ZDIR) ) * RCDZ(k)
             advch = - ( ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                       + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) )
             tmp = ( advcv + advch ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
             Sr(k,i,j) = tmp + DENS_t(k,i,j)
#ifdef HIST_TEND
             if ( lhist ) then
                advch_t(k,i,j,I_DENS) = tmp
             end if
#endif
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_sr")

       !##### momentum equation (z) #####

       ! at (x, y, z)
       ! not that z-index is added by -1
       PROFILE_START("hevi_momz_qflxhi_z")
       call ATMOS_DYN_FVM_fluxZ_XYW( qflx_hi(:,:,:,ZDIR), & ! (out)
            MOMZ, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), J33G, & ! (in)
            num_diff(:,:,:,I_MOMZ,ZDIR), & ! (in)
            CDZ, FDZ, dtrk, &
            IIS, IIE, JJS, JJE ) ! (in)
       PROFILE_STOP("hevi_momz_qflxhi_z")

       PROFILE_START("hevi_momz_qflxj")
       if ( .not. TwoD ) &
       call ATMOS_DYN_FVM_fluxJ13_XYW( qflx_J13, & ! (out)
            MOMX, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), J13G(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
            CDZ, TwoD, &
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ23_XYW( qflx_J23, & ! (out)
            MOMY, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), J23G(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
            CDZ, TwoD, &
            IIS, IIE, JJS, JJE ) ! (in)
       PROFILE_STOP("hevi_momz_qflxj")

       ! at (u, y, w)
       if ( .not. TwoD ) then
       PROFILE_START("hevi_momz_qflxhi_x")
       call ATMOS_DYN_FVM_fluxX_XYW( qflx_hi(:,:,:,XDIR), & ! (out)
            MOMX, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), & ! (in)
            num_diff(:,:,:,I_MOMZ,XDIR), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       PROFILE_STOP("hevi_momz_qflxhi_x")
       end if

       ! at (x, v, w)
       PROFILE_START("hevi_momz_qflxhi_y")
       call ATMOS_DYN_FVM_fluxY_XYW( qflx_hi(:,:,:,YDIR), & ! (out)
            MOMY, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), & ! (in)
            num_diff(:,:,:,I_MOMZ,YDIR), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       PROFILE_STOP("hevi_momz_qflxhi_y")

       !--- update momentum(z)
       PROFILE_START("hevi_sw")
       if ( TwoD ) then
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k,advcv,advch,cf,wdmp,div,tmp) &
#ifdef HIST_TEND
          !$omp shared(lhist,advcv_t,advch_t,wdmp_t,ddiv_t) &
#endif
          !$omp shared(JJS,JJE,IS,KS,KE) &
          !$omp shared(qflx_hi,qflx_J23,DDIV,MOMZ0,MOMZ_t,Sw) &
          !$omp shared(RFDZ,RCDY,FDZ,dtrk,wdamp_coef,divdmp_coef) &
          !$omp shared(MAPF,GSQRT)
          !$acc kernels
          do j = JJS, JJE
!OCL NORECURRENCE
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, qflx_hi (k  ,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi (k-1,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_J23(k  ,IS,j)        )
             call CHECK( __LINE__, qflx_J23(k-1,IS,j)        )
             call CHECK( __LINE__, qflx_hi (k  ,IS,j  ,YDIR) )
             call CHECK( __LINE__, qflx_hi (k  ,IS,j-1,YDIR) )
             call CHECK( __LINE__, DDIV(k  ,IS,j) )
             call CHECK( __LINE__, DDIV(k+1,IS,j) )
             call CHECK( __LINE__, MOMZ0(k,IS,j) )
             call CHECK( __LINE__, MOMZ_t(k,IS,j) )
#endif
             advcv = - ( qflx_hi (k,IS,j,ZDIR) - qflx_hi (k-1,IS,j  ,ZDIR) &
                       + qflx_J23(k,IS,j)      - qflx_J23(k-1,IS,j  )      ) * RFDZ(k)
             advch = - ( qflx_hi (k,IS,j,YDIR) - qflx_hi (k,  IS,j-1,YDIR) ) * RCDY(j) &
                   * MAPF(IS,j,2,I_XY)
             wdmp = - wdamp_coef(k) * MOMZ0(k,IS,j)
             div = divdmp_coef / dtrk * ( DDIV(k+1,IS,j)-DDIV(k,IS,j) ) * FDZ(k) ! divergence damping
             Sw(k,IS,j) = ( advcv + advch ) / GSQRT(k,IS,j,I_XYW) &
                        + wdmp + div + MOMZ_t(k,IS,j)
#ifdef HIST_TEND
             if ( lhist ) then
                tmp = 1.0_RP / GSQRT(k,IS,j,I_XYW)
                advcv_t(k,IS,j,I_MOMZ) = advcv * tmp
                advch_t(k,IS,j,I_MOMZ) = advch * tmp
                wdmp_t(k,IS,j) = wdmp
                ddiv_t(k,IS,j,1) = div
             endif
#endif
          enddo
          enddo
          !$acc end kernels
       else
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,advcv,advch,cf,wdmp,div,tmp) &
#ifdef HIST_TEND
          !$omp shared(lhist,advcv_t,advch_t,wdmp_t,ddiv_t) &
#endif
          !$omp shared(JJS,JJE,IIS,IIE,KS,KE) &
          !$omp shared(qflx_hi,qflx_J13,qflx_J23,DDIV,MOMZ0,MOMZ_t,Sw) &
          !$omp shared(RFDZ,RCDX,RCDY,FDZ,dtrk,wdamp_coef,divdmp_coef) &
          !$omp shared(MAPF,GSQRT)
          !$acc kernels
          do j = JJS, JJE
          do i = IIS, IIE
!OCL NORECURRENCE
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi (k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_J13(k  ,i  ,j)        )
             call CHECK( __LINE__, qflx_J13(k-1,i  ,j)        )
             call CHECK( __LINE__, qflx_J23(k  ,i  ,j)        )
             call CHECK( __LINE__, qflx_J23(k-1,i  ,j)        )
             call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_hi (k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_hi (k  ,i  ,j-1,YDIR) )
             call CHECK( __LINE__, DDIV(k  ,i,j) )
             call CHECK( __LINE__, DDIV(k+1,i,j) )
             call CHECK( __LINE__, MOMZ0(k,i,j) )
             call CHECK( __LINE__, MOMZ_t(k,i,j) )
#endif
             advcv = - ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) &
                       + qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j  )      &
                       + qflx_J23(k,i,j)      - qflx_J23(k-1,i  ,j  )      ) * RFDZ(k)
             advch = - ( ( qflx_hi (k,i,j,XDIR) - qflx_hi (k,i-1,j  ,XDIR) ) * RCDX(i) &
                       + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k,i  ,j-1,YDIR) ) * RCDY(j) ) &
                   * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY)
             wdmp = - wdamp_coef(k) * MOMZ0(k,i,j)
             div = divdmp_coef / dtrk * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) * FDZ(k) ! divergence damping
             Sw(k,i,j) = ( advcv + advch ) / GSQRT(k,i,j,I_XYW) &
                       + wdmp + div + MOMZ_t(k,i,j)
#ifdef HIST_TEND
             if ( lhist ) then
                tmp = 1.0_RP / GSQRT(k,i,j,I_XYW)
                advcv_t(k,i,j,I_MOMZ) = advcv * tmp
                advch_t(k,i,j,I_MOMZ) = advch * tmp
                wdmp_t(k,i,j) = wdmp
                ddiv_t(k,i,j,1) = div
             endif
#endif
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
       PROFILE_STOP("hevi_sw")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif


       !##### Thermodynamic Equation #####

       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JJS,JJE,IIS,IIE,KS,KE,JHALO,IHALO,RHOT,DENS,POTT)
       !$acc kernels
       do j = JJS-JHALO, JJE+JHALO
       do i = IIS-IHALO, IIE+IHALO
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, RHOT(k,i,j) )
          call CHECK( __LINE__, DENS(k,i,j) )
#endif
          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (x, y, w)
       call ATMOS_DYN_FVM_fluxZ_XYZ( tflx_hi(:,:,:,ZDIR), & ! (out)
            mflx_hi(:,:,:,ZDIR), POTT, GSQRT(:,:,:,I_XYW), & ! (in)
            num_diff(:,:,:,I_RHOT,ZDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, y, z)
       if ( .not. TwoD ) &
       call ATMOS_DYN_FVM_fluxX_XYZ( tflx_hi(:,:,:,XDIR), & ! (out)
            mflx_hi(:,:,:,XDIR), POTT, GSQRT(:,:,:,I_UYZ), & ! (in)
            num_diff(:,:,:,I_RHOT,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, v, z)
       call ATMOS_DYN_FVM_fluxY_XYZ( tflx_hi(:,:,:,YDIR), & ! (out)
            mflx_hi(:,:,:,YDIR), POTT, GSQRT(:,:,:,I_XVZ), & ! (in)
            num_diff(:,:,:,I_RHOT,YDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)


       PROFILE_START("hevi_st")
       if ( TwoD ) then
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k,advcv,advch,tmp) &
#ifdef HIST_TEND
          !$omp shared(lhist,advcv_t,advch_t) &
#endif
          !$omp shared(JJS,JJE,IS,KS,KE) &
          !$omp shared(tflx_hi,RHOT_t,St,RCDZ,RCDY) &
          !$omp shared(MAPF,GSQRT)
          !$acc kernels
          do j = JJS, JJE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, tflx_hi(k  ,IS,j  ,ZDIR) )
             call CHECK( __LINE__, tflx_hi(k-1,IS,j  ,ZDIR) )
             call CHECK( __LINE__, tflx_hi(k  ,IS,j  ,YDIR) )
             call CHECK( __LINE__, tflx_hi(k  ,IS,j-1,YDIR) )
             call CHECK( __LINE__, RHOT_t(k,IS,j) )
#endif
             advcv = - ( tflx_hi(k,IS,j,ZDIR) - tflx_hi(k-1,IS,j  ,ZDIR) ) * RCDZ(k)
             advch = - ( tflx_hi(k,IS,j,YDIR) - tflx_hi(k  ,IS,j-1,YDIR) ) * RCDY(j)
             tmp = ( advcv + advch ) * MAPF(IS,j,2,I_XY) / GSQRT(k,IS,j,I_XYZ)
             St(k,IS,j) = tmp + RHOT_t(k,IS,j)
#ifdef HIST_TEND
             if ( lhist ) then
                advch_t(k,IS,j,I_RHOT) = tmp
             endif
#endif
          enddo
          enddo
          !$acc end kernels
       else
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,advcv,advch,tmp) &
#ifdef HIST_TEND
          !$omp shared(lhist,advcv_t,advch_t) &
#endif
          !$omp shared(JJS,JJE,IIS,IIE,KS,KE) &
          !$omp shared(tflx_hi,RHOT_t,St,RCDZ,RCDX,RCDY) &
          !$omp shared(MAPF,GSQRT)
          !$acc kernels
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, tflx_hi(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, tflx_hi(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, tflx_hi(k  ,i  ,j-1,YDIR) )
             call CHECK( __LINE__, RHOT_t(k,i,j) )
#endif
             advcv = -   ( tflx_hi(k,i,j,ZDIR) - tflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k)
             advch = - ( ( tflx_hi(k,i,j,XDIR) - tflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                       + ( tflx_hi(k,i,j,YDIR) - tflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) )
             tmp = ( advcv + advch ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
             St(k,i,j) = tmp + RHOT_t(k,i,j)
#ifdef HIST_TEND
             if ( lhist ) then
                advch_t(k,i,j,I_RHOT) = tmp
             endif
#endif
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_st")

       ! implicit solver

       PROFILE_START("hevi_solver")

       call PROF_rapstart("DYN_HEVI", 3)
       NVTX_PUSH("DYN_HEVI")

#ifdef HEVI_FISSION
       ! Note: async(0) is a queue distinct from the synchronous (null) queue,
       !       so there is no implicit ordering between them. Every dependency
       !       between the split kernels below is enforced by an explicit wait.
#endif

       fact = dtrk**2 * J33G

       !OCL INDEPENDENT
!OCL PREFETCH_SEQUENTIAL(SOFT)
#ifdef HEVI_FISSION
       !$omp parallel do default(shared) OMP_SCHEDULE_ &
       !$omp private(k,i,j)
#else
#ifndef __GFORTRAN__
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(k,i,j,ii,l,A,B,pg,advcv,tmp) &
       !$omp private(PT,Ci,Co,F1,F2,F3) &
#ifdef HIST_TEND
       !$omp shared(lhist,pg_t,advcv_t) &
#endif
#ifdef DEBUG
       !$omp shared(RHOT) &
#endif
#if defined DEBUG || defined QUICKDEBUG
       !$omp shared(UNDEF) &
#endif
       !$omp shared(JJS,JJE,IIS,IIE,KA,KMAX,KS,KE) &
       !$omp shared(mflx_hi,tflx_hi,MOMZ_RK,MOMZ0) &
       !$omp shared(DENS_RK,RHOT_RK,DENS0,RHOT0,DENS,MOMZ,POTT,DPRES) &
       !$omp shared(GRAV,dtrk,fact,REF_dens,Sr,Sw,St,RT2P) &
       !$omp shared(ATMOS_DYN_FVM_flux_valueW_Z) &
       !$omp shared(MAPF,GSQRT,J33G,CDZ,RCDZ,RFDZ)
#else
       !$omp parallel do default(shared) private(i,j,k,ii,l) OMP_SCHEDULE_ &
       !$omp private(PT,Ci,Co,F1,F2,F3) &
       !$omp private(A,B,pg,advcv,tmp)
#endif
#endif
       !$acc parallel async(0)
       !$acc loop collapse(2)
#ifndef HEVI_FISSION
       !$acc private(PT,Ci,Co,F1,F2,F3,work)
#endif
       do j = JJS, JJE
#if LSIZE == 1
       do i = IIS, IIE
#else
       do ii = IIS, IIE, LSIZE

#ifndef HEVI_FISSION
#if defined DEBUG || defined QUICKDEBUG
    PT(:,:) = UNDEF
    Ci(:,:) = UNDEF
    Co(:,:) = UNDEF

    F1(:,:) = UNDEF
    F2(:,:) = UNDEF
    F3(:,:) = 0.0_RP
#endif
#endif

          do l = 1, LSIZE
             i = ii + l - 1
             if ( i > IIE ) exit
#endif

             call ATMOS_DYN_FVM_flux_valueW_Z( PT(:,l), & ! (out)
                  MOMZ(:,i,j), POTT(:,i,j), GSQRT(:,i,j,I_XYZ), & ! (in)
                  CDZ )
#ifdef HEVI_FISSION
       enddo ! i
       enddo ! j
       !$acc end parallel

       !$omp parallel do default(shared) OMP_SCHEDULE_ private(k,i,j,A,B,tmp)
       !$acc parallel async(0)
       !$acc loop collapse(3)
       do j = JJS, JJE
       do i = IIS, IIE
#endif
#ifdef _OPENACC
             do k = KS, KE-1
                tmp = fact / GSQRT(k,i,j,I_XYW)
                B = GRAV * tmp / ( CDZ(k+1) + CDZ(k) )
                tmp = tmp * RFDZ(k) * J33G
                A0 = RCDZ(k  ) * RT2P(k  ,i,j) / GSQRT(k  ,i,j,I_XYZ) * tmp
                A1 = RCDZ(k+1) * RT2P(k+1,i,j) / GSQRT(k+1,i,j,I_XYZ) * tmp
                if ( k < KE-1 ) &
                F1(k,l) =        - ( PT(k+1,l) *   A1      + B )
                F2(k,l) = 1.0_RP + ( PT(k  ,l) * ( A1+A0 )     )
                if ( k > KS ) &
                F3(k,l) =        - ( PT(k-1,l) *      A0   - B )
             end do
#else
             do k = KS, KE
                A(k) = RCDZ(k) * RT2P(k,i,j) * J33G / GSQRT(k,i,j,I_XYZ)
             enddo

             ! Note: F3(KS,l) (the sub-diagonal of the first row) and F1(KE-1,l)
             !       (the super-diagonal of the last row) are deliberately left
             !       unset. They are outside the tridiagonal system and none of
             !       the solvers currently used reads them (see the note in
             !       MATRIX_SOLVER_tridiagonal_1D_CR in scale_matrix.F90).
             !       Zeroing them costs a measurable amount of time, so it is
             !       not done. If the solver implementation is changed, check
             !       whether the new one requires them to be zero.
             tmp = fact / GSQRT(KS,i,j,I_XYW)
             B = GRAV / ( CDZ(KS+1) + CDZ(KS) )
             F1(KS,l) =        - ( PT(KS+1,l) * RFDZ(KS) *   A(KS+1)         + B ) * tmp
             F2(KS,l) = 1.0_RP + ( PT(KS  ,l) * RFDZ(KS) * ( A(KS+1)+A(KS) )     ) * tmp
             do k = KS+1, KE-2
                tmp = fact / GSQRT(k,i,j,I_XYW)
                B = GRAV / ( CDZ(k+1) + CDZ(k) )
                F1(k,l) =        - ( PT(k+1,l) * RFDZ(k) *   A(k+1)        + B ) * tmp
                F2(k,l) = 1.0_RP + ( PT(k  ,l) * RFDZ(k) * ( A(k+1)+A(k) )     ) * tmp
                F3(k,l) =        - ( PT(k-1,l) * RFDZ(k) *          A(k)   - B ) * tmp
             enddo
             tmp = fact / GSQRT(KE-1,i,j,I_XYW)
             B = GRAV / ( CDZ(KE) + CDZ(KE-1) )
             F2(KE-1,l) = 1.0_RP + ( PT(KE-1,l) * RFDZ(KE-1) * ( A(KE)+A(KE-1) )    ) * tmp
             F3(KE-1,l) =        - ( PT(KE-2,l) * RFDZ(KE-1) *         A(KE-1)  - B ) * tmp
#endif
#ifdef HEVI_FISSION
          enddo ! i
          enddo ! j
          !$acc end parallel

          !$omp parallel do default(shared) OMP_SCHEDULE_ &
          !$omp private(k,i,j,pg)
          !$acc parallel async(1)
          !$acc loop collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#endif
             do k = KS, KE-1
                ! use not density at the half level but mean density between CZ(k) and CZ(k+1)
                pg = - ( DPRES(k+1,i,j) + RT2P(k+1,i,j)*dtrk*St(k+1,i,j) &
                       - DPRES(k  ,i,j) - RT2P(k  ,i,j)*dtrk*St(k  ,i,j) ) &
                       * RFDZ(k) * J33G / GSQRT(k,i,j,I_XYW) &
                     - GRAV * 0.5_RP &
                       * ( ( DENS(k+1,i,j) - REF_dens(k+1,i,j) + Sr(k+1,i,j) * dtrk ) &
                         + ( DENS(k  ,i,j) - REF_dens(k  ,i,j) + Sr(k  ,i,j) * dtrk ) )
                Ci(k,l) = MOMZ(k,i,j) + dtrk * ( pg + Sw(k,i,j) )
#ifdef HIST_TEND
                if ( lhist ) pg_t(k,i,j,1) = pg
#endif
             enddo
#ifdef HEVI_FISSION
       enddo ! i
       enddo ! j
       !$acc end parallel
       ! wait for F1_work/F2_work/F3_work (async(0)) and Ci_work (async(1))
       !$acc wait

       ! Only the current block (IIS:IIE, JJS:JJE) of the work arrays has been
       ! filled, so the solver must be restricted to it. The lower bounds of the
       ! actual arguments are IS/JS, but they are re-mapped to 1 in the dummy
       ! arguments, hence the IS-1 / JS-1 shift.
       ! Note: the cuSPARSE branch of the solver ignores this range and always
       !       solves the whole plane, which is why cache blocking is rejected
       !       at the top of this subroutine when USE_CUDALIB is set.
       call MATRIX_SOLVER_tridiagonal( KMAX-1, 1,         KMAX-1,   &
                                       IE-IS+1, IIS-IS+1, IIE-IS+1, &
                                       JE-JS+1, JJS-JS+1, JJE-JS+1, &
                                       F1_work(:,:,:), F2_work(:,:,:), F3_work(:,:,:), & ! (in)
                                       Ci_work(:,:,:),                   & ! (in)
                                       Co_work(:,:,:)                    ) ! (out)

       ! wait for Co_work
       !$acc wait

       !$omp parallel do default(shared) OMP_SCHEDULE_ &
       !$omp private(k,i,j,tmp)
       !$acc parallel async(0)
       !$acc loop collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE

#else

#if LSIZE == 1
             call MATRIX_SOLVER_tridiagonal_1D_CR( KMAX-1, 1, KMAX-1, &
#ifdef _OPENACC
                                                   work(:,:), &
#endif
                                                   F1(:,1), F2(:,1), F3(:,1), & ! (in)
                                                   Ci(:,1),                   & ! (in)
                                                   Co(:,1)                    ) ! (out)
#else
          end do

          call MATRIX_SOLVER_tridiagonal( KMAX-1, 1, KMAX-1, &
                                          F1(:,:), F2(:,:), F3(:,:), & ! (in)
                                          Ci(:,:),                   & ! (in)
                                          Co(:,:)                    ) ! (out)

          do l = 1, LSIZE
             i = ii + l - 1
             if ( i > IIE ) exit
#endif

#endif

!OCL NORECURRENCE
             do k = KS, KE-1
#ifdef DEBUG_HEVI2HEVE
                ! for debug (change to explicit integration)
                Co(k,l) = MOMZ(k,i,j)
                tmp = J33G * MOMZ(k,i,j) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) )
                mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) + tmp
                tflx_hi(k,i,j,ZDIR) = tflx_hi(k,i,j,ZDIR) + tmp * PT(k,l)
                ! use not density at the half level but mean density between CZ(k) and CZ(k+1)
                MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                     + dtrk*( &
                     - J33G * ( DPRES(k+1,i,j)-DPRES(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,i_XYW) &
                     - GRAV * 0.5_RP * ( (DENS(k,i,j)-REF_dens(k,i,j)) + (DENS(k+1,i,j)-REF_dens(k+1,i,j)) ) &
                     + Sw(k,i,j) )
#else
                ! z-flux
                tmp = J33G * Co(k,l) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) )
                mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) + tmp
                tflx_hi(k,i,j,ZDIR) = tflx_hi(k,i,j,ZDIR) + tmp * PT(k,l)
                ! z-momentum
                MOMZ_RK(k,i,j) = MOMZ0(k,i,j) + ( Co(k,l) - MOMZ(k,i,j) )
#endif
             enddo
#ifdef HEVI_FISSION
          enddo ! i
          enddo ! j
          !$acc end parallel

          !$omp parallel do default(shared) OMP_SCHEDULE_ &
          !$omp private(i,j)
          !$acc parallel async(0)
          !$acc loop collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#endif
             MOMZ_RK(KS-1,i,j) = 0.0_RP
             MOMZ_RK(KE  ,i,j) = 0.0_RP
#ifdef HEVI_FISSION
          enddo ! i
          enddo ! j
          !$acc end parallel

          !$omp parallel do default(shared) OMP_SCHEDULE_ &
          !$omp private(k,i,j,advcv)
          !$acc parallel async(1)
          !$acc loop collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#endif
             ! density and rho*theta
!OCL NORECURRENCE
             do k = KS+1, KE-1
                advcv = - ( Co(k,l)         - Co(k-1,l) ) &
                      * J33G * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
                DENS_RK(k,i,j) = DENS0(k,i,j) + dtrk * ( advcv + Sr(k,i,j) )
#ifdef HIST_TEND
                if ( lhist ) advcv_t(k,i,j,I_DENS) = advcv
#endif
#ifdef HEVI_FISSION
             enddo ! k
          enddo ! i
          enddo ! j
          !$acc end parallel

          !$omp parallel do default(shared) OMP_SCHEDULE_ &
          !$omp private(k,i,j,advcv)
          !$acc parallel async(2)
          !$acc loop collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
!OCL NORECURRENCE
             do k = KS+1, KE-1
#endif
                advcv = - ( Co(k,l) * PT(k,l) - Co(k-1,l) * PT(k-1,l) ) &
                      * J33G * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
                RHOT_RK(k,i,j) = RHOT0(k,i,j) + dtrk * ( advcv + St(k,i,j) )
#ifdef HIST_TEND
                if ( lhist ) advcv_t(k,i,j,I_RHOT) = advcv
#endif
             enddo
#ifdef HEVI_FISSION
          enddo ! i
          enddo ! j
          !$acc end parallel

          !$omp parallel do default(shared) OMP_SCHEDULE_ &
          !$omp private(i,j,advcv)
          !$acc parallel async(1)
          !$acc loop collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#endif
             advcv = - Co(KS,l)          * J33G * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ) ! Co at KS-1 is 0
             DENS_RK(KS,i,j) = DENS0(KS,i,j) + dtrk * ( advcv + Sr(KS,i,j) )
#ifdef HIST_TEND
             if ( lhist ) advcv_t(KS,i,j,I_DENS) = advcv
#endif
             advcv = Co(KE-1,l)            * J33G * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ) ! Co at KE is 0
             DENS_RK(KE,i,j) = DENS0(KE,i,j) + dtrk * ( advcv + Sr(KE,i,j) )
#ifdef HIST_TEND
             if ( lhist ) advcv_t(KE,i,j,I_DENS) = advcv
#endif
#ifdef HEVI_FISSION
          enddo ! i
          enddo ! j
          !$acc end parallel

          !$omp parallel do default(shared) OMP_SCHEDULE_ &
          !$omp private(i,j,advcv)
          !$acc parallel async(2)
          !$acc loop collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#endif
             advcv = - Co(KS,l) * PT(KS,l) * J33G * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ) ! Co at KS-1 is 0
             RHOT_RK(KS,i,j) = RHOT0(KS,i,j) + dtrk * ( advcv + St(KS,i,j) )
#ifdef HIST_TEND
             if ( lhist ) advcv_t(KS,i,j,I_RHOT) = advcv
#endif
             advcv = Co(KE-1,l) * PT(KE-1,l) * J33G * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ) ! Co at KE is 0
             RHOT_RK(KE,i,j) = RHOT0(KE,i,j) + dtrk * ( advcv + St(KE,i,j) )
#ifdef HIST_TEND
             if ( lhist ) advcv_t(KE,i,j,I_RHOT) = advcv
#endif

#ifdef DEBUG
             call check_equation( &
                  Co(:,l), &
                  DENS(:,i,j), MOMZ(:,i,j), RHOT(:,i,j), DPRES(:,i,j), &
                  REF_dens(:,i,j), &
                  Sr(:,i,j), Sw(:,i,j), St(:,i,j), &
                  J33G, GSQRT(:,i,j,:), &
                  RT2P(:,i,j), &
                  dtrk, i, j )
#endif

#if LSIZE > 1
          end do
#endif

       enddo
       enddo
       !$acc end parallel
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifdef HEVI_FISSION
       ! wait for DENS_RK (async(1)) and RHOT_RK (async(2))
       !$acc wait
#endif

       NVTX_POP()
       call PROF_rapend("DYN_HEVI", 3)

       PROFILE_STOP("hevi_solver")


       !##### momentum equation (x) #####

       PROFILE_START("hevi_momx")
       ! at (u, y, w)
       call ATMOS_DYN_FVM_fluxZ_UYZ( qflx_hi(:,:,:,ZDIR), & ! (out)
            MOMZ, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_UYW), J33G, & ! (in)
            num_diff(:,:,:,I_MOMX,ZDIR), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       if ( .not. TwoD ) &
       call ATMOS_DYN_FVM_fluxJ13_UYZ( qflx_J13, & ! (out)
            MOMX, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_UYZ), J13G(:,:,:,I_UYW), MAPF(:,:,:,I_UY), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ23_UYZ( qflx_J23, & ! (out)
            MOMY, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_UYZ), J23G(:,:,:,I_UYW), MAPF(:,:,:,I_UY), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, y, z)
       ! note that x-index is added by -1
       if ( .not. TwoD ) &
       call ATMOS_DYN_FVM_fluxX_UYZ( qflx_hi(:,:,:,XDIR), & ! (out)
            MOMX, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
            num_diff(:,:,:,I_MOMX,XDIR), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, v, z)
       call ATMOS_DYN_FVM_fluxY_UYZ( qflx_hi(:,:,:,YDIR), & ! (out)
            MOMY, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_UVZ), MAPF(:,:,1,I_UV), & ! (in)
            num_diff(:,:,:,I_MOMX,YDIR), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       !--- update momentum(x)
       if ( TwoD ) then
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k,advch,advcv,cf,tmp) &
#ifdef HIST_TEND
          !$omp shared(lhist,advch_t,advcv_t,pg_t,cf_t,ddiv_t) &
#endif
          !$omp shared(JJS,JJE,IS,KS,KE) &
          !$omp shared(MOMX_RK,DENS,MOMX,MOMY,MOMX0,MOMX_t) &
          !$omp shared(qflx_hi,qflx_J23) &
          !$omp shared(RCDZ,RCDY,CDZ) &
          !$omp shared(MAPF,GSQRT,I_UY,I_UYZ) &
          !$omp shared(dtrk,CORIOLI)
          !$acc kernels
          do j = JJS, JJE
!OCL NORECURRENCE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_hi(k  ,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi(k-1,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,IS,j  ,YDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,IS,j-1,YDIR) )
             call CHECK( __LINE__, CORIOLI(IS,j) )
             call CHECK( __LINE__, MOMY(k,IS,j  ) )
             call CHECK( __LINE__, MOMY(k,IS,j-1) )
             call CHECK( __LINE__, MOMX0(k,IS,j) )
#endif
             advcv = - ( qflx_hi (k,IS,j,ZDIR) - qflx_hi (k-1,IS,j  ,ZDIR) &
                       + qflx_J23(k,IS,j)      - qflx_J23(k-1,IS,j       ) ) * RCDZ(k)
             advch = - ( qflx_hi (k,IS,j,YDIR) - qflx_hi (k  ,IS,j-1,YDIR) ) * RCDY(j) &
                   * MAPF(IS,j,2,I_UY)
             cf = 0.5_RP * CORIOLI(IS,j) * ( MOMY(k,IS,j) + MOMY(k,IS,j-1) )
             MOMX_RK(k,IS,j) = MOMX0(k,IS,j) &
                  + dtrk * ( ( advcv + advch ) / GSQRT(k,IS,j,I_UYZ) + cf + MOMX_t(k,IS,j) )
#ifdef HIST_TEND
             if ( lhist ) then
                tmp = 1.0_RP / GSQRT(k,IS,j,I_UYZ)
                advcv_t(k,IS,j,I_MOMX) = advcv * tmp
                advch_t(k,IS,j,I_MOMX) = advch * tmp
                pg_t(k,IS,j,2) = 0.0_RP
                cf_t(k,IS,j,1) = cf
                ddiv_t(k,IS,j,2) = 0.0_RP
             endif
#endif
          enddo
          enddo
          !$acc end kernels
       else
          iee = min(IIE,IEH)
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,advch,advcv,advc,pg,cf,momy_u,div,f2h1,f2h2,f2h1m,f2h2m,dpresm,dpresk,dpresp,tmp) &
#ifdef HEVI_FISSION
          !$omp shared(pg_work) &
#endif
#ifdef HIST_TEND
          !$omp shared(lhist,advch_t,advcv_t,pg_t,cf_t,ddiv_t) &
#endif
          !$omp shared(JJS,JJE,IIS,iee,KS,KE) &
          !$omp shared(MOMX_RK,DPRES,DENS,MOMX,MOMY,DDIV,MOMX0,MOMX_t) &
          !$omp shared(qflx_hi,qflx_J13,qflx_J23) &
          !$omp shared(RCDZ,RCDY,RFDX,CDZ,FDX) &
          !$omp shared(MAPF,GSQRT,J13G,I_UY,I_UV,I_UYW,I_UYZ) &
          !$omp shared(dtrk,CORIOLI,divdmp_coef)
#ifdef HEVI_FISSION
          !$acc kernels async(0)
#else
          !$acc kernels
#endif
          do j = JJS, JJE
          do i = IIS, iee
!OCL NORECURRENCE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
             call CHECK( __LINE__, DPRES(k,i+1,j) )
             call CHECK( __LINE__, DPRES(k,i  ,j) )
             call CHECK( __LINE__, CORIOLI(i  ,j) )
             call CHECK( __LINE__, CORIOLI(i+1,j) )
             call CHECK( __LINE__, MOMY(k,i  ,j  ) )
             call CHECK( __LINE__, MOMY(k,i+1,j  ) )
             call CHECK( __LINE__, MOMY(k,i  ,j-1) )
             call CHECK( __LINE__, MOMY(k,i+1,j-1) )
             call CHECK( __LINE__, DDIV(k,i+1,j) )
             call CHECK( __LINE__, DDIV(k,i  ,j) )
             call CHECK( __LINE__, MOMX0(k,i,j) )
#endif
             f2h1 = F2H(k,I_UYZ)
             f2h2 = 1.0_RP - f2h1
             f2h1m = F2H(k-1,I_UYZ)
             f2h2m = 1.0_RP - f2h1m
             dpresm = 0.5_RP * ( DPRES(k-1,i+1,j) + DPRES(k-1,i,j) ) ! [x,y,z->u,y,w]
             dpresk = 0.5_RP * ( DPRES(k  ,i+1,j) + DPRES(k  ,i,j) ) ! [x,y,z->u,y,w]
             dpresp = 0.5_RP * ( DPRES(k+1,i+1,j) + DPRES(k+1,i,j) ) ! [x,y,z->u,y,w]
             pg = ( ( GSQRT(k,i+1,j,I_XYZ) * DPRES(k,i+1,j) & ! [x,y,z]
                    - GSQRT(k,i  ,j,I_XYZ) * DPRES(k,i  ,j) & ! [x,y,z]
                    ) * RFDX(i) &
                  + ( J13G(k  ,i,j,I_UYW) &
                    * ( f2h1  * dpresp + f2h2  * dpresk ) & ! [u,y,z->u,y,w]
                    - J13G(k-1,i,j,I_UYW) &
                    * ( f2h1m * dpresk + f2h2m * dpresm ) & ! [u,y,z->u,y,w]
                    ) * RCDZ(k) ) &
                  * MAPF(i,j,1,I_UY)
#ifdef HEVI_FISSION
             pg_work(k,i,j) = pg
          enddo
          enddo
          enddo
          !$acc end kernels

          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,momy_u,cf)
          !$acc kernels async(1)
          do j = JJS, JJE
          do i = IIS, iee
          do k = KS, KE
#endif
             momy_u = ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) * 0.25_RP
             cf = 0.5_RP * ( CORIOLI(i+1,j  )+CORIOLI(i,j  ) ) * momy_u &
                + MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) &
                * ( momy_u      * ( 1.0_RP/MAPF(i+1,j,2,I_XY) - 1.0_RP/MAPF(i,j  ,2,I_XY) ) * RFDX(i) &
                  - MOMX(k,i,j) * ( 1.0_RP/MAPF(i  ,j,1,I_UV) - 1.0_RP/MAPF(i,j-1,1,I_UV) ) * RCDY(j) ) &
                * 2.0_RP * momy_u / ( DENS(k,i+1,j) + DENS(k,i,j) ) ! metric term
#ifdef HEVI_FISSION
             cf_work(k,i,j) = cf
          enddo
          enddo
          enddo
          !$acc end kernels

          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,advcv,advch,advc,tmp)
          !$acc kernels async(2)
          do j = JJS, JJE
          do i = IIS, iee
          do k = KS, KE
#endif
             advcv = -   ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) &
                         + qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j       ) &
                         + qflx_J23(k,i,j)      - qflx_J23(k-1,i  ,j       ) ) * RCDZ(k)
             advch = - ( ( qflx_hi (k,i,j,XDIR) - qflx_hi (k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                       + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                   * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY)
             advc = advcv + advch
#ifdef HIST_TEND
             if ( lhist ) then
                tmp = 1.0_RP / GSQRT(k,i,j,I_UYZ)
                advcv_t(k,i,j,I_MOMX) = advcv * tmp
                advch_t(k,i,j,I_MOMX) = advch * tmp
             end if
#endif
#ifdef HEVI_FISSION
             advc_work(k,i,j) = advc
          enddo
          enddo
          enddo
          !$acc end kernels

          if ( divdmp_coef > 0.0_RP ) then
          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,div)
          !$acc kernels async(3)
          do j = JJS, JJE
          do i = IIS, iee
          do k = KS, KE
#endif
             div = divdmp_coef / dtrk * ( DDIV(k,i+1,j)/MAPF(i+1,j,2,I_XY) - DDIV(k,i,j)/MAPF(i,j,1,I_XY) ) &
                 * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) * FDX(i) ! divergence damping
#ifdef HEVI_FISSION
             div_work(k,i,j) = div
          enddo
          enddo
          enddo
          !$acc end kernels
          end if

          ! wait for pg_work (async(0)), cf_work (async(1)), advc_work (async(2)), and div_work (async(3))
          !$acc wait

          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,pg,cf,advc,div)
          !$acc kernels
          do j = JJS, JJE
          do i = IIS, iee
          do k = KS, KE
             pg = pg_work(k,i,j)
             cf = cf_work(k,i,j)
             advc = advc_work(k,i,j)
             div = div_work(k,i,j)
#endif
             MOMX_RK(k,i,j) = MOMX0(k,i,j) &
                  + dtrk * ( ( advc - pg ) / GSQRT(k,i,j,I_UYZ) + cf + div + MOMX_t(k,i,j) )
#ifdef HIST_TEND
             if ( lhist ) then
                pg_t(k,i,j,2) = - pg / GSQRT(k,i,j,I_UYZ)
                cf_t(k,i,j,1) = cf
                ddiv_t(k,i,j,2) = div
             endif
#endif
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
       PROFILE_STOP("hevi_momx")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (y) #####
       PROFILE_START("hevi_momy")
       ! at (x, v, w)
       call ATMOS_DYN_FVM_fluxZ_XVZ( qflx_hi(:,:,:,ZDIR), & ! (out)
            MOMZ, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_XVW), J33G, & ! (in)
            num_diff(:,:,:,I_MOMY,ZDIR), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       if ( .not. TwoD ) &
       call ATMOS_DYN_FVM_fluxJ13_XVZ( qflx_J13, & ! (out)
            MOMX, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_XVZ), J13G(:,:,:,I_XVW), MAPF(:,:,:,I_XV), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ23_XVZ( qflx_J23, & ! (out)
            MOMY, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_XVZ), J23G(:,:,:,I_XVW), MAPF(:,:,:,I_XV), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, v, z)
       if ( .not. TwoD ) &
       call ATMOS_DYN_FVM_fluxX_XVZ( qflx_hi(:,:,:,XDIR), & ! (out)
            MOMX, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_UV), & ! (in)
            num_diff(:,:,:,I_MOMY,XDIR), & ! (in)
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, y, z)
       ! note that y-index is added by -1
       call ATMOS_DYN_FVM_fluxY_XVZ( qflx_hi(:,:,:,YDIR), & ! (out)
            MOMY, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
            num_diff(:,:,:,I_MOMY,YDIR), & ! (in
            CDZ, TwoD, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       !--- update momentum(y)
       if ( TwoD ) then
          i = IS
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k,advch,advcv,pg,cf,div,f2h1,f2h2,f2h1m,f2h2m,dpresm,dpresk,dpresp,tmp) &
#ifdef HIST_TEND
          !$omp shared(lhist,advch_t,advcv_t,pg_t,cf_t,ddiv_t) &
#endif
          !$omp shared(JJS,JJE,JEH,IS,i,KS,KE) &
          !$omp shared(MOMY_RK,DPRES,DENS,MOMX,MOMY,DDIV,MOMY0,MOMY_t) &
          !$omp shared(qflx_hi,qflx_J23) &
          !$omp shared(RCDZ,RFDY,CDZ,FDY) &
          !$omp shared(MAPF,GSQRT,J23G) &
          !$omp shared(dtrk,CORIOLI,divdmp_coef)
          !$acc kernels
          do j = JJS, min(JJE,JEH)
!OCL NORECURRENCE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_hi(k  ,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi(k-1,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,IS,j  ,YDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,IS,j-1,YDIR) )
             call CHECK( __LINE__, DPRES(k,IS,j  ) )
             call CHECK( __LINE__, DPRES(k,IS,j+1) )
             call CHECK( __LINE__, CORIOLI(IS,j  ) )
             call CHECK( __LINE__, CORIOLI(IS,j+1) )
             call CHECK( __LINE__, MOMX(k,IS,j  ) )
             call CHECK( __LINE__, MOMX(k,IS,j+1) )
             call CHECK( __LINE__, DDIV(k,IS,j+1) )
             call CHECK( __LINE__, DDIV(k,IS,j  ) )
             call CHECK( __LINE__, MOMY_t(k,IS,j) )
             call CHECK( __LINE__, MOMY0(k,IS,j) )
#endif
             advcv = - ( qflx_hi (k,IS,j,ZDIR) - qflx_hi (k-1,IS,j  ,ZDIR) &
                       + qflx_J23(k,IS,j)      - qflx_J23(k-1,IS,j  )      ) * RCDZ(k)
             advch = - ( qflx_hi (k,IS,j,YDIR) - qflx_hi (k  ,IS,j-1,YDIR) ) * RFDY(j) &
                     * MAPF(IS,j,2,I_XV)
             f2h1  = F2H(k  ,I_XVZ)
             f2h2  = 1.0_RP - f2h1
             f2h1m = F2H(k-1,I_XVZ)
             f2h2m = 1.0_RP - f2h1m
             dpresm = 0.5_RP * ( DPRES(k-1,IS,j+1) + DPRES(k-1,IS,j) ) ! [x,y,z->x,v,w]
             dpresk = 0.5_RP * ( DPRES(k  ,IS,j+1) + DPRES(k  ,IS,j) ) ! [x,y,z->x,v,w]
             dpresp = 0.5_RP * ( DPRES(k+1,IS,j+1) + DPRES(k+1,IS,j) ) ! [x,y,z->x,v,w]
             pg = ( ( GSQRT(k,IS,j+1,I_XYZ) * DPRES(k,IS,j+1) & ! [x,y,z]
                    - GSQRT(k,IS,j  ,I_XYZ) * DPRES(k,IS,j  ) & ! [x,y,z]
                    ) * RFDY(j) &
                  + ( J23G(k  ,IS,j,I_XVW) &
                    * ( f2h1  * dpresp + f2h2  * dpresk ) & ! [x,v,z->x,v,w]
                    - J23G(k-1,IS,j,I_XVW) &
                    * ( f2h1m * dpresk + f2h2m * dpresm ) & ! [x,v,z->x,v,w]
                    ) * RCDZ(k) ) &
                  * MAPF(IS,j,2,I_XV)
             cf = - 0.25_RP * ( CORIOLI(  IS,j+1)+CORIOLI(  IS,j) ) & ! [x,y,z->x,v,z]
                            * ( MOMX   (k,IS,j+1)+MOMX   (k,IS,j) )
             div = divdmp_coef / dtrk * ( DDIV(k,IS,j+1) - DDIV(k,IS,j) ) &
                 * MAPF(IS,j,2,I_XV) * FDY(j) ! divergence damping
             MOMY_RK(k,IS,j) = MOMY0(k,IS,j) &
                            + dtrk * ( ( advcv + advch - pg ) / GSQRT(k,IS,j,I_XVZ) + cf + div + MOMY_t(k,IS,j) )
#ifdef HIST_TEND
             if ( lhist ) then
                tmp = 1.0_RP / GSQRT(k,IS,j,I_XVZ)
                advcv_t(k,IS,j,I_MOMY) = advcv * tmp
                advch_t(k,IS,j,I_MOMY) = advch * tmp
                pg_t(k,IS,j,3) = - pg * tmp
                cf_t(k,IS,j,2) = cf
                ddiv_t(k,IS,j,3) = div
             endif
#endif
          enddo
          enddo
          !$acc end kernels
       else
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,advch,advcv,advc,pg,cf,momx_v,div,f2h1,f2h2,f2h1m,f2h2m,dpresm,dpresk,dpresp,tmp) &
#ifdef HEVI_FISSION
          !$omp shared(pg_work) &
#endif
#ifdef HIST_TEND
          !$omp shared(lhist,advch_t,advcv_t,pg_t,cf_t,ddiv_t) &
#endif
          !$omp shared(JJS,JJE,JEH,IIS,IIE,KS,KE) &
          !$omp shared(MOMY_RK,DPRES,DENS,MOMX,MOMY,DDIV,MOMY0,MOMY_t) &
          !$omp shared(qflx_hi,qflx_J13,qflx_J23) &
          !$omp shared(RCDZ,RCDX,RFDY,CDZ,FDY) &
          !$omp shared(MAPF,GSQRT,J23G,I_UV) &
          !$omp shared(dtrk,CORIOLI,divdmp_coef)
#ifdef HEVI_FISSION
          !$acc kernels async(0)
#else
          !$acc kernels
#endif
          do j = JJS, min(JJE,JEH)
          do i = IIS, IIE
!OCL NORECURRENCE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
             call CHECK( __LINE__, DPRES(k,i,j  ) )
             call CHECK( __LINE__, DPRES(k,i,j+1) )
             call CHECK( __LINE__, CORIOLI(i,j  ) )
             call CHECK( __LINE__, CORIOLI(i,j+1) )
             call CHECK( __LINE__, MOMX(k,i  ,j  ) )
             call CHECK( __LINE__, MOMX(k,i  ,j+1) )
             call CHECK( __LINE__, MOMX(k,i-1,j  ) )
             call CHECK( __LINE__, MOMX(k,i-1,j+1) )
             call CHECK( __LINE__, DDIV(k,i,j+1) )
             call CHECK( __LINE__, DDIV(k,i,j  ) )
             call CHECK( __LINE__, MOMY_t(k,i,j) )
             call CHECK( __LINE__, MOMY0(k,i,j) )
#endif
             f2h1  = F2H(k  ,I_XVZ)
             f2h2  = 1.0_RP - f2h1
             f2h1m = F2H(k-1,I_XVZ)
             f2h2m = 1.0_RP - f2h1m
             dpresm = 0.5_RP * ( DPRES(k-1,i,j+1) + DPRES(k-1,i,j) ) ! [x,y,z->x,v,w]
             dpresk = 0.5_RP * ( DPRES(k  ,i,j+1) + DPRES(k  ,i,j) ) ! [x,y,z->x,v,w]
             dpresp = 0.5_RP * ( DPRES(k+1,i,j+1) + DPRES(k+1,i,j) ) ! [x,y,z->x,v,w]
             pg = ( ( GSQRT(k,i,j+1,I_XYZ) * DPRES(k,i,j+1) & ! [x,y,z]
                    - GSQRT(k,i,j  ,I_XYZ) * DPRES(k,i,j  ) & ! [x,y,z]
                    ) * RFDY(j) &
                  + ( J23G(k  ,i,j,I_XVW) &
                    * ( f2h1  * dpresp + f2h2  * dpresk ) & ! [x,v,z->x,v,w]
                    - J23G(k-1,i,j,I_XVW) &
                    * ( f2h1m * dpresk + f2h2m * dpresm ) & ! [x,v,z->x,v,w]
                    ) * RCDZ(k) ) &
                  * MAPF(i,j,2,I_XV)
#ifdef HEVI_FISSION
             pg_work(k,i,j) = pg
          enddo
          enddo
          enddo
          !$acc end kernels

          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,cf,momx_v)
          !$acc kernels async(1)
          do j = JJS, min(JJE,JEH)
          do i = IIS, IIE
          do k = KS, KE
#endif
             momx_v = ( MOMX(k,i,j) + MOMX(k,i-1,j) + MOMX(k,i,j+1) + MOMX(k,i-1,j+1) ) * 0.25_RP
             cf = - 0.5_RP * ( CORIOLI(i  ,j+1)+CORIOLI(i  ,j) ) * momx_v &
                  - MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) &
                  * ( MOMY(k,i,j) * ( 1.0_RP/MAPF(i,j  ,2,I_UV) - 1.0_RP/MAPF(i-1,j,2,I_UV) ) * RCDX(i) &
                    - momx_v      * ( 1.0_RP/MAPF(i,j+1,1,I_XY) - 1.0_RP/MAPF(i  ,j,1,I_XY) ) * RFDY(j) ) &
                  * 2.0_RP * momx_v/ ( DENS(k,i,j+1) + DENS(k,i,j) ) ! metoric term
#ifdef HEVI_FISSION
             cf_work(k,i,j) = cf
          enddo
          enddo
          enddo
          !$acc end kernels

          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,advcv,advch,advc,tmp)
          !$acc kernels async(2)
          do j = JJS, min(JJE,JEH)
          do i = IIS, IIE
          do k = KS, KE
#endif
             advcv = -   ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) &
                         + qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j  )      &
                         + qflx_J23(k,i,j)      - qflx_J23(k-1,i  ,j  )      ) * RCDZ(k)
             advch = - ( ( qflx_hi (k,i,j,XDIR) - qflx_hi (k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                       + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k  ,i  ,j-1,YDIR) ) * RFDY(j) ) &
                     * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV)
             advc = advcv + advch
#ifdef HIST_TEND
             if ( lhist ) then
                tmp = 1.0_RP / GSQRT(k,i,j,I_XVZ)
                advcv_t(k,i,j,I_MOMY) = advcv * tmp
                advch_t(k,i,j,I_MOMY) = advch * tmp
             end if
#endif
#ifdef HEVI_FISSION
             advc_work(k,i,j) = advc
          enddo
          enddo
          enddo
          !$acc end kernels

          if ( divdmp_coef > 0.0_RP ) then
          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,pg,div)
          !$acc kernels async(3)
          do j = JJS, min(JJE,JEH)
          do i = IIS, IIE
          do k = KS, KE
#endif
             div = divdmp_coef / dtrk * ( DDIV(k,i,j+1)/MAPF(i,j+1,1,I_XY) - DDIV(k,i,j)/MAPF(i,j,1,I_XY) ) &
                 * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) * FDY(j) ! divergence damping
#ifdef HEVI_FISSION
             div_work(k,i,j) = div
          enddo
          enddo
          enddo
          !$acc end kernels
          end if

          ! wait for pg_work (async(0)), cf_work (async(1)), advc_work (async(2)), and div_work (async(3))
          !$acc wait

          !$omp parallel do default(shared) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,advc,div,pg,cf)
          !$acc kernels
          do j = JJS, min(JJE,JEH)
          do i = IIS, IIE
          do k = KS, KE
             pg = pg_work(k,i,j)
             cf = cf_work(k,i,j)
             advc = advc_work(k,i,j)
             div = div_work(k,i,j)
#endif
             MOMY_RK(k,i,j) = MOMY0(k,i,j) &
                            + dtrk * ( ( advc - pg ) / GSQRT(k,i,j,I_XVZ) + cf + div + MOMY_t(k,i,j) )
#ifdef HIST_TEND
             if ( lhist ) then
                pg_t(k,i,j,3) = - pg / GSQRT(k,i,j,I_XVZ)
                cf_t(k,i,j,2) = cf
                ddiv_t(k,i,j,3) = div
             endif
#endif
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
       PROFILE_STOP("hevi_momy")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

#ifdef PROFILE_FIPP
       call fipp_stop()
#endif

#ifdef HIST_TEND
    if ( lhist ) then
       call FILE_HISTORY_in(advcv_t(:,:,:,I_DENS), 'DENS_t_advcv', 'tendency of density    (vert. advection) (w/ HIST_TEND)',    'kg/m3/s'   )
       call FILE_HISTORY_in(advcv_t(:,:,:,I_MOMZ), 'MOMZ_t_advcv', 'tendency of momentum z (vert. advection) (w/ HIST_TEND)',    'kg/m2/s2', dim_type='ZHXY')
       call FILE_HISTORY_in(advcv_t(:,:,:,I_MOMX), 'MOMX_t_advcv', 'tendency of momentum x (vert. advection) (w/ HIST_TEND)',    'kg/m2/s2', dim_type='ZXHY')
       call FILE_HISTORY_in(advcv_t(:,:,:,I_MOMY), 'MOMY_t_advcv', 'tendency of momentum y (vert. advection) (w/ HIST_TEND)',    'kg/m2/s2', dim_type='ZXYH')
       call FILE_HISTORY_in(advcv_t(:,:,:,I_RHOT), 'RHOT_t_advcv', 'tendency of rho*theta  (vert. advection) (w/ HIST_TEND)',    'K kg/m3/s' )

       call FILE_HISTORY_in(advch_t(:,:,:,I_DENS), 'DENS_t_advch', 'tendency of density    (horiz. advection) (w/ HIST_TEND)',   'kg/m3/s'   )
       call FILE_HISTORY_in(advch_t(:,:,:,I_MOMZ), 'MOMZ_t_advch', 'tendency of momentum z (horiz. advection) (w/ HIST_TEND)',   'kg/m2/s2', dim_type='ZHXY')
       call FILE_HISTORY_in(advch_t(:,:,:,I_MOMX), 'MOMX_t_advch', 'tendency of momentum x (horiz. advection) (w/ HIST_TEND)',   'kg/m2/s2', dim_type='ZXHY')
       call FILE_HISTORY_in(advch_t(:,:,:,I_MOMY), 'MOMY_t_advch', 'tendency of momentum y (horiz. advection) (w/ HIST_TEND)',   'kg/m2/s2', dim_type='ZXYH')
       call FILE_HISTORY_in(advch_t(:,:,:,I_RHOT), 'RHOT_t_advch', 'tendency of rho*theta  (horiz. advection) (w/ HIST_TEND)',   'K kg/m3/s' )

       call FILE_HISTORY_in(pg_t   (:,:,:,1),      'MOMZ_t_pg',    'tendency of momentum z (pressure gradient) (w/ HIST_TEND)',  'kg/m2/s2', dim_type='ZHXY')
       if ( .not. TwoD ) &
       call FILE_HISTORY_in(pg_t   (:,:,:,2),      'MOMX_t_pg',    'tendency of momentum x (pressure gradient) (w/ HIST_TEND)',  'kg/m2/s2', dim_type='ZXHY')
       call FILE_HISTORY_in(pg_t   (:,:,:,3),      'MOMY_t_pg',    'tendency of momentum y (pressure gradient) (w/ HIST_TEND)',  'kg/m2/s2', dim_type='ZXYH')

       call FILE_HISTORY_in(wdmp_t (:,:,:),        'MOMZ_t_wdamp', 'tendency of momentum z (Rayleigh damping) (w/ HIST_TEND)',   'kg/m2/s2', dim_type='ZHXY')

       call FILE_HISTORY_in(ddiv_t (:,:,:,1),      'MOMZ_t_ddiv',  'tendency of momentum z (divergence damping) (w/ HIST_TEND)', 'kg/m2/s2', dim_type='ZHXY')
       if ( .not. TwoD ) &
       call FILE_HISTORY_in(ddiv_t (:,:,:,2),      'MOMX_t_ddiv',  'tendency of momentum x (divergence damping) (w/ HIST_TEND)', 'kg/m2/s2', dim_type='ZXHY')
       call FILE_HISTORY_in(ddiv_t (:,:,:,3),      'MOMY_t_ddiv',  'tendency of momentum y (divergence damping) (w/ HIST_TEND)', 'kg/m2/s2', dim_type='ZXYH')

       call FILE_HISTORY_in(cf_t   (:,:,:,1),      'MOMX_t_cf',    'tendency of momentum x (coliolis force) (w/ HIST_TEND)',     'kg/m2/s2', dim_type='ZXHY')
       call FILE_HISTORY_in(cf_t   (:,:,:,2),      'MOMY_t_cf',    'tendency of momentum y (coliolis force) (w/ HIST_TEND)',     'kg/m2/s2', dim_type='ZXYH')
    endif
#endif

    !$acc end data

#ifdef HEVI_FISSION
    deallocate( PT_work )
    deallocate( Ci_work )
    deallocate( Co_work )
    deallocate( F1_work )
    deallocate( F2_work )
    deallocate( F3_work )
#endif

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_hevi

#ifdef HEVI_FISSION
#undef PT
#undef Ci
#undef Co
#undef F1
#undef F2
#undef F3
#endif

#ifdef DEBUG
!OCL SERIAL
  subroutine check_equation( &
       VECT, &
       DENS, MOMZ, RHOT, DPRES, &
       REF_dens, &
       Sr, Sw, St, &
       J33G, G, &
       RT2P, &
       dt, i, j )
    use scale_const, only: &
         EPS => CONST_EPS, &
         GRAV => CONST_GRAV
    use scale_prc, only: &
         PRC_abort
    use scale_atmos_grid_cartesC, only: &
         RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
         RFDZ => ATMOS_GRID_CARTESC_RFDZ
    implicit none
    real(RP), intent(in) :: VECT(KS:KE-1)
    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: MOMZ(KA)
    real(RP), intent(in) :: RHOT(KA)
    real(RP), intent(in) :: DPRES(KA)
    real(RP), intent(in) :: REF_dens(KA)
    real(RP), intent(in) :: Sr(KA)
    real(RP), intent(in) :: Sw(KA)
    real(RP), intent(in) :: St(KA)
    real(RP), intent(in) :: J33G
    real(RP), intent(in) :: G(KA,8)
    real(RP), intent(in) :: RT2P(KA)
    real(RP), intent(in) :: dt
    integer , intent(in) :: i
    integer , intent(in) :: j

    real(RP), parameter :: small = 1e-6_RP

    real(RP) :: MOMZ_N(KA)
    real(RP) :: DENS_N(KA)
    real(RP) :: RHOT_N(KA)
    real(RP) :: DPRES_N(KA)

    real(RP) :: POTT(KA)
    real(RP) :: PT(KA)

    real(RP) :: error, lhs, rhs
    integer :: k


    do k = KS, KE-1
       MOMZ_N(k) = VECT(k)
    enddo
    MOMZ_N(:KS-1) = 0.0_RP
    MOMZ_N(KE:) = 0.0_RP

    ! density
    do k = KS+1, KE-1
       DENS_N(k) = DENS(k) &
            + dt * ( - J33G * ( MOMZ_N(k) - MOMZ_N(k-1) ) * RCDZ(k) / G(k,I_XYZ) + Sr(k) )
    enddo
    DENS_N(KS) = DENS(KS) &
         + dt * ( - J33G * MOMZ_N(KS) * RCDZ(KS) / G(KS,I_XYZ) + Sr(KS) )
    DENS_N(KE) = DENS(KE) &
         + dt * ( J33G * MOMZ_N(KE-1) * RCDZ(KE) / G(KE,I_XYZ) + Sr(KE) )

    ! rho*theta
    do k = KS, KE
       POTT(k) = RHOT(k) / DENS(k)
    enddo
    do k = KS+1, KE-2
       PT(k) = ( 7.0_RP * ( POTT(k+1) + POTT(k  ) ) &
                 -        ( POTT(k+2) + POTT(k-1) ) ) / 12.0_RP
    enddo
    PT(KS-1) = 0.0_RP
    PT(KS  ) = ( POTT(KS+1) + POTT(KS  ) ) * 0.5_RP
    PT(KE-1) = ( POTT(KE  ) + POTT(KE-1) ) * 0.5_RP
    PT(KE  ) = 0.0_RP
    do k = KS+1, KE-1
       RHOT_N(k) = RHOT(k) &
            + dt * ( - J33G * ( MOMZ_N(k)*PT(k) - MOMZ_N(k-1)*PT(k-1) ) * RCDZ(k) / G(k,I_XYZ) &
                     + St(k) )
    enddo
    RHOT_N(KS) = RHOT(KS) &
         + dt * ( - J33G * MOMZ_N(KS)*PT(KS) * RCDZ(KS) / G(KS,I_XYZ) + St(KS) )
    RHOT_N(KE) = RHOT(KE) &
         + dt * ( J33G * MOMZ_N(KE-1)*PT(KE-1) * RCDZ(KE) / G(KE-1,I_XYZ) + St(KE) )


    do k = KS, KE
       DPRES_N(k) = DPRES(k) + RT2P(k) * ( RHOT_N(k) - RHOT(k) )
    enddo

    do k = KS, KE
       lhs = ( DENS_N(k) - DENS(k) ) / dt
       rhs = - J33G * ( MOMZ_N(k) - MOMZ_N(k-1) ) * RCDZ(k) / G(k,I_XYZ) + Sr(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       endif
       if ( abs(error) > small ) then
          LOG_ERROR("check_equation",*)"DENS error", k, i, j, error, lhs, rhs
          LOG_ERROR_CONT(*)eps
          call PRC_abort
       endif
    enddo

    do k = KS, KE-1
       lhs = ( MOMZ_N(k) - MOMZ(k) ) / dt
       rhs = - J33G * ( DPRES_N(k+1) - DPRES_N(k) ) * RFDZ(k) / G(k,I_XYW) &
             - GRAV * ( ( DENS_N(k+1) - REF_dens(k+1) ) + ( DENS_N(k) - REF_dens(k) ) ) * 0.5_RP &
             + Sw(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       endif
       if ( abs(error) > small ) then
          LOG_ERROR("check_equation",*)"MOMZ error", k, i, j, error, lhs, rhs
          LOG_ERROR_CONT(*) MOMZ_N(k), MOMZ(k), dt
          LOG_ERROR_CONT(*) - J33G * ( DPRES(k+1) - DPRES(k) ) * RFDZ(k) / G(k,I_XYW) &
             - GRAV * ( ( DENS(k+1) - REF_dens(k+1) ) + ( DENS(k) - REF_dens(k) ) ) * 0.5_RP &
             + Sw(k)
          call PRC_abort
       endif
    enddo

    do k = KS, KE
       lhs = ( RHOT_N(k) - RHOT(k) ) / dt
       rhs = - J33G * ( MOMZ_N(k)*PT(k) - MOMZ_N(k-1)*PT(k-1) ) * RCDZ(k) / G(k,I_XYZ) + St(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       endif
       if ( abs(error) > small ) then
          LOG_ERROR("check_equation",*)"RHOT error", k, i, j, error, lhs, rhs
          call PRC_abort
       endif
    enddo

    return
  end subroutine check_equation
#endif

end module scale_atmos_dyn_tstep_short_fvm_hevi
