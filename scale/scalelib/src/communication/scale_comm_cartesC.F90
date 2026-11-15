!-------------------------------------------------------------------------------
!> module COMMUNICATION
!!
!! @par Description
!!          MPI Communication module for Cartesian C-grid
!!          The array given to COMM_vars, COMM_vars8, COMM_wait, COMM_vars_init and
!!          COMM_vars8_init must always be contiguous: the N/S and corner exchanges use the
!!          consecutive memory from its first elements. This is checked only in debug and
!!          quick-debug builds (SCALE_DEBUG=T, SCALE_QUICKDEBUG=T), not in a normal build.
!!          An array registered with COMM_vars_init or COMM_vars8_init is bound to persistent
!!          requests: it must stay allocated (and, in OpenACC builds, present on the device) until
!!          COMM_vars_finalize or COMM_finalize releases the persistent variables of its grid.
!!          Persistent 8-direction variables exchanged together may be declared as a group with
!!          COMM_vars_group_init (see there for the rule of the consecutive posts).
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
! is_contiguous is not in gfortran before version 9
#if defined(QUICKDEBUG) && !( defined(__GFORTRAN__) && __GNUC__ < 9 )
#define COMM_CHECK_CONTIGUOUS_ENABLED
#define COMM_CHECK_CONTIGUOUS(name,var) call COMM_check_contiguous( name, is_contiguous(var) )
#else
#define COMM_CHECK_CONTIGUOUS(name,var)
#endif
module scale_comm_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use, intrinsic :: iso_fortran_env, only: int64
#ifdef NO_MPI08
  use mpi
#else
  use mpi_f08
#endif
  use iso_c_binding
  use scale_precision
  use scale_io
  use scale_prof
  use scale_tracer
#ifdef _OPENACC
  use openacc
#endif

  use scale_prc_cartesC, only: &
     PRC_next, &
     PRC_W,    &
     PRC_N,    &
     PRC_E,    &
     PRC_S,    &
     PRC_NW,   &
     PRC_NE,   &
     PRC_SW,   &
     PRC_SE,   &
     PRC_HAS_W, &
     PRC_HAS_N, &
     PRC_HAS_E, &
     PRC_HAS_S
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: COMM_setup
  public :: COMM_regist
  public :: COMM_finalize
  public :: COMM_vars_init
  public :: COMM_vars8_init
  public :: COMM_vars_finalize
  public :: COMM_vars
  public :: COMM_vars8
  public :: COMM_wait
  public :: COMM_gather
  public :: COMM_bcast
  public :: COMM_get_vsize_max
  public :: COMM_get_pc_capacity
  public :: COMM_get_pc_info
  public :: COMM_vars_group_init

  interface COMM_address
     module procedure COMM_address_2D
     module procedure COMM_address_3D
  end interface COMM_address

  interface COMM_vars
     module procedure COMM_vars_2D
     module procedure COMM_vars_3D
  end interface COMM_vars

  interface COMM_vars8
     module procedure COMM_vars8_2D
     module procedure COMM_vars8_3D
  end interface COMM_vars8

  interface COMM_wait
     module procedure COMM_wait_2D
     module procedure COMM_wait_3D
  end interface COMM_WAIT

  interface COMM_gather
     module procedure COMM_gather_2D
     module procedure COMM_gather_3D
  end interface COMM_gather

  interface COMM_bcast
     module procedure COMM_bcast_SCR_SP
     module procedure COMM_bcast_SCR_DP
     module procedure COMM_bcast_1D_SP
     module procedure COMM_bcast_1D_DP
     module procedure COMM_bcast_2D_SP
     module procedure COMM_bcast_2D_DP
     module procedure COMM_bcast_3D_SP
     module procedure COMM_bcast_3D_DP
     module procedure COMM_bcast_4D_SP
     module procedure COMM_bcast_4D_DP
     module procedure COMM_bcast_INT_SCR
     module procedure COMM_bcast_INT_1D
     module procedure COMM_bcast_INT_2D
     module procedure COMM_bcast_LOGICAL_SCR
     module procedure COMM_bcast_LOGICAL_1D
     module procedure COMM_bcast_CHARACTER
  end interface COMM_bcast

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: COMM_datatype   !< datatype of variable
  integer, public :: COMM_world      !< communication world ID

  integer, public, parameter :: COMM_PC_GROUP_NONE  = 0 !< state of COMM_get_pc_info: the variable is not in a group
  integer, public, parameter :: COMM_PC_GROUP_CHECK = 1 !< state of COMM_get_pc_info: in a group, only checked (exchanged alone)
  integer, public, parameter :: COMM_PC_GROUP_BATCH = 2 !< state of COMM_get_pc_info: in a group exchanged in batch

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private              :: COMM_vsize_max      !< # limit of communication variables at once
  integer,  private              :: COMM_vsize_max_pc   !< initial capacity of the variables for MPI PC (extended automatically)
  integer,  private              :: COMM_tag_ub         !< upper bound of the message tags (MPI_TAG_UB)
  ! COMM_tag_width must exceed the largest direction offset; the offsets are literals in the exchange routines
  integer,  private, parameter   :: COMM_tag_width = 100 !< # of tags per variable

  logical,  private              :: COMM_IsAllPeriodic  !< periodic boundary condition?

  integer,  private              :: COMM_datasize                !< size of one element [byte]

  logical,  private              :: COMM_USE_MPI_PC       = .true.  !< MPI persistent communication
#ifdef __FUJITSU
  logical,  private              :: COMM_USE_MPI_ONESIDED = .true.  !< MPI one-sided communication
#else
  logical,  private              :: COMM_USE_MPI_ONESIDED = .false. !< MPI one-sided communication
#endif
  !> layout 2 of the persistent 8-direction exchange (issue #3594): the N/S halos are exchanged
  !! with whole rows (i = 1:IA) from and to the variable, the W/E halos with one packed message
  !! in each direction and the corners with one packed block from each diagonal neighbor
  !! (16 requests at most, instead of 12 JHALO + 4). Not used in a grid with PRC_TwoD and IHALO > 0.
  !! The default is set in COMM_setup: .true. with OpenACC, where it shortens the exchange, and
  !! .false. otherwise, where it makes no difference (issue #3615)
  logical,  private              :: COMM_USE_PC_LAYOUT2
  !> batch of the persistent exchanges of a group (COMM_vars_group_init, issue #3594): the members
  !! are packed into the buffer of the group, and one message is sent in each direction for all of
  !! them when the last member is posted (16 requests at most for a group). It requires
  !! COMM_USE_PC_LAYOUT2. Not used in 2-D (PRC_TwoD), where the N/S halos, the only ones exchanged,
  !! are sent from and to the variables without packing: the groups are checked only.
  !! The default is set in COMM_setup: .true. with OpenACC, and .false. otherwise, where the
  !! packing of the members costs more than the messages saved (issue #3615)
  logical,  private              :: COMM_USE_BATCH

  ! directions of the messages of a group, by the direction in which they travel: a message to
  ! the neighbor in N is received from the neighbor in S, and so on
  integer,  private, parameter   :: G_N  = 1
  integer,  private, parameter   :: G_S  = 2
  integer,  private, parameter   :: G_W  = 3
  integer,  private, parameter   :: G_E  = 4
  integer,  private, parameter   :: G_NW = 5
  integer,  private, parameter   :: G_NE = 6
  integer,  private, parameter   :: G_SW = 7
  integer,  private, parameter   :: G_SE = 8
  integer,  private, parameter   :: G_NDIR = 8           !< number of the directions (the tags of a group)
  integer,  private, parameter   :: G_NREQ_MAX = 2 * G_NDIR !< number of the requests of a group at most

#ifdef NO_MPI08
  integer,            private :: COMM_datatype_t
  integer,            private :: COMM_world_t
#else
  type(MPI_Datatype), private :: COMM_datatype_t
  type(MPI_Comm),     private :: COMM_world_t
#endif

  !> group of persistent variables (COMM_vars_group_init). The slot of a member is its position in
  !! the declaration. The buffers and requests exist only in the batch state; the persistent
  !! requests are bound to the buffers, which are never moved (a table of groups extended keeps them)
  type group_t
     character(len=H_MID) :: name = ''
     integer              :: state   = COMM_PC_GROUP_NONE
     integer              :: nmem    = 0       !< number of the members
     integer              :: nposted = 0       !< members posted in the round in progress
     integer              :: nwaited = 0       !< members waited in the round in progress
     logical              :: arrived = .false. !< the messages of the round have arrived (batch)
     integer              :: nreq    = 0       !< number of the requests (batch)
     integer              :: nsend   = 0       !< # of elements sent in a round (batch)
     integer              :: nbuf    = 0       !< size of each buffer (batch)
     integer, allocatable :: seq(:)            !< sequential numbers of the members, by slot
     real(RP), pointer, contiguous :: sendbuf(:) => null() !< buffer of the messages sent (batch)
     real(RP), pointer, contiguous :: recvbuf(:) => null() !< buffer of the messages received (batch)
#ifdef NO_MPI08
     integer                        :: req(G_NREQ_MAX)
#else
     type(MPI_Request)              :: req(G_NREQ_MAX)
#endif
  end type group_t

  type ginfo_t
     integer              :: KA
     integer              :: IA, IS, IE, IHALO
     integer              :: JA, JS, JE, JHALO
     integer              :: nreq_max             !< # limit of communication request at once
     integer              :: size2D_NS4           !< 2D data size (W/E    HALO, 4/8-direction comm.)
     integer              :: size2D_NS8           !< 2D data size (N/S    HALO,   4-direction comm.)
     integer              :: size2D_WE            !< 2D data size (N/S    HALO,   8-direction comm.)
     integer              :: size2D_4C            !< 2D data size (corner HALO,   8-direction comm.)
     integer              :: vars_num = 0         !< numbers of variables for persistent comm.
     logical              :: pc_layout2 = .false. !< the persistent 8-direction variables use the layout 2 (COMM_USE_PC_LAYOUT2)
     integer              :: size_L2 = 0          !< size of the packing buffer of the layout 2 for a vid
     real(RP),    pointer, contiguous :: recvpack_L2(:,:) => null() !< packing buffer of the layout 2 (receive, from W, E and the corners)
     real(RP),    pointer, contiguous :: sendpack_L2(:,:) => null() !< packing buffer of the layout 2 (send,    to   W, E and the corners)
     real(RP),    pointer, contiguous :: recvpack_WE2P(:,:,:) !< packing packet (receive, from W and E)
     real(RP),    pointer, contiguous :: sendpack_P2WE(:,:,:) !< packing packet (send,    to   W and E)
     type(c_ptr), allocatable :: recvbuf_WE(:)    !< receive buffer for MPI_Put (from W and E)
     type(c_ptr), allocatable :: recvbuf_NS(:)    !< receive buffer for MPI_Put (from N and S)
     logical,     allocatable :: put_corner(:)    !< the corners are put in the exchange (8-direction)
     integer,     allocatable :: req_cnt (:)      !< request ID of each MPI send/recv
     integer,     allocatable :: send_cnt(:)      !< # of elements sent by each MPI send/recv
     integer,     allocatable :: preq_cnt (:)     !< request ID of each MPI PC
     integer,     allocatable :: psend_cnt(:)     !< # of elements sent by each MPI PC
     integer,     allocatable :: packid(:)        !< ID of pack
     integer(c_intptr_t), allocatable :: pc_addr(:) !< address of the array registered for MPI PC
     integer,     allocatable :: pc_ndir(:)       !< 4 or 8: the exchange registered for MPI PC (COMM_vars(8)_init)
     integer,     allocatable :: pc_group(:)      !< group of the persistent variable (0: none)
     integer,     allocatable :: pc_slot(:)       !< slot of the persistent variable in its group
     integer,     allocatable :: pc_gstat(:)      !< in the round of its group: 0 not posted, 1 posted, 2 waited
     integer                  :: grp_num = 0      !< number of the groups
     type(group_t), allocatable :: grp(:)         !< groups of persistent variables
     logical                  :: comm_grp_dup = .false. !< comm_grp is made (COMM_USE_BATCH)
#ifdef NO_MPI08
     integer                        :: comm           !< communicator of the grid for non-persistent communication
     integer                        :: comm_pc        !< communicator of the grid for MPI PC
     integer                        :: comm_grp       !< communicator of the grid for the groups in batch
     integer,           allocatable :: req_list(:,:)  !< request ID set of each variables
     integer,           allocatable :: preq_list(:,:) !< request ID set of each variables for MPI PC
     integer,           allocatable :: win_packWE(:)  !< window ID for MPI onesided
     integer,           allocatable :: win_packNS(:)  !< window ID for MPI onesided
#else
     type(MPI_Comm)                 :: comm           !< communicator of the grid for non-persistent communication
     type(MPI_Comm)                 :: comm_pc        !< communicator of the grid for MPI PC
     type(MPI_Comm)                 :: comm_grp       !< communicator of the grid for the groups in batch
     type(MPI_Request), allocatable :: req_list(:,:)  !< request ID set of each variables
     type(MPI_Request), allocatable :: preq_list(:,:) !< request ID set of each variables for MPI PC
     type(MPI_Win),     allocatable :: win_packWE(:)  !< window ID for MPI onesided
     type(MPI_Win),     allocatable :: win_packNS(:)  !< window ID for MPI onesided
#endif
     integer,       allocatable :: packbuf_owner(:)   !< ID of the exchange using the packing buffer (0: not in use)
     integer,       allocatable :: packbuf_rank(:)    !< rank of the array of that exchange (2 or 3)
     integer(c_intptr_t), allocatable :: packbuf_addr(:) !< address of the array of that exchange
  end type ginfo_t

  integer, private, parameter  :: COMM_gid_max = 20
  integer, private             :: COMM_gid = 0
  type(ginfo_t), private       :: ginfo(COMM_gid_max)

#ifdef NO_MPI08
  integer,         private     :: group_packWE !< MPI_Group for pack
  integer,         private     :: group_packNS !< MPI_Group for vars
#else
  type(MPI_Group), private     :: group_packWE !< MPI_Group for pack
  type(MPI_Group), private     :: group_packNS !< MPI_Group for vars
#endif
  logical, private             :: group_packWE_created = .false.
  logical, private             :: group_packNS_created = .false.

  logical, private             :: initialized = .false.

  !> number of the groups whose posts are not complete (the first member is posted, the last is not),
  !! on all the grids: COMM_wait of a member of such a group is not allowed then (COMM_group_wait_check)
  integer, private             :: COMM_grp_posting = 0

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine COMM_setup
    use scale_prc, only: &
       PRC_abort, &
       PRC_LOCAL_COMM_WORLD
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    namelist / PARAM_COMM_CARTESC / &
       COMM_vsize_max, &
       COMM_vsize_max_pc, &
       COMM_USE_MPI_PC, &
       COMM_USE_MPI_ONESIDED, &
       COMM_USE_PC_LAYOUT2, &
       COMM_USE_BATCH

    integer         :: ranks(8)
#ifdef NO_MPI08
    integer         :: group
#else
    type(MPI_Group) :: group
#endif

    integer(kind=MPI_ADDRESS_KIND) :: tag_ub
    logical :: flag

    integer :: n, m
    integer :: ierr
    !---------------------------------------------------------------------------

    if ( initialized ) return

    LOG_NEWLINE
    LOG_INFO("COMM_setup",*) 'Setup'

    COMM_vsize_max = max( 10 + QA*2, 25 )
    COMM_vsize_max_pc = 50 + QA*2

    ! set the namelist items to the default values here (those of COMM_USE_MPI_PC and
    ! COMM_USE_MPI_ONESIDED are the same as in the declaration, which implies SAVE), because this
    ! routine is called repeatedly in the bulk-job iterations
    COMM_USE_MPI_PC = .true.
#ifdef __FUJITSU
    COMM_USE_MPI_ONESIDED = .true.
#else
    COMM_USE_MPI_ONESIDED = .false.
#endif
#ifdef _OPENACC
    COMM_USE_MPI_ONESIDED = .false.
    COMM_USE_PC_LAYOUT2   = .true.
    COMM_USE_BATCH        = .true.
#else
    COMM_USE_PC_LAYOUT2   = .false.
    COMM_USE_BATCH        = .false.
#endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COMM_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("COMM_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("COMM_setup",*) 'Not appropriate names in namelist PARAM_COMM_CARTESC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_COMM_CARTESC)

    ! by the namelist only, so that a configuration is rejected also where persistent communication
    ! is not used (e.g., one-sided communication, the default of some systems)
    if ( COMM_USE_BATCH .and. .not. COMM_USE_PC_LAYOUT2 ) then
       LOG_ERROR("COMM_setup",*) 'COMM_USE_BATCH = .true. requires COMM_USE_PC_LAYOUT2 = .true.'
#ifdef _OPENACC
       LOG_ERROR_CONT(*) 'COMM_USE_BATCH is .true. by default with OpenACC: set it to .false. together with COMM_USE_PC_LAYOUT2.'
#endif
       call PRC_abort
    end if

    if ( PRC_HAS_N .AND. PRC_HAS_S .AND. PRC_HAS_W .AND. PRC_HAS_E ) then
       COMM_IsAllPeriodic = .true.
    else
       COMM_IsAllPeriodic = .false.
    endif

    COMM_world = PRC_LOCAL_COMM_WORLD

    COMM_datasize = storage_size(0.0_RP) / 8

    if ( RP == kind(0.D0) ) then
       COMM_datatype_t = MPI_DOUBLE_PRECISION
    elseif( RP == kind(0.0) ) then
       COMM_datatype_t = MPI_REAL
    else
       LOG_ERROR("COMM_setup",*) 'precision is not supportd'
       call PRC_abort
    endif

#ifdef NO_MPI08
    COMM_world_t = COMM_world
    COMM_datatype = COMM_datatype_t
#else
    COMM_world_t%MPI_VAL = COMM_world
    COMM_datatype = COMM_datatype_t%MPI_VAL
#endif

    COMM_gid = 0

#ifdef _OPENACC
    if ( COMM_USE_MPI_ONESIDED ) then
       LOG_WARN("COMM_setup",*) "Open MPI does not support one-sided APIs with CUDA-aware UCX"
    end if
#endif

    if ( COMM_USE_MPI_ONESIDED ) then

       COMM_USE_MPI_PC = .false.

       call MPI_Comm_group( COMM_world_t, group, ierr )

       n = 0
       if ( PRC_HAS_S ) then
          n = 1
          ranks(n) = PRC_next(PRC_S)
       end if
       if ( PRC_HAS_N ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_N) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_N)
          end if
       end if
       if ( PRC_HAS_N .and. PRC_HAS_W ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_NW) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_NW)
          end if
       else if ( PRC_HAS_N ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_N) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_N)
          end if
       else if ( PRC_HAS_W ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_W) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_W)
          end if
       end if
       if ( PRC_HAS_N .and. PRC_HAS_E ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_NE) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_NE)
          end if
       else if ( PRC_HAS_N ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_N) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_N)
          end if
       else if ( PRC_HAS_E ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_E) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_E)
          end if
       end if
       if ( PRC_HAS_S .and. PRC_HAS_W ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_SW) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_SW)
          end if
       else if ( PRC_HAS_S ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_S) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_S)
          end if
       else if ( PRC_HAS_W ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_W) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_W)
          end if
       end if
       if ( PRC_HAS_S .and. PRC_HAS_E ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_SE) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_SE)
          end if
       else if ( PRC_HAS_S ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_S) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_S)
          end if
       else if ( PRC_HAS_E ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_E) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_E)
          end if
       end if
       if ( n > 0 ) then
          call MPI_Group_incl( group, n, ranks, group_packNS, ierr )
          group_packNS_created = .true.
       else
          group_packNS = MPI_GROUP_EMPTY
          group_packNS_created = .false.
       end if

       n = 0
       if ( .not. PRC_TwoD ) then
          if ( PRC_HAS_W ) then
             n = 1
             ranks(n) = PRC_next(PRC_W)
          end if
          if ( PRC_HAS_E ) then
             if ( n == 0 .or. ranks(1) .ne. PRC_next(PRC_E) ) then
                n = n + 1
                ranks(n) = PRC_next(PRC_E)
             end if
          end if
       end if
       if ( n > 0 ) then
          call MPI_Group_incl( group, n, ranks, group_packWE, ierr )
          group_packWE_created = .true.
       else
          group_packWE = MPI_GROUP_EMPTY
          group_packWE_created = .false.
       end if

       call MPI_Group_free( group, ierr )
    end if

    ! MPI_TAG_UB is an attribute of MPI_COMM_WORLD
    call MPI_Comm_get_attr( MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, flag, ierr )
    COMM_tag_ub = 32767 ! the minimum in the MPI standard
    if ( flag ) COMM_tag_ub = int( min( tag_ub, int( huge(0), MPI_ADDRESS_KIND ) ) )

    LOG_NEWLINE
    LOG_INFO("COMM_setup",*) 'Communication information'
    LOG_INFO_CONT(*)         'Maximum number of vars for one communication: ', COMM_vsize_max
    LOG_INFO_CONT(*)         'All side is periodic?                       : ', COMM_IsAllPeriodic
    LOG_INFO_CONT(*)         'Upper bound of the message tags             : ', COMM_tag_ub
    LOG_INFO_CONT(*)         'Layout 2 of the persistent exchange?        : ', COMM_USE_MPI_PC .and. COMM_USE_PC_LAYOUT2
    LOG_INFO_CONT(*)         'Batch of the groups of persistent variables?: ', COMM_USE_MPI_PC .and. COMM_USE_BATCH &
                                                                               .and. .not. PRC_TwoD
    if ( COMM_USE_MPI_PC .and. COMM_USE_BATCH .and. PRC_TwoD ) then
       LOG_INFO_CONT(*)      'COMM_USE_BATCH is not used in 2-D: the groups are checked only'
    end if

    initialized = .true.

    return
  end subroutine COMM_setup

  !-----------------------------------------------------------------------------
  !> Regist grid
  !! This is collective over the local communicator: all the processes must call it in the same order.
  subroutine COMM_regist( &
       KA, IA, JA, IHALO, JHALO, &
       gid )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    integer, intent(in)  :: KA, IA, JA, IHALO, JHALO
    integer, intent(out) :: gid

    integer :: IMAX, JMAX
    integer :: nreq_NS, nreq_WE, nreq_4C

#ifdef NO_MPI08
    integer        :: win_info
#else
    type(MPI_Info) :: win_info
#endif

    integer(kind=MPI_ADDRESS_KIND) :: size
    integer(int64) :: datasize_k !< size of one column [byte]

    integer :: ierr
    integer :: n

    if ( .not. initialized ) then
       LOG_ERROR("COMM_regist",*) 'COMM_setup must be called before calling COMM_regist'
       call PRC_abort
    end if

    COMM_gid = COMM_gid + 1
    if ( COMM_gid > COMM_gid_max ) then
       LOG_ERROR("COMM_regist",*) 'number of registed grid size exceeds the limit'
       call PRC_abort
    end if
    gid = COMM_gid

    if ( IA < IHALO * 3 ) then
       LOG_ERROR("COMM_regist",*) 'IA must be >= IHALO * 3'
       call PRC_abort
    end if
    if ( JA < JHALO * 3 ) then
       LOG_ERROR("COMM_regist",*) 'JA must be >= JHALO * 3'
       call PRC_abort
    end if
    if ( .not. COMM_USE_MPI_ONESIDED .and. JHALO > 10 ) then ! the direction offsets of the message tags overlap otherwise
       LOG_ERROR("COMM_regist",*) 'JHALO must be <= 10'
       call PRC_abort
    end if

    IMAX = IA - IHALO * 2
    JMAX = JA - JHALO * 2

    ginfo(gid)%KA    = KA
    ginfo(gid)%IA    = IA
    ginfo(gid)%IS    = IHALO + 1
    ginfo(gid)%IE    = IA - IHALO
    ginfo(gid)%IHALO = IHALO
    ginfo(gid)%JA    = JA
    ginfo(gid)%JS    = JHALO + 1
    ginfo(gid)%JE    = JA - JHALO
    ginfo(gid)%JHALO = JHALO

    nreq_NS  = 2 * JHALO !--- send x JHALO, recv x JHALO
    nreq_WE  = 2         !--- send x 1    , recv x 1
    nreq_4C  = 2 * JHALO !--- send x JHALO, recv x JHALO

    if ( COMM_USE_MPI_PC ) then
       ginfo(gid)%nreq_MAX = 2 * nreq_NS + 2 * nreq_WE + 4 * nreq_4C + 1
    else
       ginfo(gid)%nreq_MAX = 2 * nreq_NS + 2 * nreq_WE + 4 * nreq_4C
    end if

    ginfo(gid)%size2D_NS4 = IA   * JHALO
    ginfo(gid)%size2D_NS8 = IMAX
    ginfo(gid)%size2D_WE  = JMAX * IHALO
    ginfo(gid)%size2D_4C  =        IHALO

    allocate( ginfo(gid)%sendpack_P2WE(ginfo(gid)%size2D_WE * KA, 2, COMM_vsize_max) )
    !$acc enter data create(ginfo(gid)%sendpack_P2WE)

    allocate( ginfo(gid)%packbuf_owner(COMM_vsize_max) )
    allocate( ginfo(gid)%packbuf_rank (COMM_vsize_max) )
    allocate( ginfo(gid)%packbuf_addr (COMM_vsize_max) )
    ginfo(gid)%packbuf_owner(:) = 0
    ginfo(gid)%packbuf_rank (:) = 0
    ginfo(gid)%packbuf_addr (:) = 0

    if ( COMM_USE_MPI_ONESIDED ) then

       allocate( ginfo(gid)%recvbuf_WE(COMM_vsize_max) )
       allocate( ginfo(gid)%recvbuf_NS(COMM_vsize_max) )
       allocate( ginfo(gid)%put_corner(COMM_vsize_max) )
       ginfo(gid)%put_corner(:) = .false.

       allocate( ginfo(gid)%win_packWE(COMM_vsize_max) )
       allocate( ginfo(gid)%win_packNS(COMM_vsize_max) )

       call MPI_Info_create(win_info, ierr)
       call MPI_Info_set(win_info, "no_locks", "true", ierr)
       call MPI_Info_set(win_info, "same_size", "true", ierr)
       call MPI_Info_set(win_info, "same_disp_unit", "true", ierr)

       do n = 1, COMM_vsize_max
          size = int( ginfo(gid)%size2D_WE, kind=MPI_ADDRESS_KIND ) * KA * 2 * COMM_datasize
#ifdef _OPENACC
          block
            real(RP), pointer, contiguous :: pack(:)
            call MPI_Alloc_mem(size, MPI_INFO_NULL, ginfo(gid)%recvbuf_WE(n), ierr)
            call C_F_pointer(ginfo(gid)%recvbuf_WE(n), pack, (/ ginfo(gid)%size2D_WE*KA*2 /))
            !$acc enter data create(pack)
            !$acc host_data use_device(pack)
            call MPI_Win_create(pack, size, COMM_datasize, &
                                win_info, COMM_world_t, &
                                ginfo(gid)%win_packWE(n), ierr)
            !$acc end host_data
          end block
#else
          call MPI_Win_allocate(size, COMM_datasize, &
                                win_info, COMM_world_t, &
                                ginfo(gid)%recvbuf_WE(n), ginfo(gid)%win_packWE(n), ierr)
#endif
          size = int( ginfo(gid)%size2D_NS4, kind=MPI_ADDRESS_KIND ) * KA * 2 * COMM_datasize
#ifdef _OPENACC
          block
            real(RP), pointer, contiguous :: pack(:)
            call MPI_Alloc_mem(size, MPI_INFO_NULL, ginfo(gid)%recvbuf_NS(n), ierr)
            call C_F_pointer(ginfo(gid)%recvbuf_NS(n), pack, (/ ginfo(gid)%size2D_NS4*KA*2 /))
            !$acc enter data create(pack)
            !$acc host_data use_device(pack)
            call MPI_Win_create(pack, size, COMM_datasize, &
                                win_info, COMM_world_t, &
                                ginfo(gid)%win_packNS(n), ierr)
            !$acc end host_data
          end block
#else
          call MPI_Win_allocate(size, COMM_datasize, &
                                win_info, COMM_world_t, &
                                ginfo(gid)%recvbuf_NS(n), ginfo(gid)%win_packNS(n), ierr)
#endif
       end do

       call MPI_Info_free(win_info, ierr)

       do n = 1, COMM_vsize_max
          call MPI_Win_post( group_packWE, MPI_MODE_NOSTORE, ginfo(gid)%win_packWE(n), ierr )
          call MPI_Win_post( group_packNS, MPI_MODE_NOSTORE, ginfo(gid)%win_packNS(n), ierr )
       end do

       ginfo(gid)%vars_num = 0

    else

       allocate( ginfo(gid)%recvpack_WE2P(ginfo(gid)%size2D_WE * KA, 2, COMM_vsize_max) )
       !$acc enter data create(ginfo(gid)%recvpack_WE2P)

       allocate( ginfo(gid)%req_cnt (                     COMM_vsize_max) )
       allocate( ginfo(gid)%send_cnt(                     COMM_vsize_max) )
       allocate( ginfo(gid)%req_list(ginfo(gid)%nreq_MAX, COMM_vsize_max) )
       ginfo(gid)%req_cnt (:)   = -1
       ginfo(gid)%send_cnt(:)   = 0
       ginfo(gid)%req_list(:,:) = MPI_REQUEST_NULL

       ! each grid has its own communicators, which separate the message tags between the
       ! grids and between persistent and non-persistent communication
       call MPI_Comm_dup( COMM_world_t, ginfo(gid)%comm, ierr )

       if ( COMM_USE_MPI_PC ) then
          call MPI_Comm_dup( COMM_world_t, ginfo(gid)%comm_pc, ierr )

          ginfo(gid)%vars_num = 0
          call COMM_grow_pc( gid ) ! allocate the arrays for the persistent variables

          ! the groups: the table is allocated at the first declaration and freed in COMM_finalize;
          ! the communicator of the groups in batch is made and freed on the same condition, which
          ! only COMM_setup changes
          if ( allocated( ginfo(gid)%grp ) ) then
             LOG_ERROR("COMM_regist",*) 'the table of the groups is not freed: grid', gid
             call PRC_abort
          end if
          ginfo(gid)%grp_num = 0
          if ( COMM_USE_BATCH ) then
             if ( ginfo(gid)%comm_grp_dup ) then
                LOG_ERROR("COMM_regist",*) 'the communicator of the groups is not freed: grid', gid
                call PRC_abort
             end if
             call MPI_Comm_dup( COMM_world_t, ginfo(gid)%comm_grp, ierr )
             ginfo(gid)%comm_grp_dup = .true.
          end if
       end if

    end if

    ! the layout 2 of the persistent 8-direction exchange. In 2-D (PRC_TwoD) the exchange in x is
    ! skipped, which the whole rows in y would not keep (IHALO > 0): the current layout is used.
    ! The buffers are allocated (with size 0 where they are not used) and freed in COMM_finalize
    ! on the same condition, which only COMM_setup changes
    ginfo(gid)%pc_layout2 = COMM_USE_MPI_PC .and. COMM_USE_PC_LAYOUT2 .and. .not. ( PRC_TwoD .and. IHALO > 0 )
    ginfo(gid)%size_L2 = 0
    if ( COMM_USE_MPI_PC .and. COMM_USE_PC_LAYOUT2 ) then
       if ( associated( ginfo(gid)%sendpack_L2 ) .or. associated( ginfo(gid)%recvpack_L2 ) ) then
          LOG_ERROR("COMM_regist",*) 'the packing buffers of the layout 2 are not freed: grid', gid
          call PRC_abort
       end if
       ! W and E (JA rows at most), and the four corners
       if ( ginfo(gid)%pc_layout2 ) ginfo(gid)%size_L2 = KA * IHALO * ( 2 * JA + 4 * JHALO )
       allocate( ginfo(gid)%sendpack_L2(ginfo(gid)%size_L2, COMM_vsize_max) )
       allocate( ginfo(gid)%recvpack_L2(ginfo(gid)%size_L2, COMM_vsize_max) )
       !$acc enter data create(ginfo(gid)%sendpack_L2, ginfo(gid)%recvpack_L2)

    end if


    datasize_k = int( COMM_datasize, kind=int64 ) * KA

    LOG_NEWLINE
    LOG_INFO("COMM_regist",*) 'Register grid: id=', gid
    LOG_INFO_CONT(*)          'Data size of var (3D,including halo) [byte] : ', datasize_k*IA*JA
    LOG_INFO_CONT(*)          'Data size of halo                    [byte] : ', datasize_k*(2*IA*JHALO+2*JMAX*IHALO)
    LOG_INFO_CONT(*)          'Ratio of halo against the whole 3D grid     : ', real(2*IA*JHALO+2*JMAX*IHALO) / real(IA*JA)
    LOG_INFO_CONT(*)          'Layout 2 of the persistent 8-dir. exchange? : ', ginfo(gid)%pc_layout2

    return
  end subroutine COMM_regist

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine COMM_finalize
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer :: gid
    integer :: i, ierr
    !---------------------------------------------------------------------------

    ! an exchange posted but never waited: its requests or windows are still active
    do gid = 1, COMM_gid
       if ( any( ginfo(gid)%packbuf_owner(:) /= 0 ) ) then
          LOG_ERROR("COMM_finalize",*) 'an exchange was posted but never waited (COMM_wait) on grid', gid
          do i = 1, COMM_vsize_max
             if ( ginfo(gid)%packbuf_owner(i) /= 0 ) then
                LOG_ERROR_CONT(*) 'vid of the packing buffer, ID posted: ', i, ginfo(gid)%packbuf_owner(i)
             end if
          end do
          call PRC_abort
       end if
    end do

    do gid = 1, COMM_gid

       if ( COMM_USE_MPI_ONESIDED ) then

          do i = 1, COMM_vsize_max
             call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(i), ierr )
             call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(i), ierr )
         end do

          do i = 1, COMM_vsize_max
             call MPI_Win_complete( ginfo(gid)%win_packWE(i), ierr )
             call MPI_Win_complete( ginfo(gid)%win_packNS(i), ierr )
          end do

          do i = 1, COMM_vsize_max
             call MPI_Win_wait( ginfo(gid)%win_packWE(i), ierr )
             call MPI_Win_wait( ginfo(gid)%win_packNS(i), ierr )
          end do

          do i = 1, COMM_vsize_max
             call MPI_Win_free(ginfo(gid)%win_packWE(i), ierr)
             call MPI_Win_free(ginfo(gid)%win_packNS(i), ierr)
#ifdef _OPENACC
             block
               real(RP), pointer, contiguous :: pack(:)
               integer :: KA
               KA = ginfo(gid)%KA
               call C_F_pointer( ginfo(gid)%recvbuf_WE(i), pack, (/ginfo(gid)%size2D_WE*KA*2/) )
               !$acc exit data delete(pack)
               call MPI_Free_mem(pack, ierr)
               call C_F_pointer( ginfo(gid)%recvbuf_NS(i), pack, (/ginfo(gid)%size2D_NS4*KA*2/) )
               !$acc exit data delete(pack)
               call MPI_Free_mem(pack, ierr)
             end block
#endif
          end do

          ginfo(gid)%vars_num = 0

          deallocate( ginfo(gid)%win_packWE )
          deallocate( ginfo(gid)%win_packNS )

          deallocate( ginfo(gid)%recvbuf_WE )
          deallocate( ginfo(gid)%recvbuf_NS )
          deallocate( ginfo(gid)%put_corner )

       else

          if ( COMM_USE_MPI_PC ) then

             call COMM_pc_free( "COMM_finalize", gid ) ! persistent variables not released by COMM_vars_finalize
             deallocate( ginfo(gid)%preq_cnt )
             deallocate( ginfo(gid)%psend_cnt )
             deallocate( ginfo(gid)%preq_list )
             deallocate( ginfo(gid)%packid )
             deallocate( ginfo(gid)%pc_addr )
             deallocate( ginfo(gid)%pc_ndir )
             deallocate( ginfo(gid)%pc_group )
             deallocate( ginfo(gid)%pc_slot )
             deallocate( ginfo(gid)%pc_gstat )
             ! the buffers and the requests of the groups are freed in COMM_pc_free
             if ( allocated( ginfo(gid)%grp ) ) then
                do i = 1, size( ginfo(gid)%grp )
                   if ( associated( ginfo(gid)%grp(i)%sendbuf ) .or. associated( ginfo(gid)%grp(i)%recvbuf ) ) then
                      LOG_ERROR("COMM_finalize",*) 'the buffers of a group are not freed: grid, group', gid, i
                      call PRC_abort
                   end if
                end do
                deallocate( ginfo(gid)%grp )
             end if
             ginfo(gid)%grp_num = 0

             call MPI_Comm_free( ginfo(gid)%comm_pc, ierr )

             if ( COMM_USE_BATCH ) then ! the condition in COMM_regist
                call MPI_Comm_free( ginfo(gid)%comm_grp, ierr )
                ginfo(gid)%comm_grp_dup = .false.
             end if

             if ( COMM_USE_PC_LAYOUT2 ) then ! the condition of the allocation in COMM_regist
                !$acc exit data delete(ginfo(gid)%sendpack_L2, ginfo(gid)%recvpack_L2)
                deallocate( ginfo(gid)%sendpack_L2 )
                deallocate( ginfo(gid)%recvpack_L2 )
                nullify( ginfo(gid)%sendpack_L2, ginfo(gid)%recvpack_L2 )
             end if
             ginfo(gid)%pc_layout2 = .false.
             ginfo(gid)%size_L2    = 0

          end if

          call MPI_Comm_free( ginfo(gid)%comm, ierr )

          deallocate( ginfo(gid)%req_cnt )
          deallocate( ginfo(gid)%send_cnt )
          deallocate( ginfo(gid)%req_list )

          !$acc exit data delete(ginfo(gid)%recvpack_WE2P)
          deallocate( ginfo(gid)%recvpack_WE2P )

       end if

       !$acc exit data delete(ginfo(gid)%sendpack_P2WE)
       deallocate( ginfo(gid)%sendpack_P2WE )
       deallocate( ginfo(gid)%packbuf_owner )
       deallocate( ginfo(gid)%packbuf_rank  )
       deallocate( ginfo(gid)%packbuf_addr  )

    end do

    if ( COMM_USE_MPI_ONESIDED ) then
       if ( group_packWE_created ) then
          call MPI_Group_free(group_packWE, ierr)
          group_packWE_created = .false.
       end if
       if ( group_packNS_created ) then
          call MPI_Group_free(group_packNS, ierr)
          group_packNS_created = .false.
       end if
    end if


    COMM_gid = 0
    COMM_grp_posting = 0

    initialized = .false.

    return
  end subroutine COMM_finalize

  !-----------------------------------------------------------------------------
  !> Register a variable for the persistent 4-direction exchange (COMM_vars)
  !! The vid given here selects the packing buffer of the variable, which is shared with the
  !! other persistent variables registered with that vid and the non-persistent exchanges of
  !! that vid: no other exchange with the same vid may be in progress on the grid at the same time.
  !! This is checked in all builds: the buffer is taken by the ID given when an exchange is posted
  !! (COMM_vars, COMM_vars8) and released only by COMM_wait with the same ID, so a post that is not
  !! followed by COMM_wait aborts at the next post with that vid, and a wait with another ID or
  !! without a post aborts. 2D and 3D exchanges share the buffer of a vid.
  !! On return, vid is the ID of the persistent variable (> COMM_vsize_max), or unchanged if the
  !! variable is not registered (persistent communication is not used, or, in OpenACC builds,
  !! the array is not on the device).
  !! The persistent requests are bound to var: it must stay allocated (and, in OpenACC builds,
  !! present on the device) until COMM_vars_finalize or COMM_finalize releases the persistent
  !! variables of the grid. The code that owns the grid calls COMM_vars_finalize before the
  !! components that registered variables free their arrays.
  subroutine COMM_vars_init( &
       varname, &
       var,     &
       vid,     &
       gid      )
    implicit none

    character(len=*), intent(in)    :: varname    !< variable name
    real(RP),         intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID

    integer,          intent(in), optional :: gid
    !---------------------------------------------------------------------------

    call COMM_pc_regist( "COMM_vars_init", varname, var, vid, gid, .false. )

    return
  end subroutine COMM_vars_init

  !-----------------------------------------------------------------------------
  !> Register a variable for the persistent 8-direction exchange (COMM_vars8); see COMM_vars_init
  subroutine COMM_vars8_init( &
       varname, &
       var,     &
       vid,     &
       gid      )
    implicit none

    character(len=*), intent(in)    :: varname    !< variable name
    real(RP),         intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID

    integer,          intent(in), optional :: gid
    !---------------------------------------------------------------------------

    call COMM_pc_regist( "COMM_vars8_init", varname, var, vid, gid, .true. )

    return
  end subroutine COMM_vars8_init

  !-----------------------------------------------------------------------------
  !> Release the persistent variables (COMM_vars_init, COMM_vars8_init) of a grid, or of all the
  !! grids if gid is omitted: their persistent requests are freed, so the registered arrays may be
  !! freed afterwards. The IDs returned by COMM_vars_init and COMM_vars8_init become invalid, and
  !! the next registration on the grid starts again from the first ID. No persistent exchange may
  !! be in progress on the grid. All the processes must call it at the same point.
  !! It does nothing for a grid without persistent variables (e.g., at the second call), and
  !! nothing at all before COMM_setup or after COMM_finalize.
  subroutine COMM_vars_finalize( gid )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in), optional :: gid

    integer :: g
    !---------------------------------------------------------------------------

    if ( COMM_gid == 0 ) return

    if ( present(gid) ) then
       if ( gid < 1 .or. gid > COMM_gid ) then
          LOG_ERROR("COMM_vars_finalize",*) 'gid is invalid', gid, COMM_gid
          call PRC_abort
       end if
       call COMM_pc_free( "COMM_vars_finalize", gid )
    else
       do g = 1, COMM_gid
          call COMM_pc_free( "COMM_vars_finalize", g )
       end do
    end if

    return
  end subroutine COMM_vars_finalize

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_3D(var, vid, gid)
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< atmospheric 3D variable to communication

    integer,  intent(in)    :: vid        !< request ID

    integer,  intent(in), optional :: gid

    integer :: gid_
    integer :: pid
    integer(c_intptr_t) :: addr
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR("COMM_vars_3D",*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    COMM_CHECK_CONTIGUOUS("COMM_vars_3D",var)
    addr = COMM_address( var )
    call COMM_check_vid( "COMM_vars_3D", vid, gid_, .true., 4, addr, pid )
    call COMM_packbuf_acquire( "COMM_vars_3D", gid_, pid, vid, 3, addr )

    if ( vid > COMM_vsize_max ) then
       call PROF_rapstart('COMM_vars_pers', 2)
       call vars_3D_mpi_pc(var, gid_, vid-COMM_vsize_max, pid)
       call PROF_rapend  ('COMM_vars_pers', 2)
    else
       call PROF_rapstart('COMM_vars', 2)
       !$acc wait
       if ( COMM_USE_MPI_ONESIDED ) then
          call vars_3D_mpi_onesided(var, gid_, vid)
       else
          call vars_3D_mpi(var, gid_, vid)
       end if
       call PROF_rapend  ('COMM_vars', 2)
    end if

    return
  end subroutine COMM_vars_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_3D(var, vid, gid)
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(inout) :: var(:,:,:)

    integer,  intent(in)    :: vid

    integer,  intent(in), optional :: gid

    integer :: gid_
    integer :: pid
    integer(c_intptr_t) :: addr
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR("COMM_vars8_3D",*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    COMM_CHECK_CONTIGUOUS("COMM_vars8_3D",var)
    addr = COMM_address( var )
    call COMM_check_vid( "COMM_vars8_3D", vid, gid_, .true., 8, addr, pid )
    call COMM_packbuf_acquire( "COMM_vars8_3D", gid_, pid, vid, 3, addr )

    if ( vid > COMM_vsize_max ) then
       call PROF_rapstart('COMM_vars_pers', 2)
       call vars_3D_mpi_pc(var, gid_, vid-COMM_vsize_max, pid)
       call PROF_rapend  ('COMM_vars_pers', 2)
    else
       call PROF_rapstart('COMM_vars', 2)
       !$acc wait
       if ( COMM_USE_MPI_ONESIDED ) then
          call vars8_3D_mpi_onesided(var, gid_, vid)
       else
          call vars8_3D_mpi(var, gid_, vid)
       end if
       call PROF_rapend  ('COMM_vars', 2)
    end if

    return
  end subroutine COMM_vars8_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_3D(var, vid, FILL_BND, gid)
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(inout) :: var(:,:,:)

    integer, intent(in)     :: vid

    logical, intent(in), optional :: FILL_BND
    integer, intent(in), optional :: gid

    logical :: FILL_BND_
    integer :: gid_
    integer :: pid
    integer(c_intptr_t) :: addr
    !---------------------------------------------------------------------------

    FILL_BND_ = .true.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR("COMM_wait_3D",*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    COMM_CHECK_CONTIGUOUS("COMM_wait_3D",var)
    addr = COMM_address( var )
    call COMM_check_vid( "COMM_wait_3D", vid, gid_, .true., 0, addr, pid )
    if ( COMM_grp_posting > 0 ) call COMM_group_wait_check( "COMM_wait_3D", vid, gid_ )
    call COMM_packbuf_release( "COMM_wait_3D", gid_, pid, vid, 3, addr )

    if ( vid > COMM_vsize_max ) then
       call PROF_rapstart('COMM_wait_pers', 2)
       call wait_3D_mpi_pc(var, gid_, vid-COMM_vsize_max, pid)
    else
       call PROF_rapstart('COMM_wait', 2)
       if ( COMM_USE_MPI_ONESIDED ) then
          call wait_3D_mpi_onesided(var, gid_, vid)
       else
          call wait_3D_mpi(var, gid_, vid)
       end if
    end if

    ! copy inner data to boundary
    if ( .NOT. COMM_IsAllPeriodic ) then
       if ( FILL_BND_ ) then
          ! no barrier: this item is not called in the processes that have all the neighbors
          call PROF_rapstart('COMM_copy_bnd', 3, disable_barrier=.true.)
          call copy_boundary_3D(var, gid_)
          call PROF_rapend  ('COMM_copy_bnd', 3, disable_barrier=.true.)
       end if
    end if

    if ( vid > COMM_vsize_max ) then
       call PROF_rapend  ('COMM_wait_pers', 2)
    else
       call PROF_rapend  ('COMM_wait', 2)
    end if

    return
  end subroutine COMM_wait_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_2D(var, vid, gid)
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(inout) :: var(:,:)

    integer,  intent(in)    :: vid

    integer,  intent(in), optional :: gid

    integer :: gid_
    integer :: pid
    integer(c_intptr_t) :: addr
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR("COMM_vars_2D",*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    COMM_CHECK_CONTIGUOUS("COMM_vars_2D",var)
    addr = COMM_address( var )
    call COMM_check_vid( "COMM_vars_2D", vid, gid_, .false., 4, addr, pid )
    call COMM_packbuf_acquire( "COMM_vars_2D", gid_, pid, vid, 2, addr )

    call PROF_rapstart('COMM_vars', 2)
    !$acc wait
    if ( COMM_USE_MPI_ONESIDED ) then
       call vars_2D_mpi_onesided(var, gid_, vid)
    else
       call vars_2D_mpi(var, gid_, vid)
    end if
    call PROF_rapend  ('COMM_vars', 2)

    return
  end subroutine COMM_vars_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_2D(var, vid, gid)
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(inout) :: var(:,:)

    integer,  intent(in)    :: vid

    integer,  intent(in), optional :: gid

    integer :: gid_
    integer :: pid
    integer(c_intptr_t) :: addr
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR("COMM_vars8_2D",*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    COMM_CHECK_CONTIGUOUS("COMM_vars8_2D",var)
    addr = COMM_address( var )
    call COMM_check_vid( "COMM_vars8_2D", vid, gid_, .false., 8, addr, pid )
    call COMM_packbuf_acquire( "COMM_vars8_2D", gid_, pid, vid, 2, addr )

    call PROF_rapstart('COMM_vars', 2)
    !$acc wait
    if ( COMM_USE_MPI_ONESIDED ) then
       call vars8_2D_mpi_onesided(var, gid_, vid)
    else
       call vars8_2D_mpi(var, gid_, vid)
    end if
    call PROF_rapend  ('COMM_vars', 2)

    return
  end subroutine COMM_vars8_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_2D(var, vid, FILL_BND, gid)
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(inout) :: var(:,:)

    integer,  intent(in)    :: vid

    logical,  intent(in), optional :: FILL_BND
    integer,  intent(in), optional :: gid

    logical :: FILL_BND_
    integer :: gid_
    integer :: pid
    integer(c_intptr_t) :: addr
    !---------------------------------------------------------------------------

    FILL_BND_ = .true.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR("COMM_wait_2D",*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    COMM_CHECK_CONTIGUOUS("COMM_wait_2D",var)
    addr = COMM_address( var )
    call COMM_check_vid( "COMM_wait_2D", vid, gid_, .false., 0, addr, pid )
    if ( COMM_grp_posting > 0 ) call COMM_group_wait_check( "COMM_wait_2D", vid, gid_ )
    call COMM_packbuf_release( "COMM_wait_2D", gid_, pid, vid, 2, addr )

    call PROF_rapstart('COMM_wait', 2)
    if ( COMM_USE_MPI_ONESIDED ) then
       call wait_2D_mpi_onesided(var, gid_, vid)
    else
       call wait_2D_mpi(var, gid_, vid)
    end if

    if( .NOT. COMM_IsAllPeriodic ) then
       if ( FILL_BND_ ) then
          ! no barrier: this item is not called in the processes that have all the neighbors
          call PROF_rapstart('COMM_copy_bnd', 3, disable_barrier=.true.)
          call copy_boundary_2D(var, gid_)
          call PROF_rapend  ('COMM_copy_bnd', 3, disable_barrier=.true.)
       end if
    end if
    call PROF_rapend  ('COMM_wait', 2)

    return
  end subroutine COMM_wait_2D

  !-----------------------------------------------------------------------------
  !> Get data from whole process value in 2D field
  subroutine COMM_gather_2D( &
       IA, JA, &
       send, &
       recv  )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)  :: IA, JA      !< dimension size
    real(RP), intent(in)  :: send(IA,JA) !< send buffer

    real(RP), intent(out) :: recv(:,:,:) !< receive buffer (IA,JA,nprcs)

    integer :: sendcounts, recvcounts
    integer :: ierr
    !---------------------------------------------------------------------------

    sendcounts = IA * JA
    recvcounts = IA * JA

    !$acc host_data use_device(send, recv) if(acc_is_present(send))
    call PROF_rapstart('COMM_Gather', 2)
    call MPI_GATHER( send(:,:),       &
                     sendcounts,      &
                     COMM_datatype_t, &
                     recv(:,:,:),     &
                     recvcounts,      &
                     COMM_datatype_t, &
                     PRC_masterrank,  &
                     COMM_world_t,    &
                     ierr             )
    call PROF_rapend  ('COMM_Gather', 2)
    !$acc end host_data

    return
  end subroutine COMM_gather_2D

  !-----------------------------------------------------------------------------
  !> Get data from whole process value in 3D field
  subroutine COMM_gather_3D( &
       KA, IA, JA, &
       send, &
       recv  )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)  :: KA, IA, JA     !< dimension size
    real(RP), intent(in)  :: send(KA,IA,JA) !< send buffer

    real(RP), intent(out) :: recv(:,:,:,:) !< receive buffer(KA,IA,JA,nprcs)

    integer :: sendcounts, recvcounts
    integer :: ierr
    !---------------------------------------------------------------------------

    sendcounts = KA * IA * JA
    recvcounts = KA * IA * JA

    !$acc host_data use_device(send, recv) if(acc_is_present(send))
    call PROF_rapstart('COMM_Gather', 2)
    call MPI_GATHER( send(:,:,:),     &
                     sendcounts,      &
                     COMM_datatype_t, &
                     recv(:,:,:,:),   &
                     recvcounts,      &
                     COMM_datatype_t, &
                     PRC_masterrank,  &
                     COMM_world_t,    &
                     ierr             )
    call PROF_rapend  ('COMM_Gather', 2)
    !$acc end host_data

    return
  end subroutine COMM_gather_3D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in scalar field
  subroutine COMM_bcast_SCR_SP( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    real(SP), intent(inout) :: var  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_SCR_SP
  subroutine COMM_bcast_SCR_DP( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    real(DP), intent(inout) :: var  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_SCR_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 1D field
  subroutine COMM_bcast_1D_SP( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA       !< dimension size

    real(SP), intent(inout) :: var(IA)  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:),         &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_1D_SP
  subroutine COMM_bcast_1D_DP( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA       !< dimension size

    real(DP), intent(inout) :: var(IA)  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:),         &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_1D_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 2D field
  subroutine COMM_bcast_2D_SP( IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA, JA     !< dimension size

    real(SP), intent(inout) :: var(IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:),       &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_2D_SP
  subroutine COMM_bcast_2D_DP( IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA, JA     !< dimension size

    real(DP), intent(inout) :: var(IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:),       &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_2D_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 3D field
  subroutine COMM_bcast_3D_SP( KA, IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA    !< dimension size

    real(SP), intent(inout) :: var(KA,IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = KA * IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:,:),     &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_3D_SP
  subroutine COMM_bcast_3D_DP( KA, IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA    !< dimension size

    real(DP), intent(inout) :: var(KA,IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = KA * IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:,:),     &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_3D_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 4D field
  subroutine COMM_bcast_4D_SP( KA, IA, JA, NT, var )
    use scale_prc, only: &
       PRC_abort, &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA, NT   !< dimension size

    real(SP), intent(inout) :: var(KA,IA,JA,NT) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = KA * IA * JA * NT
    if ( KA>0 .AND. IA>0 .AND. JA>0 .AND. NT>0 .AND. &
         counts < 0 ) then
       LOG_ERROR("COMM_bcast_4D",*) 'counts overflow'
       call PRC_abort
    end if

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:,:,:),   &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_4D_SP
  subroutine COMM_bcast_4D_DP( KA, IA, JA, NT, var )
    use scale_prc, only: &
       PRC_abort, &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA, NT   !< dimension size

    real(DP), intent(inout) :: var(KA,IA,JA,NT) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = KA * IA * JA * NT
    if ( KA>0 .AND. IA>0 .AND. JA>0 .AND. NT>0 .AND. &
         counts < 0 ) then
       LOG_ERROR("COMM_bcast_4D",*) 'counts overflow'
       call PRC_abort
    end if

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:,:,:),   &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_4D_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in scalar (integer)
  subroutine COMM_bcast_INT_SCR( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer, intent(inout) :: var !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_INTEGER,    &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_INT_SCR

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 1D field (integer)
  subroutine COMM_bcast_INT_1D( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer, intent(in)    :: IA      !< dimension size
    integer, intent(inout) :: var(IA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    call MPI_BCAST( var(:),         &
                    counts,         &
                    MPI_INTEGER,    &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_INT_1D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 2D field (integer)
  subroutine COMM_bcast_INT_2D( IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer, intent(in)    :: IA, JA     !< dimension size

    integer, intent(inout) :: var(IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:),       &
                    counts,         &
                    MPI_INTEGER,    &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_INT_2D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in scalar (logical)
  subroutine COMM_bcast_LOGICAL_SCR( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    logical, intent(inout) :: var   !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_LOGICAL,    &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_LOGICAL_SCR

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 1D (logical)
  subroutine COMM_bcast_LOGICAL_1D( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer, intent(in)    :: IA      !< dimension size
    logical, intent(inout) :: var(IA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:),         &
                    counts,         &
                    MPI_LOGICAL,    &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_LOGICAL_1D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in character
  subroutine COMM_bcast_CHARACTER( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    character(len=*), intent(inout) :: var   !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = len(var)

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_CHARACTER,  &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_CHARACTER

!-------------------------------------------------------------------------------
! private routines
!-------------------------------------------------------------------------------
  !> Base of the message tags of a variable: the tag of a message is id * 100 + the offset of
  !> the direction (0 to 59 for the corners and N/S, 60 for E, 70 for W in 8 directions, 1 to 4 in
  !> 4 directions). The whole width of 100 is reserved for a variable and checked against
  !> MPI_TAG_UB. id is vid in a non-persistent exchange and the sequential number for a persistent
  !> variable. The grid and the kind of exchange are identified by the communicator.
  pure function COMM_tag_base( id ) result( tag )
    implicit none
    integer, intent(in) :: id
    integer :: tag

    tag = id * COMM_tag_width

    return
  end function COMM_tag_base

#ifdef COMM_CHECK_CONTIGUOUS_ENABLED
  !-----------------------------------------------------------------------------
  !> Stop if the array given to the communication is not contiguous
  subroutine COMM_check_contiguous( routine_name, contiguous )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: routine_name
    logical,          intent(in) :: contiguous
    !---------------------------------------------------------------------------

    if ( .not. contiguous ) then
       LOG_ERROR(routine_name,*) 'the array must be contiguous: copy it to a work array, if it is a section with a stride'
       call PRC_abort
    end if

    return
  end subroutine COMM_check_contiguous
#endif

  !-----------------------------------------------------------------------------
  !> Address of the first element of an array, as an integer (0 for an empty array)
  !! (an integer, not a c_ptr, is kept: the actual argument usually lacks the target attribute)
  function COMM_address_3D( var ) result(addr)
    implicit none

    real(RP), intent(in), target :: var(:,:,:)
    integer(c_intptr_t) :: addr
    !---------------------------------------------------------------------------

    addr = 0
    if ( size(var) > 0 ) addr = transfer( c_loc(var(1,1,1)), 0_c_intptr_t )

    return
  end function COMM_address_3D

  function COMM_address_2D( var ) result(addr)
    implicit none

    real(RP), intent(in), target :: var(:,:)
    integer(c_intptr_t) :: addr
    !---------------------------------------------------------------------------

    addr = 0
    if ( size(var) > 0 ) addr = transfer( c_loc(var(1,1)), 0_c_intptr_t )

    return
  end function COMM_address_2D

  !-----------------------------------------------------------------------------
  !> Check the variable ID given to an exchange or a wait routine, for a valid gid, and return the
  !> ID of its packing buffer (pid).
  !> An ID above COMM_vsize_max is checked against the number of persistent variables initialized
  !> on the grid, against the direction of its registration (4 or 8, for an exchange), and against
  !> the address of the array registered with the ID:
  !> this detects an ID obtained on another grid and a non-persistent vid in the range of
  !> COMM_vsize_max + 1 to COMM_vsize_max + that number, unless the very array registered with
  !> that ID is given. Not detected: such a vid given with the registered array, which is then
  !> exchanged as the persistent variable, and an ID of another grid on which the same array is
  !> registered with the same sequential number, which is then exchanged as the persistent
  !> variable of the grid given.
  subroutine COMM_check_vid( routine_name, vid, gid, allow_persistent, ndir, addr, pid )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: routine_name
    integer,          intent(in) :: vid
    integer,          intent(in) :: gid
    logical,          intent(in) :: allow_persistent !< the ID of a persistent variable is allowed
    integer,          intent(in) :: ndir             !< 4 or 8: the direction of the exchange posted (0 for a wait)
    integer(c_intptr_t), intent(in) :: addr          !< address of the array given with the ID
    integer,          intent(out) :: pid             !< ID of the packing buffer
    !---------------------------------------------------------------------------

    if ( vid < 1 ) then
       LOG_ERROR(routine_name,*) 'vid is invalid', vid, COMM_vsize_max, gid
       call PRC_abort
    else if ( vid <= COMM_vsize_max ) then
       if ( .not. COMM_USE_MPI_ONESIDED ) then ! one-sided communication uses no tags
          if ( COMM_tag_base( vid ) + COMM_tag_width - 1 > COMM_tag_ub ) then
             LOG_ERROR(routine_name,*) 'message tags exceed MPI_TAG_UB: ID, MPI_TAG_UB = ', vid, COMM_tag_ub
             call PRC_abort
          end if
       end if
       pid = vid
    else
       if ( .not. allow_persistent ) then
          LOG_ERROR(routine_name,*) 'vid is invalid', vid, COMM_vsize_max, gid
          call PRC_abort
       else if ( .not. COMM_USE_MPI_PC ) then
          LOG_ERROR(routine_name,*) 'vid is invalid: persistent communication is not used', vid, COMM_vsize_max, gid
          call PRC_abort
       else if ( vid - COMM_vsize_max > ginfo(gid)%vars_num ) then
          LOG_ERROR(routine_name,*) 'vid is invalid: not initialized', vid, COMM_vsize_max, ginfo(gid)%vars_num, gid
          call PRC_abort
       else if ( addr /= ginfo(gid)%pc_addr(vid-COMM_vsize_max) ) then
          LOG_ERROR(routine_name,*) 'vid is invalid: the array is not the one registered with this ID', vid, gid
          LOG_ERROR_CONT(*) '(an ID of another grid, or a non-persistent vid above COMM_vsize_max)'
          call PRC_abort
       else if ( ndir /= 0 .and. ndir /= ginfo(gid)%pc_ndir(vid-COMM_vsize_max) ) then
          LOG_ERROR(routine_name,*) 'vid is invalid: the persistent variable is registered for the other direction', vid, gid
          LOG_ERROR_CONT(*) 'directions of the registration and of the exchange: ', &
                            ginfo(gid)%pc_ndir(vid-COMM_vsize_max), ndir
          call PRC_abort
       end if
       pid = ginfo(gid)%packid(vid-COMM_vsize_max)
    end if

    return
  end subroutine COMM_check_vid

  !-----------------------------------------------------------------------------
  !> Allocate the arrays for the persistent variables of the grid (COMM_vsize_max_pc), or extend
  !! them to twice the size. The request handles are copied as values;
  !! the persistent requests are bound to the registered arrays and the packing buffers,
  !! which are not moved.
  subroutine COMM_grow_pc( gid )
    implicit none

    integer, intent(in) :: gid

    integer,           allocatable :: packid(:)
    integer(c_intptr_t), allocatable :: pc_addr(:)
    integer,           allocatable :: pc_ndir  (:)
    integer,           allocatable :: pc_group (:)
    integer,           allocatable :: pc_slot  (:)
    integer,           allocatable :: pc_gstat (:)
    integer,           allocatable :: preq_cnt (:)
    integer,           allocatable :: psend_cnt(:)
#ifdef NO_MPI08
    integer,           allocatable :: preq_list(:,:)
#else
    type(MPI_Request), allocatable :: preq_list(:,:)
#endif

    integer :: nold, nnew
    !---------------------------------------------------------------------------

    if ( allocated( ginfo(gid)%packid ) ) then
       nold = size(ginfo(gid)%packid)
       nnew = max( nold * 2, 1 )
    else
       nold = 0
       nnew = max( COMM_vsize_max_pc, 1 )
    end if

    allocate( packid   (                      nnew) )
    allocate( pc_addr  (                      nnew) )
    allocate( pc_ndir  (                      nnew) )
    allocate( pc_group (                      nnew) )
    allocate( pc_slot  (                      nnew) )
    allocate( pc_gstat (                      nnew) )
    allocate( preq_cnt (                      nnew) )
    allocate( psend_cnt(                      nnew) )
    allocate( preq_list(ginfo(gid)%nreq_MAX+1,nnew) )

    if ( nold > 0 ) then
       packid   (     1:nold) = ginfo(gid)%packid   (     1:nold)
       pc_addr  (     1:nold) = ginfo(gid)%pc_addr  (     1:nold)
       pc_ndir  (     1:nold) = ginfo(gid)%pc_ndir  (     1:nold)
       pc_group (     1:nold) = ginfo(gid)%pc_group (     1:nold)
       pc_slot  (     1:nold) = ginfo(gid)%pc_slot  (     1:nold)
       pc_gstat (     1:nold) = ginfo(gid)%pc_gstat (     1:nold)
       preq_cnt (     1:nold) = ginfo(gid)%preq_cnt (     1:nold)
       psend_cnt(     1:nold) = ginfo(gid)%psend_cnt(     1:nold)
       preq_list(:,   1:nold) = ginfo(gid)%preq_list(:,   1:nold)
    end if

    ! initial values (no request: COMM_pc_free frees only the requests made)
    packid   (  nold+1:nnew) = 0
    pc_addr  (  nold+1:nnew) = 0
    pc_ndir  (  nold+1:nnew) = 0
    pc_group (  nold+1:nnew) = 0
    pc_slot  (  nold+1:nnew) = 0
    pc_gstat (  nold+1:nnew) = 0
    preq_cnt (  nold+1:nnew) = -1
    psend_cnt(  nold+1:nnew) = 0
    preq_list(:,nold+1:nnew) = MPI_REQUEST_NULL

    call move_alloc( packid,    ginfo(gid)%packid    )
    call move_alloc( pc_addr,   ginfo(gid)%pc_addr   )
    call move_alloc( pc_ndir,   ginfo(gid)%pc_ndir   )
    call move_alloc( pc_group,  ginfo(gid)%pc_group  )
    call move_alloc( pc_slot,   ginfo(gid)%pc_slot   )
    call move_alloc( pc_gstat,  ginfo(gid)%pc_gstat  )
    call move_alloc( preq_cnt,  ginfo(gid)%preq_cnt  )
    call move_alloc( psend_cnt, ginfo(gid)%psend_cnt )
    call move_alloc( preq_list, ginfo(gid)%preq_list )

    if ( nold > 0 ) then
       LOG_INFO("COMM_grow_pc",*) 'The capacity of the variables for MPI PC is extended (grid ID, size): ', gid, nnew
    end if

    return
  end subroutine COMM_grow_pc

  !-----------------------------------------------------------------------------
  !> Take the packing buffer pid for the exchange posted with the ID vid (pid itself, or
  !! COMM_vsize_max + the sequential number of a persistent variable registered with pid) of an
  !! array of rank rank (2 or 3) at address addr
  subroutine COMM_packbuf_acquire( routine_name, gid, pid, vid, rank, addr )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: routine_name
    integer,          intent(in) :: gid
    integer,          intent(in) :: pid  !< ID of the packing buffer
    integer,          intent(in) :: vid  !< ID given to the exchange
    integer,          intent(in) :: rank !< rank of the array (2 or 3)
    integer(c_intptr_t), intent(in) :: addr !< address of the array
    !---------------------------------------------------------------------------

    if ( ginfo(gid)%packbuf_owner(pid) /= 0 ) then
       LOG_ERROR(routine_name,*) 'packing buffer of vid is already in use on grid: ', pid, gid
       LOG_ERROR_CONT(*) 'another exchange with the same vid is in progress'
       LOG_ERROR_CONT(*) 'ID of the exchange in progress, ID posted: ', ginfo(gid)%packbuf_owner(pid), vid
       call PRC_abort
    end if
    ginfo(gid)%packbuf_owner(pid) = vid
    ginfo(gid)%packbuf_rank (pid) = rank
    ginfo(gid)%packbuf_addr (pid) = addr

    return
  end subroutine COMM_packbuf_acquire

  !-----------------------------------------------------------------------------
  !> Release the packing buffer pid at the wait with the ID vid of an array of rank rank at address
  !! addr: the exchange in progress with the buffer must be the one posted with the same ID and array
  subroutine COMM_packbuf_release( routine_name, gid, pid, vid, rank, addr )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: routine_name
    integer,          intent(in) :: gid
    integer,          intent(in) :: pid  !< ID of the packing buffer
    integer,          intent(in) :: vid  !< ID given to the wait
    integer,          intent(in) :: rank !< rank of the array (2 or 3)
    integer(c_intptr_t), intent(in) :: addr !< address of the array
    !---------------------------------------------------------------------------

    if ( ginfo(gid)%packbuf_owner(pid) == 0 ) then
       LOG_ERROR(routine_name,*) 'wait without a matching post: no exchange with this ID is in progress on grid', vid, gid
       call PRC_abort
    else if ( ginfo(gid)%packbuf_owner(pid) /= vid ) then
       LOG_ERROR(routine_name,*) 'packing buffer of vid is in use by another ID: ', pid, gid
       LOG_ERROR_CONT(*) 'the wait does not match the exchange in progress'
       LOG_ERROR_CONT(*) 'ID of the exchange in progress, ID waited: ', ginfo(gid)%packbuf_owner(pid), vid
       call PRC_abort
    else if ( ginfo(gid)%packbuf_rank(pid) /= rank ) then
       LOG_ERROR(routine_name,*) 'the rank of the array waited differs from that of the exchange in progress (2D/3D)', pid, gid
       LOG_ERROR_CONT(*) 'rank of the array posted, rank of the array waited: ', ginfo(gid)%packbuf_rank(pid), rank
       call PRC_abort
    else if ( ginfo(gid)%packbuf_addr(pid) /= addr ) then
       LOG_ERROR(routine_name,*) 'the array waited is not the array posted with this ID', vid, gid
       call PRC_abort
    end if
    ginfo(gid)%packbuf_owner(pid) = 0
    ginfo(gid)%packbuf_rank (pid) = 0
    ginfo(gid)%packbuf_addr (pid) = 0

    return
  end subroutine COMM_packbuf_release

  !-----------------------------------------------------------------------------
  !> COMM_vsize_max: the IDs of non-persistent exchanges (and of the packing buffers) are 1 to it
  function COMM_get_vsize_max() result(vsize_max)
    implicit none
    integer :: vsize_max
    !---------------------------------------------------------------------------

    vsize_max = COMM_vsize_max

    return
  end function COMM_get_vsize_max

  !-----------------------------------------------------------------------------
  !> Current capacity of the persistent variables of the grid (0 without persistent communication)
  function COMM_get_pc_capacity( gid ) result(capacity)
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: gid
    integer :: capacity
    !---------------------------------------------------------------------------

    if ( gid < 1 .or. gid > COMM_gid ) then
       LOG_ERROR("COMM_get_pc_capacity",*) 'gid is invalid', gid, COMM_gid
       call PRC_abort
    end if

    capacity = 0
    if ( allocated( ginfo(gid)%packid ) ) capacity = size( ginfo(gid)%packid )

    return
  end function COMM_get_pc_capacity

  !-----------------------------------------------------------------------------
  !> Information of the persistent exchange of a variable, for tests and measurements: nreq is the
  !! number of the persistent requests of the variable, started by an exchange of the variable
  !! unless it is in a group in batch (0 for the vid of a non-persistent exchange); state is the
  !! state of its group (COMM_PC_GROUP_NONE: not in a group, or in a group that is not checked;
  !! COMM_PC_GROUP_CHECK: in a group, exchanged alone; COMM_PC_GROUP_BATCH: in a group exchanged in
  !! batch); nreq_group is the number of the requests of its group in batch, started once in a
  !! round of the group (0 otherwise), and nmem_group the number of the members of its group (0
  !! if it is not in a group). It does not change the communication.
  subroutine COMM_get_pc_info( &
       vid,        &
       nreq,       &
       state,      &
       nreq_group, &
       nmem_group, &
       gid         )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in)  :: vid   !< ID of the exchange (COMM_vars_init, COMM_vars8_init, or a vid)
    integer, intent(out) :: nreq  !< number of the persistent requests
    integer, intent(out) :: state !< state of the group of the variable
    integer, intent(out), optional :: nreq_group !< number of the requests of the group in batch
    integer, intent(out), optional :: nmem_group !< number of the members of the group
    integer, intent(in),  optional :: gid

    integer :: gid_
    integer :: ig
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR("COMM_get_pc_info",*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    nreq  = 0
    state = COMM_PC_GROUP_NONE

    if ( vid < 1 ) then
       LOG_ERROR("COMM_get_pc_info",*) 'vid is invalid', vid, COMM_vsize_max, gid_
       call PRC_abort
    else if ( vid > COMM_vsize_max ) then
       if ( .not. COMM_USE_MPI_PC ) then
          LOG_ERROR("COMM_get_pc_info",*) 'vid is invalid: persistent communication is not used', vid, COMM_vsize_max, gid_
          call PRC_abort
       else if ( vid - COMM_vsize_max > ginfo(gid_)%vars_num ) then
          LOG_ERROR("COMM_get_pc_info",*) 'vid is invalid: not initialized', vid, COMM_vsize_max, ginfo(gid_)%vars_num, gid_
          call PRC_abort
       end if
       nreq = ginfo(gid_)%preq_cnt(vid-COMM_vsize_max)
       ig = ginfo(gid_)%pc_group(vid-COMM_vsize_max)
       if ( ig > 0 ) then
          state = ginfo(gid_)%grp(ig)%state
          if ( present(nreq_group) ) nreq_group = ginfo(gid_)%grp(ig)%nreq
          if ( present(nmem_group) ) nmem_group = ginfo(gid_)%grp(ig)%nmem
          return
       end if
    end if
    if ( present(nreq_group) ) nreq_group = 0
    if ( present(nmem_group) ) nmem_group = 0

    return
  end subroutine COMM_get_pc_info

  !-----------------------------------------------------------------------------
  !> Declare a group of persistent variables of the 8-direction exchange (COMM_vars8_init) on a
  !! grid, which are exchanged together. This is collective over the local communicator: all the
  !! processes must call it in the same order, with the same IDs (it aborts otherwise).
  !! The rule of the consecutive posts: between the post (COMM_vars8) of the first member of the
  !! group in a round and that of the last member,
  !!   - no member of the group may be waited (it aborts), and
  !!   - no operation that waits for other processes (a COMM_wait of another exchange, a collective
  !!     communication, a barrier, the completion of a message) may be placed where the processes
  !!     differ (it may hang; this is not detected): the messages of a batch are started only at
  !!     the post of the last member. The same calls at the same place in all the processes are
  !!     allowed, e.g., an exchange of a variable outside the group that is posted and waited
  !!     between the posts of the members in all the processes.
  !! A member posted must not be changed until it is waited: in a batch it is packed by a kernel in
  !! the queue, which is waited only at the post of the last member.
  !! All the members are posted in a round before they are waited; the order of the posts and of
  !! the waits may differ among the processes. A member waited may not be posted again before all
  !! the members are waited (it aborts).
  !! The state of the group (see COMM_get_pc_info):
  !!   not checked (COMM_PC_GROUP_NONE): persistent communication is not used, or a member is not a
  !!     persistent variable (e.g., the array is not on the device in OpenACC builds), or no member
  !!   checked (COMM_PC_GROUP_CHECK): the rule is checked, and the members are exchanged alone
  !!     (COMM_USE_BATCH = .false., or one member, or in 2-D (PRC_TwoD))
  !!   batch (COMM_PC_GROUP_BATCH): the rule is checked, the members are packed into the buffer of
  !!     the group, and one message is sent in each direction for all the members when the last
  !!     member is posted (COMM_USE_BATCH = .true., not in 2-D)
  !! The group lives as long as its members: COMM_vars_finalize and COMM_finalize release it. A
  !! variable may be in one group only, and must not be in an exchange in progress (posted and not
  !! waited) at the declaration (it aborts). The IDs of the variables do not change.
  !! In 2-D (PRC_TwoD), a member that is not exchanged in 2-D must not be in the group.
  subroutine COMM_vars_group_init( &
       name, &
       vids, &
       gid   )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    character(len=*), intent(in) :: name    !< name of the group
    integer,          intent(in) :: vids(:) !< IDs of the members (COMM_vars8_init)

    integer,          intent(in), optional :: gid

    integer :: gid_
    integer :: n, m, m2, seq, ig, state
    integer :: hdr(2), hdr_min(2), hdr_max(2)
    integer, allocatable :: key(:), key_min(:), key_max(:)
    integer(int64) :: bytes
    integer :: ierr
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR("COMM_vars_group_init",*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    if ( .not. COMM_USE_MPI_PC ) then ! the same in all the processes
       LOG_INFO("COMM_vars_group_init",*) 'Group ', trim(name), ' (grid ID = ', gid_, &
                                          '): not checked (persistent communication is not used)'
       return
    end if

    n = size(vids)

    ! the IDs
    do m = 1, n
       if ( vids(m) < 1 .or. vids(m) > COMM_vsize_max + ginfo(gid_)%vars_num ) then
          LOG_ERROR("COMM_vars_group_init",*) 'vid of a member is invalid: ', trim(name), vids(m), gid_
          call PRC_abort
       else if ( vids(m) > COMM_vsize_max ) then
          if ( ginfo(gid_)%pc_ndir(vids(m)-COMM_vsize_max) /= 8 ) then
             LOG_ERROR("COMM_vars_group_init",*) 'a member is registered for the 4-direction exchange: ', trim(name), vids(m), gid_
             call PRC_abort
          end if
       end if
    end do

    ! the definition must be the same in all the processes
    hdr(1) = n
    hdr(2) = gid_
    call MPI_Allreduce( hdr, hdr_min, 2, MPI_INTEGER, MPI_MIN, COMM_world_t, ierr )
    call MPI_Allreduce( hdr, hdr_max, 2, MPI_INTEGER, MPI_MAX, COMM_world_t, ierr )
    if ( any( hdr_min(:) /= hdr_max(:) ) ) then
       LOG_ERROR("COMM_vars_group_init",*) 'the group differs among the processes: ', trim(name)
       LOG_ERROR_CONT(*) 'number of the members (min, max): ', hdr_min(1), hdr_max(1)
       LOG_ERROR_CONT(*) 'grid ID (min, max): ', hdr_min(2), hdr_max(2)
       call PRC_abort
    end if
    if ( n > 0 ) then
       ! the IDs (the sequential numbers, and whether persistent), and the directions
       allocate( key(2*n), key_min(2*n), key_max(2*n) )
       do m = 1, n
          key(m) = vids(m)
          key(n+m) = 0
          if ( vids(m) > COMM_vsize_max ) key(n+m) = ginfo(gid_)%pc_ndir(vids(m)-COMM_vsize_max)
       end do
       call MPI_Allreduce( key, key_min, 2*n, MPI_INTEGER, MPI_MIN, COMM_world_t, ierr )
       call MPI_Allreduce( key, key_max, 2*n, MPI_INTEGER, MPI_MAX, COMM_world_t, ierr )
       if ( any( key_min(:) /= key_max(:) ) ) then
          LOG_ERROR("COMM_vars_group_init",*) 'the IDs of the members differ among the processes: ', trim(name)
          do m = 1, n
             LOG_ERROR_CONT(*) 'member, ID (min, max), direction (min, max): ', m, key_min(m), key_max(m), &
                                                                                key_min(n+m), key_max(n+m)
          end do
          call PRC_abort
       end if
       deallocate( key, key_min, key_max )
    end if

    ! not checked: the same in all the processes from here
    if ( n == 0 ) then
       LOG_INFO("COMM_vars_group_init",*) 'Group ', trim(name), ' (grid ID = ', gid_, '): not checked (no member)'
       return
    end if
    if ( any( vids(:) <= COMM_vsize_max ) ) then
       LOG_INFO("COMM_vars_group_init",*) 'Group ', trim(name), ' (grid ID = ', gid_, &
                                          '): not checked (a member is not a persistent variable)'
       return
    end if

    ! local checks
    do m = 1, n
       seq = vids(m) - COMM_vsize_max
       if ( ginfo(gid_)%pc_group(seq) /= 0 ) then
          LOG_ERROR("COMM_vars_group_init",*) 'a member is already in another group: ', trim(name), vids(m), gid_
          LOG_ERROR_CONT(*) 'the other group: ', trim(ginfo(gid_)%grp(ginfo(gid_)%pc_group(seq))%name)
          call PRC_abort
       end if
       ! posted and not waited: its wait would be taken as a wait in a round of the group
       if ( ginfo(gid_)%packbuf_owner( ginfo(gid_)%packid(seq) ) == vids(m) ) then
          LOG_ERROR("COMM_vars_group_init",*) 'a member is in an exchange in progress (posted and not waited): ', &
                                              trim(name), vids(m), gid_
          call PRC_abort
       end if
       do m2 = 1, m - 1
          if ( vids(m2) == vids(m) ) then
             LOG_ERROR("COMM_vars_group_init",*) 'a member is given twice: ', trim(name), vids(m), gid_
             call PRC_abort
          else if ( ginfo(gid_)%packid(vids(m2)-COMM_vsize_max) == ginfo(gid_)%packid(seq) ) then
             LOG_ERROR("COMM_vars_group_init",*) 'two members share the packing buffer of vid: ', trim(name), &
                                                 ginfo(gid_)%packid(seq), gid_
             LOG_ERROR_CONT(*) 'IDs of the members: ', vids(m2), vids(m)
             call PRC_abort
          end if
       end do
    end do

    ig = ginfo(gid_)%grp_num + 1
    if ( ( ig + 1 ) * G_NDIR - 1 > COMM_tag_ub ) then
       LOG_ERROR("COMM_vars_group_init",*) 'message tags of the groups exceed MPI_TAG_UB: number of the groups, MPI_TAG_UB = ', &
                                           ig, COMM_tag_ub
       call PRC_abort
    end if
    if ( .not. allocated( ginfo(gid_)%grp ) ) then
       call COMM_grow_grp( gid_ )
    else if ( ig > size( ginfo(gid_)%grp ) ) then
       call COMM_grow_grp( gid_ )
    end if
    ginfo(gid_)%grp_num = ig

    ! not in batch in 2-D: only the N/S halos are exchanged, from and to the variables without
    ! packing, and the packing of the batch costs more than the messages saved (issue #3615)
    if ( COMM_USE_BATCH .and. n >= 2 .and. .not. PRC_TwoD ) then
       state = COMM_PC_GROUP_BATCH
    else
       state = COMM_PC_GROUP_CHECK
    end if

    ginfo(gid_)%grp(ig)%name    = name
    ginfo(gid_)%grp(ig)%state   = state
    ginfo(gid_)%grp(ig)%nmem    = n
    ginfo(gid_)%grp(ig)%nposted = 0
    ginfo(gid_)%grp(ig)%nwaited = 0
    ginfo(gid_)%grp(ig)%arrived = .false.
    ginfo(gid_)%grp(ig)%nreq    = 0
    ginfo(gid_)%grp(ig)%nsend   = 0
    ginfo(gid_)%grp(ig)%nbuf    = 0
    ginfo(gid_)%grp(ig)%req(:)  = MPI_REQUEST_NULL
    allocate( ginfo(gid_)%grp(ig)%seq(n) )
    do m = 1, n
       seq = vids(m) - COMM_vsize_max
       ginfo(gid_)%grp(ig)%seq(m) = seq
       ginfo(gid_)%pc_group(seq) = ig
       ginfo(gid_)%pc_slot (seq) = m
       ginfo(gid_)%pc_gstat(seq) = 0
    end do

    if ( state == COMM_PC_GROUP_BATCH ) then
       call COMM_group_init_batch( gid_, ig )
       bytes = int( ginfo(gid_)%grp(ig)%nbuf, kind=int64 ) * 2 * COMM_datasize
       LOG_INFO("COMM_vars_group_init",*) 'Group ', trim(name), ' (grid ID = ', gid_, '): batch, members = ', n, &
                                          ', requests = ', ginfo(gid_)%grp(ig)%nreq, ', buffers [byte] = ', bytes
    else
       LOG_INFO("COMM_vars_group_init",*) 'Group ', trim(name), ' (grid ID = ', gid_, '): checked, members = ', n
    end if
    LOG_INFO_CONT(*) 'IDs of the members: ', vids(:)

    return
  end subroutine COMM_vars_group_init

  !-----------------------------------------------------------------------------
  !> Allocate the table of the groups of the grid, or extend it to twice the size. The request
  !! handles are copied as values, and the buffers are not moved (pointer components)
  subroutine COMM_grow_grp( gid )
    implicit none

    integer, intent(in) :: gid

    type(group_t), allocatable :: grp(:)
    integer :: nold, n
    !---------------------------------------------------------------------------

    if ( allocated( ginfo(gid)%grp ) ) then
       nold = size( ginfo(gid)%grp )
       allocate( grp(nold*2) )
       grp(1:nold) = ginfo(gid)%grp(1:nold)
       call move_alloc( grp, ginfo(gid)%grp )
       LOG_INFO("COMM_grow_grp",*) 'The capacity of the groups is extended (grid ID, size): ', gid, nold*2
    else
       nold = 0
       allocate( ginfo(gid)%grp(4) )
    end if
    do n = nold + 1, size( ginfo(gid)%grp )
       ginfo(gid)%grp(n)%req(:) = MPI_REQUEST_NULL
    end do

    return
  end subroutine COMM_grow_grp

  !-----------------------------------------------------------------------------
  !> Elements of a member in the message of each direction of a group (by the direction in which
  !! the message travels: see G_N etc.), offsets of the sections of the directions in the buffers,
  !! and the size of each buffer. A section has the slots of the members one after another, so
  !! that the message of a direction is contiguous. N/S: whole rows (i = IS:IE in 2-D);
  !! W/E: the rows of L2_range_WE; corners: IHALO x JHALO blocks. No W/E and corners with
  !! PRC_TwoD or IHALO = 0 (see pack8_sections). The sections have the largest sizes, independent
  !! of the neighbors
  subroutine G_layout( gid, nmem, elem, off, nbuf )
    implicit none

    integer, intent(in)  :: gid
    integer, intent(in)  :: nmem
    integer, intent(out) :: elem(G_NDIR)
    integer, intent(out) :: off (G_NDIR)
    integer, intent(out) :: nbuf

    integer :: KA, IA, JA, IHALO, JHALO
    integer :: i0, i1, jw0, jw1
    logical :: with_x, has_N, has_S, has_W, has_E
    integer :: emax(G_NDIR)
    integer :: d
    !---------------------------------------------------------------------------

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    JA    = ginfo(gid)%JA
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    call pack8_sections( gid, with_x, has_N, has_S, has_W, has_E, i0, i1 )

    elem(G_N) = KA * ( i1 - i0 + 1 ) * JHALO
    elem(G_S) = elem(G_N)
    emax(G_N) = KA * IA * JHALO
    emax(G_S) = emax(G_N)

    if ( with_x ) then
       call L2_range_WE( gid, jw0, jw1 )
       elem(G_W) = KA * IHALO * ( jw1 - jw0 + 1 )
       emax(G_W) = KA * IHALO * JA
       elem(G_NW) = KA * IHALO * JHALO
       emax(G_NW) = elem(G_NW)
    else
       elem(G_W)  = 0
       emax(G_W)  = 0
       elem(G_NW) = 0
       emax(G_NW) = 0
    end if
    elem(G_E) = elem(G_W)
    emax(G_E) = emax(G_W)
    elem(G_NE:G_SE) = elem(G_NW)
    emax(G_NE:G_SE) = emax(G_NW)

    off(1) = 0
    do d = 2, G_NDIR
       off(d) = off(d-1) + nmem * emax(d-1)
    end do
    nbuf = off(G_NDIR) + nmem * emax(G_NDIR)

    return
  end subroutine G_layout

  !-----------------------------------------------------------------------------
  !> Sections of the packed 8-direction exchange of a grid (the layout 2 and the groups in batch),
  !! for both the requests and the packing: with_x is false with PRC_TwoD or IHALO = 0 (no W/E
  !! halos and no corners), has_* tell the neighbors with which the messages are exchanged (the
  !! corners with the diagonal neighbors, where the neighbors in x and in y exist), and i0:i1 are
  !! the columns of the whole rows in y (the inner columns only in 2-D)
  subroutine pack8_sections( gid, with_x, has_N, has_S, has_W, has_E, i0, i1 )
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    integer, intent(in)  :: gid
    logical, intent(out) :: with_x
    logical, intent(out) :: has_N, has_S, has_W, has_E
    integer, intent(out) :: i0, i1
    !---------------------------------------------------------------------------

    with_x = ( .not. PRC_TwoD ) .and. ginfo(gid)%IHALO > 0
    has_N  = PRC_HAS_N
    has_S  = PRC_HAS_S
    has_W  = with_x .and. PRC_HAS_W
    has_E  = with_x .and. PRC_HAS_E

    if ( PRC_TwoD ) then
       i0 = ginfo(gid)%IS
       i1 = ginfo(gid)%IE
    else
       i0 = 1
       i1 = ginfo(gid)%IA
    end if

    return
  end subroutine pack8_sections

  !-----------------------------------------------------------------------------
  !> Offsets of the sections of a member of a group in batch (slot m) in the buffer of the group,
  !! by the neighbor of each direction (see pack8_3D and unpack8_3D): to the neighbor for the
  !! sending buffer, and from the neighbor for the receiving one (the section of the message in
  !! the opposite direction)
  subroutine G_offsets( gid, nmem, m, recv, off )
    implicit none
    integer, intent(in)  :: gid
    integer, intent(in)  :: nmem      !< number of the members
    integer, intent(in)  :: m         !< slot of the member
    logical, intent(in)  :: recv      !< offsets in the receiving buffer
    integer, intent(out) :: off(G_NDIR)

    integer :: elem(G_NDIR), offd(G_NDIR), nbuf
    integer :: d
    !---------------------------------------------------------------------------

    call G_layout( gid, nmem, elem, offd, nbuf )
    do d = 1, G_NDIR
       offd(d) = offd(d) + ( m - 1 ) * elem(d)
    end do

    if ( recv ) then
       off(G_N ) = offd(G_S )
       off(G_S ) = offd(G_N )
       off(G_W ) = offd(G_E )
       off(G_E ) = offd(G_W )
       off(G_NW) = offd(G_SE)
       off(G_NE) = offd(G_SW)
       off(G_SW) = offd(G_NE)
       off(G_SE) = offd(G_NW)
    else
       off(:) = offd(:)
    end if

    return
  end subroutine G_offsets

  !-----------------------------------------------------------------------------
  !> Buffers and persistent requests of a group in batch: one message in each direction with a
  !! neighbor, for all the members (16 requests at most). The tags are those of the directions of
  !! the group, on the communicator of the groups: a neighbor in more than one direction (e.g., the
  !! process itself) is told apart. The buffers must not be associated (they are freed in
  !! COMM_group_free)
  subroutine COMM_group_init_batch( gid, ig )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in) :: gid
    integer, intent(in) :: ig

    integer :: elem(G_NDIR), off(G_NDIR)
    integer :: n, nbuf, tag, ireq, nsend
    integer :: i0, i1
    logical :: with_x, has_N, has_S, has_W, has_E
    logical :: flag
    integer :: i, ierr

    real(RP), pointer, contiguous :: ptr(:)
    !---------------------------------------------------------------------------

    if ( associated( ginfo(gid)%grp(ig)%sendbuf ) .or. associated( ginfo(gid)%grp(ig)%recvbuf ) ) then
       LOG_ERROR("COMM_group_init_batch",*) 'the buffers of the group are not freed: grid, group', gid, ig
       call PRC_abort
    end if

    n = ginfo(gid)%grp(ig)%nmem
    call G_layout( gid, n, elem, off, nbuf )
    ginfo(gid)%grp(ig)%nbuf = nbuf

    allocate( ginfo(gid)%grp(ig)%sendbuf(nbuf) )
    allocate( ginfo(gid)%grp(ig)%recvbuf(nbuf) )
    ginfo(gid)%grp(ig)%sendbuf(:) = 0.0_RP
    ginfo(gid)%grp(ig)%recvbuf(:) = 0.0_RP
    ptr => ginfo(gid)%grp(ig)%sendbuf
    !$acc enter data copyin(ptr)
    ptr => ginfo(gid)%grp(ig)%recvbuf
    !$acc enter data copyin(ptr)

    ! the same sections as in the packing (pack8_3D, unpack8_3D)
    call pack8_sections( gid, with_x, has_N, has_S, has_W, has_E, i0, i1 )

    tag   = ig * G_NDIR - 1 ! + direction
    ireq  = 0
    nsend = 0

    !--- receive, from the neighbor opposite to the direction of the message
    ptr => ginfo(gid)%grp(ig)%recvbuf
    !$acc host_data use_device(ptr)
    if ( has_S ) then
       ireq = ireq + 1
       call MPI_RECV_INIT( ptr(off(G_N)+1), n*elem(G_N), COMM_datatype_t, PRC_next(PRC_S), tag+G_N, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
    end if
    if ( has_N ) then
       ireq = ireq + 1
       call MPI_RECV_INIT( ptr(off(G_S)+1), n*elem(G_S), COMM_datatype_t, PRC_next(PRC_N), tag+G_S, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
    end if
    if ( has_E ) then
       ireq = ireq + 1
       call MPI_RECV_INIT( ptr(off(G_W)+1), n*elem(G_W), COMM_datatype_t, PRC_next(PRC_E), tag+G_W, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
    end if
    if ( has_W ) then
       ireq = ireq + 1
       call MPI_RECV_INIT( ptr(off(G_E)+1), n*elem(G_E), COMM_datatype_t, PRC_next(PRC_W), tag+G_E, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
    end if
    if ( has_S .and. has_E ) then
       ireq = ireq + 1
       call MPI_RECV_INIT( ptr(off(G_NW)+1), n*elem(G_NW), COMM_datatype_t, PRC_next(PRC_SE), tag+G_NW, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
    end if
    if ( has_S .and. has_W ) then
       ireq = ireq + 1
       call MPI_RECV_INIT( ptr(off(G_NE)+1), n*elem(G_NE), COMM_datatype_t, PRC_next(PRC_SW), tag+G_NE, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
    end if
    if ( has_N .and. has_E ) then
       ireq = ireq + 1
       call MPI_RECV_INIT( ptr(off(G_SW)+1), n*elem(G_SW), COMM_datatype_t, PRC_next(PRC_NE), tag+G_SW, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
    end if
    if ( has_N .and. has_W ) then
       ireq = ireq + 1
       call MPI_RECV_INIT( ptr(off(G_SE)+1), n*elem(G_SE), COMM_datatype_t, PRC_next(PRC_NW), tag+G_SE, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
    end if
    !$acc end host_data

    !--- send
    ptr => ginfo(gid)%grp(ig)%sendbuf
    !$acc host_data use_device(ptr)
    if ( has_N ) then
       ireq = ireq + 1
       call MPI_SEND_INIT( ptr(off(G_N)+1), n*elem(G_N), COMM_datatype_t, PRC_next(PRC_N), tag+G_N, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
       nsend = nsend + n*elem(G_N)
    end if
    if ( has_S ) then
       ireq = ireq + 1
       call MPI_SEND_INIT( ptr(off(G_S)+1), n*elem(G_S), COMM_datatype_t, PRC_next(PRC_S), tag+G_S, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
       nsend = nsend + n*elem(G_S)
    end if
    if ( has_W ) then
       ireq = ireq + 1
       call MPI_SEND_INIT( ptr(off(G_W)+1), n*elem(G_W), COMM_datatype_t, PRC_next(PRC_W), tag+G_W, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
       nsend = nsend + n*elem(G_W)
    end if
    if ( has_E ) then
       ireq = ireq + 1
       call MPI_SEND_INIT( ptr(off(G_E)+1), n*elem(G_E), COMM_datatype_t, PRC_next(PRC_E), tag+G_E, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
       nsend = nsend + n*elem(G_E)
    end if
    if ( has_N .and. has_W ) then
       ireq = ireq + 1
       call MPI_SEND_INIT( ptr(off(G_NW)+1), n*elem(G_NW), COMM_datatype_t, PRC_next(PRC_NW), tag+G_NW, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
       nsend = nsend + n*elem(G_NW)
    end if
    if ( has_N .and. has_E ) then
       ireq = ireq + 1
       call MPI_SEND_INIT( ptr(off(G_NE)+1), n*elem(G_NE), COMM_datatype_t, PRC_next(PRC_NE), tag+G_NE, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
       nsend = nsend + n*elem(G_NE)
    end if
    if ( has_S .and. has_W ) then
       ireq = ireq + 1
       call MPI_SEND_INIT( ptr(off(G_SW)+1), n*elem(G_SW), COMM_datatype_t, PRC_next(PRC_SW), tag+G_SW, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
       nsend = nsend + n*elem(G_SW)
    end if
    if ( has_S .and. has_E ) then
       ireq = ireq + 1
       call MPI_SEND_INIT( ptr(off(G_SE)+1), n*elem(G_SE), COMM_datatype_t, PRC_next(PRC_SE), tag+G_SE, &
                           ginfo(gid)%comm_grp, ginfo(gid)%grp(ig)%req(ireq), ierr )
       nsend = nsend + n*elem(G_SE)
    end if
    !$acc end host_data

    ginfo(gid)%grp(ig)%nreq  = ireq
    ginfo(gid)%grp(ig)%nsend = nsend

    ! to finish initial processes of MPI
    do i = 1, 32
       call MPI_TESTALL( ireq, ginfo(gid)%grp(ig)%req(1:ireq), flag, MPI_STATUSES_IGNORE, ierr )
    enddo

    return
  end subroutine COMM_group_init_batch

  !-----------------------------------------------------------------------------
  !> Free the requests and the buffers of a group, and reset its entry (COMM_pc_free)
  subroutine COMM_group_free( gid, ig )
    implicit none

    integer, intent(in) :: gid
    integer, intent(in) :: ig

    real(RP), pointer, contiguous :: ptr(:)
    integer :: i, ierr
    !---------------------------------------------------------------------------

    do i = 1, ginfo(gid)%grp(ig)%nreq
       if ( ginfo(gid)%grp(ig)%req(i) /= MPI_REQUEST_NULL ) then
          call MPI_REQUEST_FREE( ginfo(gid)%grp(ig)%req(i), ierr )
       end if
    end do
    ginfo(gid)%grp(ig)%req(:) = MPI_REQUEST_NULL

    if ( associated( ginfo(gid)%grp(ig)%sendbuf ) ) then
       ptr => ginfo(gid)%grp(ig)%sendbuf
       !$acc exit data delete(ptr)
       deallocate( ginfo(gid)%grp(ig)%sendbuf )
    end if
    if ( associated( ginfo(gid)%grp(ig)%recvbuf ) ) then
       ptr => ginfo(gid)%grp(ig)%recvbuf
       !$acc exit data delete(ptr)
       deallocate( ginfo(gid)%grp(ig)%recvbuf )
    end if
    nullify( ginfo(gid)%grp(ig)%sendbuf, ginfo(gid)%grp(ig)%recvbuf )

    if ( allocated( ginfo(gid)%grp(ig)%seq ) ) deallocate( ginfo(gid)%grp(ig)%seq )
    ginfo(gid)%grp(ig)%name    = ''
    ginfo(gid)%grp(ig)%state   = COMM_PC_GROUP_NONE
    ginfo(gid)%grp(ig)%nmem    = 0
    ginfo(gid)%grp(ig)%nposted = 0
    ginfo(gid)%grp(ig)%nwaited = 0
    ginfo(gid)%grp(ig)%arrived = .false.
    ginfo(gid)%grp(ig)%nreq    = 0
    ginfo(gid)%grp(ig)%nsend   = 0
    ginfo(gid)%grp(ig)%nbuf    = 0

    return
  end subroutine COMM_group_free

  !-----------------------------------------------------------------------------
  !> A member of a group is posted (vars_3D_mpi_pc): the round of the group is checked and counted
  subroutine COMM_group_post( gid, seq, ig )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in) :: gid
    integer, intent(in) :: seq !< sequential number of the member
    integer, intent(in) :: ig  !< group of the member
    !---------------------------------------------------------------------------

    if ( ginfo(gid)%pc_gstat(seq) /= 0 ) then
       LOG_ERROR("COMM_vars8_3D",*) 'a member of the group is posted again before all the members are waited: ', &
                                    trim(ginfo(gid)%grp(ig)%name)
       LOG_ERROR_CONT(*) 'grid, ID posted, members posted, members waited: ', gid, seq + COMM_vsize_max, &
                         ginfo(gid)%grp(ig)%nposted, ginfo(gid)%grp(ig)%nwaited
       call PRC_abort
    end if

    ginfo(gid)%pc_gstat(seq) = 1
    ginfo(gid)%grp(ig)%nposted = ginfo(gid)%grp(ig)%nposted + 1
    if ( ginfo(gid)%grp(ig)%nmem > 1 ) then
       if ( ginfo(gid)%grp(ig)%nposted == 1 ) COMM_grp_posting = COMM_grp_posting + 1
       if ( ginfo(gid)%grp(ig)%nposted == ginfo(gid)%grp(ig)%nmem ) COMM_grp_posting = COMM_grp_posting - 1
    end if

    return
  end subroutine COMM_group_post

  !-----------------------------------------------------------------------------
  !> A member of a group is waited (wait_3D_mpi_pc): the round ends when all the members are waited
  subroutine COMM_group_waited( gid, seq, ig )
    implicit none

    integer, intent(in) :: gid
    integer, intent(in) :: seq !< sequential number of the member
    integer, intent(in) :: ig  !< group of the member

    integer :: m
    !---------------------------------------------------------------------------

    ginfo(gid)%pc_gstat(seq) = 2
    ginfo(gid)%grp(ig)%nwaited = ginfo(gid)%grp(ig)%nwaited + 1
    if ( ginfo(gid)%grp(ig)%nwaited == ginfo(gid)%grp(ig)%nmem ) then
       do m = 1, ginfo(gid)%grp(ig)%nmem
          ginfo(gid)%pc_gstat(ginfo(gid)%grp(ig)%seq(m)) = 0
       end do
       ginfo(gid)%grp(ig)%nposted = 0
       ginfo(gid)%grp(ig)%nwaited = 0
       ginfo(gid)%grp(ig)%arrived = .false.
    end if

    return
  end subroutine COMM_group_waited

  !-----------------------------------------------------------------------------
  !> COMM_wait while the posts of a group are not complete (COMM_grp_posting > 0): abort if the
  !! exchange waited is that of a member of such a group, with the group and the members not posted
  subroutine COMM_group_wait_check( routine_name, vid, gid )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: routine_name
    integer,          intent(in) :: vid
    integer,          intent(in) :: gid

    integer :: seq, ig, g, m
    !---------------------------------------------------------------------------

    ig = 0
    if ( vid > COMM_vsize_max .and. COMM_USE_MPI_PC ) then ! vid is checked (COMM_check_vid)
       seq = vid - COMM_vsize_max
       ig = ginfo(gid)%pc_group(seq)
       if ( ig > 0 ) then
          if ( ginfo(gid)%grp(ig)%nposted == 0 .or. ginfo(gid)%grp(ig)%nposted == ginfo(gid)%grp(ig)%nmem ) ig = 0
       end if
    end if

    ! an exchange outside the groups with posts not complete may be waited (at the same place in
    ! all the processes: see COMM_vars_group_init)
    if ( ig == 0 ) return

    LOG_ERROR(routine_name,*) 'a member of the group is waited before all its members are posted: ', &
                              trim(ginfo(gid)%grp(ig)%name)
    LOG_ERROR_CONT(*) 'grid, ID waited: ', gid, vid
    do m = 1, ginfo(gid)%grp(ig)%nmem
       if ( ginfo(gid)%pc_gstat(ginfo(gid)%grp(ig)%seq(m)) == 0 ) then
          LOG_ERROR_CONT(*) 'ID of a member not posted: ', ginfo(gid)%grp(ig)%seq(m) + COMM_vsize_max
       end if
    end do
    do g = 1, COMM_gid
       do m = 1, ginfo(g)%grp_num
          if ( ginfo(g)%grp(m)%nposted > 0 .and. ginfo(g)%grp(m)%nposted < ginfo(g)%grp(m)%nmem ) then
             LOG_ERROR_CONT(*) 'group with posts not complete (name, grid, members posted, members): ', &
                               trim(ginfo(g)%grp(m)%name), g, ginfo(g)%grp(m)%nposted, ginfo(g)%grp(m)%nmem
          end if
       end do
    end do
    call PRC_abort

    return
  end subroutine COMM_group_wait_check

  !-----------------------------------------------------------------------------
  !> Register a variable for persistent communication (COMM_vars_init, COMM_vars8_init)
  subroutine COMM_pc_regist( routine_name, varname, var, vid, gid, dir8 )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in)    :: routine_name
    character(len=*), intent(in)    :: varname    !< variable name
    real(RP),         intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID
    integer,          intent(in), optional :: gid
    logical,          intent(in)    :: dir8       !< 8-direction exchange (4-direction if .false.)

    integer :: gid_
    integer :: vars_id
    !---------------------------------------------------------------------------

    COMM_CHECK_CONTIGUOUS(routine_name,var)

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ < 1 .or. gid_ > COMM_gid ) then
       LOG_ERROR(routine_name,*) 'gid is invalid', gid_, COMM_gid
       call PRC_abort
    end if

    if ( vid < 1 .or. vid > COMM_vsize_max ) then
       LOG_ERROR(routine_name,*) 'vid is invalid', vid, COMM_vsize_max, gid_
       call PRC_abort
    end if

    if ( .not. COMM_USE_MPI_PC ) return
#ifdef _OPENACC
    if ( .not. acc_is_present(var) ) return
#endif

    call PROF_rapstart('COMM_init_pers', 2)

    ginfo(gid_)%vars_num = ginfo(gid_)%vars_num + 1
    if ( ginfo(gid_)%vars_num > size(ginfo(gid_)%packid) ) call COMM_grow_pc( gid_ )

    vars_id = ginfo(gid_)%vars_num
    ginfo(gid_)%packid(vars_id) = vid
    ginfo(gid_)%pc_addr(vars_id) = COMM_address( var )
    if ( dir8 ) then
       ginfo(gid_)%pc_ndir(vars_id) = 8
    else
       ginfo(gid_)%pc_ndir(vars_id) = 4
    end if

    if ( COMM_tag_base( vars_id ) + COMM_tag_width - 1 > COMM_tag_ub ) then
       LOG_ERROR(routine_name,*) 'message tags exceed MPI_TAG_UB: sequential number, MPI_TAG_UB = ', vars_id, COMM_tag_ub
       call PRC_abort
    end if

    if ( dir8 .and. ginfo(gid_)%pc_layout2 ) then
       call vars8_init_mpi_pc_L2(var, gid_, vars_id, vid)
    else if ( dir8 ) then
       call vars8_init_mpi_pc(var, gid_, vars_id, vid)
    else
       call vars_init_mpi_pc(var, gid_, vars_id, vid)
    end if

    vid = vars_id + COMM_vsize_max

    LOG_INFO(routine_name,'(1x,A,I3.3,A,I0,2A)') 'Initialize variable (grid ID = ', gid_, '): ID = ', vid, &
                                                                               ', name = ', trim(varname)

    call PROF_rapend  ('COMM_init_pers', 2)

    return
  end subroutine COMM_pc_regist

  !-----------------------------------------------------------------------------
  !> Free the persistent requests of all the persistent variables of the grid and reset their
  !! slots (COMM_vars_finalize, COMM_finalize). The tables keep their capacity.
  subroutine COMM_pc_free( routine_name, gid )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: routine_name
    integer,          intent(in) :: gid

    integer :: nvars
    integer :: i, n, ierr
    !---------------------------------------------------------------------------

    if ( .not. COMM_USE_MPI_PC ) return
    nvars = ginfo(gid)%vars_num
    if ( nvars == 0 ) return

    ! a round of a group in progress (some members are posted, and not all are waited). A member
    ! posted and not waited is also found by the check below; this one names the group
    do n = 1, ginfo(gid)%grp_num
       if ( ginfo(gid)%grp(n)%nposted > 0 ) then
          LOG_ERROR(routine_name,*) 'a round of the group is in progress at the release: ', trim(ginfo(gid)%grp(n)%name)
          LOG_ERROR_CONT(*) 'grid, members posted, members waited: ', gid, ginfo(gid)%grp(n)%nposted, ginfo(gid)%grp(n)%nwaited
          call PRC_abort
       end if
    end do

    ! a persistent exchange posted and not waited: its requests are still active
    if ( any( ginfo(gid)%packbuf_owner(:) > COMM_vsize_max ) ) then
       LOG_ERROR(routine_name,*) 'a persistent exchange is still active at the release: grid', gid
       LOG_ERROR_CONT(*) '(posted without COMM_wait before ', trim(routine_name), ')'
       do i = 1, COMM_vsize_max
          if ( ginfo(gid)%packbuf_owner(i) > COMM_vsize_max ) then
             LOG_ERROR_CONT(*) 'vid of the packing buffer, ID posted: ', i, ginfo(gid)%packbuf_owner(i)
          end if
       end do
       call PRC_abort
    end if

    ! the groups: their lifetime is that of the members
    do n = 1, ginfo(gid)%grp_num
       call COMM_group_free( gid, n )
    end do
    if ( ginfo(gid)%grp_num > 0 ) then
       LOG_INFO(routine_name,*) 'Release the groups (grid ID, number): ', gid, ginfo(gid)%grp_num
    end if
    ginfo(gid)%grp_num = 0

    do n = 1, nvars
       do i = 1, ginfo(gid)%nreq_MAX+1
          if ( ginfo(gid)%preq_list(i,n) /= MPI_REQUEST_NULL ) then
             call MPI_REQUEST_FREE( ginfo(gid)%preq_list(i,n), ierr )
          end if
       end do
    end do

    ginfo(gid)%packid   (  1:nvars) = 0
    ginfo(gid)%pc_addr  (  1:nvars) = 0
    ginfo(gid)%pc_ndir  (  1:nvars) = 0
    ginfo(gid)%pc_group (  1:nvars) = 0
    ginfo(gid)%pc_slot  (  1:nvars) = 0
    ginfo(gid)%pc_gstat (  1:nvars) = 0
    ginfo(gid)%preq_cnt (  1:nvars) = -1
    ginfo(gid)%psend_cnt(  1:nvars) = 0
    ginfo(gid)%preq_list(:,1:nvars) = MPI_REQUEST_NULL
    ginfo(gid)%vars_num = 0

    LOG_INFO(routine_name,*) 'Release the persistent variables (grid ID, number): ', gid, nvars

    return
  end subroutine COMM_pc_free

  !-----------------------------------------------------------------------------
  subroutine vars_init_mpi_pc(var, gid, seqid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in) :: gid
    integer,  intent(in) :: seqid
    integer,  intent(in) :: vid

    integer :: ireq, tag, ierr
    integer :: nsend
    logical :: flag

    integer :: KA
    integer :: JA, JS, JE, JHALO

    integer :: nreq
    integer :: i

#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
#endif

    tag  = COMM_tag_base( seqid )
    ireq = 1
    nsend = 0

    KA    = ginfo(gid)%KA
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$acc host_data use_device(var)

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype_t,                 &
                        MPI_PROC_NULL, tag, ginfo(gid)%comm_pc,                 &
                        ginfo(gid)%preq_list(ginfo(gid)%nreq_max+1,seqid), ierr )

    !--- From 4-Direction HALO communicate
    ! From S
    if ( PRC_HAS_S ) then
       call MPI_RECV_INIT( var(:,:,1:JS-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,                        &
                           PRC_next(PRC_S), tag+1, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
       ireq = ireq + 1
    end if
    ! From N
    if ( PRC_HAS_N ) then
       call MPI_RECV_INIT( var(:,:,JE+1:JA), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,                       &
                           PRC_next(PRC_N), tag+2, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
       ireq = ireq + 1
    end if
    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
       !$acc host_data use_device(ptr)
#endif

       ! From E
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_RECV_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#endif
                              PRC_next(PRC_E), tag+3, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
       end if
       ! From W
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_RECV_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#endif
                              PRC_next(PRC_W), tag+4, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
       end if
       !$acc end host_data
    end if

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc host_data use_device(ptr)
#endif
       ! To W HALO
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_SEND_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#endif
                              PRC_next(PRC_W), tag+3, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          ireq = ireq + 1
       end if
       ! To E HALO
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_SEND_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#endif
                              PRC_next(PRC_E), tag+4, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          ireq = ireq + 1
       end if
       !$acc end host_data
    end if
    ! To N HALO
    if ( PRC_HAS_N ) then
       call MPI_SEND_INIT( var(:,:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,                 &
                           PRC_next(PRC_N), tag+1, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4*KA
       ireq = ireq + 1
    end if
    ! To S HALO
    if ( PRC_HAS_S ) then
       call MPI_SEND_INIT( var(:,:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,                 &
                           PRC_next(PRC_S), tag+2, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4*KA
       ireq = ireq + 1
    end if

    ginfo(gid)%preq_cnt(seqid) = ireq - 1
    ginfo(gid)%psend_cnt(seqid) = nsend

    ! to finish initial processes of MPI
    nreq = ginfo(gid)%preq_cnt(seqid)
    do i = 1, 32
       call MPI_TESTALL( nreq, ginfo(gid)%preq_list(1:nreq,seqid), &
                         flag, MPI_STATUSES_IGNORE, ierr         )
    enddo

    !$acc end host_data

    return
  end subroutine vars_init_mpi_pc

  subroutine vars8_init_mpi_pc(var, gid, seqid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in) :: gid
    integer,  intent(in) :: seqid
    integer,  intent(in) :: vid

    integer :: ireq, tag, tagc
    integer :: ierr
    integer :: nsend
    logical :: flag

    integer :: KA
    integer :: IS, IE, IHALO
    integer :: JA, JS, JE, JHALO

    integer :: nreq
    integer :: i, j

#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
#endif

    KA    = ginfo(gid)%KA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    tag  = COMM_tag_base( seqid )
    ireq = 1
    nsend = 0

    !$acc host_data use_device(var)

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype_t,                 &
                        MPI_PROC_NULL, tag, ginfo(gid)%comm_pc,                 &
                        ginfo(gid)%preq_list(ginfo(gid)%nreq_max+1,seqid), ierr )


    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          tagc = 0
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                               &
                                 PRC_next(PRC_SE), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From SW
          tagc = 10
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                   &
                                  PRC_next(PRC_SW), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NE
          tagc = 20
          do j = JE+1, JA
             call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                               &
                                 PRC_next(PRC_NE), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NW
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                  &
                                 PRC_next(PRC_NW), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From E
          tagc = 60
#ifdef _OPENACC
          ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
          !$acc host_data use_device(ptr)
          call MPI_RECV_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
          call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,          &
#endif
                              PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
          ! From W
          tagc = 70
#ifdef _OPENACC
          call MPI_RECV_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
          call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,          &
#endif
                              PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          !$acc end host_data
          ireq = ireq + 1
       end if
       ! From S
       tagc = 40
       do j = 1, JS-1
          call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                               &
                              PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JA
          call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                               &
                              PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

       !--- To 8-Direction HALO communicate
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                           &
                          PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + ginfo(gid)%size2D_NS8*KA
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                           &
                          PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + ginfo(gid)%size2D_NS8*KA
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       if ( .not. PRC_TwoD ) then
          ! To W HALO
          tagc = 60
#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr)
          call MPI_SEND_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,          &
#endif
                              PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          ireq = ireq + 1
          ! To E HALO
          tagc = 70
#ifdef _OPENACC
          call MPI_SEND_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,          &
#endif
                              PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          !$acc end host_data
          ireq = ireq + 1
          ! To NW HALO
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                 PRC_next(PRC_NW), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To NE HALO
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                 PRC_next(PRC_NE), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SW HALO
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                 PRC_next(PRC_SW), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SE HALO
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                 PRC_next(PRC_SE), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       end if

    else ! non-periodic condition

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                               &
                                    PRC_next(PRC_SE), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                              &
                                    PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                              &
                                    PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From SW
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                  &
                                    PRC_next(PRC_SW), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                    PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                    PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NE
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                               &
                                    PRC_next(PRC_NE), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                              &
                                    PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                              &
                                    PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NW
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                  &
                                    PRC_next(PRC_NW), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                    PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                    PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
#ifdef _OPENACC
          ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
          !$acc host_data use_device(ptr)
#endif
          ! From E
          if ( PRC_HAS_E ) then
             tagc = 60
#ifdef _OPENACC
             call MPI_RECV_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
             call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,          &
#endif
                                 PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             ireq = ireq + 1
          endif
          ! From W
          if ( PRC_HAS_W ) then
             tagc = 70
#ifdef _OPENACC
             call MPI_RECV_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
             call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,          &
#endif
                                 PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             ireq = ireq + 1
          endif
          !$acc end host_data
       end if
       ! From S
       if ( PRC_HAS_S ) then
          tagc = 40
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                           &
                             PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_HAS_N ) then
          tagc = 50
          do j = JE+1, JA
             call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                           &
                             PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

       !--- To 8-Direction HALO communicate
       ! To N HALO
       if ( PRC_HAS_N ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                           &
                             PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             nsend = nsend + ginfo(gid)%size2D_NS8*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                           &
                             PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             nsend = nsend + ginfo(gid)%size2D_NS8*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr)
#endif
          ! To W HALO
          if ( PRC_HAS_W ) then
             tagc = 60
#ifdef _OPENACC
             call MPI_SEND_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
             call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,          &
#endif
                                 PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             nsend = nsend + ginfo(gid)%size2D_WE*KA
             ireq = ireq + 1
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             tagc = 70
#ifdef _OPENACC
             call MPI_SEND_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
             call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,          &
#endif
                                 PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
             nsend = nsend + ginfo(gid)%size2D_WE*KA
             ireq = ireq + 1
          endif
          !$acc end host_data
          ! To NW HALO
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                    PRC_next(PRC_NW), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                    PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                &
                                    PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To NE HALO
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_NE), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                              &
                                    PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                    PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SW HALO
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                    PRC_next(PRC_SW), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                 &
                                    PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                                &
                                    PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SE HALO
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_SE), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                              &
                                    PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                    PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
       end if

    endif

    ginfo(gid)%preq_cnt(seqid) = ireq - 1
    ginfo(gid)%psend_cnt(seqid) = nsend

    ! to finish initial processes of MPI
    nreq = ginfo(gid)%preq_cnt(seqid)
    do i = 1, 32
       call MPI_TESTALL( nreq, ginfo(gid)%preq_list(1:nreq,seqid), &
                         flag, MPI_STATUSES_IGNORE, ierr         )
    enddo

    !$acc end host_data

    return
  end subroutine vars8_init_mpi_pc

  !-----------------------------------------------------------------------------
  !> Range in y of the W/E halos exchanged in the layout 2: the inner rows on the side with a
  !! neighbor in y (its halo rows come with the whole rows from that neighbor and the corners from
  !! the diagonal one), and the halo rows too on the side without it. The neighbors in x have the
  !! same neighbors in y, so the range is the same in the sender and the receiver
  subroutine L2_range_WE( gid, jw0, jw1 )
    implicit none
    integer, intent(in)  :: gid
    integer, intent(out) :: jw0, jw1

    jw0 = 1
    if ( PRC_HAS_S ) jw0 = ginfo(gid)%JS
    jw1 = ginfo(gid)%JA
    if ( PRC_HAS_N ) jw1 = ginfo(gid)%JE

    return
  end subroutine L2_range_WE

  !-----------------------------------------------------------------------------
  !> Offsets in the packing buffer of the layout 2 of a vid: W (from W, to W), E, and the corners
  !! NW, NE, SW, SE (from or to the diagonal neighbor in that direction)
  subroutine L2_offsets( gid, off_E, off_NW, off_NE, off_SW, off_SE )
    implicit none
    integer, intent(in)  :: gid
    integer, intent(out) :: off_E  !< offset of the section of E (that of W is 0)
    integer, intent(out) :: off_NW, off_NE, off_SW, off_SE !< offsets of the sections of the corners

    integer :: KA, IHALO, JA, JHALO

    KA    = ginfo(gid)%KA
    IHALO = ginfo(gid)%IHALO
    JA    = ginfo(gid)%JA
    JHALO = ginfo(gid)%JHALO

    off_E  = KA * IHALO * JA
    off_NW = 2 * KA * IHALO * JA
    off_NE = off_NW + KA * IHALO * JHALO
    off_SW = off_NE + KA * IHALO * JHALO
    off_SE = off_SW + KA * IHALO * JHALO

    return
  end subroutine L2_offsets

  !-----------------------------------------------------------------------------
  !> Persistent requests of the 8-direction exchange in the layout 2 (COMM_USE_PC_LAYOUT2):
  !! N/S: the whole JHALO rows (i = 1:IA) from and to the variable, without packing;
  !! W/E: one packed message (rows in L2_range_WE); corners: one packed IHALO x JHALO block, only
  !! from and to the diagonal neighbors (where the neighbors in x and in y exist).
  !! The tags are those of the current layout without the row offset, one for each direction of
  !! the message, so that a neighbor in more than one direction (e.g., the process itself) is
  !! told apart. The W/E halos and the corners are unpacked after all the requests complete
  !! (wait_3D_mpi_pc): the N/S messages also write the corners, which the corner blocks overwrite.
  subroutine vars8_init_mpi_pc_L2(var, gid, seqid, pid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in) :: gid
    integer,  intent(in) :: seqid
    integer,  intent(in) :: pid   !< ID of the packing buffer

    integer :: ireq, tag
    integer :: ierr
    integer :: nsend
    logical :: flag

    integer :: KA
    integer :: IHALO
    integer :: JA, JS, JE, JHALO
    integer :: jw0, jw1
    integer :: n_NS, n_WE, n_C
    integer :: off_E, off_NW, off_NE, off_SW, off_SE

    integer :: nreq
    integer :: i

    real(RP), pointer, contiguous :: ptr(:)

    KA    = ginfo(gid)%KA
    IHALO = ginfo(gid)%IHALO
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    call L2_range_WE( gid, jw0, jw1 )
    call L2_offsets( gid, off_E, off_NW, off_NE, off_SW, off_SE )
    n_NS = ginfo(gid)%size2D_NS4 * KA       ! IA * JHALO * KA
    n_WE = KA * IHALO * ( jw1 - jw0 + 1 )
    n_C  = KA * IHALO * JHALO

    tag  = COMM_tag_base( seqid )
    ireq = 1
    nsend = 0

    !$acc host_data use_device(var)

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype_t,                 &
                        MPI_PROC_NULL, tag, ginfo(gid)%comm_pc,                 &
                        ginfo(gid)%preq_list(ginfo(gid)%nreq_max+1,seqid), ierr )

    !--- receive
    ! From S
    if ( PRC_HAS_S ) then
       call MPI_RECV_INIT( var(1,1,1), n_NS, COMM_datatype_t,                                                &
                           PRC_next(PRC_S), tag+40, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
       ireq = ireq + 1
    end if
    ! From N
    if ( PRC_HAS_N ) then
       call MPI_RECV_INIT( var(1,1,JE+1), n_NS, COMM_datatype_t,                                             &
                           PRC_next(PRC_N), tag+50, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
       ireq = ireq + 1
    end if

    !$acc end host_data

    if ( IHALO > 0 ) then
       ptr => ginfo(gid)%recvpack_L2(:,pid)
       !$acc host_data use_device(ptr)
       ! From E
       if ( PRC_HAS_E ) then
          call MPI_RECV_INIT( ptr(off_E+1), n_WE, COMM_datatype_t,                                              &
                              PRC_next(PRC_E), tag+60, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
       end if
       ! From W
       if ( PRC_HAS_W ) then
          call MPI_RECV_INIT( ptr(1), n_WE, COMM_datatype_t,                                                    &
                              PRC_next(PRC_W), tag+70, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
       end if
       ! From NW
       if ( PRC_HAS_N .and. PRC_HAS_W ) then
          call MPI_RECV_INIT( ptr(off_NW+1), n_C, COMM_datatype_t,                                            &
                              PRC_next(PRC_NW), tag+30, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
       end if
       ! From NE
       if ( PRC_HAS_N .and. PRC_HAS_E ) then
          call MPI_RECV_INIT( ptr(off_NE+1), n_C, COMM_datatype_t,                                            &
                              PRC_next(PRC_NE), tag+20, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
       end if
       ! From SW
       if ( PRC_HAS_S .and. PRC_HAS_W ) then
          call MPI_RECV_INIT( ptr(off_SW+1), n_C, COMM_datatype_t,                                            &
                              PRC_next(PRC_SW), tag+10, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
       end if
       ! From SE
       if ( PRC_HAS_S .and. PRC_HAS_E ) then
          call MPI_RECV_INIT( ptr(off_SE+1), n_C, COMM_datatype_t,                                            &
                              PRC_next(PRC_SE), tag+0,  ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          ireq = ireq + 1
       end if
       !$acc end host_data
    end if

    !--- send
    !$acc host_data use_device(var)
    ! To N
    if ( PRC_HAS_N ) then
       call MPI_SEND_INIT( var(1,1,JE-JHALO+1), n_NS, COMM_datatype_t,                                       &
                           PRC_next(PRC_N), tag+40, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
       nsend = nsend + n_NS
       ireq = ireq + 1
    end if
    ! To S
    if ( PRC_HAS_S ) then
       call MPI_SEND_INIT( var(1,1,JS), n_NS, COMM_datatype_t,                                               &
                           PRC_next(PRC_S), tag+50, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
       nsend = nsend + n_NS
       ireq = ireq + 1
    end if
    !$acc end host_data

    if ( IHALO > 0 ) then
       ptr => ginfo(gid)%sendpack_L2(:,pid)
       !$acc host_data use_device(ptr)
       ! To W
       if ( PRC_HAS_W ) then
          call MPI_SEND_INIT( ptr(1), n_WE, COMM_datatype_t,                                                    &
                              PRC_next(PRC_W), tag+60, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + n_WE
          ireq = ireq + 1
       end if
       ! To E
       if ( PRC_HAS_E ) then
          call MPI_SEND_INIT( ptr(off_E+1), n_WE, COMM_datatype_t,                                              &
                              PRC_next(PRC_E), tag+70, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + n_WE
          ireq = ireq + 1
       end if
       ! To NW
       if ( PRC_HAS_N .and. PRC_HAS_W ) then
          call MPI_SEND_INIT( ptr(off_NW+1), n_C, COMM_datatype_t,                                            &
                              PRC_next(PRC_NW), tag+0,  ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + n_C
          ireq = ireq + 1
       end if
       ! To NE
       if ( PRC_HAS_N .and. PRC_HAS_E ) then
          call MPI_SEND_INIT( ptr(off_NE+1), n_C, COMM_datatype_t,                                            &
                              PRC_next(PRC_NE), tag+10, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + n_C
          ireq = ireq + 1
       end if
       ! To SW
       if ( PRC_HAS_S .and. PRC_HAS_W ) then
          call MPI_SEND_INIT( ptr(off_SW+1), n_C, COMM_datatype_t,                                            &
                              PRC_next(PRC_SW), tag+20, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + n_C
          ireq = ireq + 1
       end if
       ! To SE
       if ( PRC_HAS_S .and. PRC_HAS_E ) then
          call MPI_SEND_INIT( ptr(off_SE+1), n_C, COMM_datatype_t,                                            &
                              PRC_next(PRC_SE), tag+30, ginfo(gid)%comm_pc, ginfo(gid)%preq_list(ireq,seqid), ierr )
          nsend = nsend + n_C
          ireq = ireq + 1
       end if
       !$acc end host_data
    end if

    ginfo(gid)%preq_cnt(seqid) = ireq - 1
    ginfo(gid)%psend_cnt(seqid) = nsend

    ! to finish initial processes of MPI
    nreq = ginfo(gid)%preq_cnt(seqid)
    do i = 1, 32
       call MPI_TESTALL( nreq, ginfo(gid)%preq_list(1:nreq,seqid), &
                         flag, MPI_STATUSES_IGNORE, ierr         )
    enddo

    return
  end subroutine vars8_init_mpi_pc_L2

  subroutine vars_3D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< atmospheric 3D variable to communication
    integer,  intent(in)    :: gid        !< grid ID
    integer,  intent(in)    :: vid        !< request ID


    integer :: ireq, tag

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer :: ierr
    integer :: nsend
#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

    tag  = COMM_tag_base( vid )
    ireq = 1
    nsend = 0

    call PROF_rapstart('COMM_mpi_post', 3)

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    !$acc host_data use_device(var) if(flag_device)

    !--- From 4-Direction HALO communicate
    ! From S
    if ( PRC_HAS_S ) then
       call MPI_IRECV( var(:,:,1:JS-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,                  &
                       PRC_next(PRC_S), tag+1, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif
    ! From N
    if ( PRC_HAS_N ) then
       call MPI_IRECV( var(:,:,JE+1:JA), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,                 &
                       PRC_next(PRC_N), tag+2, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif
    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
       !$acc host_data use_device(ptr) if(flag_device)
#endif
       ! From E
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_E), tag+3, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_W), tag+4, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       !$acc end host_data
    end if

    !$acc end host_data

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then
       call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
       call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, vid)
       call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)
    end if

    !$acc host_data use_device(var) if(flag_device)

    ! To N HALO
    if ( PRC_HAS_N ) then
       call MPI_ISEND( var(:,:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,           &
                       PRC_next(PRC_N), tag+1, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4*KA
       ireq = ireq + 1
    endif
    ! To S HALO
    if ( PRC_HAS_S ) then
       call MPI_ISEND( var(:,:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,           &
                       PRC_next(PRC_S), tag+2, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4*KA
       ireq = ireq + 1
    endif

    !$acc end host_data

    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc host_data use_device(ptr) if(flag_device)
#endif
       ! To W HALO
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_W), tag+3, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          ireq = ireq + 1
       endif
       ! To E HALO
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_E), tag+4, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          ireq = ireq + 1
       endif

       !$acc end host_data
    end if

    call PROF_rapend  ('COMM_mpi_post', 3)

    ginfo(gid)%send_cnt(vid) = nsend ! counted in COMM_wait

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars_3D_mpi

  subroutine vars_3D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< atmospheric 3D variable to communication
    integer,  intent(in)    :: gid        !< grid ID
    integer,  intent(in)    :: vid        !< request ID

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer(kind=MPI_ADDRESS_KIND) :: disp

    integer :: ierr
    integer :: nsend
#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
#endif
    !---------------------------------------------------------------------------

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    !$acc data copyin(var)

    nsend = 0
    ginfo(gid)%put_corner(vid) = .false.
    call PROF_rapstart('COMM_mpi_post', 3)
    call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(vid), ierr )

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then
       call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
       call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, vid)
       call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)
    end if

    !$acc host_data use_device(var)

    ! To N HALO
    if ( PRC_HAS_N ) then
       disp = 0
       call MPI_PUT( var(1,1,JE-JHALO+1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t, &
                     PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS4*KA, COMM_datatype_t, &
                     ginfo(gid)%win_packNS(vid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4*KA
    endif
    ! To S HALO
    if ( PRC_HAS_S ) then
       disp = KA * IA * JHALO
       call MPI_PUT( var(1,1,JS), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t, &
                     PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS4*KA, COMM_datatype_t, &
                     ginfo(gid)%win_packNS(vid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4*KA
    endif

    !$acc end host_data

    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc host_data use_device(ptr)
#endif

       ! To W HALO
       if ( PRC_HAS_W ) then
          disp = int( ginfo(gid)%size2D_WE, kind=MPI_ADDRESS_KIND ) * KA
#ifdef _OPENACC
          call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                        PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
       endif
       ! To E HALO
       if ( PRC_HAS_E ) then
          disp = 0
#ifdef _OPENACC
          call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                        PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
       endif

       !$acc end host_data
    end if

    call MPI_Win_complete( ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_complete( ginfo(gid)%win_packNS(vid), ierr )
    call PROF_rapend  ('COMM_mpi_post', 3)

    ! the data is transferred until MPI_Win_complete
    call count_bytes( 'COMM_vars', 'COMM_mpi_post', nsend )

    !$acc end data

    return
  end subroutine vars_3D_mpi_onesided

  subroutine vars8_3D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: ireq, tag, tagc

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer :: ierr
    integer :: nsend
    integer :: j
#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

    tag  = COMM_tag_base( vid )
    ireq = 1
    nsend = 0

    call PROF_rapstart('COMM_mpi_post', 3)

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !$acc host_data use_device(var) if(flag_device)

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          tagc = 0
          do j = 1, JS-1
             call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                             PRC_next(PRC_SE), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From SW
          tagc = 10
          do j = 1, JS-1
             call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                            &
                             PRC_next(PRC_SW), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NE
          tagc = 20
          do j = JE+1, JA
             call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                             PRC_next(PRC_NE), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NW
          tagc = 30
          do j = JE+1, JA
             call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                            &
                             PRC_next(PRC_NW), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
#ifdef _OPENACC
          ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
          !$acc host_data use_device(ptr) if(flag_device)
#endif
          ! From E
          tagc = 60
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#endif
                          PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! From W
          tagc = 70
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#endif
                          PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          !$acc end host_data
       end if
       ! From S
       tagc = 40
       do j = 1, JS-1
          call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                         &
                          PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JA
          call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                         &
                          PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

       !--- To 8-Direction HALO communicate
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                         &
                          PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_NS8*KA
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                         &
                          PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_NS8*KA
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
          call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)
          call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

          !$acc host_data use_device(var) if(flag_device)

          ! To NW HALO
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                             PRC_next(PRC_NW), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To NE HALO
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                   &
                             PRC_next(PRC_NE), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SW HALO
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                             PRC_next(PRC_SW), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SE HALO
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                   &
                             PRC_next(PRC_SE), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo

          !$acc end host_data

#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr) if(flag_device)
#endif
          ! To W HALO
          tagc = 60
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#endif
                          PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          ireq = ireq + 1
          ! To E HALO
          tagc = 70
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#endif
                          PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          !$acc end host_data
          ireq = ireq + 1

       end if

    else ! non-periodic condition

       !$acc host_data use_device(var) if(flag_device)

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                PRC_next(PRC_SE), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From SW
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                            &
                                PRC_next(PRC_SW), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NE
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                PRC_next(PRC_NE), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NW
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                            &
                                PRC_next(PRC_NW), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
#ifdef _OPENACC
          ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
          !$acc host_data use_device(ptr) if(flag_device)
#endif
          ! From E
          if ( PRC_HAS_E ) then
             tagc = 60
#ifdef _OPENACC
             call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
             call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#endif
                             PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! From W
          if ( PRC_HAS_W ) then
             tagc = 70
#ifdef _OPENACC
             call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
             call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#endif
                             PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          !$acc end host_data
       end if
       ! From S
       if ( PRC_HAS_S ) then
          tagc = 40
          do j = 1, JS-1
             call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                         &
                             PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_HAS_N ) then
          tagc = 50
          do j = JE+1, JA
             call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                         &
                             PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

       !--- To 8-Direction HALO communicate
       ! To N HALO
       if ( PRC_HAS_N ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                         &
                             PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_NS8*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                         &
                             PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_NS8*KA
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
          call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)
          call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

          !$acc host_data use_device(var) if(flag_device)

          ! To NW HALO
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                PRC_next(PRC_NW), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                          &
                                PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To NE HALO
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                   &
                                PRC_next(PRC_NE), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_N), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                  &
                                PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SW HALO
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                PRC_next(PRC_SW), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                          &
                                PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SE HALO
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                   &
                                PRC_next(PRC_SE), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_S), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                  &
                                PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif

          !$acc end host_data

#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr) if(flag_device)
#endif

          ! To W HALO
          if ( PRC_HAS_W ) then
             tagc = 60
#ifdef _OPENACC
             call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
             call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#endif
                             PRC_next(PRC_W), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_WE*KA
             ireq = ireq + 1
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             tagc = 70
#ifdef _OPENACC
             call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
             call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#endif
                             PRC_next(PRC_E), tag+tagc, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_WE*KA
             ireq = ireq + 1
          endif
          !$acc end host_data

       end if

    endif

    call PROF_rapend  ('COMM_mpi_post', 3)

    ginfo(gid)%send_cnt(vid) = nsend ! counted in COMM_wait

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars8_3D_mpi

  subroutine vars8_3D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer(kind=MPI_ADDRESS_KIND) :: disp

    integer :: ierr
    integer :: nsend
    integer :: j
#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
#endif
    !---------------------------------------------------------------------------

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    !$acc data copyin(var)
    nsend = 0
    ginfo(gid)%put_corner(vid) = .true.
    call PROF_rapstart('COMM_mpi_post', 3)
    call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(vid), ierr )

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !$acc host_data use_device(var)

       !--- To 8-Direction HALO communicate
       ! To N HALO
       do j = JE-JHALO+1, JE
          disp = KA * ( IHALO + IA * ( j - JE+JHALO-1 ) )
          call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                        PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packNS(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_NS8*KA
       enddo
       ! To S HALO
       do j = JS, JS+JHALO-1
          disp = KA * ( IHALO + IA * ( j - JS + JHALO ) )
          call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                        PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packNS(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_NS8*KA
       enddo

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
          call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)
          call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

          !$acc host_data use_device(var)
          ! To NW HALO
          do j = JE-JHALO+1, JE
             disp = KA * ( IE + IA * ( j - JE+JHALO-1 ) )
             call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           PRC_next(PRC_NW), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
          enddo
          ! To NE HALO
          do j = JE-JHALO+1, JE
             disp = KA * ( IA * ( j - JE+JHALO-1 ) )
             call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           PRC_next(PRC_NE), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
          enddo
          ! To SW HALO
          do j = JS, JS+JHALO-1
             disp = KA * ( IE + IA * ( j - JS + JHALO ) )
             call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           PRC_next(PRC_SW), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
          enddo
          ! To SE HALO
          do j = JS, JS+JHALO-1
             disp = KA * ( IA * ( j - JS + JHALO ) )
             call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           PRC_next(PRC_SE), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C*KA
          enddo

          !$acc end host_data

#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr)
#endif
          ! To W HALO
          disp = int( ginfo(gid)%size2D_WE, kind=MPI_ADDRESS_KIND ) * KA
#ifdef _OPENACC
          call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                        PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          ! To E HALO
          disp = 0
#ifdef _OPENACC
          call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                        PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE*KA
          !$acc end host_data

       end if

    else ! non-periodic condition

       !$acc host_data use_device(var)

       !--- To 8-Direction HALO communicate
       ! To N HALO
       if ( PRC_HAS_N ) then
          do j = JE-JHALO+1, JE
             disp = KA * ( IHALO + IA * ( j - JE+JHALO-1 ) )
             call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                           PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_NS8*KA
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          do j = JS, JS+JHALO-1
             disp = KA * ( IHALO + IA * ( j - JS + JHALO ) )
             call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                           PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_NS8*KA
          enddo
       endif

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
          call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)
          call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

          !$acc host_data use_device(var)

          ! To NW HALO
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             do j = JE-JHALO+1, JE
                disp = KA * ( IE + IA * ( j - JE+JHALO-1 ) )
                call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_NW), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          else if ( PRC_HAS_N ) then
             do j = JE-JHALO+1, JE
                disp = KA * ( IA * ( j - JE+JHALO-1 ) )
                call MPI_PUT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_N), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          else if ( PRC_HAS_W ) then
             do j = JE+1, JA
                disp = KA * ( IE + IA * ( j - JE-1 + JHALO ) )
                call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_W), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          endif
          ! To NE HALO
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             do j = JE-JHALO+1, JE
                disp = KA * ( IA * ( j - JE+JHALO-1 ) )
                call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_NE), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          else if ( PRC_HAS_N ) then
             do j = JE-JHALO+1, JE
                disp = KA * ( IE + IA * ( j - JE+JHALO-1 ) )
                call MPI_PUT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_N), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          else if ( PRC_HAS_E ) then
             do j = JE+1, JA
                disp = KA * IA * ( j - JE-1 + JHALO )
                call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_E), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          endif
          ! To SW HALO
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             do j = JS, JS+JHALO-1
                disp = KA * ( IE + IA * ( j - JS + JHALO ) )
                call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_SW), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          else if ( PRC_HAS_S ) then
             do j = JS, JS+JHALO-1
                disp = KA * ( IA * ( j - JS + JHALO ) )
                call MPI_PUT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_S), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          else if ( PRC_HAS_W ) then
             do j = 1, JS-1
                disp = KA * ( IE + IA * (j-1) )
                call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_W), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          endif
          ! To SE HALO
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             do j = JS, JS+JHALO-1
                disp = KA * ( IA * ( j - JS + JHALO ) )
                call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_SE), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          else if ( PRC_HAS_S ) then
             do j = JS, JS+JHALO-1
                disp = KA * ( IE + IA * ( j - JS + JHALO ) )
                call MPI_PUT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_S), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          else if ( PRC_HAS_E ) then
             do j = 1, JS-1
                disp = KA * ( IA * ( j - 1 ) )
                call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_E), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C*KA
             enddo
          endif

          !$acc end host_data

#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr)
#endif

          ! To W HALO
          if ( PRC_HAS_W ) then
             disp = int( ginfo(gid)%size2D_WE, kind=MPI_ADDRESS_KIND ) * KA
#ifdef _OPENACC
             call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
             call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                           PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packWE(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_WE*KA
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             disp = 0
#ifdef _OPENACC
             call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
             call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                           PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packWE(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_WE*KA
          endif
          !$acc end host_data

       end if

    endif

    call MPI_Win_complete( ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_complete( ginfo(gid)%win_packNS(vid), ierr )
    call PROF_rapend  ('COMM_mpi_post', 3)

    ! the data is transferred until MPI_Win_complete
    call count_bytes( 'COMM_vars', 'COMM_mpi_post', nsend )

    !$acc end data

    return
  end subroutine vars8_3D_mpi_onesided

  subroutine vars_2D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer :: ireq, tag
    integer :: ierr
    integer :: nsend
#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    tag  = COMM_tag_base( vid )
    ireq = 1
    nsend = 0

    call PROF_rapstart('COMM_mpi_post', 3)

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    !$acc host_data use_device(var) if(flag_device)

    !--- From 4-Direction HALO communicate
    ! From S
    if ( PRC_HAS_S ) then
       call MPI_IRECV( var(:,1:JS-1), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                       PRC_next(PRC_S), tag+1, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif
    ! From N
    if ( PRC_HAS_N ) then
       call MPI_IRECV( var(:,JE+1:JA), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                       PRC_next(PRC_N), tag+2, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif

    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
       !$acc host_data use_device(ptr) if(flag_device)
#endif
       ! From E
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t,    &
#endif
                          PRC_next(PRC_E), tag+3, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t,    &
#endif
                          PRC_next(PRC_W), tag+4, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       !$acc end host_data
    end if

    !$acc end host_data

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then

       call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
       call packWE_2D( IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, vid)
       call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc host_data use_device(ptr) if(flag_device)
#endif

       ! To W HALO communicate
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t,    &
#endif
                          PRC_next(PRC_W), tag+3, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE
          ireq = ireq + 1
       endif
       ! To E HALO communicate
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t,    &
#endif
                          PRC_next(PRC_E), tag+4, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE
          ireq = ireq + 1
       endif

       !$acc end host_data

    end if

    !$acc host_data use_device(var) if(flag_device)

    ! To N HALO communicate
    if ( PRC_HAS_N ) then
       call MPI_ISEND( var(:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
            PRC_next(PRC_N), tag+1, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4
       ireq = ireq + 1
    endif
    ! To S HALO communicate
    if ( PRC_HAS_S ) then
       call MPI_ISEND( var(:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
            PRC_next(PRC_S), tag+2, ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4
       ireq = ireq + 1
    endif

    !$acc end host_data

    call PROF_rapend  ('COMM_mpi_post', 3)

    ginfo(gid)%send_cnt(vid) = nsend ! counted in COMM_wait

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars_2D_mpi

  subroutine vars_2D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer(kind=MPI_ADDRESS_KIND) :: disp

    integer :: ierr
    integer :: nsend
#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
#endif
    !---------------------------------------------------------------------------

    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    !$acc data copyin(var)

    nsend = 0
    ginfo(gid)%put_corner(vid) = .false.
    call PROF_rapstart('COMM_mpi_post', 3)
    call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(vid), ierr )

    !--- To 4-Direction HALO communicate

    if ( .not. PRC_TwoD ) then

       call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
       call packWE_2D( IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, vid)
       call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc host_data use_device(ptr)
#endif

       ! To W HALO communicate
       if ( PRC_HAS_W ) then
          disp = int( ginfo(gid)%size2D_WE, kind=MPI_ADDRESS_KIND ) * ginfo(gid)%KA
#ifdef _OPENACC
          call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                        PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE
       endif
       ! To E HALO communicate
       if ( PRC_HAS_E ) then
          disp = 0
#ifdef _OPENACC
          call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                        PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE
       endif

       !$acc end host_data

    end if

    !$acc host_data use_device(var)

    ! To N HALO communicate
    if ( PRC_HAS_N ) then
       disp = 0
       call MPI_PUT( var(1,JE-JHALO+1), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                     PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                     ginfo(gid)%win_packNS(vid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4
    endif
    ! To S HALO communicate
    if ( PRC_HAS_S ) then
       disp = IA * JHALO
       call MPI_PUT( var(1,JS), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                     PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                     ginfo(gid)%win_packNS(vid), ierr )
       nsend = nsend + ginfo(gid)%size2D_NS4
    endif

    !$acc end host_data

    call MPI_Win_complete( ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_complete( ginfo(gid)%win_packNS(vid), ierr )
    call PROF_rapend  ('COMM_mpi_post', 3)

    ! the data is transferred until MPI_Win_complete
    call count_bytes( 'COMM_vars', 'COMM_mpi_post', nsend )

    !$acc end data

    return
  end subroutine vars_2D_mpi_onesided

  subroutine vars8_2D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer :: ireq, tag, tagc

    integer :: ierr
    integer :: nsend
    integer :: j
#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    tag  = COMM_tag_base( vid )
    ireq = 1
    nsend = 0

    call PROF_rapstart('COMM_mpi_post', 3)

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- From 8-Direction HALO communicate
        !$acc host_data use_device(var) if(flag_device)

        if ( .not. PRC_TwoD ) then
           ! From SE
           tagc = 0
           do j = 1, JS-1
              call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                              COMM_datatype_t, PRC_next(PRC_SE), tag+tagc,         &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From SW
           tagc = 10
           do j = 1, JS-1
              call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                      &
                              COMM_datatype_t, PRC_next(PRC_SW), tag+tagc,         &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From NE
           tagc = 20
           do j = JE+1, JA
              call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                              COMM_datatype_t, PRC_next(PRC_NE), tag+tagc,         &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From NW
           tagc = 30
           do j = JE+1, JA
              call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                      &
                              COMM_datatype_t, PRC_next(PRC_NW), tag+tagc,         &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
#ifdef _OPENACC
           ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
           !$acc host_data use_device(ptr) if(flag_device)
#endif
           ! From E
#ifdef _OPENACC
           call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE, &
#else
           call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE, &
#endif
                           COMM_datatype_t, PRC_next(PRC_E), tag+60,                &
                           ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr     )
           ireq = ireq + 1
           ! From W
#ifdef _OPENACC
           call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE, &
#else
           call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE, &
#endif
                           COMM_datatype_t, PRC_next(PRC_W), tag+70,                &
                           ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr     )
           ireq = ireq + 1
           !$acc end host_data
        end if
        ! From S
        tagc = 40
        do j = 1, JS-1
            call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,                    &
                            COMM_datatype_t, PRC_next(PRC_S), tag+tagc,          &
                            ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
        enddo
        ! From N
        tagc = 50
        do j = JE+1, JA
            call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,                    &
                            COMM_datatype_t, PRC_next(PRC_N), tag+tagc,          &
                            ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo


        !--- To 8-Direction HALO communicate

        ! To N HALO communicate
        tagc = 40
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,                    &
                            COMM_datatype_t, PRC_next(PRC_N), tag+tagc,          &
                            ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
            nsend = nsend + ginfo(gid)%size2D_NS8
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To S HALO communicate
        tagc = 50
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,                    &
                            COMM_datatype_t, PRC_next(PRC_S), tag+tagc,          &
                            ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
            nsend = nsend + ginfo(gid)%size2D_NS8
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        !$acc end host_data

        if ( .not. PRC_TwoD ) then

           call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
           call packWE_2D( IA, IS, IE, JA, JS, JE, &
                           IHALO, &
                           var, gid, vid)
           call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

           !$acc host_data use_device(var) if(flag_device)
#ifdef _OPENACC
           ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
           !$acc host_data use_device(ptr) if(flag_device)
#endif

           ! To W HALO communicate
#ifdef _OPENACC
           call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE, &
#else
           call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, &
#endif
                           COMM_datatype_t, PRC_next(PRC_W), tag+60,                &
                           ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr     )
           nsend = nsend + ginfo(gid)%size2D_WE
           ireq = ireq + 1

           ! To E HALO communicate
#ifdef _OPENACC
           call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE, &
#else
           call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, &
#endif
                           COMM_datatype_t, PRC_next(PRC_E), tag+70,                &
                           ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr     )
           nsend = nsend + ginfo(gid)%size2D_WE
           ireq = ireq + 1
           !$acc end host_data

           ! To NW HALO communicate
           tagc = 0
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                     &
                              COMM_datatype_t, PRC_next(PRC_NW), tag+tagc,         &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              nsend = nsend + ginfo(gid)%size2D_4C
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To NE HALO communicate
           tagc = 10
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,             &
                              COMM_datatype_t, PRC_next(PRC_NE), tag+tagc,         &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              nsend = nsend + ginfo(gid)%size2D_4C
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To SW HALO communicate
           tagc = 20
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                     &
                              COMM_datatype_t, PRC_next(PRC_SW), tag+tagc,         &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              nsend = nsend + ginfo(gid)%size2D_4C
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To SE HALO communicate
           tagc = 30
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,             &
                              COMM_datatype_t, PRC_next(PRC_SE), tag+tagc,         &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              nsend = nsend + ginfo(gid)%size2D_4C
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           !$acc end host_data

        end if

    else
    !--- non-periodic condition
        !--- From 8-Direction HALO communicate

        !$acc host_data use_device(var) if(flag_device)

        if ( .not. PRC_TwoD ) then
           ! From SE
           if ( PRC_HAS_S .AND. PRC_HAS_E ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_SE), tag+tagc,         &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_S), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_E), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From SW
           if ( PRC_HAS_S .AND. PRC_HAS_W ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                      &
                                 COMM_datatype_t, PRC_next(PRC_SW), tag+tagc,         &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                      &
                                 COMM_datatype_t, PRC_next(PRC_S), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                      &
                                 COMM_datatype_t, PRC_next(PRC_W), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From NE
           if ( PRC_HAS_N .AND. PRC_HAS_E ) then
              tagc = 20
              do j = JE+1, JE+JHALO
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_NE), tag+tagc,         &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 20
              do j = JE+1, JA
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_N), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 20
              do j = JE+1, JA
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_E), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From NW
           if ( PRC_HAS_N .AND. PRC_HAS_W ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                      &
                                 COMM_datatype_t, PRC_next(PRC_NW), tag+tagc,         &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                      &
                                 COMM_datatype_t, PRC_next(PRC_N), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                      &
                                 COMM_datatype_t, PRC_next(PRC_W), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

#ifdef _OPENACC
           ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
           !$acc host_data use_device(ptr) if(flag_device)
#endif
           ! From E
           if ( PRC_HAS_E ) then
#ifdef _OPENACC
              call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE, &
#else
              call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE, &
#endif
                              COMM_datatype_t, PRC_next(PRC_E), tag+60,                &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr     )
              ireq = ireq + 1
           endif

           ! From W
           if ( PRC_HAS_W ) then
#ifdef _OPENACC
              call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE, &
#else
              call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE, &
#endif
                              COMM_datatype_t, PRC_next(PRC_W), tag+70,                &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr     )
              ireq = ireq + 1
           endif
           !$acc end host_data

        end if

        ! From S
        if ( PRC_HAS_S ) then
           tagc = 40
           do j = 1, JS-1
              call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,                    &
                              COMM_datatype_t, PRC_next(PRC_S), tag+tagc,          &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        ! From N
        if ( PRC_HAS_N ) then
           tagc = 50
           do j = JE+1, JA
              call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,                    &
                              COMM_datatype_t, PRC_next(PRC_N), tag+tagc,          &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif


        !! RECEIVE

        ! To N HALO communicate
        if ( PRC_HAS_N ) then
           tagc = 40
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,                    &
                              COMM_datatype_t, PRC_next(PRC_N), tag+tagc,          &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              nsend = nsend + ginfo(gid)%size2D_NS8
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        ! To S HALO communicate
        if ( PRC_HAS_S ) then
           tagc = 50
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,                    &
                              COMM_datatype_t, PRC_next(PRC_S), tag+tagc,          &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
              nsend = nsend + ginfo(gid)%size2D_NS8
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        !$acc end host_data

        if ( .not. PRC_TwoD ) then

           call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
           call packWE_2D( IA, IS, IE, JA, JS, JE, &
                           IHALO, &
                           var, gid, vid)
           call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

           !$acc host_data use_device(var) if(flag_device)
#ifdef _OPENACC
           ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
           !$acc host_data use_device(ptr) if(flag_device)
#endif

           ! To W HALO communicate
           if ( PRC_HAS_W ) then
#ifdef _OPENACC
              call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE, &
#else
              call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, &
#endif
                              COMM_datatype_t, PRC_next(PRC_W), tag+60,                &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr     )
              nsend = nsend + ginfo(gid)%size2D_WE
              ireq = ireq + 1
           endif

           ! To E HALO communicate
           if ( PRC_HAS_E ) then
#ifdef _OPENACC
              call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE, &
#else
              call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, &
#endif
                              COMM_datatype_t, PRC_next(PRC_E), tag+70,                &
                              ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr     )
              nsend = nsend + ginfo(gid)%size2D_WE
              ireq = ireq + 1
           endif
           !$acc end host_data

           ! To NW HALO communicate
           if ( PRC_HAS_N .AND. PRC_HAS_W ) then
              tagc = 0
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                     &
                                 COMM_datatype_t, PRC_next(PRC_NW), tag+tagc,         &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 10
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(1,j), ginfo(gid)%size2D_4C,                      &
                                 COMM_datatype_t, PRC_next(PRC_N), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 20
              do j = JE+1, JA
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                     &
                                 COMM_datatype_t, PRC_next(PRC_W), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To NE HALO communicate
           if ( PRC_HAS_N .AND. PRC_HAS_E ) then
              tagc = 10
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,             &
                                 COMM_datatype_t, PRC_next(PRC_NE), tag+tagc,         &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 0
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_N), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,             &
                                 COMM_datatype_t, PRC_next(PRC_E), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To SW HALO communicate
           if ( PRC_HAS_S .AND. PRC_HAS_W ) then
              tagc = 20
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                     &
                                 COMM_datatype_t, PRC_next(PRC_SW), tag+tagc,         &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 30
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(1,j), ginfo(gid)%size2D_4C,                      &
                                 COMM_datatype_t, PRC_next(PRC_S), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                     &
                                 COMM_datatype_t, PRC_next(PRC_W), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To SE HALO communicate
           if ( PRC_HAS_S .AND. PRC_HAS_E ) then
              tagc = 30
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,             &
                                 COMM_datatype_t, PRC_next(PRC_SE), tag+tagc,         &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 20
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IE+1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_S), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,             &
                                 COMM_datatype_t, PRC_next(PRC_E), tag+tagc,          &
                                 ginfo(gid)%comm, ginfo(gid)%req_list(ireq,vid), ierr )
                 nsend = nsend + ginfo(gid)%size2D_4C
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           !$acc end host_data

        end if

    endif

    call PROF_rapend  ('COMM_mpi_post', 3)

    ginfo(gid)%send_cnt(vid) = nsend ! counted in COMM_wait

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars8_2D_mpi

  subroutine vars8_2D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: IA, IS, IE, IHALO
    integer :: JA, JS, JE, JHALO

    integer(kind=MPI_ADDRESS_KIND) :: disp

    integer :: ierr
    integer :: nsend
    integer :: j
#ifdef _OPENACC
    real(RP), pointer, contiguous :: ptr(:,:)
#endif
    !---------------------------------------------------------------------------

    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$acc data copyin(var)

    nsend = 0
    ginfo(gid)%put_corner(vid) = .true.
    call PROF_rapstart('COMM_mpi_post', 3)
    call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(vid), ierr )

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition

       !--- To 8-Direction HALO communicate

       !$acc host_data use_device(var)

       ! To N HALO communicate
       do j = JE-JHALO+1, JE
          disp = IHALO + IA * ( j - JE+JHALO-1 )
          call MPI_PUT( var(IS,j), ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                        PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                        ginfo(gid)%win_packNS(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_NS8
       enddo
       ! To S HALO communicate
       do j = JS, JS+JHALO-1
          disp = IHALO + IA * ( j - JS + JHALO )
          call MPI_PUT( var(IS,j), ginfo(gid)%size2D_NS8, COMM_datatype_t, &
               PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS8, COMM_datatype_t, &
               ginfo(gid)%win_packNS(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_NS8
       enddo

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
          call packWE_2D( IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)
          call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

          !$acc host_data use_device(var)
#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr)
#endif

          ! To W HALO communicate
          disp = int( ginfo(gid)%size2D_WE, kind=MPI_ADDRESS_KIND ) * ginfo(gid)%KA
#ifdef _OPENACC
          call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                        PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE
          ! To E HALO communicate
          disp = 0
#ifdef _OPENACC
          call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
               PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
               ginfo(gid)%win_packWE(vid), ierr )
          nsend = nsend + ginfo(gid)%size2D_WE
          !$acc end host_data
          ! To NW HALO communicate
          do j = JE-JHALO+1, JE
             disp = IE + IA * ( j - JE+JHALO-1 )
             call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           PRC_next(PRC_NW), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C
          enddo
          ! To NE HALO communicate
          do j = JE-JHALO+1, JE
             disp = IA * ( j - JE+JHALO-1 )
             call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           PRC_next(PRC_NE), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C
          enddo
          ! To SW HALO communicate
          do j = JS, JS+JHALO-1
             disp = IE + IA * ( j - JS + JHALO )
             call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           PRC_next(PRC_SW), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C
          enddo
          ! To SE HALO communicate
          do j = JS, JS+JHALO-1
             disp = IA * ( j - JS + JHALO )
             call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           PRC_next(PRC_SE), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_4C
          enddo

          !$acc end host_data

       end if
    else
    !--- non-periodic condition

       !$acc host_data use_device(var)

       ! To N HALO communicate
       if ( PRC_HAS_N ) then
          do j = JE-JHALO+1, JE
             disp = IHALO + IA * ( j - JE+JHALO-1 )
             call MPI_PUT( var(IS,j), ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                           PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_NS8
          enddo
       endif
       ! To S HALO communicate
       if ( PRC_HAS_S ) then
          do j = JS, JS+JHALO-1
             disp = IHALO + IA * ( j - JS + JHALO )
             call MPI_PUT( var(IS,j), ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                           PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_NS8
          enddo
       endif

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call PROF_rapend  ('COMM_mpi_post', 3, disable_barrier=.true.)
          call packWE_2D( IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)
          call PROF_rapstart('COMM_mpi_post', 3, disable_barrier=.true.)

          !$acc host_data use_device(var)
#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr)
#endif

          ! To W HALO communicate
          if ( PRC_HAS_W ) then
             disp = int( ginfo(gid)%size2D_WE, kind=MPI_ADDRESS_KIND ) * ginfo(gid)%KA
#ifdef _OPENACC
             call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
             call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                           PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                           ginfo(gid)%win_packWE(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_WE
          endif
          ! To E HALO communicate
          if ( PRC_HAS_E ) then
             disp = 0
#ifdef _OPENACC
             call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
             call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                           PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                           ginfo(gid)%win_packWE(vid), ierr )
             nsend = nsend + ginfo(gid)%size2D_WE
          endif
          !$acc end host_data
          ! To NW HALO communicate
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             do j = JE-JHALO+1, JE
                disp = IE + IA * ( j - JE+JHALO-1 )
                call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_NW), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          else if ( PRC_HAS_N ) then
             do j = JE-JHALO+1, JE
                disp = IA * ( j - JE+JHALO-1 )
                call MPI_PUT( var(1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_N), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          else if ( PRC_HAS_W ) then
             do j = JE+1, JA
                disp = IE + IA * ( j - JE-1 + JHALO )
                call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_W), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          endif
          ! To NE HALO communicate
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             do j = JE-JHALO+1, JE
                disp = IA * ( j - JE+JHALO-1 )
                call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_NE), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          else if ( PRC_HAS_N ) then
             do j = JE-JHALO+1, JE
                disp = IE + IA * ( j - JE+JHALO-1 )
                call MPI_PUT( var(IE+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_N), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          else if ( PRC_HAS_E ) then
             do j = JE+1, JA
                disp = IA * ( j - JE-1 + JHALO )
                call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_E), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          endif
          ! To SW HALO communicate
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             do j = JS, JS+JHALO-1
                disp = IE + IA * ( j - JS + JHALO )
                call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_SW), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          else if ( PRC_HAS_S ) then
             do j = JS, JS+JHALO-1
                disp = IA * ( j - JS + JHALO )
                call MPI_PUT( var(1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_S), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          else if ( PRC_HAS_W ) then
             do j = 1, JS-1
                disp = IE + IA * ( j - 1 )
                call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_W), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          endif
          ! To SE HALO communicate
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             do j = JS, JS+JHALO-1
                disp = IA * ( j - JS + JHALO )
                call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_SE), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          else if ( PRC_HAS_S ) then
             do j = JS, JS+JHALO-1
                disp = IE + IA * ( j - JS + JHALO )
                call MPI_PUT( var(IE+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_S), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          else if ( PRC_HAS_E ) then
             do j = 1, JS-1
                disp = IA * ( j - 1 )
                call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_E), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
                nsend = nsend + ginfo(gid)%size2D_4C
             enddo
          endif

          !$acc end host_data
       end if

    endif

    call MPI_Win_complete( ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_complete( ginfo(gid)%win_packNS(vid), ierr )
    call PROF_rapend  ('COMM_mpi_post', 3)

    ! the data is transferred until MPI_Win_complete
    call count_bytes( 'COMM_vars', 'COMM_mpi_post', nsend )

    !$acc end data

    return
  end subroutine vars8_2D_mpi_onesided

  subroutine vars_3D_mpi_pc(var, gid, vid, pid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid        !< sequential number of the persistent variable
    integer, intent(in)     :: pid        !< ID of the packing buffer of the variable

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO
    integer :: ig
    integer :: off(G_NDIR) ! offsets of the sections of the packed exchange

    integer :: ierr
    !---------------------------------------------------------------------------

    KA = ginfo(gid)%KA
    IA = ginfo(gid)%IA
    IS = ginfo(gid)%IS
    IE = ginfo(gid)%IE
    JA = ginfo(gid)%JA
    JS = ginfo(gid)%JS
    JE = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO

    ig = ginfo(gid)%pc_group(vid)
    if ( ig > 0 ) then
       call COMM_group_post( gid, vid, ig )
       if ( ginfo(gid)%grp(ig)%state == COMM_PC_GROUP_BATCH ) then
          ! the member is packed into its slot; all the messages of the group are started when its
          ! last member is posted. This is decided by the group, not by the layout of the grid
          ! (a grid with PRC_TwoD and IHALO > 0 does not use the layout 2)
          call G_offsets( gid, ginfo(gid)%grp(ig)%nmem, ginfo(gid)%pc_slot(vid), .false., off )
          call pack8_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, ginfo(gid)%JHALO,  &
                         var, ginfo(gid)%grp(ig)%sendbuf, ginfo(gid)%grp(ig)%nbuf, &
                         gid, .true., off )
          if ( ginfo(gid)%grp(ig)%nposted == ginfo(gid)%grp(ig)%nmem ) then
             ! the packing kernels of all the members are waited once, before the messages are started
             ! (COMM_pack of pack8_3D is the time to queue the kernel only)
             call PROF_rapstart('COMM_pack_wait', 3)
             !$acc wait
             call PROF_rapend  ('COMM_pack_wait', 3)
             call PROF_rapstart('COMM_mpi_post', 3)
             call MPI_STARTALL( ginfo(gid)%grp(ig)%nreq, ginfo(gid)%grp(ig)%req(1:ginfo(gid)%grp(ig)%nreq), ierr )
             call PROF_rapend  ('COMM_mpi_post', 3)
          end if
          return
       end if
    end if

    if ( ginfo(gid)%pc_layout2 .and. ginfo(gid)%pc_ndir(vid) == 8 ) then
       if ( IHALO > 0 ) then
          ! the W/E halos and the corners (the N/S rows are sent from the variable)
          off(G_N:G_W) = 0
          call L2_offsets( gid, off(G_E), off(G_NW), off(G_NE), off(G_SW), off(G_SE) )
          call pack8_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, ginfo(gid)%JHALO, &
                         var, ginfo(gid)%sendpack_L2(:,pid), ginfo(gid)%size_L2, &
                         gid, .false., off )
       end if
    else if ( .not. PRC_TwoD ) then
       call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, pid)
    end if

    call PROF_rapstart('COMM_pack_wait', 3)
    !$acc wait
    call PROF_rapend  ('COMM_pack_wait', 3)
    call PROF_rapstart('COMM_mpi_post', 3)
    call MPI_STARTALL(ginfo(gid)%preq_cnt(vid), ginfo(gid)%preq_list(1:ginfo(gid)%preq_cnt(vid),vid), ierr)
    call PROF_rapend  ('COMM_mpi_post', 3)

    return
  end subroutine vars_3D_mpi_pc

  subroutine wait_3D_mpi(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- wait packets
    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_WAITALL( ginfo(gid)%req_cnt (vid),                           &
                      ginfo(gid)%req_list(1:ginfo(gid)%req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                &
                      ierr                                                )
    call PROF_rapend  ('COMM_mpi_wait', 3)

    call count_bytes( 'COMM_wait', 'COMM_mpi_wait', ginfo(gid)%send_cnt(vid) )
    if ( .not. PRC_TwoD ) then
       KA = ginfo(gid)%KA
       IA = ginfo(gid)%IA
       IS = ginfo(gid)%IS
       IE = ginfo(gid)%IE
       JA = ginfo(gid)%JA
       JS = ginfo(gid)%JS
       JE = ginfo(gid)%JE
       IHALO = ginfo(gid)%IHALO
       call unpackWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, ginfo(gid)%recvpack_WE2P(:,:,vid) )
    end if

    return
  end subroutine wait_3D_mpi

  subroutine wait_3D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    real(RP), pointer, contiguous :: pack(:)

    integer :: ierr
    !---------------------------------------------------------------------------

    KA = ginfo(gid)%KA
    IA = ginfo(gid)%IA
    IS = ginfo(gid)%IS
    IE = ginfo(gid)%IE
    JA = ginfo(gid)%JA
    JS = ginfo(gid)%JS
    JE = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_Win_wait( ginfo(gid)%win_packWE(vid), ierr )
    call PROF_rapend  ('COMM_mpi_wait', 3)
    if ( .not. PRC_TwoD ) then
       call C_F_pointer( ginfo(gid)%recvbuf_WE(vid), pack, (/ginfo(gid)%size2D_WE*KA*2/) )
       call unpackWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, pack )
    end if

    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_Win_wait( ginfo(gid)%win_packNS(vid), ierr )
    call PROF_rapend  ('COMM_mpi_wait', 3)
    call C_F_pointer( ginfo(gid)%recvbuf_NS(vid), pack, (/ginfo(gid)%size2D_NS4*KA*2/) )
    call unpackNS_3D( KA, IA, IS, IE, JA, JS, JE, &
                      JHALO, ginfo(gid)%put_corner(vid), &
                      var, pack )

    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_Win_post( group_packWE, MPI_MODE_NOSTORE, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_post( group_packNS, MPI_MODE_NOSTORE, ginfo(gid)%win_packNS(vid), ierr )
    call PROF_rapend  ('COMM_mpi_wait', 3)

    return
  end subroutine wait_3D_mpi_onesided

  subroutine wait_2D_mpi(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- wait packets
    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_WAITALL( ginfo(gid)%req_cnt(vid),                            &
                      ginfo(gid)%req_list(1:ginfo(gid)%req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                &
                      ierr                                                )
    call PROF_rapend  ('COMM_mpi_wait', 3)

    call count_bytes( 'COMM_wait', 'COMM_mpi_wait', ginfo(gid)%send_cnt(vid) )
    if ( .not. PRC_TwoD ) then
       KA = ginfo(gid)%KA
       IA = ginfo(gid)%IA
       IS = ginfo(gid)%IS
       IE = ginfo(gid)%IE
       JA = ginfo(gid)%JA
       JS = ginfo(gid)%JS
       JE = ginfo(gid)%JE
       IHALO = ginfo(gid)%IHALO
       call unpackWE_2D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, ginfo(gid)%recvpack_WE2P(:,:,vid) )
    end if

    return
  end subroutine wait_2D_mpi

  subroutine wait_2D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    real(RP), pointer, contiguous :: pack(:)

    integer :: ierr
    !---------------------------------------------------------------------------

    KA = ginfo(gid)%KA
    IA = ginfo(gid)%IA
    IS = ginfo(gid)%IS
    IE = ginfo(gid)%IE
    JA = ginfo(gid)%JA
    JS = ginfo(gid)%JS
    JE = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_Win_wait( ginfo(gid)%win_packWE(vid), ierr )
    call PROF_rapend  ('COMM_mpi_wait', 3)
    if ( .not. PRC_TwoD ) then
       call C_F_pointer( ginfo(gid)%recvbuf_WE(vid), pack, (/ginfo(gid)%size2D_WE*KA*2/)  )
       call unpackWE_2D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, pack )
    end if

    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_Win_wait( ginfo(gid)%win_packNS(vid), ierr )
    call PROF_rapend  ('COMM_mpi_wait', 3)
    call C_F_pointer( ginfo(gid)%recvbuf_NS(vid), pack, (/ginfo(gid)%size2D_NS4*KA*2/)  )
    call unpackNS_2D( IA, IS, IE, JA, JS, JE, &
                      JHALO, ginfo(gid)%put_corner(vid), &
                      var, pack )

    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_Win_post( group_packWE, MPI_MODE_NOSTORE, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_post( group_packNS, MPI_MODE_NOSTORE, ginfo(gid)%win_packNS(vid), ierr )
    call PROF_rapend  ('COMM_mpi_wait', 3)

    return
  end subroutine wait_2D_mpi_onesided

  subroutine wait_3D_mpi_pc(var, gid, vid, pid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid        !< sequential number of the persistent variable
    integer,  intent(in)    :: pid        !< ID of the packing buffer of the variable

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO
    integer :: ig
    integer :: off(G_NDIR) ! offsets of the sections of the packed exchange

    integer :: ierr
    !---------------------------------------------------------------------------

    ig = ginfo(gid)%pc_group(vid)
    if ( ig > 0 ) then
       ! the member is waited in a round of its group with all the members posted (the checks of
       ! COMM_wait and of the declaration ensure it; the messages of a batch are not waited otherwise)
       if ( ginfo(gid)%pc_gstat(vid) /= 1 .or. ginfo(gid)%grp(ig)%nposted /= ginfo(gid)%grp(ig)%nmem ) then
          LOG_ERROR("COMM_wait_3D",*) 'the round of the group is not in the state of a wait: ', trim(ginfo(gid)%grp(ig)%name)
          LOG_ERROR_CONT(*) 'grid, ID waited, state of the member, members posted, members: ', gid, vid + COMM_vsize_max, &
                            ginfo(gid)%pc_gstat(vid), ginfo(gid)%grp(ig)%nposted, ginfo(gid)%grp(ig)%nmem
          call PRC_abort
       end if
       if ( ginfo(gid)%grp(ig)%state == COMM_PC_GROUP_BATCH ) then
          ! all the messages of the group are waited at the first wait of the round
          if ( .not. ginfo(gid)%grp(ig)%arrived ) then
             call PROF_rapstart('COMM_mpi_wait', 3)
             call MPI_WAITALL( ginfo(gid)%grp(ig)%nreq,                          &
                               ginfo(gid)%grp(ig)%req(1:ginfo(gid)%grp(ig)%nreq), &
                               MPI_STATUSES_IGNORE,                              &
                               ierr                                              )
             call PROF_rapend  ('COMM_mpi_wait', 3)
             call count_bytes( 'COMM_wait_pers', 'COMM_mpi_wait', ginfo(gid)%grp(ig)%nsend )
             ginfo(gid)%grp(ig)%arrived = .true.
          end if
          call G_offsets( gid, ginfo(gid)%grp(ig)%nmem, ginfo(gid)%pc_slot(vid), .true., off )
          call unpack8_3D( ginfo(gid)%KA, ginfo(gid)%IA, ginfo(gid)%IS, ginfo(gid)%IE, &
                           ginfo(gid)%JA, ginfo(gid)%JS, ginfo(gid)%JE,                &
                           ginfo(gid)%IHALO, ginfo(gid)%JHALO,                         &
                           var, ginfo(gid)%grp(ig)%recvbuf, ginfo(gid)%grp(ig)%nbuf,  &
                           gid, .true., off                                           )
          call COMM_group_waited( gid, vid, ig )
          return
       end if
       call COMM_group_waited( gid, vid, ig )
    end if

    !--- wait packets
    call PROF_rapstart('COMM_mpi_wait', 3)
    call MPI_WAITALL( ginfo(gid)%preq_cnt (vid),                            &
                      ginfo(gid)%preq_list(1:ginfo(gid)%preq_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                  &
                      ierr                                                  )
    call PROF_rapend  ('COMM_mpi_wait', 3)

    call count_bytes( 'COMM_wait_pers', 'COMM_mpi_wait', ginfo(gid)%psend_cnt(vid) )
    KA = ginfo(gid)%KA
    IA = ginfo(gid)%IA
    IS = ginfo(gid)%IS
    IE = ginfo(gid)%IE
    JA = ginfo(gid)%JA
    JS = ginfo(gid)%JS
    JE = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    if ( ginfo(gid)%pc_layout2 .and. ginfo(gid)%pc_ndir(vid) == 8 ) then
       ! all the requests have completed: the N/S halos (whole rows) are in var
       if ( IHALO > 0 ) then
          off(G_N:G_W) = 0
          call L2_offsets( gid, off(G_E), off(G_NW), off(G_NE), off(G_SW), off(G_SE) )
          call unpack8_3D( KA, IA, IS, IE, JA, JS, JE, &
                           IHALO, ginfo(gid)%JHALO, &
                           var, ginfo(gid)%recvpack_L2(:,pid), ginfo(gid)%size_L2, &
                           gid, .false., off )
       end if
    else if ( .not. PRC_TwoD ) then
       call unpackWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, ginfo(gid)%recvpack_WE2P(:,:,pid) )
    end if

    return
  end subroutine wait_3D_mpi_pc

  !-----------------------------------------------------------------------------
  !> Count the size of the sent data
  !! This tree does not have the data size statistics of PROF (PROF_count_bytes, PROF_rapbytes) of
  !! the upstream SCALE, so this routine does nothing.
  subroutine count_bytes( rapname2, rapname3, nsend )
    implicit none

    character(len=*), intent(in) :: rapname2 !< name of the level-2 item in which the data is transferred
    character(len=*), intent(in) :: rapname3 !< name of the level-3 item in which the data is transferred
    integer,          intent(in) :: nsend    !< # of the sent elements
    !---------------------------------------------------------------------------

    return
  end subroutine count_bytes

  subroutine packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                        IHALO, &
                        var, gid, vid)
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO
    real(RP), intent(in) :: var(KA,IA,JA)
    integer,  intent(in) :: gid
    integer,  intent(in) :: vid

    integer :: k, i, j, n

#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:,:)
    ptr => ginfo(gid)%sendpack_P2WE
#endif

    !$acc data copyin(var) if(acc_is_present(var))

    call PROF_rapstart('COMM_pack', 3)

    if ( PRC_HAS_W ) then
       !--- packing packets to West
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       !$acc parallel if(acc_is_present(var)) async
       !$acc loop collapse(2) gang
       do j = JS, JE
       do i = IS, IS+IHALO-1
       !$acc loop independent vector
       do k = 1, KA
          n = (j-JS) * KA * IHALO &
            + (i-IS) * KA         &
            + k
#ifdef _OPENACC
          ptr(n,1,vid) = var(k,i,j)
#else
          ginfo(gid)%sendpack_P2WE(n,1,vid) = var(k,i,j)
#endif
       enddo
       enddo
       enddo
       !$acc end parallel
    end if

    if ( PRC_HAS_E ) then
       !--- packing packets to East
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       !$acc parallel if(acc_is_present(var)) async
       !$acc loop collapse(2) gang
       do j = JS, JE
       do i = IE-IHALO+1, IE
       !$acc loop independent vector
       do k = 1, KA
          n = (j-JS)         * KA * IHALO &
            + (i-IE+IHALO-1) * KA         &
            + k
#ifdef _OPENACC
          ptr(n,2,vid) = var(k,i,j)
#else
          ginfo(gid)%sendpack_P2WE(n,2,vid) = var(k,i,j)
#endif
       enddo
       enddo
       enddo
       !$acc end parallel
    end if

    !$acc wait

    call PROF_rapend('COMM_pack', 3)

    !$acc end data

    return
  end subroutine packWE_3D

  !-----------------------------------------------------------------------------
  !> Pack the halos of the packed 8-direction exchange (the layout 2 and the groups in batch) into
  !! the sections of a buffer, in one kernel for all the directions. off gives the offset of the
  !! section of each direction (G_N etc.: to the neighbor in that direction). The rows in y (N/S)
  !! are packed only with with_NS (the batch; the layout 2 sends them from the variable), the W/E
  !! halos are the rows of L2_range_WE, and the corners are IHALO x JHALO blocks to the diagonal
  !! neighbors. The sections are those of pack8_sections, as in the requests.
  !! The kernel is left in the queue: the caller waits for it before the messages are started (one
  !! wait for all the members of a batch). The variable and the buffer are on the device (a
  !! persistent variable, and the buffers of the layout 2 and of a group), where the data region
  !! below copies nothing
  subroutine pack8_3D( KA, IA, IS, IE, JA, JS, JE, &
                       IHALO, JHALO, &
                       var, buf, nbuf, gid, with_NS, off )
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO, JHALO
    real(RP), intent(in)    :: var(KA,IA,JA)
    integer,  intent(in)    :: nbuf
    real(RP), intent(inout) :: buf(nbuf)
    integer,  intent(in)    :: gid
    logical,  intent(in)    :: with_NS     !< the rows in y are packed
    integer,  intent(in)    :: off(G_NDIR) !< offsets of the sections (to the neighbor)

    integer :: oN, oS, oW, oE, oNW, oNE, oSW, oSE
    integer :: i0, i1, ni, jw0, jw1
    logical :: with_x, has_N, has_S, has_W, has_E, do_N, do_S
    integer :: ngang
    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! local copies used in the compute region
    oN  = off(G_N )
    oS  = off(G_S )
    oW  = off(G_W )
    oE  = off(G_E )
    oNW = off(G_NW)
    oNE = off(G_NE)
    oSW = off(G_SW)
    oSE = off(G_SE)

    call pack8_sections( gid, with_x, has_N, has_S, has_W, has_E, i0, i1 )
    ni = i1 - i0 + 1
    do_N = with_NS .and. has_N
    do_S = with_NS .and. has_S
    call L2_range_WE( gid, jw0, jw1 )

    ! the largest loop
    ngang = 1
    if ( with_NS ) ngang = max( ngang, JHALO * ni )
    if ( with_x  ) ngang = max( ngang, ( jw1 - jw0 + 1 ) * IHALO, JHALO * IHALO )

    !$acc data copyin(var) copy(buf) if(acc_is_present(var))

    call PROF_rapstart('COMM_pack', 3)

    ! one kernel for all the directions: the blocks write disjoint ranges.
    ! The number of gangs is given explicitly: without it, the compiler takes it from one of the loops
    ! (e.g., a corner block, JHALO*IHALO gangs), and the larger loops run on a few gangs
    !$acc parallel num_gangs(ngang) if(acc_is_present(var)) async
    if ( do_N ) then
       !--- to N
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JE-JHALO+1, JE
       do i = i0, i1
       !$acc loop independent vector
       do k = 1, KA
          buf( oN + ( ( j - JE + JHALO - 1 ) * ni + i - i0 ) * KA + k ) = var(k,i,j)
       enddo
       enddo
       enddo
    end if

    if ( do_S ) then
       !--- to S
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JS, JS+JHALO-1
       do i = i0, i1
       !$acc loop independent vector
       do k = 1, KA
          buf( oS + ( ( j - JS ) * ni + i - i0 ) * KA + k ) = var(k,i,j)
       enddo
       enddo
       enddo
    end if

    if ( has_W ) then
       !--- to W
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = jw0, jw1
       do i = IS, IS+IHALO-1
       !$acc loop independent vector
       do k = 1, KA
          buf( oW + ( ( j - jw0 ) * IHALO + i - IS ) * KA + k ) = var(k,i,j)
       enddo
       enddo
       enddo
    end if

    if ( has_E ) then
       !--- to E
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = jw0, jw1
       do i = IE-IHALO+1, IE
       !$acc loop independent vector
       do k = 1, KA
          buf( oE + ( ( j - jw0 ) * IHALO + i - IE + IHALO - 1 ) * KA + k ) = var(k,i,j)
       enddo
       enddo
       enddo
    end if

    if ( has_N .and. has_W ) then
       !--- to NW
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JE-JHALO+1, JE
       do i = IS, IS+IHALO-1
       !$acc loop independent vector
       do k = 1, KA
          buf( oNW + ( ( j - JE + JHALO - 1 ) * IHALO + i - IS ) * KA + k ) = var(k,i,j)
       enddo
       enddo
       enddo
    end if

    if ( has_N .and. has_E ) then
       !--- to NE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JE-JHALO+1, JE
       do i = IE-IHALO+1, IE
       !$acc loop independent vector
       do k = 1, KA
          buf( oNE + ( ( j - JE + JHALO - 1 ) * IHALO + i - IE + IHALO - 1 ) * KA + k ) = var(k,i,j)
       enddo
       enddo
       enddo
    end if

    if ( has_S .and. has_W ) then
       !--- to SW
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JS, JS+JHALO-1
       do i = IS, IS+IHALO-1
       !$acc loop independent vector
       do k = 1, KA
          buf( oSW + ( ( j - JS ) * IHALO + i - IS ) * KA + k ) = var(k,i,j)
       enddo
       enddo
       enddo
    end if

    if ( has_S .and. has_E ) then
       !--- to SE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JS, JS+JHALO-1
       do i = IE-IHALO+1, IE
       !$acc loop independent vector
       do k = 1, KA
          buf( oSE + ( ( j - JS ) * IHALO + i - IE + IHALO - 1 ) * KA + k ) = var(k,i,j)
       enddo
       enddo
       enddo
    end if
    !$acc end parallel

    call PROF_rapend('COMM_pack', 3)

    !$acc end data

    return
  end subroutine pack8_3D

  !-----------------------------------------------------------------------------
  !> Unpack the halos of the packed 8-direction exchange (see pack8_3D) from the sections of a
  !! buffer, after all its messages have arrived, in one kernel for all the directions. off gives
  !! the offset of the section from the neighbor in each direction (G_N etc.). The ranges do not
  !! overlap: the N/S halos (only with with_NS) exclude the halo columns on the sides with a
  !! neighbor in x, the W/E halos the halo rows on the sides with a neighbor in y (L2_range_WE),
  !! and the corners come from the diagonal neighbors, which exist only where the neighbors in x
  !! and in y exist
  subroutine unpack8_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, JHALO, &
                         var, buf, nbuf, gid, with_NS, off )
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO, JHALO
    real(RP), intent(inout) :: var(KA,IA,JA)
    integer,  intent(in)    :: nbuf
    real(RP), intent(in)    :: buf(nbuf)
    integer,  intent(in)    :: gid
    logical,  intent(in)    :: with_NS     !< the rows in y are unpacked
    integer,  intent(in)    :: off(G_NDIR) !< offsets of the sections (from the neighbor)

    integer :: oN, oS, oW, oE, oNW, oNE, oSW, oSE
    integer :: i0, i1, ni, iu0, iu1, jw0, jw1
    logical :: with_x, has_N, has_S, has_W, has_E, do_N, do_S
    integer :: ngang
    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! local copies used in the compute region
    oN  = off(G_N )
    oS  = off(G_S )
    oW  = off(G_W )
    oE  = off(G_E )
    oNW = off(G_NW)
    oNE = off(G_NE)
    oSW = off(G_SW)
    oSE = off(G_SE)

    ! columns of the rows in y in the message (i0:i1) and those unpacked (iu0:iu1)
    call pack8_sections( gid, with_x, has_N, has_S, has_W, has_E, i0, i1 )
    ni  = i1 - i0 + 1
    iu0 = i0
    iu1 = i1
    if ( has_W ) iu0 = IS
    if ( has_E ) iu1 = IE
    do_N = with_NS .and. has_N
    do_S = with_NS .and. has_S
    call L2_range_WE( gid, jw0, jw1 )

    ! the largest loop
    ngang = 1
    if ( with_NS ) ngang = max( ngang, JHALO * ( iu1 - iu0 + 1 ) )
    if ( with_x  ) ngang = max( ngang, ( jw1 - jw0 + 1 ) * IHALO, JHALO * IHALO )

    !$acc data copy(var) copyin(buf) if(acc_is_present(var))

    call PROF_rapstart('COMM_unpack', 3)

    ! one kernel for all the directions: the blocks write disjoint ranges.
    ! The number of gangs is given explicitly (the largest loop; see pack8_3D)
    !$acc parallel num_gangs(ngang) if(acc_is_present(var)) async
    if ( do_S ) then
       !--- from S
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = 1, JS-1
       do i = iu0, iu1
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf( oS + ( ( j - 1 ) * ni + i - i0 ) * KA + k )
       enddo
       enddo
       enddo
    end if

    if ( do_N ) then
       !--- from N
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JE+1, JA
       do i = iu0, iu1
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf( oN + ( ( j - JE - 1 ) * ni + i - i0 ) * KA + k )
       enddo
       enddo
       enddo
    end if

    if ( has_W ) then
       !--- from W
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = jw0, jw1
       do i = 1, IS-1
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf( oW + ( ( j - jw0 ) * IHALO + i - 1 ) * KA + k )
       enddo
       enddo
       enddo
    end if

    if ( has_E ) then
       !--- from E
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = jw0, jw1
       do i = IE+1, IA
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf( oE + ( ( j - jw0 ) * IHALO + i - IE - 1 ) * KA + k )
       enddo
       enddo
       enddo
    end if

    if ( has_N .and. has_W ) then
       !--- from NW
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JE+1, JA
       do i = 1, IS-1
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf( oNW + ( ( j - JE - 1 ) * IHALO + i - 1 ) * KA + k )
       enddo
       enddo
       enddo
    end if

    if ( has_N .and. has_E ) then
       !--- from NE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = JE+1, JA
       do i = IE+1, IA
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf( oNE + ( ( j - JE - 1 ) * IHALO + i - IE - 1 ) * KA + k )
       enddo
       enddo
       enddo
    end if

    if ( has_S .and. has_W ) then
       !--- from SW
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = 1, JS-1
       do i = 1, IS-1
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf( oSW + ( ( j - 1 ) * IHALO + i - 1 ) * KA + k )
       enddo
       enddo
       enddo
    end if

    if ( has_S .and. has_E ) then
       !--- from SE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc loop collapse(2) gang
       do j = 1, JS-1
       do i = IE+1, IA
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf( oSE + ( ( j - 1 ) * IHALO + i - IE - 1 ) * KA + k )
       enddo
       enddo
       enddo
    end if
    !$acc end parallel

    !$acc wait

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpack8_3D

  subroutine packWE_2D( IA, IS, IE, JA, JS, JE, &
                        IHALO, &
                        var, gid, vid)
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO
    real(RP), intent(in) :: var(IA,JA)
    integer,  intent(in) :: vid
    integer,  intent(in) :: gid

    integer :: i, j, n

#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:,:)
    ptr => ginfo(gid)%sendpack_P2WE
#endif
    !$acc data copyin(var) if(acc_is_present(var))

    call PROF_rapstart('COMM_pack', 3)

    if ( PRC_HAS_W ) then
       !--- To 4-Direction HALO communicate
       !--- packing packets to West
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_
       !$acc kernels if(acc_is_present(var)) async
       !$acc loop independent
       do j = JS, JE
       !$acc loop independent
       do i = IS, IS+IHALO-1
          n = (j-JS) * IHALO &
            + (i-IS) + 1
#ifdef _OPENACC
          ptr(n,1,vid) = var(i,j)
#else
          ginfo(gid)%sendpack_P2WE(n,1,vid) = var(i,j)
#endif
       enddo
       enddo
       !$acc end kernels
    end if

    if ( PRC_HAS_E ) then
       !--- packing packets to East
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_
       !$acc kernels if(acc_is_present(var)) async
       !$acc loop independent
       do j = JS, JE
       !$acc loop independent
       do i = IE-IHALO+1, IE
          n = (j-JS)         * IHALO &
            + (i-IE+IHALO-1) + 1
#ifdef _OPENACC
          ptr(n,2,vid) = var(i,j)
#else
          ginfo(gid)%sendpack_P2WE(n,2,vid) = var(i,j)
#endif
       enddo
       enddo
       !$acc end kernels
    end if

    !$acc wait

    call PROF_rapend('COMM_pack', 3)

    !$acc end data

    return
  end subroutine packWE_2D

  subroutine unpackWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, buf )
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO
    real(RP), intent(inout) :: var(KA,IA,JA)
    real(RP), intent(in)    :: buf(KA,IHALO,JS:JE,2)

    integer :: i, j, k
    !---------------------------------------------------------------------------

    !$acc data copy(var) copyin(buf) if(acc_is_present(var))

    call PROF_rapstart('COMM_unpack', 3)

    if ( PRC_HAS_E ) then
       !--- unpacking packets from East
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc parallel if(acc_is_present(var)) async
       !$acc loop collapse(2) gang
       do j = JS, JE
       do i = IE+1, IA
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf(k,i-IE,j,2)
       enddo
       enddo
       enddo
       !$acc end parallel
    end if

    if ( PRC_HAS_W ) then
       !--- unpacking packets from West
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc parallel if(acc_is_present(var)) async
       !$acc loop collapse(2) gang
       do j = JS, JE
       do i = 1, IS-1
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf(k,i,j,1)
       enddo
       enddo
       enddo
       !$acc end parallel
    end if

    !$acc wait

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpackWE_3D

  subroutine unpackWE_2D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, buf )
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO
    real(RP), intent(inout) :: var(IA,JA)
    real(RP), intent(in)    :: buf(IHALO,JS:JE,KA,2)

    integer :: i, j
    !---------------------------------------------------------------------------

    !$acc data copy(var) copyin(buf) if(acc_is_present(var))

    call PROF_rapstart('COMM_unpack', 3)

    if( PRC_HAS_E ) then
        !--- unpacking packets from East
        !$omp parallel do private(i,j) OMP_SCHEDULE_
       !$acc kernels if(acc_is_present(var)) async
        do j = JS, JE
        do i = IE+1, IE+IHALO
           var(i,j) = buf(i-IE,j,1,2)
        enddo
        enddo
       !$acc end kernels
     end if

     if( PRC_HAS_W ) then
        !--- unpacking packets from West
        !$omp parallel do private(i,j) OMP_SCHEDULE_
       !$acc kernels if(acc_is_present(var)) async
        do j = JS, JE
        do i = IS-IHALO, IS-1
           var(i,j) = buf(i,j,1,1)
        enddo
        enddo
       !$acc end kernels
    end if

    !$acc wait

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpackWE_2D

  subroutine unpackNS_3D( KA, IA, IS, IE, JA, JS, JE, &
                          JHALO, put_corner, &
                          var, buf )
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: JHALO
    logical, intent(in) :: put_corner !< the corners are in the buffer (not in the 4-direction exchange)
    real(RP), intent(inout) :: var(KA,IA,JA)
    real(RP), intent(in)    :: buf(KA,IA,JHALO,2)

    integer :: i, j, k
    !---------------------------------------------------------------------------

    !$acc data copy(var) copyin(buf)

    call PROF_rapstart('COMM_unpack', 3)

    if ( PRC_HAS_S ) then
       !--- unpacking packets from S, SW, and SE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc kernels async
       do j = 1, JS-1
       do i = 1, IA
       do k = 1, KA
          var(k,i,j) = buf(k,i,j,1)
       enddo
       enddo
       enddo
       !$acc end kernels
    else if ( put_corner ) then
       if ( PRC_HAS_W ) then
          !--- unpacking packets from SW
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !$acc kernels async
          do j = 1, JS-1
          do i = 1, IS-1
          do k = 1, KA
             var(k,i,j) = buf(k,i,j,1)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
       if ( PRC_HAS_E ) then
          !--- unpacking packets from SE
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !$acc kernels async
          do j = 1, JS-1
          do i = IE+1, IA
          do k = 1, KA
             var(k,i,j) = buf(k,i,j,1)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
    end if

    if ( PRC_HAS_N ) then
       !--- unpacking packets from N, NW, and NE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc kernels async
       do j = JE+1, JA
       do i = 1, IA
       do k = 1, KA
          var(k,i,j) = buf(k,i,j-JE,2)
       enddo
       enddo
       enddo
       !$acc end kernels
    else if ( put_corner ) then
       if ( PRC_HAS_W ) then
          !--- unpacking packets from NW
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !$acc kernels async
          do j = JE+1, JA
          do i = 1, IS-1
          do k = 1, KA
             var(k,i,j) = buf(k,i,j-JE,2)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
       if ( PRC_HAS_E ) then
          !--- unpacking packets from NE
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !$acc kernels async
          do j = JE+1, JA
          do i = IE+1, IA
          do k = 1, KA
             var(k,i,j) = buf(k,i,j-JE,2)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
    end if

    !$acc wait

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpackNS_3D

  subroutine unpackNS_2D( IA, IS, IE, JA, JS, JE, &
                          JHALO, put_corner, &
                          var, buf )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: JHALO
    logical, intent(in) :: put_corner !< the corners are in the buffer (not in the 4-direction exchange)
    real(RP), intent(inout) :: var(IA,JA)
    real(RP), intent(in)    :: buf(IA,JHALO,2)

    integer :: i, j
    !---------------------------------------------------------------------------

    !$acc data copy(var) copyin(buf)

    call PROF_rapstart('COMM_unpack', 3)

    if ( PRC_HAS_S ) then
       !--- unpacking packets from S, SW, and SE
       !$omp parallel do private(i,j) OMP_SCHEDULE_
       !$acc kernels async
       do j = 1, JS-1
       do i = 1, IA
          var(i,j) = buf(i,j,1)
       enddo
       enddo
       !$acc end kernels
    else if ( put_corner ) then
       if ( PRC_HAS_W ) then
          !--- unpacking packets from SW
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          !$acc kernels async
          do j = 1, JS-1
          do i = 1, IS-1
             var(i,j) = buf(i,j,1)
          enddo
          enddo
          !$acc end kernels
       end if
       if ( PRC_HAS_E ) then
          !--- unpacking packets from SE
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          !$acc kernels async
          do j = 1, JS-1
          do i = IE+1, IA
             var(i,j) = buf(i,j,1)
          enddo
          enddo
          !$acc end kernels
       end if
    end if

    if ( PRC_HAS_N ) then
       !--- unpacking packets from N, NW, and NE
       !$omp parallel do private(i,j) OMP_SCHEDULE_
       !$acc kernels async
       do j = JE+1, JA
       do i = 1, IA
          var(i,j) = buf(i,j-JE,2)
       enddo
       enddo
       !$acc end kernels
    else if ( put_corner ) then
       if ( PRC_HAS_W ) then
          !--- unpacking packets from NW
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          !$acc kernels async
          do j = JE+1, JA
          do i = 1, IS-1
             var(i,j) = buf(i,j-JE,2)
          enddo
          enddo
          !$acc end kernels
       end if
       if ( PRC_HAS_E ) then
          !--- unpacking packets from NE
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          !$acc kernels async
          do j = JE+1, JA
          do i = IE+1, IA
             var(i,j) = buf(i,j-JE,2)
          enddo
          enddo
          !$acc end kernels
       end if
    end if

    !$acc wait

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpackNS_2D

  subroutine copy_boundary_3D(var, gid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid

    integer :: KA
    integer :: IS, IE, IHALO
    integer :: JS, JE, JHALO

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data copy(var)

    KA    = ginfo(gid)%KA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$omp parallel

    !--- copy inner data to HALO(North)
    if ( .NOT. PRC_HAS_N ) then
       !$acc kernels async
       do j = JE+1, JE+JHALO
       !$omp do
       do i = IS, IE
       do k = 1, KA
          var(k,i,j) = var(k,i,JE)
       enddo
       enddo
       !$omp end do nowait
       enddo
       !$acc end kernels
    endif

    !--- copy inner data to HALO(South)
    if ( .NOT. PRC_HAS_S ) then
       !$acc kernels async
       !$acc loop independent
       do j = JS-JHALO, JS-1
       !$omp do
       do i = IS, IE
       do k = 1, KA
          var(k,i,j) = var(k,i,JS)
       enddo
       enddo
       !$omp end do nowait
       enddo
       !$acc end kernels
    endif

    if ( .not. PRC_TwoD ) then

       !--- copy inner data to HALO(East)
       if ( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          !$omp do
          do j = JS, JE
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
          !$omp end do nowait
          !$acc end kernels
       end if

       !--- copy inner data to HALO(West)
       if ( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          !$omp do
          do j = JS, JE
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(:,i,j) = var(:,IS,j)
          enddo
          enddo
          !$omp end do nowait
          !$acc end kernels
       end if

       !--- copy inner data to HALO(NorthWest)
       if ( .NOT. PRC_HAS_N .AND. &
            .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          !$acc loop independent
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,JE)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_N ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,i,JE)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          !$acc loop independent
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(SouthWest)
       if ( .NOT. PRC_HAS_S .AND. &
            .NOT. PRC_HAS_W ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          !$acc loop independent
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,JS)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_S ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,i,JS)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          !$acc loop independent
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(NorthEast)
       if ( .NOT. PRC_HAS_N .AND. &
            .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,JE)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_N ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,i,JE)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(SouthEast)
       if ( .NOT. PRC_HAS_S .AND. &
            .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,JS)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_S ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,i,JS)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

    end if

    !$omp end parallel

    !$acc wait

    !$acc end data

    return
  end subroutine copy_boundary_3D

  subroutine copy_boundary_2D(var, gid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid

    integer :: IS, IE, IHALO
    integer :: JS, JE, JHALO

    integer :: i, j
    !---------------------------------------------------------------------------

    !$acc data copy(var)

    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$omp parallel

    !--- copy inner data to HALO(North)
    if( .NOT. PRC_HAS_N ) then
       !$acc kernels async
       do j = JE+1, JE+JHALO
       !$omp do
       do i = IS, IE
          var(i,j) = var(i,JE)
       enddo
       !$omp end do nowait
       enddo
       !$acc end kernels
    endif

    !--- copy inner data to HALO(South)
    if( .NOT. PRC_HAS_S ) then
       !$acc kernels async
       !$acc loop independent
       do j = JS-JHALO, JS-1
       !$omp do
       do i = IS, IE
          var(i,j) = var(i,JS)
       enddo
       !$omp end do nowait
       enddo
       !$acc end kernels
    endif

    if ( .not. PRC_TwoD ) then

       if( .NOT. PRC_HAS_E ) then
          !$omp do
          !$acc kernels async
          do j = JS, JE
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       endif

       if( .NOT. PRC_HAS_W ) then
          !$omp do
          !$acc kernels async
          do j = JS, JE
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       endif

       !--- copy inner data to HALO(NorthWest)
       if( .NOT. PRC_HAS_N .AND. .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,JE)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_N ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
             var(i,j) = var(i,JE)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(SouthWest)
       if( .NOT. PRC_HAS_S .AND. .NOT. PRC_HAS_W ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,JS)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_S ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          do i = IS-IHALO, IS-1
             var(i,j) = var(i,JS)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(NorthEast)
       if( .NOT. PRC_HAS_N .AND. .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,JE)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_N ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(i,JE)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(SouthEast)
       if( .NOT. PRC_HAS_S .AND. .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,JS)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_S ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(i,JS)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
          !$acc end kernels
       endif

    end if

    !$omp end parallel

    !$acc wait

    !$acc end data

    return
  end subroutine copy_boundary_2D

end module scale_comm_cartesC
