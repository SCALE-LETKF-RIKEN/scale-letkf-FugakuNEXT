!-------------------------------------------------------------------------------
!> Standalone test of the order of the halo exchanges in COMM_cartesC (issue #3594)
!!
!! The exchanges that are in progress at the same time may be posted and waited in an order that
!! differs among the processes: every exchange is started when it is posted, so the result does not
!! depend on the order. This test fixes that freedom of the current code, so that a change of the
!! communication (e.g., grouped exchanges) cannot remove it unnoticed. The cases are:
!!   post order of three persistent variables of a grid (8-direction, and 4-direction)
!!   wait order of three persistent variables of a grid (8-direction, and 4-direction)
!!   wait order (and post order) of persistent variables on two grids
!!   a variable posted and waited while two others are in progress (the form of the numerical
!!   filter), with the inner variable on the same grid or on the other grid
!!   a variable posted and waited between the posts of two variables in some processes and after
!!   them in the others, on the same or the other grid (this form hangs if a group of exchanges
!!   is started only when its last variable is posted)
!!   persistent (4- and 8-direction) and non-persistent exchanges in progress at the same time
!!   persistent and non-persistent (3-D and 2-D) exchanges with the same packing buffer, one after
!!   another
!!   post and wait order of all the 8-direction persistent variables of the two grids
!! The order depends on the position of the process. An order of three or more variables is
!! reversed where x is odd and rotated by y, so that all the neighbors (in x, in y and diagonal)
!! use different orders in all the layouts of make run-comm-order (a rotation needs at least as
!! many variables as processes in y, or the periodic neighbors in y may have the same order). An
!! order of two variables has only two forms: it is reversed where x + y is odd, so the diagonal
!! neighbors have the same order; so do the orders of case 4 below.
!!
!! Each case is run twice (the persistent requests are started again) with other values. Before an
!! exchange, the inner cells have values of the variable and the global position, and the halos
!! have sentinel values. The result is compared on every cell, halos and corners included, with an
!! exchange of one variable at a time in the same order in all the processes (the reference), and
!! the halos whose source cell is in the domain are compared with the value of that cell (the
!! reference too). Between the reference and the case, all the variables are exchanged with other
!! values, so that no value of the reference is left in the buffers (a buffer not packed again
!! would be sent with them unnoticed). All the cases are run without FILL_BND, so that the halos at
!! the non-periodic boundaries keep the values left by the exchange, and again with FILL_BND.
!!
!! The result of each variable after a case (in each round, without and with FILL_BND) does not
!! depend on the case: the result with the first communication configuration is kept, and the
!! results of all the cases with the other configurations are compared with it on every cell
!! (e.g., a new option of the communication against the current code).
!! The values of the inner cells are integers below 2**24 (exact in single precision), unique while
!! the global numbers of cells in x and in y and KA (of grid 2) are below 100; keep the layouts so
!! (run_comm_large.conf, whose messages exceed the eager limit of Open MPI, is 96 x 96 x 64). The
!! sentinel values exceed 2**24, so they are rounded in single precision and tell the halos apart
!! less well (the test does not fail by that).
!! With PROF_MPI_BARRIER = .true. (PARAM_PROF), set SKIP_ASYMMETRIC = .true. (PARAM_TEST_COMM_ORDER,
!! in the first file): the posts and waits differ among the processes in case 4, which hangs with
!! the barriers of the current code, so it is skipped (run_comm_barrier.conf).
!!
!! The number of the persistent requests of each variable (COMM_get_pc_info) is checked against
!! the layout of the exchange that PARAM_COMM_CARTESC selects (COMM_USE_PC_LAYOUT2), so that a
!! configuration cannot fall back to the current layout unnoticed.
!!
!! Groups (COMM_vars_group_init): after the cases above (phase A), the persistent variables are
!! released and registered again, and groups are declared (phase B): a group of five variables of
!! grid 1 (four in 2-D, without the one that stands for a variable not exchanged in 2-D), a group
!! of two posted after it (the form of PROG), a group of one, and a group of three of grid 2; a
!! group without members and one with a non-persistent vid are not checked. Groups of one of q4
!! registered again extend the table of the groups after the first two, which are exchanged
!! afterwards. The variables registered after the groups extend the table of the persistent
!! variables (with a small capacity). The IDs must be those
!! of the same registrations without groups, and the state of each group (COMM_get_pc_info) that
!! of the configuration: checked, or in batch with COMM_USE_BATCH (the number of the requests of
!! the group from the neighbors). The cases of phase B post all the members of the groups in an
!! order that depends on the process, before they are waited in another such order: with the
!! posts of the groups on the two grids interleaved, with a group posted after the other, with
!! persistent and non-persistent exchanges of other variables, with an exchange of another
!! variable between the waits, with exchanges of other variables between the posts (at the same
!! place in all the processes), and in two rounds. The reference exchange of one variable at a
!! time is not allowed for a member of a group: the result is compared with the result kept from
!! phase A of the first configuration and with the source cells, and the values are left in the
!! buffers by a round of all the variables posted and waited in the same order in all the
!! processes. The setups again below declare a group, and declare another one with the variables
!! registered again without the release.
!!
!! The first argument is the configuration file of the grid and the process layout (run_comm_*.conf).
!! The other arguments are configuration files of the communication (PARAM_COMM_CARTESC): all the
!! cases are run for each of them in turn (COMM_finalize and COMM_setup in between), and each must
!! use persistent communication. Without them, the namelist in the first file is used. With two or
!! more of them, the results with the first one are compared with those with the others, and all
!! of them are then set up again twice in turn with a single exchange and a round of a group (a
!! buffer or a communicator not freed at COMM_finalize aborts the next COMM_regist or the next
!! declaration that allocates it).
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
program test_comm_order
  use scale
  use scale_prc, only: &
     PRC_myrank
  use scale_prc_cartesC, only: &
     PRC_CARTESC_setup, &
     PRC_CARTESC_finalize, &
     PRC_2Drank, &
     PRC_NUM_X, &
     PRC_NUM_Y, &
     PRC_PERIODIC_X, &
     PRC_PERIODIC_Y, &
     PRC_TwoD, &
     PRC_HAS_W, &
     PRC_HAS_N, &
     PRC_HAS_E, &
     PRC_HAS_S
  use scale_atmos_grid_cartesC_index
  use scale_comm_cartesC, only: &
     COMM_setup, &
     COMM_regist, &
     COMM_finalize, &
     COMM_vars_init, &
     COMM_vars8_init, &
     COMM_vars_finalize, &
     COMM_vars, &
     COMM_vars8, &
     COMM_wait, &
     COMM_get_vsize_max, &
     COMM_get_pc_info, &
     COMM_vars_group_init, &
     COMM_PC_GROUP_NONE, &
     COMM_PC_GROUP_CHECK, &
     COMM_PC_GROUP_BATCH
  use dc_test, only: &
     AssertEqual
  implicit none

  type grid_t
     integer :: gid
     integer :: KA, IA, JA
     integer :: IHALO, JHALO
     integer :: IS, IE, JS, JE
  end type grid_t

  type field_t
     character(len=2)      :: name
     integer               :: grid         !< index of the grid (1 or 2)
     integer               :: ndir         !< 8 or 4
     logical               :: persistent
     integer               :: pid          !< vid given at the registration (packing buffer)
     integer               :: vid          !< ID of the exchange
     real(RP), pointer     :: v(:,:,:) => null()
     real(RP), allocatable :: ref(:,:,:)
     real(RP), allocatable :: expect(:,:,:)
  end type field_t

  type saved_t
     real(RP), allocatable :: v(:,:,:)
  end type saved_t

  ! fields; the packing buffers (pid) of the fields of a grid are all different, so that any of
  ! them may be in progress at the same time
  integer, parameter :: V_A1 = 1, V_B1 = 2,  V_C1 = 3  ! 8-direction persistent, grid 1
  integer, parameter :: V_D1 = 4, V_E1 = 5,  V_F1 = 6  ! 4-direction persistent, grid 1
  integer, parameter :: V_Z1 = 7                       ! 8-direction persistent, grid 1
  integer, parameter :: V_N1 = 8                       ! 8-direction non-persistent, grid 1
  integer, parameter :: V_A2 = 9, V_B2 = 10, V_Z2 = 11 ! 8-direction persistent, grid 2
  integer, parameter :: V_M1 = 12                      ! 8-direction non-persistent, grid 1, packing buffer of a1
  integer, parameter :: V_Q1 = 13, V_Q2 = 14           ! 8-direction persistent, grid 1
  integer, parameter :: V_Q3 = 15, V_Q4 = 16           ! 8-direction persistent, grid 1
  integer, parameter :: NF = 16
  integer, parameter :: OP_2D = NF + 1                 ! (in a case) 2-D exchanges with the packing buffer of a1
  integer, parameter :: NREENTRY = 2                   ! rounds of the setups again with all the configurations
  integer, parameter :: NROUND = 2
  integer, parameter :: NREQ_GROUP_MAX = 16            ! requests of a group in batch at most

  real(RP), allocatable, target :: a1(:,:,:), b1(:,:,:), c1(:,:,:), d1(:,:,:), e1(:,:,:), f1(:,:,:)
  real(RP), allocatable, target :: z1(:,:,:), n1(:,:,:), a2(:,:,:), b2(:,:,:), z2(:,:,:), m1(:,:,:)
  real(RP), allocatable, target :: q1(:,:,:), q2(:,:,:), q3(:,:,:), q4(:,:,:)
  real(RP), allocatable         :: s1(:,:) ! 2-D, grid 1

  type(grid_t)  :: grd(2)
  type(field_t) :: fld(NF)

  ! results with the first communication configuration: (field, round, FILL_BND)
  type(saved_t) :: saved(NF,NROUND,0:1)
  integer       :: ncompare     !< number of the results compared with saved
  logical       :: fill_bnd     !< FILL_BND of COMM_wait

  ! groups of phase B (fields): ga, gp and g1 on grid 1, gb on grid 2; lone: a field of ga in 3-D,
  ! which is not in it in 2-D (a variable not exchanged in 2-D), and is then exchanged alone
  integer, allocatable :: ga(:), lone(:)
  integer, parameter   :: gp(2) = (/ V_Z1, V_Q3 /)
  integer, parameter   :: g1(1) = (/ V_Q4 /)
  integer, parameter   :: gb(3) = (/ V_A2, V_B2, V_Z2 /)
  logical :: batch              !< COMM_USE_BATCH of the configuration
  logical :: batch_nonpc        !< COMM_USE_NONPC_BATCH of the configuration
  logical :: unbound = .false.  !< phase C: the members of the groups are IDs not bound to the arrays
  logical :: skip_case4         !< skip case 4 (SKIP_ASYMMETRIC of PARAM_TEST_COMM_ORDER)

  character(len=H_LONG) :: comm_conf
  integer :: nconf, iconf, irep
  integer :: f
  integer :: ix, iy
  logical :: odd_x, odd_xy
  logical :: has_neighbor

  call SCALE_init( "Unit test of the order of the halo exchanges" )
  call PRC_CARTESC_setup
  call ATMOS_GRID_CARTESC_INDEX_setup
  call read_test_conf

  ! order of the process (see perm)
  ix = PRC_2Drank(PRC_myrank,1)
  iy = PRC_2Drank(PRC_myrank,2)
  odd_x  = mod( ix, 2 ) == 1
  odd_xy = mod( ix + iy, 2 ) == 1

  ! the halos have source cells in the domain (the neighbor may be the process itself)
  has_neighbor = PRC_HAS_N .or. PRC_HAS_S .or. ( IHALO > 0 .and. ( PRC_HAS_W .or. PRC_HAS_E ) )

  ! grid 1: the grid of the model; grid 2: other halos and vertical size, the same inner cells.
  ! In 2-D (PRC_TwoD, IHALO = 0), grid 2 has IHALO = 1, which is not exchanged (the only value that
  ! COMM_regist accepts with IMAX = 1: IA >= IHALO * 3)
  call set_grid( grd(1), KA, IHALO, JHALO )
  call set_grid( grd(2), KA + 2, IHALO + 1, JHALO + 1 )

  ! the values of the inner cells are unique (see inner_value)
  if ( IMAX * PRC_NUM_X >= 100 .or. JMAX * PRC_NUM_Y >= 100 .or. grd(2)%KA >= 100 ) then
     call AssertEqual( "the global numbers of cells and KA are below 100", .true., .false. )
  end if

  call alloc( a1, 1 )
  call alloc( b1, 1 )
  call alloc( c1, 1 )
  call alloc( d1, 1 )
  call alloc( e1, 1 )
  call alloc( f1, 1 )
  call alloc( z1, 1 )
  call alloc( n1, 1 )
  call alloc( a2, 2 )
  call alloc( b2, 2 )
  call alloc( z2, 2 )
  call alloc( m1, 1 )
  call alloc( q1, 1 )
  call alloc( q2, 1 )
  call alloc( q3, 1 )
  call alloc( q4, 1 )
  allocate( s1(grd(1)%IA,grd(1)%JA) )
  s1(:,:) = 0.0_RP

  call set_field( V_A1, "a1", a1, 1, 8, .true.,  1 )
  call set_field( V_B1, "b1", b1, 1, 8, .true.,  2 )
  call set_field( V_C1, "c1", c1, 1, 8, .true.,  3 )
  call set_field( V_D1, "d1", d1, 1, 4, .true.,  4 )
  call set_field( V_E1, "e1", e1, 1, 4, .true.,  5 )
  call set_field( V_F1, "f1", f1, 1, 4, .true.,  6 )
  call set_field( V_Z1, "z1", z1, 1, 8, .true.,  7 )
  call set_field( V_N1, "n1", n1, 1, 8, .false., 8 )
  call set_field( V_A2, "a2", a2, 2, 8, .true.,  1 )
  call set_field( V_B2, "b2", b2, 2, 8, .true.,  2 )
  call set_field( V_Z2, "z2", z2, 2, 8, .true.,  3 )
  call set_field( V_M1, "m1", m1, 1, 8, .false., 1 )
  call set_field( V_Q1, "q1", q1, 1, 8, .true.,  9 )
  call set_field( V_Q2, "q2", q2, 1, 8, .true., 10 )
  call set_field( V_Q3, "q3", q3, 1, 8, .true., 11 )
  call set_field( V_Q4, "q4", q4, 1, 8, .true., 12 )

  if ( PRC_TwoD ) then
     allocate( ga(4), lone(1) )
     ga(:) = (/ V_A1, V_B1, V_Q1, V_Q2 /)
     lone(:) = (/ V_C1 /)
  else
     allocate( ga(5), lone(0) )
     ga(:) = (/ V_A1, V_B1, V_C1, V_Q1, V_Q2 /)
  end if

  ! the variables must be on the device to be registered in the OpenACC build, and they must
  ! outlive the registration
  !$acc enter data copyin(a1, b1, c1, d1, e1, f1, z1, n1, a2, b2, z2, m1, q1, q2, q3, q4, s1)

  nconf = command_argument_count() - 1
  ncompare = 0
  do iconf = 1, max( nconf, 1 )
     call setup_comm( (/ ( .true., f = 1, NF ) /) )

     ! phase A
     fill_bnd = .false.
     call run_cases
     fill_bnd = .true.
     call run_cases

     ! phase B: groups
     call COMM_vars_finalize
     call setup_groups
     fill_bnd = .false.
     call run_group_cases
     fill_bnd = .true.
     call run_group_cases

     ! phase C: the same groups, with IDs that are not bound to the arrays (issue #3624). The IDs
     ! of the first two members of the first group are swapped after the declaration (in all the
     ! processes): an ID is given with an array other than that of its registration
     call COMM_vars_finalize
     unbound = .true.
     call setup_groups
     f = fld(ga(1))%vid
     fld(ga(1))%vid = fld(ga(2))%vid
     fld(ga(2))%vid = f
     fill_bnd = .false.
     call run_group_cases
     fill_bnd = .true.
     call run_group_cases
     unbound = .false.

     call COMM_vars_finalize
     call COMM_finalize
  end do

  ! the comparison of the configurations must not be empty
  if ( nconf > 1 ) then
     LOG_INFO("test_comm_order",*) 'results compared with the first configuration: ', ncompare
     if ( ncompare == 0 ) then
        call AssertEqual( "number of the results compared with the first configuration", 1, 0 )
     end if
  end if

  ! the configurations are set up again in turn (e.g., the layout 2 and the current one), with an
  ! exchange of a variable of each kind and a round of a group: the buffers and the communicators
  ! must be freed at each COMM_finalize
  if ( nconf > 1 ) then
     fill_bnd = .true.
     do irep = 1, NREENTRY
     do iconf = 1, nconf
        call setup_comm( (/ ( f == V_A1 .or. f == V_B1 .or. f == V_D1 .or. f == V_A2, f = 1, NF ) /) )
        call run_single( (/ V_A1, V_D1, V_A2 /) )
        call run_group_single
        call COMM_vars_finalize
        call COMM_finalize
     end do
     end do
     if ( PRC_myrank == 0 ) write(*,*) 'test_comm_order: setups again with all the configurations OK'
  end if

  !$acc exit data delete(a1, b1, c1, d1, e1, f1, z1, n1, a2, b2, z2, m1, q1, q2, q3, q4, s1)

  call PRC_CARTESC_finalize
  call SCALE_finalize

contains

  !-----------------------------------------------------------------------------
  subroutine run_cases
    integer :: p3(3), w3(3), p4(4), w4(4), p2(2), w2(2)

    ! 1. post order depends on the process, common wait order
    p3 = perm( (/ V_A1, V_B1, V_C1 /) )
    call run_case( "post order, 8-direction", (/ p3, -V_A1, -V_B1, -V_C1 /) )
    p3 = perm( (/ V_D1, V_E1, V_F1 /) )
    call run_case( "post order, 4-direction", (/ p3, -V_D1, -V_E1, -V_F1 /) )

    ! 2. common post order, wait order depends on the process
    w3 = perm( (/ V_A1, V_B1, V_C1 /) )
    call run_case( "wait order, 8-direction", (/ V_A1, V_B1, V_C1, -w3 /) )
    w3 = perm( (/ V_D1, V_E1, V_F1 /) )
    call run_case( "wait order, 4-direction", (/ V_D1, V_E1, V_F1, -w3 /) )
    w4 = perm( (/ V_A1, V_B1, V_A2, V_B2 /) )
    call run_case( "wait order, two grids", (/ V_A1, V_B1, V_A2, V_B2, -w4 /) )
    p4 = perm( (/ V_A1, V_A2, V_B1, V_B2 /) )
    w4 = perm( (/ V_B2, V_A1, V_B1, V_A2 /) )
    call run_case( "post and wait order, two grids", (/ p4, -w4 /) )

    ! 3. nested: a and b are posted, c is posted and waited, a and b are waited (numerical filter)
    p2 = perm( (/ V_A1, V_B1 /) )
    w2 = perm( (/ V_B1, V_A1 /) )
    call run_case( "nested, same grid", (/ p2, V_C1, -V_C1, -w2 /) )
    call run_case( "nested, other grid", (/ p2, V_A2, -V_A2, -w2 /) )

    ! 4. z is posted and waited between the posts of a and b in some processes, and after the
    !    wait of a in the others (skipped with the barriers of PROF_MPI_BARRIER, with which it hangs)
    if ( skip_case4 ) then
       if ( PRC_myrank == 0 ) write(*,*) 'test_comm_order: post/wait between posts skipped (SKIP_ASYMMETRIC)'
    else if ( odd_xy ) then
       call run_case( "post/wait between posts, same grid",  (/ V_A1, V_Z1, -V_Z1, V_B1, -V_A1, -V_B1 /) )
       call run_case( "post/wait between posts, other grid", (/ V_A1, V_Z2, -V_Z2, V_B1, -V_A1, -V_B1 /) )
    else
       call run_case( "post/wait between posts, same grid",  (/ V_A1, V_B1, -V_A1, V_Z1, -V_Z1, -V_B1 /) )
       call run_case( "post/wait between posts, other grid", (/ V_A1, V_B1, -V_A1, V_Z2, -V_Z2, -V_B1 /) )
    end if

    ! 5. persistent (8- and 4-direction) and non-persistent exchanges at the same time
    p4 = perm( (/ V_A1, V_N1, V_D1, V_A2 /) )
    w4 = perm( (/ V_D1, V_A2, V_N1, V_A1 /) )
    call run_case( "persistent and non-persistent", (/ p4, -w4 /) )

    ! 6. persistent and non-persistent (3-D, and 2-D in OP_2D) exchanges with the same packing
    !    buffer, one after another (they cannot be in progress at the same time)
    call run_case( "persistent and non-persistent, same buffer", &
                   (/ V_A1, -V_A1, V_M1, -V_M1, OP_2D, V_A1, -V_A1, OP_2D, V_M1, -V_M1, V_A1, -V_A1 /) )

    ! 7. all the 8-direction persistent variables of the two grids (the results kept for phase B)
    call run_case( "post and wait order, all the 8-direction persistent", &
                   (/ perm( (/ V_A1, V_B1, V_C1, V_Z1, V_Q1, V_Q2, V_Q3, V_Q4, V_A2, V_B2, V_Z2 /) ), &
                     -perm( (/ V_Z2, V_Q4, V_A2, V_Q2, V_C1, V_B2, V_Z1, V_A1, V_Q3, V_B1, V_Q1 /) ) /) )

  end subroutine run_cases

  !-----------------------------------------------------------------------------
  !> COMM_setup with the communication configuration iconf (the first configuration file if none
  !! is given), the grids, and the persistent variables of the fields in sel; the number of the
  !! persistent requests of each of them is checked
  subroutine setup_comm( sel )
    logical, intent(in) :: sel(NF)
    logical :: layout2
    integer :: f

    if ( nconf > 0 ) then
       ! COMM_setup reads its namelist from IO_FID_CONF
       call get_command_argument( iconf + 1, comm_conf )
       close( IO_FID_CONF )
       open( IO_FID_CONF, file=trim(comm_conf), status="old", form="formatted", action="read" )
    else
       comm_conf = "(the first configuration file)"
    end if
    LOG_INFO("test_comm_order",*) 'communication configuration: ', trim(comm_conf)

    call COMM_setup
    call COMM_regist( grd(1)%KA, grd(1)%IA, grd(1)%JA, grd(1)%IHALO, grd(1)%JHALO, grd(1)%gid )
    call COMM_regist( grd(2)%KA, grd(2)%IA, grd(2)%JA, grd(2)%IHALO, grd(2)%JHALO, grd(2)%gid )

    do f = 1, NF
       fld(f)%vid = fld(f)%pid
       if ( .not. ( fld(f)%persistent .and. sel(f) ) ) cycle
       if ( fld(f)%ndir == 8 ) then
          call COMM_vars8_init( fld(f)%name, fld(f)%v, fld(f)%vid, gid=grd(fld(f)%grid)%gid )
       else
          call COMM_vars_init ( fld(f)%name, fld(f)%v, fld(f)%vid, gid=grd(fld(f)%grid)%gid )
       end if
       ! without a persistent ID the exchange silently falls back to the non-persistent one
       if ( fld(f)%vid <= COMM_get_vsize_max() ) then
          call AssertEqual( "persistent ID of "//fld(f)%name//" ("//trim(comm_conf)//")", .true., .false. )
       end if
    end do

    call read_layout2( layout2 )
    call check_pc_info( sel, layout2 )
  end subroutine setup_comm

  !-----------------------------------------------------------------------------
  !> the layout 2 of the persistent exchange is selected in PARAM_COMM_CARTESC of IO_FID_CONF (the
  !! namelist of COMM_setup, with the same items and default values), and so is the batch of the
  !! groups (batch)
  subroutine read_layout2( layout2 )
    logical, intent(out) :: layout2

    integer :: COMM_vsize_max, COMM_vsize_max_pc
    logical :: COMM_USE_MPI_PC, COMM_USE_MPI_ONESIDED, COMM_USE_PC_LAYOUT2, COMM_USE_BATCH, COMM_USE_NONPC_LAYOUT2, COMM_USE_NONPC_BATCH
    namelist / PARAM_COMM_CARTESC / &
       COMM_vsize_max, &
       COMM_vsize_max_pc, &
       COMM_USE_MPI_PC, &
       COMM_USE_MPI_ONESIDED, &
       COMM_USE_PC_LAYOUT2, &
       COMM_USE_BATCH, &
       COMM_USE_NONPC_LAYOUT2, &
       COMM_USE_NONPC_BATCH
    integer :: ierr

    COMM_USE_MPI_PC = .true.
    COMM_USE_MPI_ONESIDED = .false.
#ifdef _OPENACC
    COMM_USE_PC_LAYOUT2 = .true.
    COMM_USE_BATCH = .true.
    COMM_USE_NONPC_LAYOUT2 = .true.
    COMM_USE_NONPC_BATCH = .true.
#else
    COMM_USE_PC_LAYOUT2 = .false.
    COMM_USE_BATCH = .false.
    COMM_USE_NONPC_LAYOUT2 = .false.
    COMM_USE_NONPC_BATCH = .false.
#endif
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_COMM_CARTESC, iostat=ierr )
    if ( ierr > 0 ) then
       call AssertEqual( "PARAM_COMM_CARTESC read by the test ("//trim(comm_conf)//")", 0, ierr )
    end if
    layout2 = COMM_USE_PC_LAYOUT2 .and. COMM_USE_MPI_PC .and. .not. COMM_USE_MPI_ONESIDED
    batch   = COMM_USE_BATCH .and. layout2 .and. .not. PRC_TwoD ! the batch is not used in 2-D
    batch_nonpc = COMM_USE_NONPC_BATCH
  end subroutine read_layout2

  !-----------------------------------------------------------------------------
  !> the namelist of the test in the first configuration file (before it is replaced by setup_comm)
  subroutine read_test_conf
    logical :: SKIP_ASYMMETRIC
    namelist / PARAM_TEST_COMM_ORDER / &
       SKIP_ASYMMETRIC
    integer :: ierr

    SKIP_ASYMMETRIC = .false.
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_TEST_COMM_ORDER, iostat=ierr )
    if ( ierr > 0 ) then
       call AssertEqual( "PARAM_TEST_COMM_ORDER read by the test", 0, ierr )
    end if
    skip_case4 = SKIP_ASYMMETRIC
    LOG_INFO("test_comm_order",*) 'skip the cases with posts and waits that differ among the processes: ', skip_case4
  end subroutine read_test_conf

  !-----------------------------------------------------------------------------
  !> phase B: the persistent variables are registered again (the members of the groups first), the
  !! groups are declared, and the other persistent variables are registered after them. The IDs
  !! must be those of the same registrations without groups, and the states of the groups those
  !! of the configuration
  subroutine setup_groups
    integer :: nreg(2) !< number of the persistent variables registered on each grid
    integer, allocatable :: none(:)
    integer :: n, k

    nreg(:) = 0
    do n = 1, size(ga)
       call regist_field( ga(n), nreg )
    end do
    do n = 1, size(lone)
       call regist_field( lone(n), nreg )
    end do
    call regist_field( V_Z1, nreg )
    call regist_field( V_Q3, nreg )
    call regist_field( V_Q4, nreg )
    call regist_field( V_A2, nreg )
    call regist_field( V_B2, nreg )
    call regist_field( V_Z2, nreg )

    ! not checked: no member, a non-persistent vid (before the group of one with the same variable)
    allocate( none(0) )
    call COMM_vars_group_init( "test-none", none, gid=grd(1)%gid )
    call COMM_vars_group_init( "test-nonpersistent", (/ fld(V_Q4)%vid, fld(V_N1)%pid /), gid=grd(1)%gid )
    call COMM_vars_group_init( "test-ga", fld(ga)%vid, gid=grd(1)%gid )
    call COMM_vars_group_init( "test-gp", fld(gp)%vid, gid=grd(1)%gid )
    ! groups of one with q4 registered again: the table of the groups (4 at first) is extended
    ! with ga and gp in it, which are exchanged afterwards (their requests and buffers are kept)
    do k = 1, 3
       call COMM_vars_group_init( "test-extra", fld(V_Q4:V_Q4)%vid, gid=grd(1)%gid )
       call regist_field( V_Q4, nreg )
    end do
    call COMM_vars_group_init( "test-g1", fld(g1)%vid, gid=grd(1)%gid )
    call COMM_vars_group_init( "test-gb", fld(gb)%vid, gid=grd(2)%gid )

    ! after the groups (the table of the persistent variables is extended with a small capacity)
    call regist_field( V_D1, nreg )
    call regist_field( V_E1, nreg )
    call regist_field( V_F1, nreg )
    do n = 1, NF
       if ( .not. fld(n)%persistent ) fld(n)%vid = fld(n)%pid
    end do

    call check_groups
  end subroutine setup_groups

  !> register a persistent field again: the ID must follow the registrations on its grid
  subroutine regist_field( f, nreg )
    integer, intent(in)    :: f
    integer, intent(inout) :: nreg(2)

    fld(f)%vid = fld(f)%pid
    if ( is_unbound( f ) ) then
       ! an ID that is not bound to the array (issue #3624)
       call COMM_vars8_init( fld(f)%name, fld(f)%vid, gid=grd(fld(f)%grid)%gid )
    else if ( fld(f)%ndir == 8 ) then
       call COMM_vars8_init( fld(f)%name, fld(f)%v, fld(f)%vid, gid=grd(fld(f)%grid)%gid )
    else
       call COMM_vars_init ( fld(f)%name, fld(f)%v, fld(f)%vid, gid=grd(fld(f)%grid)%gid )
    end if
    nreg(fld(f)%grid) = nreg(fld(f)%grid) + 1
    call AssertEqual( "ID of "//fld(f)%name//" registered with the groups ("//trim(comm_conf)//")", &
                      COMM_get_vsize_max() + nreg(fld(f)%grid), fld(f)%vid )
  end subroutine regist_field

  !> the field is registered with an ID not bound to its array: the members of the groups of two
  !! or more in phase C
  logical function is_unbound( f )
    integer, intent(in) :: f

    is_unbound = unbound .and. group_size( f ) >= 2
  end function is_unbound

  !> number of the members of the group of a field in phase B (0 if not in a group)
  function group_size( f ) result(n)
    integer, intent(in) :: f
    integer :: n

    n = 0
    if ( any( ga(:) == f ) ) n = size(ga)
    if ( any( gp(:) == f ) ) n = size(gp)
    if ( any( g1(:) == f ) ) n = size(g1)
    if ( any( gb(:) == f ) ) n = size(gb)
  end function group_size

  !-----------------------------------------------------------------------------
  !> the state of the group of each persistent field (COMM_get_pc_info): checked, or in batch with
  !! COMM_USE_BATCH if it has two or more members, with the number of the requests of the group
  !! from the neighbors; the requests of the variables themselves are unchanged
  subroutine check_groups
    integer :: f, n, nreq, state, nreq_group, nmem_group, state_ex, nreq_ex, nbatch
    character(len=H_LONG) :: label

    nbatch = 0
    do f = 1, NF
       if ( .not. fld(f)%persistent ) cycle
       label = trim(fld(f)%name)//" ("//trim(comm_conf)//")"
       call COMM_get_pc_info( fld(f)%vid, nreq=nreq, state=state, nreq_group=nreq_group, nmem_group=nmem_group, &
                             gid=grd(fld(f)%grid)%gid )
       n = group_size( f )
       if ( n == 0 ) then
          state_ex = COMM_PC_GROUP_NONE
       else if ( batch .and. n >= 2 .and. ( batch_nonpc .or. .not. is_unbound( f ) ) ) then
          state_ex = COMM_PC_GROUP_BATCH
       else
          state_ex = COMM_PC_GROUP_CHECK
       end if
       nreq_ex = 0
       if ( state_ex == COMM_PC_GROUP_BATCH ) then
          nreq_ex = expected_group_nreq( fld(f)%grid )
          nbatch = nbatch + 1
       end if
       LOG_INFO("test_comm_order",*) 'group of ', trim(label), ': state, members, requests of the group = ', &
                                     state, nmem_group, nreq_group
       call AssertEqual( "group of "//trim(label), state_ex, state )
       call AssertEqual( "members of the group of "//trim(label), n, nmem_group )
       call AssertEqual( "requests of the group of "//trim(label), nreq_ex, nreq_group )
       call AssertEqual( "requests of the group of "//trim(label)//" at most", .true., nreq_group <= NREQ_GROUP_MAX )
       if ( is_unbound( f ) ) then
          call AssertEqual( "number of the persistent requests of "//trim(label)//" (not bound to the array)", 0, nreq )
       else
          call AssertEqual( "number of the persistent requests of "//trim(label)//" with the groups", &
                            expected_nreq( f, read_layout2_value() ), nreq )
       end if
    end do
    if ( batch .and. ( batch_nonpc .or. .not. unbound ) .and. nbatch == 0 ) then
       call AssertEqual( "variables in a group in batch ("//trim(comm_conf)//")", 1, 0 )
    end if
  end subroutine check_groups

  logical function read_layout2_value()
    call read_layout2( read_layout2_value )
  end function read_layout2_value

  !> number of the requests of a group in batch: one message to and from each neighbor in y, and
  !! in x and diagonal (not with PRC_TwoD or IHALO = 0)
  function expected_group_nreq( ig ) result(n)
    integer, intent(in) :: ig
    integer :: n

    integer :: nNS, nWE, ndiag

    nNS = count( (/ PRC_HAS_N, PRC_HAS_S /) )
    nWE = count( (/ PRC_HAS_W, PRC_HAS_E /) )
    ndiag = count( (/ PRC_HAS_N .and. PRC_HAS_W, PRC_HAS_N .and. PRC_HAS_E, &
                      PRC_HAS_S .and. PRC_HAS_W, PRC_HAS_S .and. PRC_HAS_E /) )
    n = 2 * nNS
    if ( .not. PRC_TwoD .and. grd(ig)%IHALO > 0 ) n = n + 2 * ( nWE + ndiag )
  end function expected_group_nreq

  !-----------------------------------------------------------------------------
  !> the cases of phase B: all the members of a group are posted before any of them is waited
  subroutine run_group_cases
    integer, allocatable :: pa(:), wa(:)

    ! posts and waits of a group in orders that depend on the process (issue #3594 (1))
    call run_gcase( "group: post and wait order", (/ perm(ga), lone, -perm(ga(size(ga):1:-1)), -lone /) )

    ! the groups of the two grids, their posts interleaved
    call run_gcase( "groups on two grids", &
                    (/ perm( (/ V_A1, V_A2, V_B1, V_B2, V_C1, V_Z2, V_Q1, V_Q2 /) ), &
                      -perm( (/ V_Z2, V_Q2, V_A1, V_B2, V_C1, V_A2, V_Q1, V_B1 /) ) /) )

    ! a group posted after the other (the form of PROG), the waits of both interleaved
    pa = (/ V_Q3, ga, lone, V_Z1 /)
    call run_gcase( "group posted after a group", (/ perm(ga), lone, perm(gp), -perm(pa) /) )

    ! with a group of one, a variable alone (4-direction) and a non-persistent exchange
    pa = (/ ga, lone, V_Q4, V_D1, V_N1 /)
    wa = (/ V_N1, V_D1, V_Q4, lone, ga(size(ga):1:-1) /)
    call run_gcase( "group with other exchanges", (/ perm(pa), -perm(wa) /) )

    ! an exchange of another variable between the waits of a group (the round is in progress, but
    ! the posts are complete)
    wa = perm( ga(size(ga):1:-1) )
    call run_gcase( "exchange between the waits of a group", &
                    (/ perm(ga), lone, -wa(1), V_D1, -V_D1, -wa(2:), -lone /) )

    ! exchanges of other variables posted and waited between the posts of a group, at the same
    ! place in all the processes (a 4-direction persistent variable, a non-persistent one, the
    ! 8-direction persistent variables of another group of the grid, 2-D exchanges with the packing
    ! buffer of the non-persistent one, and the group of grid 2); the first member posted depends on
    ! the process
    pa = perm(ga)
    call run_gcase( "exchanges between the posts of a group", &
                    (/ pa(1), V_D1, -V_D1, V_N1, -V_N1, perm(gp), -perm(gp), OP_2D, perm(gb), -perm(gb), pa(2:), lone, &
                       -perm(ga), -lone /) )

    ! two rounds of a group in a case
    call run_gcase( "two rounds of a group", &
                    (/ perm(ga), lone, -perm(ga), -lone, perm(ga(size(ga):1:-1)), lone, -perm(ga(size(ga):1:-1)), -lone /) )
  end subroutine run_group_cases

  !-----------------------------------------------------------------------------
  !> run a sequence of posts and waits of phase B (see run_case). The result is compared with the
  !! source cells and with the result kept from phase A of the first configuration (which also
  !! gives the halos without a source cell); the reference exchange of one variable at a time is
  !! not made (it is not allowed for a member of a group)
  subroutine run_gcase( cname, ops )
    character(len=*), intent(in) :: cname
    integer,          intent(in) :: ops(:)

    logical :: used(NF)
    integer :: n, r, f, ncheck, nhalo, ib

    used(:) = .false.
    do n = 1, size(ops)
       if ( abs(ops(n)) <= NF ) used(abs(ops(n))) = .true.
    end do
    ib = merge( 1, 0, fill_bnd )

    nhalo = 0
    do r = 1, NROUND
       do f = 1, NF
          if ( .not. used(f) ) cycle
          if ( .not. allocated( saved(f,r,ib)%v ) ) then
             call AssertEqual( cname//": "//fld(f)%name//" (result of phase A kept)", .true., .false. )
          end if
          fld(f)%ref(:,:,:) = saved(f,r,ib)%v(:,:,:)
          call make_expect( f, r, ncheck )
          nhalo = nhalo + ncheck
       end do

       ! other values in the buffers: all the variables are posted, and then waited, in the same
       ! order in all the processes (a round of each group)
       call fill_all( r + NROUND, used )
       do f = 1, NF
          if ( used(f) ) call post( f )
       end do
       do f = 1, NF
          if ( used(f) ) call wait( f )
       end do
       call to_host

       ! the case
       call fill_all( r, used )
       do n = 1, size(ops)
          if ( ops(n) == OP_2D ) then
             ! with the packing buffer of n1 (that of a1 may be in use: a1 is a member of ga)
             call exchange_2d( r, fld(V_N1)%pid )
          else if ( ops(n) > 0 ) then
             call post( ops(n) )
          else
             call wait( -ops(n) )
          end if
       end do
       call to_host
       do f = 1, NF
          if ( .not. used(f) ) cycle
          call check( cname//": "//fld(f)%name//" (source cells)", fld(f)%expect, fld(f)%v )
          call compare_saved( cname, f, r )
       end do
    end do

    if ( has_neighbor .and. nhalo == 0 ) then
       call AssertEqual( cname//": number of halo cells compared with the source", 1, 0 )
    end if

    LOG_INFO("test_comm_order",*) trim(cname), ': OK (FILL_BND = ', fill_bnd, &
                                  ', halo cells compared with the source: ', nhalo, ')'
    if ( PRC_myrank == 0 ) write(*,*) 'test_comm_order: ', trim(cname), ' OK (FILL_BND = ', fill_bnd, ')'
  end subroutine run_gcase

  !-----------------------------------------------------------------------------
  !> in the setups again: a group of a1 and b1 (in the order of irep), a round checked with the
  !! source cells; then a1 and b1 are registered again without the release (the form of the setups
  !! again of the dynamics tests), and a group of the new IDs is declared and exchanged
  subroutine run_group_single
    integer :: rep, n, f, ncheck, nreq, state, nreq_group, nmem_group, state_ex
    integer :: mem(2), ord(2)

    if ( mod( irep, 2 ) == 1 ) then
       mem(:) = (/ V_A1, V_B1 /)
    else
       mem(:) = (/ V_B1, V_A1 /)
    end if
    state_ex = COMM_PC_GROUP_CHECK
    if ( batch ) state_ex = COMM_PC_GROUP_BATCH

    do rep = 1, 2
       if ( rep == 2 ) then
          do n = 1, 2
             f = mem(n)
             fld(f)%vid = fld(f)%pid
             call COMM_vars8_init( fld(f)%name, fld(f)%v, fld(f)%vid, gid=grd(fld(f)%grid)%gid )
          end do
       end if
       call COMM_vars_group_init( "test-again", fld(mem)%vid, gid=grd(1)%gid )
       do n = 1, 2
          f = mem(n)
          call COMM_get_pc_info( fld(f)%vid, nreq=nreq, state=state, nreq_group=nreq_group, nmem_group=nmem_group, &
                                gid=grd(1)%gid )
          call AssertEqual( "group of "//fld(f)%name//", setup again ("//trim(comm_conf)//")", state_ex, state )
       end do

       call fill_all( 1 )
       ord = perm( mem )
       do n = 1, 2
          call post( ord(n) )
       end do
       do n = 2, 1, -1
          call wait( ord(n) )
       end do
       call to_host
       do n = 1, 2
          f = mem(n)
          fld(f)%ref(:,:,:) = fld(f)%v(:,:,:)
          call make_expect( f, 1, ncheck )
          call check( "group, setup again ("//trim(comm_conf)//"): "//fld(f)%name//" (source cells)", fld(f)%expect, fld(f)%v )
       end do
    end do
  end subroutine run_group_single

  !-----------------------------------------------------------------------------
  !> the number of the persistent requests of each variable registered (fields in sel), from the
  !! neighbors (see expected_nreq), and no group; 0 for a non-persistent vid. With the layout 2, it
  !! must differ from the current layout for some variable where a process has a neighbor.
  subroutine check_pc_info( sel, layout2 )
    logical, intent(in) :: sel(NF)
    logical, intent(in) :: layout2

    integer :: f, nreq, state, ndiffer
    character(len=H_LONG) :: label

    ndiffer = 0
    do f = 1, NF
       label = trim(fld(f)%name)//" ("//trim(comm_conf)//")"
       if ( fld(f)%persistent .and. sel(f) ) then
          call COMM_get_pc_info( fld(f)%vid, nreq=nreq, state=state, gid=grd(fld(f)%grid)%gid )
          LOG_INFO("test_comm_order",*) 'number of the persistent requests of ', trim(label), ': ', nreq
          call AssertEqual( "number of the persistent requests of "//trim(label), expected_nreq( f, layout2 ), nreq )
          call AssertEqual( "group of "//trim(label), COMM_PC_GROUP_NONE, state )
          if ( nreq /= expected_nreq( f, .false. ) ) ndiffer = ndiffer + 1
       else
          call COMM_get_pc_info( fld(f)%pid, nreq=nreq, state=state, gid=grd(fld(f)%grid)%gid )
          call AssertEqual( "number of the persistent requests of a vid, "//trim(label), 0, nreq )
       end if
    end do

    if ( layout2 .and. has_neighbor .and. any( sel(:) .and. fld(:)%ndir == 8 .and. fld(:)%persistent ) ) then
       if ( ndiffer == 0 ) then
          call AssertEqual( "variables in the layout 2 ("//trim(comm_conf)//")", 1, 0 )
       end if
    end if
  end subroutine check_pc_info

  !-----------------------------------------------------------------------------
  !> number of the persistent requests of a variable (sends and receives), from the neighbors.
  !! Layout 2 (8-direction, not in a grid with PRC_TwoD and IHALO > 0): one message to and from
  !! each neighbor in y, in x (IHALO > 0) and diagonal (IHALO > 0, where the neighbors in x and in
  !! y exist). Current 8-direction: JHALO messages to and from each neighbor in y, one to and from
  !! each neighbor in x, and JHALO for each corner where a neighbor in x or in y exists (not in
  !! 2-D). Current 4-direction: one to and from each neighbor in y and in x (not in 2-D).
  function expected_nreq( f, layout2 ) result(n)
    integer, intent(in) :: f
    logical, intent(in) :: layout2
    integer :: n

    integer :: nNS, nWE, ndiag, nany

    nNS = count( (/ PRC_HAS_N, PRC_HAS_S /) )
    nWE = count( (/ PRC_HAS_W, PRC_HAS_E /) )
    ndiag = count( (/ PRC_HAS_N .and. PRC_HAS_W, PRC_HAS_N .and. PRC_HAS_E, &
                      PRC_HAS_S .and. PRC_HAS_W, PRC_HAS_S .and. PRC_HAS_E /) )
    nany  = count( (/ PRC_HAS_N .or.  PRC_HAS_W, PRC_HAS_N .or.  PRC_HAS_E, &
                      PRC_HAS_S .or.  PRC_HAS_W, PRC_HAS_S .or.  PRC_HAS_E /) )
    associate( g => grd(fld(f)%grid) )
    if ( fld(f)%ndir == 4 ) then
       n = 2 * nNS
       if ( .not. PRC_TwoD ) n = n + 2 * nWE
    else if ( layout2 .and. .not. ( PRC_TwoD .and. g%IHALO > 0 ) ) then
       n = 2 * nNS
       if ( g%IHALO > 0 ) n = n + 2 * ( nWE + ndiag )
    else
       n = 2 * g%JHALO * nNS
       if ( .not. PRC_TwoD ) n = n + 2 * nWE + 2 * g%JHALO * nany
    end if
    end associate
  end function expected_nreq

  !-----------------------------------------------------------------------------
  !> a single exchange of each field in fs, one after another, checked with the source cells
  subroutine run_single( fs )
    integer, intent(in) :: fs(:)
    integer :: n, f, ncheck

    call fill_all( 1 )
    do n = 1, size(fs)
       call post( fs(n) )
       call wait( fs(n) )
    end do
    call to_host
    do n = 1, size(fs)
       f = fs(n)
       fld(f)%ref(:,:,:) = fld(f)%v(:,:,:)
       call make_expect( f, 1, ncheck )
       call check( "setup again ("//trim(comm_conf)//"): "//fld(f)%name//" (source cells)", fld(f)%expect, fld(f)%v )
    end do
  end subroutine run_single

  !-----------------------------------------------------------------------------
  !> 2-D exchanges (8- and 4-direction) of s1 with the packing buffer of a1 (OP_2D in a case),
  !! checked with the source cells
  subroutine exchange_2d( r, pid )
    integer, intent(in) :: r
    integer, intent(in), optional :: pid !< packing buffer (default: that of a1)
    integer :: ndir, pid_

    pid_ = fld(V_A1)%pid
    if ( present(pid) ) pid_ = pid

    do ndir = 8, 4, -4
       call fill_2d( r )
       !$acc update device(s1)
       if ( ndir == 8 ) then
          call COMM_vars8( s1, pid_, gid=grd(1)%gid )
       else
          call COMM_vars ( s1, pid_, gid=grd(1)%gid )
       end if
       call COMM_wait( s1, pid_, fill_bnd, gid=grd(1)%gid )
       !$acc wait
       !$acc update host(s1)
       call check_2d( r, ndir )
    end do
  end subroutine exchange_2d

  subroutine fill_2d( r )
    integer, intent(in) :: r
    integer  :: i, j
    real(RP) :: sgn

    sgn = real( 3 - 2 * r, RP )
    associate( g => grd(1) )
    do j = 1, g%JA
    do i = 1, g%IA
       if ( i >= g%IS .and. i <= g%IE .and. j >= g%JS .and. j <= g%JE ) then
          s1(i,j) = sgn * inner_value( OP_2D, gi(g,i), gj(g,j), 1 )
       else
          s1(i,j) = sgn * ( 2.0E7_RP + real(OP_2D,RP) * 1.0E5_RP + real(PRC_myrank,RP) * 1.0E4_RP &
                          + real(i,RP) * 100.0_RP + real(j,RP) )
       end if
    end do
    end do
    end associate
  end subroutine fill_2d

  !> the halos whose source cell is in the domain (as in make_expect)
  subroutine check_2d( r, ndir )
    integer, intent(in) :: r, ndir
    integer  :: i, j, ig, jg, nxg, nyg
    logical  :: inner_i, inner_j, defined
    real(RP) :: sgn, val
    character(len=H_SHORT) :: label

    write(label,'(a,i1,a)') '2-D, ', ndir, '-direction'
    sgn = real( 3 - 2 * r, RP )
    associate( g => grd(1) )
    nxg = ( g%IE - g%IS + 1 ) * PRC_NUM_X
    nyg = ( g%JE - g%JS + 1 ) * PRC_NUM_Y
    do j = 1, g%JA
    do i = 1, g%IA
       inner_i = i >= g%IS .and. i <= g%IE
       inner_j = j >= g%JS .and. j <= g%JE
       ig = gi(g,i)
       jg = gj(g,j)
       defined = .true.
       if ( ig < 1 .or. ig > nxg ) then
          if ( PRC_PERIODIC_X .and. .not. PRC_TwoD ) then
             ig = modulo( ig - 1, nxg ) + 1
          else
             defined = .false.
          end if
       end if
       if ( jg < 1 .or. jg > nyg ) then
          if ( PRC_PERIODIC_Y ) then
             jg = modulo( jg - 1, nyg ) + 1
          else
             defined = .false.
          end if
       end if
       if ( ndir == 4 .and. .not. inner_i .and. .not. inner_j ) defined = .false.
       if ( .not. defined ) cycle
       val = sgn * inner_value( OP_2D, ig, jg, 1 )
       if ( s1(i,j) /= val ) then
          call AssertEqual( "persistent and non-persistent, same buffer: s1 ("//trim(label)//", source cells)", &
                            val, s1(i,j) )
       end if
    end do
    end do
    end associate
  end subroutine check_2d

  !-----------------------------------------------------------------------------
  !> run a sequence of posts (positive field index) and waits (negative field index)
  subroutine run_case( cname, ops )
    character(len=*), intent(in) :: cname
    integer,          intent(in) :: ops(:)

    logical :: used(NF)
    integer :: n, r, f, ncheck, nhalo

    used(:) = .false.
    do n = 1, size(ops)
       if ( abs(ops(n)) <= NF ) used(abs(ops(n))) = .true.
    end do

    nhalo = 0
    do r = 1, NROUND
       ! reference: one variable at a time, in the same order in all the processes
       call fill_all( r, used )
       do f = 1, NF
          if ( .not. used(f) ) cycle
          call post( f )
          call wait( f )
       end do
       call to_host
       do f = 1, NF
          if ( .not. used(f) ) cycle
          fld(f)%ref(:,:,:) = fld(f)%v(:,:,:)
          call make_expect( f, r, ncheck )
          nhalo = nhalo + ncheck
          call check( cname//" (reference): "//fld(f)%name, fld(f)%expect, fld(f)%v )
       end do

       ! other values in the buffers
       call fill_all( r + NROUND, used )
       do f = 1, NF
          if ( .not. used(f) ) cycle
          call post( f )
          call wait( f )
       end do
       call to_host

       ! the case
       call fill_all( r, used )
       do n = 1, size(ops)
          if ( ops(n) == OP_2D ) then
             call exchange_2d( r )
          else if ( ops(n) > 0 ) then
             call post( ops(n) )
          else
             call wait( -ops(n) )
          end if
       end do
       call to_host
       do f = 1, NF
          if ( .not. used(f) ) cycle
          call check( cname//": "//fld(f)%name//" (source cells)", fld(f)%expect, fld(f)%v )
          call check( cname//": "//fld(f)%name//" (reference)",    fld(f)%ref,    fld(f)%v )
          call compare_saved( cname, f, r )
       end do
    end do

    ! the comparison with the source cells must not be empty
    if ( has_neighbor .and. nhalo == 0 ) then
       call AssertEqual( cname//": number of halo cells compared with the source", 1, 0 )
    end if

    LOG_INFO("test_comm_order",*) trim(cname), ': OK (FILL_BND = ', fill_bnd, &
                                  ', halo cells compared with the source: ', nhalo, ')'
    if ( PRC_myrank == 0 ) write(*,*) 'test_comm_order: ', trim(cname), ' OK (FILL_BND = ', fill_bnd, ')'

  end subroutine run_case

  !-----------------------------------------------------------------------------
  subroutine post( f )
    integer, intent(in) :: f

    if ( fld(f)%ndir == 8 ) then
       call COMM_vars8( fld(f)%v, fld(f)%vid, gid=grd(fld(f)%grid)%gid )
    else
       call COMM_vars ( fld(f)%v, fld(f)%vid, gid=grd(fld(f)%grid)%gid )
    end if
  end subroutine post

  subroutine wait( f )
    integer, intent(in) :: f

    call COMM_wait( fld(f)%v, fld(f)%vid, fill_bnd, gid=grd(fld(f)%grid)%gid )
  end subroutine wait

  !-----------------------------------------------------------------------------
  !> keep the result of the first communication configuration (the first case that uses the
  !! field), and compare the others with it
  subroutine compare_saved( cname, f, r )
    character(len=*), intent(in) :: cname
    integer,          intent(in) :: f, r

    integer :: ib

    ib = merge( 1, 0, fill_bnd )
    if ( .not. allocated( saved(f,r,ib)%v ) ) then
       if ( iconf > 1 ) then ! every result is kept with the first configuration
          call AssertEqual( cname//": "//fld(f)%name//" (result of the first configuration kept)", &
                            .true., .false. )
       end if
       allocate( saved(f,r,ib)%v, source=fld(f)%v )
    else
       call check( cname//": "//fld(f)%name//" (first configuration, "//trim(comm_conf)//")", &
                   saved(f,r,ib)%v, fld(f)%v )
       if ( iconf > 1 ) ncompare = ncompare + 1
    end if
  end subroutine compare_saved

  !-----------------------------------------------------------------------------
  !> compare the whole array; AssertEqual (which aborts with the first difference) only on a mismatch
  subroutine check( message, answer, val )
    character(len=*), intent(in) :: message
    real(RP),         intent(in) :: answer(:,:,:)
    real(RP),         intent(in) :: val(:,:,:)

    if ( any( answer(:,:,:) /= val(:,:,:) ) ) then
       call AssertEqual( message, answer, val )
    end if
  end subroutine check

  !-----------------------------------------------------------------------------
  !> the inner cells: unique integer values (exact in single precision) of the variable and the
  !! global position; the halos: sentinel values of the variable, the process and the position.
  !! The sign changes with the round; r > NROUND gives the values of round r - NROUND halved, which
  !! differ from those of all the rounds (exactly in single precision too). Only the fields in used
  !! if it is given
  subroutine fill_all( r, used )
    integer, intent(in) :: r
    logical, intent(in), optional :: used(NF)
    integer  :: f, i, j, k
    real(RP) :: sgn

    if ( r > NROUND ) then
       sgn = real( 3 - 2 * ( r - NROUND ), RP ) * 0.5_RP ! 0.5, -0.5
    else
       sgn = real( 3 - 2 * r, RP ) ! 1, -1
    end if
    do f = 1, NF
       if ( present(used) ) then
          if ( .not. used(f) ) cycle
       end if
       associate( g => grd(fld(f)%grid), v => fld(f)%v )
       do j = 1, g%JA
       do i = 1, g%IA
       do k = 1, g%KA
          if ( i >= g%IS .and. i <= g%IE .and. j >= g%JS .and. j <= g%JE ) then
             v(k,i,j) = sgn * inner_value( f, gi(g,i), gj(g,j), k )
          else
             v(k,i,j) = sgn * ( 2.0E7_RP + real(f,RP) * 1.0E5_RP + real(PRC_myrank,RP) * 1.0E4_RP &
                              + real(i,RP) * 100.0_RP + real(j,RP) + real(k,RP) * 0.01_RP )
          end if
       end do
       end do
       end do
       end associate
    end do
    !$acc update device(a1, b1, c1, d1, e1, f1, z1, n1, a2, b2, z2, m1, q1, q2, q3, q4)
  end subroutine fill_all

  subroutine to_host
    !$acc wait
    !$acc update host(a1, b1, c1, d1, e1, f1, z1, n1, a2, b2, z2, m1, q1, q2, q3, q4)
  end subroutine to_host

  !> unique while ig, jg and k are below 100
  function inner_value( f, ig, jg, k ) result(val)
    integer, intent(in) :: f, ig, jg, k
    real(RP) :: val

    val = real(f,RP) * 1.0E6_RP + real(jg,RP) * 1.0E4_RP + real(ig,RP) * 100.0_RP + real(k,RP)
  end function inner_value

  !> global x index of the local index i (out of 1:IMAX*PRC_NUM_X in the halos)
  function gi( g, i )
    type(grid_t), intent(in) :: g
    integer,      intent(in) :: i
    integer :: gi

    gi = PRC_2Drank(PRC_myrank,1) * ( g%IE - g%IS + 1 ) + i - g%IS + 1
  end function gi

  function gj( g, j )
    type(grid_t), intent(in) :: g
    integer,      intent(in) :: j
    integer :: gj

    gj = PRC_2Drank(PRC_myrank,2) * ( g%JE - g%JS + 1 ) + j - g%JS + 1
  end function gj

  !-----------------------------------------------------------------------------
  !> expected values: the value of the source cell where it is in the domain (with periodicity),
  !! the reference elsewhere (halos at the non-periodic boundaries, the corners of the 4-direction
  !! exchange, which come from the halos of the neighbors in y)
  subroutine make_expect( f, r, ncheck )
    integer, intent(in)  :: f, r
    integer, intent(out) :: ncheck

    integer  :: i, j, k, ig, jg, nxg, nyg
    logical  :: inner_i, inner_j, defined
    real(RP) :: sgn

    sgn = real( 3 - 2 * r, RP )
    ncheck = 0
    associate( g => grd(fld(f)%grid), ex => fld(f)%expect )
    ex(:,:,:) = fld(f)%ref(:,:,:)
    nxg = ( g%IE - g%IS + 1 ) * PRC_NUM_X
    nyg = ( g%JE - g%JS + 1 ) * PRC_NUM_Y
    do j = 1, g%JA
    do i = 1, g%IA
       inner_i = i >= g%IS .and. i <= g%IE
       inner_j = j >= g%JS .and. j <= g%JE
       ig = gi(g,i)
       jg = gj(g,j)
       defined = .true.
       if ( ig < 1 .or. ig > nxg ) then
          if ( PRC_PERIODIC_X .and. .not. PRC_TwoD ) then
             ig = modulo( ig - 1, nxg ) + 1
          else
             defined = .false.
          end if
       end if
       if ( jg < 1 .or. jg > nyg ) then
          if ( PRC_PERIODIC_Y ) then
             jg = modulo( jg - 1, nyg ) + 1
          else
             defined = .false.
          end if
       end if
       if ( fld(f)%ndir == 4 .and. .not. inner_i .and. .not. inner_j ) defined = .false.
       if ( .not. defined ) cycle
       if ( .not. ( inner_i .and. inner_j ) ) ncheck = ncheck + 1
       do k = 1, g%KA
          ex(k,i,j) = sgn * inner_value( f, ig, jg, k )
       end do
    end do
    end do
    end associate
  end subroutine make_expect

  !-----------------------------------------------------------------------------
  subroutine set_grid( g, KA_, IHALO_, JHALO_ )
    type(grid_t), intent(out) :: g
    integer,      intent(in)  :: KA_, IHALO_, JHALO_

    g%KA    = KA_
    g%IHALO = IHALO_
    g%JHALO = JHALO_
    g%IA    = IMAX + IHALO_ * 2
    g%JA    = JMAX + JHALO_ * 2
    g%IS    = IHALO_ + 1
    g%IE    = IHALO_ + IMAX
    g%JS    = JHALO_ + 1
    g%JE    = JHALO_ + JMAX
    g%gid   = 0
  end subroutine set_grid

  subroutine alloc( v, ig )
    real(RP), allocatable, intent(inout) :: v(:,:,:)
    integer,               intent(in)    :: ig

    allocate( v(grd(ig)%KA,grd(ig)%IA,grd(ig)%JA) )
    v(:,:,:) = 0.0_RP
  end subroutine alloc

  subroutine set_field( f, name, v, ig, ndir, persistent, pid )
    integer,          intent(in)            :: f
    character(len=*), intent(in)            :: name
    real(RP),         intent(inout), target :: v(:,:,:)
    integer,          intent(in)            :: ig, ndir, pid
    logical,          intent(in)            :: persistent

    fld(f)%name       = name
    fld(f)%v          => v
    fld(f)%grid       = ig
    fld(f)%ndir       = ndir
    fld(f)%persistent = persistent
    fld(f)%pid        = pid
    fld(f)%vid        = pid
    allocate( fld(f)%ref   (grd(ig)%KA,grd(ig)%IA,grd(ig)%JA) )
    allocate( fld(f)%expect(grd(ig)%KA,grd(ig)%IA,grd(ig)%JA) )
  end subroutine set_field

  !-----------------------------------------------------------------------------
  !> order of the fields in this process. More than two fields: rotated by y and reversed where x is
  !! odd, so that all the neighbors have different orders (a reversed order is not a rotation).
  !! Two fields: reversed where x + y is odd (the diagonal neighbors have the same order).
  function perm( a ) result(b)
    integer, intent(in) :: a(:)
    integer :: b(size(a))
    integer :: n, m

    n = size(a)
    if ( n > 2 ) then
       do m = 1, n
          b(m) = a( mod( m - 1 + iy, n ) + 1 )
       end do
       if ( odd_x ) b(:) = b(n:1:-1)
    else
       b(:) = a(:)
       if ( odd_xy ) b(:) = b(n:1:-1)
    end if
  end function perm

end program test_comm_order
