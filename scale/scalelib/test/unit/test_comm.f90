module test_comm

  !-----------------------------------------------------------------------------
  use scale_precision
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  public :: test_comm_run

  real(RP) :: data_P
  real(RP) :: data_W
  real(RP) :: data_N
  real(RP) :: data_E
  real(RP) :: data_S
  real(RP) :: data_NW
  real(RP) :: data_NE
  real(RP) :: data_SW
  real(RP) :: data_SE

  character(len=1), parameter :: dir_name(2) = (/ "x", "y" /)

  ! The persistent requests are bound to the memory of the registered array: each test releases
  ! the persistent variables of its grids (COMM_vars_finalize) before it frees the arrays.

contains

  subroutine test_comm_run
    use scale_atmos_grid_cartesC_index
    use scale_prc, only: &
       PRC_myrank
    use scale_prc_cartesC, only: &
       PRC_next, &
       PRC_W, &
       PRC_N, &
       PRC_E, &
       PRC_S, &
       PRC_NW, &
       PRC_NE, &
       PRC_SW, &
       PRC_SE, &
       PRC_HAS_W, &
       PRC_HAS_N, &
       PRC_HAS_E, &
       PRC_HAS_S
    use scale_comm_cartesC, only: &
       COMM_regist
    implicit none

    integer, parameter :: KA2 = 10, IA2 = 12, IHALO2 = 4, JA2 = 9, JHALO2 = 3
    integer :: IS2, IE2, JS2, JE2
    integer :: gid
    integer :: gid_prev
    !===========================================================================

    data_P = real(PRC_myrank, RP)
    if ( PRC_HAS_W ) then
       data_W = real(PRC_next(PRC_W), RP)
    else
       data_W = data_P
    end if
    if ( PRC_HAS_N ) then
       data_N = real(PRC_next(PRC_N), RP)
    else
       data_N = data_P
    end if
    if ( PRC_HAS_E ) then
       data_E = real(PRC_next(PRC_E), RP)
    else
       data_E = data_P
    end if
    if ( PRC_HAS_S ) then
       data_S = real(PRC_next(PRC_S), RP)
    else
       data_S = data_P
    end if
    if ( PRC_HAS_N ) then
       if ( PRC_HAS_W ) then
          data_NW = real(PRC_next(PRC_NW), RP)
       else
          data_NW = real(PRC_next(PRC_N), RP)
       end if
       if ( PRC_HAS_E ) then
          data_NE = real(PRC_next(PRC_NE), RP)
       else
          data_NE = real(PRC_next(PRC_N), RP)
       end if
    else
       if ( PRC_HAS_W ) then
          data_NW = real(PRC_next(PRC_W), RP)
       else
          data_NW = data_P
       end if
       if ( PRC_HAS_E ) then
          data_NE = real(PRC_next(PRC_E), RP)
       else
          data_NE = data_P
       end if
    end if
    if ( PRC_HAS_S ) then
       if ( PRC_HAS_W ) then
          data_SW = real(PRC_next(PRC_SW), RP)
       else
          data_SW = real(PRC_next(PRC_S), RP)
       end if
       if ( PRC_HAS_E ) then
          data_SE = real(PRC_next(PRC_SE), RP)
       else
          data_SE = real(PRC_next(PRC_S), RP)
       end if
    else
       if ( PRC_HAS_W ) then
          data_SW = real(PRC_next(PRC_W), RP)
       else
          data_SW = data_P
       end if
       if ( PRC_HAS_E ) then
          data_SE = real(PRC_next(PRC_E), RP)
       else
          data_SE = data_P
       end if
    end if

    gid = 1
    call test_vars( KA, IA, IS, IE, JA, JS, JE, gid )

    call test_vars8( KA, IA, IS, IE, JA, JS, JE, gid )

    IS2 = IHALO2 + 1
    IE2 = IA2 - IHALO2
    JS2 = JHALO2 + 1
    JE2 = JA2 - JHALO2

    call COMM_regist( KA2, IA2, JA2, IHALO2, JHALO2, gid )

    call test_vars( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid )

    call test_vars8( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid )

    gid_prev = gid
    call COMM_regist( KA2, IA2, JA2, IHALO2, JHALO2, gid )
    call test_vars8_two_grids( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid_prev, gid )
    call test_vars_two_grids ( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid_prev, gid )

    ! Each test releases the persistent variables of its grids at its end, so every grid has no
    ! persistent variables at the start of the next test.

    ! persistent variables with the same sequential number on two grids
    gid_prev = gid
    call COMM_regist( KA2, IA2, JA2, IHALO2, JHALO2, gid )
    call test_vars8_pc_two_grids( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid_prev, gid )
    call test_vars8_valid_ids   ( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid_prev, gid )

    ! a persistent variable and a non-persistent one on the same grid.
    ! The 8-direction test is on grid 1, where the tags of vars8_3D_mpi used to collide
    call test_vars8_pc_same_grid( KA, IA, IS, IE, JA, JS, JE, 1 )

    ! a persistent variable (4-direction) and a non-persistent one on the same grid:
    ! the sequential number of the persistent variable is 1
    call COMM_regist( KA2, IA2, JA2, IHALO2, JHALO2, gid )
    call test_vars_pc_same_grid( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid )

    ! more persistent variables than the initial capacity of a new grid
    call COMM_regist( KA2, IA2, JA2, IHALO2, JHALO2, gid )
    call test_vars8_pc_grow( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid )

    ! release and register again
    call COMM_regist( KA2, IA2, JA2, IHALO2, JHALO2, gid )
    call test_vars8_pc_reregister( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid )

  end subroutine test_comm_run
!=============================================================================


  subroutine test_vars( KA, IA, IS, IE, JA, JS, JE, gid )
    use scale_prc, only: &
         PRC_myrank
    use scale_comm_cartesC, only: &
         COMM_vars_init, &
         COMM_vars_finalize, &
         COMM_vars, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid

    real(RP) :: data3d(KA,IA,JA), data2d(IA,JA)
    real(RP), allocatable :: pc3d(:,:,:)
    integer :: vid

    write(*,*) "Test vars"


    ! MPI
    data2d(:,:) = data_P
    call COMM_vars( data2d, 1, gid=gid )
    call COMM_wait( data2d, 1, .false., gid=gid )
    call check_vars2d( IA, IS, IE, JA, JS, JE, data2d )

    data3d(:,:,:) = data_P
    call COMM_vars( data3d, 1, gid=gid )
    call COMM_wait( data3d, 1, .false., gid=gid )
    call check_vars3d( KA, IA, IS, IE, JA, JS, JE, data3d )

    ! MPI PC
!!$    data2d(:,:) = data_P
!!$    vid = 1
!!$    call COMM_vars_init( 'testdata2d', data2d, vid, gid=gid )
!!$    call COMM_vars( data2d, vid, gid=gid )
!!$    call COMM_wait( data2d, vid, .false., gid=gid )
!!$    call check_vars2d( IA, IS, IE, JA, JS, JE, data2d )

    allocate( pc3d(KA,IA,JA) )
    !$acc enter data create(pc3d)
    pc3d(:,:,:) = data_P
    !$acc update device(pc3d)
    vid = 2
    call COMM_vars_init( 'testdata3d', pc3d, vid, gid=gid )
    call COMM_vars( pc3d, vid, gid=gid )
    call COMM_wait( pc3d, vid, .false., gid=gid )
    !$acc update host(pc3d)
    call check_vars3d( KA, IA, IS, IE, JA, JS, JE, pc3d )

    call COMM_vars_finalize( gid )
    !$acc exit data delete(pc3d)
    deallocate( pc3d )

  end subroutine test_vars

  subroutine test_vars8( KA, IA, IS, IE, JA, JS, JE, gid )
    use scale_prc, only: &
         PRC_myrank
    use scale_comm_cartesC, only: &
         COMM_vars8_init, &
         COMM_vars_finalize, &
         COMM_vars8, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid

    integer, parameter :: nmax = 2
!    integer, parameter :: nmax = 5
    real(RP) :: data3d(KA,IA,JA,nmax), data2d(IA,JA,nmax)
    real(RP), allocatable :: pc4d(:,:,:,:)
    integer :: vid(nmax)
    integer :: n

    character(len=22) :: title

    write(*,*) "Test vars8"

    write(title(12:22),'(a3,i2.2,a5,i1)') " @ ", PRC_myrank, " gid=", gid

    ! MPI
    data2d(:,:,:) = data_P
    do n = 1, nmax
       write(title(1:11),'(a9,i2.2)') "MPI1(2D)-", n
       call COMM_vars8( data2d(:,:,n), 1, gid=gid )
       call COMM_wait ( data2d(:,:,n), 1, .false., gid=gid )
       call check_vars82d( IA, IS, IE, JA, JS, JE, data2d(:,:,n), title )
    end do

    data3d(:,:,:,:) = data_P
    do n = 1, nmax
       write(title(1:11),'(a9,i2.2)') "MPI1(3D)-", n
       call COMM_vars8( data3d(:,:,:,n), 1, gid=gid )
       call COMM_wait ( data3d(:,:,:,n), 1, .false., gid=gid )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, data3d(:,:,:,n), title )
    end do

    data2d(:,:,:) = data_P
    do n = 1, nmax
       call COMM_vars8( data2d(:,:,n), n, gid=gid )
    end do
    do n = 1, nmax
       write(title(1:11),'(a9,i2.2)') "MPI2(2D)-", n
       call COMM_wait ( data2d(:,:,n), n, .false., gid=gid )
       call check_vars82d( IA, IS, IE, JA, JS, JE, data2d(:,:,n), title )
    end do

    data3d(:,:,:,:) = data_P
    do n = 1, nmax
       call COMM_vars8( data3d(:,:,:,n), n, gid=gid )
    end do
    do n = 1, nmax
       write(title(1:11),'(a9,i2.2)') "MPI2(3D)-", n
       call COMM_wait ( data3d(:,:,:,n), n, .false., gid=gid )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, data3d(:,:,:,n), title )
    end do

    ! MPI PC
!!$    data2d(:,:,:) = data_P
!!$    vid(:) = 1
!!$    do n = 1, nmax
!!$       call COMM_vars8_init( 'PC1-testdata2d', data2d(:,:,n), vid(n), gid=gid )
!!$    end do
!!$    do n = 1, nmax
!!$       write(title(1:11),'(a8,i2.2,a1)') "PC1(2D)-",n," "
!!$       call COMM_vars8( data2d(:,:,n), vid(n), gid=gid )
!!$       call COMM_wait ( data2d(:,:,n), vid(n), .false., gid=gid )
!!$       call check_vars82d( IA, IS, IE, JA, JS, JE, data2d(:,:,n), title )
!!$    end do

    allocate( pc4d(KA,IA,JA,nmax) )
    !$acc enter data create(pc4d)

    pc4d(:,:,:,:) = data_P
    !$acc update device(pc4d)
    vid(:) = 1
    do n = 1, nmax
       call COMM_vars8_init( 'PC1-testdata', pc4d(:,:,:,n), vid(n), gid=gid )
    end do
    do n = 1, nmax
       write(title(1:11),'(a8,i2.2,a1)') "PC1(3D)-",n," "
       call COMM_vars8( pc4d(:,:,:,n), vid(n), gid=gid )
       call COMM_wait ( pc4d(:,:,:,n), vid(n), .false., gid=gid )
       !$acc update host(pc4d(:,:,:,n))
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, pc4d(:,:,:,n), title )
    end do

!!$    data2d(:,:,:) = data_P
!!$    do n = 1, nmax
!!$       vid(n) = n
!!$       call COMM_vars8_init( 'PC2-testdata2d', data2d(:,:,n), vid(n), gid=gid )
!!$    end do
!!$    do n = 1, nmax
!!$       call COMM_vars8( data2d(:,:,n), vid(n), gid=gid )
!!$    end do
!!$    do n = 1, nmax
!!$       write(title(1:11),'(a8,i02,a1)') "PC2(2D)-",n," "
!!$       call COMM_wait ( data2d(:,:,n), vid(n), .false., gid=gid )
!!$       call check_vars82d( IA, IS, IE, JA, JS, JE, data2d(:,:,n), title )
!!$    end do

    pc4d(:,:,:,:) = data_P
    !$acc update device(pc4d)
    do n = 1, nmax
       vid(n) = n
       call COMM_vars8_init( 'PC2-testdata3d', pc4d(:,:,:,n), vid(n), gid=gid )
    end do
    do n = 1, nmax
       call COMM_vars8( pc4d(:,:,:,n), vid(n), gid=gid )
    end do
    do n = 1, nmax
       write(title(1:11),'(a8,i02,a1)') "PC2(3D)-",n," "
       call COMM_wait ( pc4d(:,:,:,n), vid(n), .false., gid=gid )
       !$acc update host(pc4d(:,:,:,n))
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, pc4d(:,:,:,n), title )
    end do

    call COMM_vars_finalize( gid )
    !$acc exit data delete(pc4d)
    deallocate( pc4d )

  end subroutine test_vars8

  !> The post order of the exchanges of the tests below: only the tags and the communicators
  !> can separate the messages when the neighbors post in the opposite order. The order is
  !> flipped by the position of the process in the x (dir = 1) or y (dir = 2) direction, since
  !> the neighbors in the other direction post in the same order.
  function post_first( dir )
    use scale_prc, only: &
         PRC_myrank
    use scale_prc_cartesC, only: &
         PRC_2Drank
    implicit none
    integer, intent(in) :: dir
    logical :: post_first

    post_first = mod( PRC_2Drank(PRC_myrank,dir), 2 ) == 0

  end function post_first

  !> Exchange the same vid on two grids of the same shape at the same time,
  !> posted in a rank-dependent order: only the communicators of the grids can separate
  !> the messages, since the tags do not include the grid.
  subroutine test_vars8_two_grids( KA, IA, IS, IE, JA, JS, JE, gid1, gid2 )
    use scale_comm_cartesC, only: &
         COMM_vars8, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid1, gid2

    real(RP), parameter :: offset = 100.0_RP
    integer,  parameter :: vid = 1
    real(RP) :: a2d(IA,JA),    b2d(IA,JA)
    real(RP) :: a3d(KA,IA,JA), b3d(KA,IA,JA)
    integer :: dir

    write(*,*) "Test vars8 on two grids"

    do dir = 1, 2

       ! control: the former tags of vars8_2D_mpi included the grid ID, so this did not fail with them
       a2d(:,:) = data_P
       b2d(:,:) = data_P + offset
       if ( post_first(dir) ) then
          call COMM_vars8( a2d, vid, gid=gid1 )
          call COMM_vars8( b2d, vid, gid=gid2 )
       else
          call COMM_vars8( b2d, vid, gid=gid2 )
          call COMM_vars8( a2d, vid, gid=gid1 )
       end if
       call COMM_wait( a2d, vid, .false., gid=gid1 )
       call COMM_wait( b2d, vid, .false., gid=gid2 )
       call check_vars82d( IA, IS, IE, JA, JS, JE, a2d,        "TwoGrids(2D)-1 "//dir_name(dir) )
       call check_vars82d( IA, IS, IE, JA, JS, JE, b2d-offset, "TwoGrids(2D)-2 "//dir_name(dir) )

       ! regression test of #3385: the former tags of vars8_3D_mpi did not include the grid ID
       a3d(:,:,:) = data_P
       b3d(:,:,:) = data_P + offset
       if ( post_first(dir) ) then
          call COMM_vars8( a3d, vid, gid=gid1 )
          call COMM_vars8( b3d, vid, gid=gid2 )
       else
          call COMM_vars8( b3d, vid, gid=gid2 )
          call COMM_vars8( a3d, vid, gid=gid1 )
       end if
       call COMM_wait( a3d, vid, .false., gid=gid1 )
       call COMM_wait( b3d, vid, .false., gid=gid2 )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, a3d,        "TwoGrids(3D)-1 "//dir_name(dir) )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, b3d-offset, "TwoGrids(3D)-2 "//dir_name(dir) )

    end do

  end subroutine test_vars8_two_grids

  !> Exchange the same vid in 4 directions on two grids of the same shape at the same time,
  !> posted in a rank-dependent order: only the communicators of the grids can separate the
  !> messages. This guards their separation; the tags with the grid ID also separated these
  !> messages, so this is not a regression test of a former defect.
  subroutine test_vars_two_grids( KA, IA, IS, IE, JA, JS, JE, gid1, gid2 )
    use scale_comm_cartesC, only: &
         COMM_vars, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid1, gid2

    real(RP), parameter :: offset = 100.0_RP
    integer,  parameter :: vid = 1
    real(RP) :: a2d(IA,JA),    b2d(IA,JA)
    real(RP) :: a3d(KA,IA,JA), b3d(KA,IA,JA)
    integer :: dir

    write(*,*) "Test vars on two grids"

    do dir = 1, 2

       a2d(:,:) = data_P
       b2d(:,:) = data_P + offset
       if ( post_first(dir) ) then
          call COMM_vars( a2d, vid, gid=gid1 )
          call COMM_vars( b2d, vid, gid=gid2 )
       else
          call COMM_vars( b2d, vid, gid=gid2 )
          call COMM_vars( a2d, vid, gid=gid1 )
       end if
       call COMM_wait( a2d, vid, .false., gid=gid1 )
       call COMM_wait( b2d, vid, .false., gid=gid2 )
       write(*,*) "TwoGrids4(2D)-1 ", dir_name(dir)
       call check_vars2d( IA, IS, IE, JA, JS, JE, a2d )
       write(*,*) "TwoGrids4(2D)-2 ", dir_name(dir)
       call check_vars2d( IA, IS, IE, JA, JS, JE, b2d-offset )

       a3d(:,:,:) = data_P
       b3d(:,:,:) = data_P + offset
       if ( post_first(dir) ) then
          call COMM_vars( a3d, vid, gid=gid1 )
          call COMM_vars( b3d, vid, gid=gid2 )
       else
          call COMM_vars( b3d, vid, gid=gid2 )
          call COMM_vars( a3d, vid, gid=gid1 )
       end if
       call COMM_wait( a3d, vid, .false., gid=gid1 )
       call COMM_wait( b3d, vid, .false., gid=gid2 )
       write(*,*) "TwoGrids4(3D)-1 ", dir_name(dir)
       call check_vars3d( KA, IA, IS, IE, JA, JS, JE, a3d )
       write(*,*) "TwoGrids4(3D)-2 ", dir_name(dir)
       call check_vars3d( KA, IA, IS, IE, JA, JS, JE, b3d-offset )

    end do

  end subroutine test_vars_two_grids

  !> Exchange persistent variables with the same sequential number on two grids at the
  !> same time, posted in a rank-dependent order: only the communicators of the grids can
  !> separate the messages. This guards their separation; the tags with the grid ID also
  !> separated these messages, so this is not a regression test of a former defect.
  !> The two grids must have the same number of persistent variables.
  subroutine test_vars8_pc_two_grids( KA, IA, IS, IE, JA, JS, JE, gid1, gid2 )
    use dc_test, only: &
         AssertEqual
    use scale_comm_cartesC, only: &
         COMM_vars8_init, &
         COMM_vars_finalize, &
         COMM_vars8, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid1, gid2

    real(RP), parameter :: offset = 100.0_RP
    real(RP), allocatable :: a3d(:,:,:)
    real(RP), allocatable :: b3d(:,:,:)
    integer :: vid1, vid2
    integer :: dir

    write(*,*) "Test vars8 of persistent variables on two grids"

    allocate( a3d(KA,IA,JA) )
    allocate( b3d(KA,IA,JA) )
    !$acc enter data create(a3d, b3d)

    vid1 = 1
    vid2 = 1
    call COMM_vars8_init( 'PC-two-grids', a3d, vid1, gid=gid1 )
    call COMM_vars8_init( 'PC-two-grids', b3d, vid2, gid=gid2 )
    if ( .not. ( pc_assigned( "TwoGridsPC-1", vid1, gid1 ) .and. pc_assigned( "TwoGridsPC-2", vid2, gid2 ) ) ) then
       write(*,*) "Skip: persistent communication is not used"
    else
       ! the returned ID is the sequential number + COMM_vsize_max
       call AssertEqual( "TwoGridsPC: the same sequential number", vid1, vid2 )

       do dir = 1, 2
          a3d(:,:,:) = data_P
          b3d(:,:,:) = data_P + offset
          !$acc update device(a3d, b3d)
          if ( post_first(dir) ) then
             call COMM_vars8( a3d, vid1, gid=gid1 )
             call COMM_vars8( b3d, vid2, gid=gid2 )
          else
             call COMM_vars8( b3d, vid2, gid=gid2 )
             call COMM_vars8( a3d, vid1, gid=gid1 )
          end if
          call COMM_wait( a3d, vid1, .false., gid=gid1 )
          call COMM_wait( b3d, vid2, .false., gid=gid2 )
          !$acc update host(a3d, b3d)
          call check_vars83d( KA, IA, IS, IE, JA, JS, JE, a3d,        "TwoGridsPC(3D)-1 "//dir_name(dir) )
          call check_vars83d( KA, IA, IS, IE, JA, JS, JE, b3d-offset, "TwoGridsPC(3D)-2 "//dir_name(dir) )
       end do
    end if

    call COMM_vars_finalize( gid1 )
    call COMM_vars_finalize( gid2 )
    !$acc exit data delete(a3d, b3d)
    deallocate( a3d )
    deallocate( b3d )

  end subroutine test_vars8_pc_two_grids

  !> Valid uses of the IDs, which must pass the checks of the packing buffers and of the arrays
  !> given with persistent IDs (the invalid uses are run by unit_comm_check): the same vid on two
  !> grids and two vids on a grid at the same time, a vid again after its wait, and persistent
  !> variables with their registered arrays and grids
  subroutine test_vars8_valid_ids( KA, IA, IS, IE, JA, JS, JE, gid1, gid2 )
    use scale_comm_cartesC, only: &
         COMM_vars8_init, &
         COMM_vars_finalize, &
         COMM_vars8, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid1, gid2

    real(RP), allocatable :: a3d(:,:,:)
    real(RP), allocatable :: b3d(:,:,:)
    real(RP) :: c3d(KA,IA,JA) ! not registered: it is exchanged on the host in the OpenACC build
    integer :: vid1, vid2

    write(*,*) "Test vars8 of valid uses of the IDs"

    allocate( a3d(KA,IA,JA) )
    allocate( b3d(KA,IA,JA) )
    !$acc enter data create(a3d, b3d)

    vid1 = 1
    call COMM_vars8_init( 'valid-ids-a', a3d, vid1, gid=gid1 )
    vid2 = 1
    call COMM_vars8_init( 'valid-ids-b', b3d, vid2, gid=gid2 )

    ! the same vid on two grids, and two vids on a grid, at the same time
    a3d(:,:,:) = data_P
    b3d(:,:,:) = data_P
    c3d(:,:,:) = data_P
    !$acc update device(a3d, b3d)
    call COMM_vars8( a3d, 1, gid=gid1 )
    call COMM_vars8( b3d, 1, gid=gid2 )
    call COMM_vars8( c3d, 2, gid=gid1 )
    call COMM_wait ( a3d, 1, .false., gid=gid1 )
    call COMM_wait ( b3d, 1, .false., gid=gid2 )
    call COMM_wait ( c3d, 2, .false., gid=gid1 )
    !$acc update host(a3d, b3d)
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, a3d, "ValidIDs(3D)-1" )
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, b3d, "ValidIDs(3D)-2" )
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, c3d, "ValidIDs(3D)-3" )

    ! the same vid again after its wait
    c3d(:,:,:) = data_P
    call COMM_vars8( c3d, 2, gid=gid1 )
    call COMM_wait ( c3d, 2, .false., gid=gid1 )
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, c3d, "ValidIDs(3D)-4" )

    ! the persistent variables with their registered arrays and grids (if IDs were assigned)
    a3d(:,:,:) = data_P
    b3d(:,:,:) = data_P
    !$acc update device(a3d, b3d)
    call COMM_vars8( a3d, vid1, gid=gid1 )
    call COMM_vars8( b3d, vid2, gid=gid2 )
    call COMM_wait ( a3d, vid1, .false., gid=gid1 )
    call COMM_wait ( b3d, vid2, .false., gid=gid2 )
    !$acc update host(a3d, b3d)
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, a3d, "ValidIDs(3D)-5" )
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, b3d, "ValidIDs(3D)-6" )

    call COMM_vars_finalize( gid1 )
    call COMM_vars_finalize( gid2 )
    !$acc exit data delete(a3d, b3d)
    deallocate( a3d )
    deallocate( b3d )

  end subroutine test_vars8_valid_ids

  !> Exchange a persistent variable (8-direction) and a non-persistent one on the same grid
  !> at the same time, posted in a rank-dependent order. The ID of the non-persistent one
  !> equals the sequential number of the persistent one, so their tags are the same: only
  !> their communicators can separate the messages (regression test of #3385).
  subroutine test_vars8_pc_same_grid( KA, IA, IS, IE, JA, JS, JE, gid )
    use dc_test, only: &
         AssertLessThan
    use scale_comm_cartesC, only: &
         COMM_vars8_init, &
         COMM_vars_finalize, &
         COMM_vars8, &
         COMM_wait, &
         COMM_get_vsize_max, &
         COMM_get_pc_capacity
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid

    real(RP), parameter :: offset = 100.0_RP
    real(RP), allocatable :: a3d(:,:,:)
    real(RP) :: b3d(KA,IA,JA)
    integer :: vsize_max
    integer :: vid1, vid2
    integer :: dir

    write(*,*) "Test vars8 of a persistent and a non-persistent variable on a grid"

    if ( COMM_get_pc_capacity( gid ) == 0 ) then
       write(*,*) "Skip: persistent communication is not used"
       return
    end if
    vsize_max = COMM_get_vsize_max()

    ! b3d is not registered: it is exchanged on the host in the OpenACC build
    allocate( a3d(KA,IA,JA) )
    !$acc enter data create(a3d)

    ! the ID given here selects the packing buffer, which must differ from that of vid2
    vid1 = vsize_max
    call COMM_vars8_init( 'PC-same-grid', a3d, vid1, gid=gid )
    if ( pc_assigned( "SameGridPC", vid1, gid ) ) then
       vid2 = vid1 - vsize_max ! the sequential number
       call AssertLessThan( "SameGridPC: the sequential number is below COMM_vsize_max", vsize_max, vid2 )

       do dir = 1, 2
          a3d(:,:,:) = data_P
          b3d(:,:,:) = data_P + offset
          !$acc update device(a3d)
          if ( post_first(dir) ) then
             call COMM_vars8( a3d, vid1, gid=gid )
             call COMM_vars8( b3d, vid2, gid=gid )
          else
             call COMM_vars8( b3d, vid2, gid=gid )
             call COMM_vars8( a3d, vid1, gid=gid )
          end if
          call COMM_wait( a3d, vid1, .false., gid=gid )
          call COMM_wait( b3d, vid2, .false., gid=gid )
          !$acc update host(a3d)
          call check_vars83d( KA, IA, IS, IE, JA, JS, JE, a3d,        "SameGridPC(3D)-1 "//dir_name(dir) )
          call check_vars83d( KA, IA, IS, IE, JA, JS, JE, b3d-offset, "SameGridPC(3D)-2 "//dir_name(dir) )
       end do
    end if

    call COMM_vars_finalize( gid )
    !$acc exit data delete(a3d)
    deallocate( a3d )

  end subroutine test_vars8_pc_same_grid

  !> Exchange a persistent variable (4-direction) and a non-persistent one on the same grid
  !> at the same time, posted in a rank-dependent order. The sequential number of the
  !> persistent variable and the ID of the non-persistent one are both 1, so their tags are
  !> the same: only their communicators can separate the messages (regression test of #3385).
  !> The grid must have no persistent variables.
  subroutine test_vars_pc_same_grid( KA, IA, IS, IE, JA, JS, JE, gid )
    use scale_comm_cartesC, only: &
         COMM_vars_init, &
         COMM_vars_finalize, &
         COMM_vars, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid

    real(RP), parameter :: offset = 100.0_RP
    integer,  parameter :: vid2 = 1
    real(RP), allocatable :: a3d(:,:,:)
    real(RP) :: b3d(KA,IA,JA)
    integer :: vid1
    integer :: dir

    write(*,*) "Test vars of a persistent and a non-persistent variable on a grid"

    ! b3d is not registered: it is exchanged on the host in the OpenACC build
    allocate( a3d(KA,IA,JA) )
    !$acc enter data create(a3d)

    ! the ID given here selects the packing buffer, which must differ from that of vid2
    vid1 = 2
    call COMM_vars_init( 'PC-same-grid', a3d, vid1, gid=gid )
    if ( .not. pc_assigned( "SameGridPC4", vid1, gid ) ) then
       write(*,*) "Skip: persistent communication is not used"
    else
       do dir = 1, 2
          a3d(:,:,:) = data_P
          b3d(:,:,:) = data_P + offset
          !$acc update device(a3d)
          if ( post_first(dir) ) then
             call COMM_vars( a3d, vid1, gid=gid )
             call COMM_vars( b3d, vid2, gid=gid )
          else
             call COMM_vars( b3d, vid2, gid=gid )
             call COMM_vars( a3d, vid1, gid=gid )
          end if
          call COMM_wait( a3d, vid1, .false., gid=gid )
          call COMM_wait( b3d, vid2, .false., gid=gid )
          !$acc update host(a3d)
          write(*,*) "SameGridPC4(3D)-1 ", dir_name(dir)
          call check_vars3d( KA, IA, IS, IE, JA, JS, JE, a3d )
          write(*,*) "SameGridPC4(3D)-2 ", dir_name(dir)
          call check_vars3d( KA, IA, IS, IE, JA, JS, JE, b3d-offset )
       end do
    end if

    call COMM_vars_finalize( gid )
    !$acc exit data delete(a3d)
    deallocate( a3d )

  end subroutine test_vars_pc_same_grid

  !> Register more persistent variables than the initial capacity (COMM_vsize_max_pc), so that
  !> the arrays for them are extended, and exchange the first and the last one.
  !> The grid must have no persistent variables.
  subroutine test_vars8_pc_grow( KA, IA, IS, IE, JA, JS, JE, gid )
    use dc_test, only: &
         AssertLessThan
    use scale_comm_cartesC, only: &
         COMM_vars8_init, &
         COMM_vars_finalize, &
         COMM_vars8, &
         COMM_wait, &
         COMM_get_pc_capacity
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid

    real(RP), allocatable :: data3d(:,:,:,:)
    integer, allocatable :: vid(:)
    integer :: capacity
    integer :: npc
    integer :: n

    write(*,*) "Test vars8 of more persistent variables than the initial capacity"

    ! more than the current (initial) capacity of the grid, 75 with run.conf
    capacity = COMM_get_pc_capacity( gid )
    npc = capacity + 5
    allocate( vid(npc) )
    allocate( data3d(KA,IA,JA,npc) )
    !$acc enter data create(data3d)

    vid(:) = 1
    do n = 1, npc
       call COMM_vars8_init( 'PC-grow', data3d(:,:,:,n), vid(n), gid=gid )
    end do
    if ( .not. pc_assigned( "GrowPC", vid(1), gid ) ) then
       write(*,*) "Skip: persistent communication is not used"
    else
       call AssertLessThan( "GrowPC: the capacity is extended", COMM_get_pc_capacity( gid ), capacity )
       data3d(:,:,:,:) = data_P
       !$acc update device(data3d)
       call COMM_vars8( data3d(:,:,:,1), vid(1), gid=gid )
       call COMM_wait ( data3d(:,:,:,1), vid(1), .false., gid=gid )
       call COMM_vars8( data3d(:,:,:,npc), vid(npc), gid=gid )
       call COMM_wait ( data3d(:,:,:,npc), vid(npc), .false., gid=gid )
       !$acc update host(data3d)
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, data3d(:,:,:,1),   "GrowPC(3D)-first" )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, data3d(:,:,:,npc), "GrowPC(3D)-last " )
    end if

    call COMM_vars_finalize( gid )
    !$acc exit data delete(data3d)
    deallocate( data3d )
    deallocate( vid )

  end subroutine test_vars8_pc_grow

  !> Release the persistent variables of a grid (COMM_vars_finalize, twice: the second call does
  !> nothing), free their arrays, and register new arrays: the registration starts again from the
  !> first ID and the capacity of the grid is kept. The new persistent variables and a
  !> non-persistent one are exchanged after the release.
  !> The grid must have no persistent variables.
  subroutine test_vars8_pc_reregister( KA, IA, IS, IE, JA, JS, JE, gid )
    use dc_test, only: &
         AssertEqual
    use scale_comm_cartesC, only: &
         COMM_vars_init, &
         COMM_vars8_init, &
         COMM_vars_finalize, &
         COMM_vars, &
         COMM_vars8, &
         COMM_wait, &
         COMM_get_pc_capacity
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid

    real(RP), allocatable :: a3d(:,:,:)
    real(RP), allocatable :: b3d(:,:,:)
    real(RP) :: c3d(KA,IA,JA) ! not registered: it is exchanged on the host in the OpenACC build
    integer :: vid1, vid2
    integer :: vid1_first, vid2_first
    integer :: capacity
    logical :: assigned

    write(*,*) "Test vars8 of persistent variables registered again after their release"

    ! first registration
    allocate( a3d(KA,IA,JA) )
    allocate( b3d(KA,IA,JA) )
    !$acc enter data create(a3d, b3d)

    vid1 = 1
    vid2 = 2
    call COMM_vars8_init( 'PC-reregist-a', a3d, vid1, gid=gid )
    call COMM_vars8_init( 'PC-reregist-b', b3d, vid2, gid=gid )
    assigned = pc_assigned( "ReregistPC-1", vid1, gid ) .and. pc_assigned( "ReregistPC-2", vid2, gid )
    vid1_first = vid1
    vid2_first = vid2
    capacity = COMM_get_pc_capacity( gid )

    a3d(:,:,:) = data_P
    b3d(:,:,:) = data_P
    !$acc update device(a3d, b3d)
    call COMM_vars8( a3d, vid1, gid=gid )
    call COMM_vars8( b3d, vid2, gid=gid )
    call COMM_wait ( a3d, vid1, .false., gid=gid )
    call COMM_wait ( b3d, vid2, .false., gid=gid )
    !$acc update host(a3d, b3d)
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, a3d, "ReregistPC(3D)-1" )
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, b3d, "ReregistPC(3D)-2" )

    ! release (the second call does nothing) and free the arrays
    call COMM_vars_finalize( gid )
    call COMM_vars_finalize( gid )
    !$acc exit data delete(a3d, b3d)
    deallocate( a3d )
    deallocate( b3d )

    ! register new arrays: a 4-direction one and an 8-direction one
    allocate( a3d(KA,IA,JA) )
    allocate( b3d(KA,IA,JA) )
    !$acc enter data create(a3d, b3d)

    vid1 = 1
    vid2 = 2
    call COMM_vars_init ( 'PC-reregist-a', a3d, vid1, gid=gid )
    call COMM_vars8_init( 'PC-reregist-b', b3d, vid2, gid=gid )
    if ( assigned ) then
       call AssertEqual( "ReregistPC: the first ID again",       vid1_first, vid1 )
       call AssertEqual( "ReregistPC: the second ID again",      vid2_first, vid2 )
       call AssertEqual( "ReregistPC: the capacity is kept",     capacity,   COMM_get_pc_capacity( gid ) )
    else
       write(*,*) "Skip the checks of the IDs: persistent communication is not used"
    end if

    a3d(:,:,:) = data_P
    b3d(:,:,:) = data_P
    c3d(:,:,:) = data_P
    !$acc update device(a3d, b3d)
    call COMM_vars ( a3d, vid1, gid=gid )
    call COMM_vars8( b3d, vid2, gid=gid )
    call COMM_vars8( c3d, 3,    gid=gid )
    call COMM_wait ( a3d, vid1, .false., gid=gid )
    call COMM_wait ( b3d, vid2, .false., gid=gid )
    call COMM_wait ( c3d, 3,    .false., gid=gid )
    !$acc update host(a3d, b3d)
    write(*,*) "ReregistPC(3D)-3"
    call check_vars3d( KA, IA, IS, IE, JA, JS, JE, a3d )
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, b3d, "ReregistPC(3D)-4" )
    call check_vars83d( KA, IA, IS, IE, JA, JS, JE, c3d, "ReregistPC(3D)-5" )

    call COMM_vars_finalize( gid )
    !$acc exit data delete(a3d, b3d)
    deallocate( a3d )
    deallocate( b3d )

  end subroutine test_vars8_pc_reregister

  !> Whether the ID vid returned by the registration on the grid gid is a persistent ID: .false. if
  !> persistent communication is not used (the test is skipped); the test stops if persistent
  !> communication is used but no persistent ID is assigned
  function pc_assigned( name, vid, gid ) result(assigned)
    use dc_test, only: &
         AssertEqual
    use scale_comm_cartesC, only: &
         COMM_get_vsize_max, &
         COMM_get_pc_capacity
    implicit none
    character(len=*), intent(in) :: name
    integer,          intent(in) :: vid
    integer,          intent(in) :: gid
    logical :: assigned

    assigned = .false.
    if ( COMM_get_pc_capacity( gid ) == 0 ) return

    call AssertEqual( name//": a persistent ID is assigned (persistent communication is used)", &
                      .true., vid > COMM_get_vsize_max() )
    assigned = .true.

  end function pc_assigned

  subroutine check_vars2d( &
       IA, IS, IE, JA, JS, JE, &
       data )
    use dc_test, only: &
         AssertEqual
    implicit none
    integer, intent(in) :: IA, IS, IE, JA, JS, JE
    real(RP), intent(in) :: data(IA,JA)

    real(RP) :: expect(IA,JA)

    integer :: i, j

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(i,j) = data_N
    end do
    end do
    call AssertEqual("North(2D)", expect(IS:IE,JE+1:JA), data(IS:IE,JE+1:JA))

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(i,j) = data_S
    end do
    end do
    call AssertEqual("South(2D)", expect(IS:IE,1:JS-1), data(IS:IE,1:JS-1))

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(i,j) = data_E
    end do
    end do
    call AssertEqual("East(2D)", expect(IE+1:IA,JS:JE), data(IE+1:IA,JS:JE))

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(i,j) = data_W
    end do
    end do
    call AssertEqual("West(2D)", expect(1:IS-1,JS:JE), data(1:IS-1,JS:JE))

  end subroutine check_vars2d

  subroutine check_vars3d( &
       KA, IA, IS, IE, JA, JS, JE, &
       data )
    use dc_test, only: &
         AssertEqual
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    real(RP), intent(in) :: data(KA,IA,JA)

    real(RP) :: expect(KA,IA,JA)

    integer :: i, j

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(:,i,j) = data_N
    end do
    end do
    call AssertEqual("North(3D)", expect(:,IS:IE,JE+1:JA), data(:,IS:IE,JE+1:JA))

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(:,i,j) = data_S
    end do
    end do
    call AssertEqual("South(3D)", expect(:,IS:IE,1:JS-1), data(:,IS:IE,1:JS-1))

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(:,i,j) = data_E
    end do
    end do
    call AssertEqual("East(3D)", expect(:,IE+1:IA,JS:JE), data(:,IE+1:IA,JS:JE))

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(:,i,j) = data_W
    end do
    end do
    call AssertEqual("West(3D)", expect(:,1:IS-1,JS:JE), data(:,1:IS-1,JS:JE))

  end subroutine check_vars3d

  subroutine check_vars82d( &
       IA, IS, IE, JA, JS, JE, &
       data, title )
    use dc_test, only: &
         AssertEqual
    implicit none
    integer, intent(in) :: IA, IS, IE, JA, JS, JE
    real(RP), intent(in) :: data(IA,JA)
    character(len=*), intent(in) :: title

    real(RP) :: expect(IA,JA)

    integer :: i, j

    do j = JS, JE
    do i = IS, IE
       expect(i,j) = data_P
    end do
    end do

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(i,j) = data_N
    end do
    end do

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(i,j) = data_S
    end do
    end do

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(i,j) = data_E
    end do
    end do

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(i,j) = data_W
    end do
    end do

    ! North-East
    do j = JE+1, JA
    do i = IE+1, IA
       expect(i,j) = data_NE
    end do
    end do

    ! North-West
    do j = JE+1, JA
    do i = 1, IS-1
       expect(i,j) = data_NW
    end do
    end do

    ! South-East
    do j = 1, JS-1
    do i = IE+1, IA
       expect(i,j) = data_SE
    end do
    end do

    ! South-West
    do j = 1, JS-1
    do i = 1, IS-1
       expect(i,j) = data_SW
    end do
    end do

    call AssertEqual(title, expect, data)

  end subroutine check_vars82d

  subroutine check_vars83d( &
       KA, IA, IS, IE, JA, JS, JE, &
       data, title )
    use dc_test, only: &
         AssertEqual
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    real(RP), intent(in) :: data(KA,IA,JA)
    character(len=*), intent(in) :: title

    real(RP) :: expect(KA,IA,JA)

    integer :: i, j

    do j = JS, JE
    do i = IS, IE
       expect(:,i,j) = data_P
    end do
    end do

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(:,i,j) = data_N
    end do
    end do

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(:,i,j) = data_S
    end do
    end do

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(:,i,j) = data_E
    end do
    end do

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(:,i,j) = data_W
    end do
    end do

    ! North-East
    do j = JE+1, JA
    do i = IE+1, IA
       expect(:,i,j) = data_NE
    end do
    end do

    ! North-West
    do j = JE+1, JA
    do i = 1, IS-1
       expect(:,i,j) = data_NW
    end do
    end do

    ! South-East
    do j = 1, JS-1
    do i = IE+1, IA
       expect(:,i,j) = data_SE
    end do
    end do

    ! South-West
    do j = 1, JS-1
    do i = 1, IS-1
       expect(:,i,j) = data_SW
    end do
    end do

    call AssertEqual(title, expect, data)

  end subroutine check_vars83d

end module test_comm




