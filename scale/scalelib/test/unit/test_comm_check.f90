!-------------------------------------------------------------------------------
!> Standalone test of the checks of the variable ID in COMM_cartesC (issue #3426)
!!
!! The checks abort the run, so each mode is run as a separate MPI job. The mode is the
!! second argument:
!!   packbuf       : the same vid is posted twice before its wait (the packing buffer is in use)
!!   wait_nopost   : a wait without a post
!!   wait_other    : a wait with a non-persistent vid for the exchange of a persistent variable
!!                   registered with that vid (the packing buffer is in use by another ID)
!!   other_grid    : a persistent ID is given with another grid that has as many persistent variables
!!   foreign_array : a persistent ID is given with an array other than the registered one
!!   post_nowait   : an exchange is posted and never waited before COMM_finalize
!!   wait_rank     : an exchange of a 3D array is waited with a 2D array (the same vid)
!!   wait_array    : an exchange of a 3D array is waited with another 3D array (the same vid)
!!   pc_same_vid   : two persistent variables registered with the same vid on a grid (they share the
!!                   packing buffer) are exchanged at the same time
!!   pc_dir        : a persistent variable registered for the 4-direction exchange (COMM_vars_init)
!!                   is posted with COMM_vars8
!!   pc_released   : a persistent ID is used after COMM_vars_finalize released the variables of its grid
!!   pc_rel_nowait : COMM_vars_finalize is called while a persistent exchange is posted and not waited
!! The groups of persistent variables (COMM_vars_group_init, issue #3594), with a group of two
!! variables on grid 1 (vid1 and vid3):
!!   grp_wait_mem  : the member posted is waited before the other member is posted
!!   grp_repost    : a member waited is posted again before the other member is waited
!!   grp_rel_round : COMM_vars_finalize is called in a round of the group (a member not waited)
!!   grp_in_flight : the group is declared while a member is posted and not waited
!!   grp_same_pid  : the members share the packing buffer (the same vid at the registration)
!!   grp_dup       : a variable is declared in two groups
!!   grp_4dir      : a member is registered for the 4-direction exchange
!!   grp_range     : the ID of a member is not that of a persistent variable of the grid
!!   grp_rank_n    : the number of the members differs among the processes
!!   grp_rank_id   : the IDs of the members differ among the processes (their order)
!!   batch_nol2    : COMM_USE_BATCH without COMM_USE_PC_LAYOUT2 (COMM_setup aborts; the mode does
!!                   nothing itself): comm_check_batch_nolayout2.conf, and comm_check_batch_nopc.conf
!!                   without persistent communication
!! The configuration file (the first argument) gives the process layout and the communication mode.
!! The valid uses of the IDs, which must pass the checks, are tested in test_comm (./unit).
!! wait_other, other_grid, foreign_array, pc_same_vid, pc_dir, pc_released and pc_rel_nowait need
!! persistent IDs and must be run with persistent communication (comm_check_pc.conf); without a
!! persistent ID they stop with a distinct error.
!! The modes of the groups need persistent IDs too; the runtime checks (grp_wait_mem, grp_repost,
!! grp_rel_round, grp_in_flight) are made in both states of the group (checked: comm_check_pc.conf; in batch:
!! comm_check_batch.conf).
!! The checks are done in every build, so this program does not depend on -DDEBUG.
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
program test_comm_check
  use scale
  use scale_prc, only: &
     PRC_abort, &
     PRC_myrank
  use scale_io, only: &
     IO_UNIVERSALRANK, &
     IO_LOCALRANK, &
     IO_JOBID, &
     IO_DOMAINID
  use scale_prc_cartesC, only: &
     PRC_CARTESC_setup, &
     PRC_CARTESC_finalize
  use scale_atmos_grid_cartesC_index
  use scale_comm_cartesC, only: &
     COMM_setup, &
     COMM_regist, &
     COMM_finalize, &
     COMM_vars_init, &
     COMM_vars8_init, &
     COMM_vars_finalize, &
     COMM_vars8, &
     COMM_wait, &
     COMM_vars_group_init
  use dc_test, only: &
     AssertEqual
  implicit none

  real(RP), allocatable :: a3d(:,:,:), b3d(:,:,:), c3d(:,:,:), d3d(:,:,:), e3d(:,:,:)
  real(RP), allocatable :: c2d(:,:)
  character(len=16) :: mode
  integer :: gid1, gid2
  integer :: vid1, vid2, vid3, vid4, vid5
  logical :: grp

  call SCALE_init( "Unit test of the checks of COMM_cartesC" )
  call get_command_argument(2, mode)
  call AssertEqual( 'comm check mode', .true., &
       mode == 'packbuf' .or. mode == 'wait_nopost' .or. mode == 'wait_other' &
       .or. mode == 'other_grid' .or. mode == 'foreign_array' .or. mode == 'post_nowait' &
       .or. mode == 'wait_rank' .or. mode == 'wait_array' .or. mode == 'pc_same_vid' .or. mode == 'pc_dir' &
       .or. mode == 'pc_released' .or. mode == 'pc_rel_nowait' &
       .or. mode == 'grp_wait_mem' .or. mode == 'grp_repost' .or. mode == 'grp_rel_round' .or. mode == 'grp_same_pid' &
       .or. mode == 'grp_dup' .or. mode == 'grp_4dir' .or. mode == 'grp_range' .or. mode == 'grp_rank_n' &
       .or. mode == 'grp_rank_id' .or. mode == 'grp_in_flight' .or. mode == 'batch_nol2' )
  grp = mode(1:4) == 'grp_'

  call PRC_CARTESC_setup
  call ATMOS_GRID_CARTESC_INDEX_setup( KMAX=4, IMAX=6, JMAX=6 )
  call COMM_setup
  call COMM_regist( KA, IA, JA, IHALO, JHALO, gid1 )
  call COMM_regist( KA, IA, JA, IHALO, JHALO, gid2 )

  allocate( a3d(KA,IA,JA), b3d(KA,IA,JA), c3d(KA,IA,JA), d3d(KA,IA,JA), e3d(KA,IA,JA), c2d(IA,JA) )
  c2d(:,:) = 4.0_RP
  a3d(:,:,:) = 1.0_RP
  b3d(:,:,:) = 2.0_RP
  c3d(:,:,:) = 3.0_RP
  d3d(:,:,:) = 5.0_RP
  e3d(:,:,:) = 6.0_RP

  ! the variables must be on the device to be registered in the OpenACC build, and they must
  ! outlive the registration, which lasts until COMM_vars_finalize or COMM_finalize
  !$acc enter data copyin(a3d, b3d, c3d, d3d, e3d)

  ! one persistent variable on each grid: the two IDs are the same
  vid1 = 1
  call COMM_vars8_init( 'check-a', a3d, vid1, gid=gid1 )
  vid2 = 1
  call COMM_vars8_init( 'check-b', b3d, vid2, gid=gid2 )
  LOG_INFO("test_comm_check",*) 'mode = ', trim(mode), ', IDs = ', vid1, vid2

  vid3 = 1
  vid4 = 4
  vid5 = 5
  if ( grp ) then
     ! a group of a3d and c3d on grid 1 (other packing buffers), and d3d outside the group
     vid3 = 3
     call COMM_vars8_init( 'check-c', c3d, vid3, gid=gid1 )
     call COMM_vars8_init( 'check-d', d3d, vid4, gid=gid1 )
     if ( mode == 'grp_same_pid' ) then
        vid5 = 1 ! the packing buffer of a3d
        call COMM_vars8_init( 'check-e', e3d, vid5, gid=gid1 )
     else if ( mode == 'grp_4dir' ) then
        call COMM_vars_init( 'check-e4', e3d, vid5, gid=gid1 )
     end if
     if ( vid1 == 1 .or. vid2 == 1 .or. vid3 == 3 .or. vid4 == 4 &
          .or. ( mode == 'grp_same_pid' .and. vid5 == 1 ) .or. ( mode == 'grp_4dir' .and. vid5 == 5 ) ) then
        LOG_ERROR("test_comm_check",*) 'test error: no persistent ID is assigned', vid1, vid2, vid3, vid4, vid5
        LOG_ERROR_CONT(*) '(this mode needs persistent communication)'
        call PRC_abort
     end if
     select case ( mode )
     case ( 'grp_same_pid', 'grp_4dir' )
        ! must abort here
        call COMM_vars_group_init( 'check-group', (/ vid1, vid5 /), gid=gid1 )
     case ( 'grp_range' )
        ! must abort here: not an ID of a persistent variable of grid 1
        call COMM_vars_group_init( 'check-group', (/ vid1, vid4 + 1 /), gid=gid1 )
     case ( 'grp_rank_n' )
        ! must abort here, in all the processes
        if ( PRC_myrank == 0 ) then
           call COMM_vars_group_init( 'check-group', (/ vid1 /), gid=gid1 )
        else
           call COMM_vars_group_init( 'check-group', (/ vid1, vid3 /), gid=gid1 )
        end if
     case ( 'grp_in_flight' )
        ! must abort here: a3d is posted and not waited
        call COMM_vars8( a3d, vid1, gid=gid1 )
        call COMM_vars_group_init( 'check-group', (/ vid1, vid3 /), gid=gid1 )
        call COMM_wait ( a3d, vid1, .false., gid=gid1 )
     case ( 'grp_rank_id' )
        ! must abort here, in all the processes
        if ( PRC_myrank == 0 ) then
           call COMM_vars_group_init( 'check-group', (/ vid3, vid1 /), gid=gid1 )
        else
           call COMM_vars_group_init( 'check-group', (/ vid1, vid3 /), gid=gid1 )
        end if
     case default
        call COMM_vars_group_init( 'check-group', (/ vid1, vid3 /), gid=gid1 )
     end select
  end if
  if ( mode == 'pc_same_vid' ) then
     ! another persistent variable with the vid of a3d on grid 1
     call COMM_vars8_init( 'check-c', c3d, vid3, gid=gid1 )
  else if ( mode == 'pc_dir' ) then
     ! a persistent variable for the 4-direction exchange on grid 1
     vid3 = 3
     call COMM_vars_init( 'check-c4', c3d, vid3, gid=gid1 )
  end if

  if ( mode == 'wait_other' .or. mode == 'other_grid' .or. mode == 'foreign_array' &
       .or. mode == 'pc_same_vid' .or. mode == 'pc_dir' .or. mode == 'pc_released' .or. mode == 'pc_rel_nowait' ) then
     if ( vid1 == 1 .or. vid2 == 1 .or. ( mode == 'pc_same_vid' .and. vid3 == 1 ) &
          .or. ( mode == 'pc_dir' .and. vid3 == 3 ) ) then
        LOG_ERROR("test_comm_check",*) 'test error: no persistent ID is assigned', vid1, vid2
        LOG_ERROR_CONT(*) '(this mode needs persistent communication)'
        call PRC_abort
     end if
  end if

  select case ( mode )
  case ( 'packbuf' )
     ! must abort at the second post: the packing buffer of vid 1 is in use
     call COMM_vars8( a3d, 1, gid=gid1 )
     call COMM_vars8( c3d, 1, gid=gid1 )
     call COMM_wait ( a3d, 1, .false., gid=gid1 )
     call COMM_wait ( c3d, 1, .false., gid=gid1 )
  case ( 'wait_nopost' )
     ! must abort: nothing is posted with vid 2
     call COMM_wait ( c3d, 2, .false., gid=gid1 )
  case ( 'wait_other' )
     ! must abort at the wait: the packing buffer of vid 1 is taken by the persistent ID vid1
     call COMM_vars8( a3d, vid1, gid=gid1 )
     call COMM_wait ( a3d, 1, .false., gid=gid1 )
  case ( 'other_grid' )
     ! must abort: vid1 is the ID of a3d on grid 1, and grid 2 has as many persistent variables
     call COMM_vars8( a3d, vid1, gid=gid2 )
     call COMM_wait ( a3d, vid1, .false., gid=gid2 )
  case ( 'wait_rank' )
     ! must abort at the wait: the exchange of vid 2 is of a 3D array
     call COMM_vars8( c3d, 2, gid=gid1 )
     call COMM_wait ( c2d, 2, .false., gid=gid1 )
  case ( 'pc_same_vid' )
     ! must abort at the second post: vid1 and vid3 share the packing buffer of vid 1
     call COMM_vars8( a3d, vid1, gid=gid1 )
     call COMM_vars8( c3d, vid3, gid=gid1 )
  case ( 'pc_dir' )
     ! must abort at the post: vid3 is registered for the 4-direction exchange
     call COMM_vars8( c3d, vid3, gid=gid1 )
     call COMM_wait ( c3d, vid3, .false., gid=gid1 )
  case ( 'wait_array' )
     ! must abort at the wait: the exchange of vid 2 is of c3d
     call COMM_vars8( c3d, 2, gid=gid1 )
     call COMM_wait ( b3d, 2, .false., gid=gid1 )
  case ( 'post_nowait' )
     ! must abort at COMM_finalize: the exchange of vid 2 is never waited
     call COMM_vars8( c3d, 2, gid=gid1 )
  case ( 'foreign_array' )
     ! must abort: vid1 is the ID of a3d, not of c3d
     call COMM_vars8( c3d, vid1, gid=gid1 )
     call COMM_wait ( c3d, vid1, .false., gid=gid1 )
  case ( 'pc_released' )
     ! must abort at the post: the persistent variables of grid 1 are released
     call COMM_vars_finalize( gid1 )
     call COMM_vars8( a3d, vid1, gid=gid1 )
     call COMM_wait ( a3d, vid1, .false., gid=gid1 )
  case ( 'pc_rel_nowait' )
     ! must abort at the release: the persistent exchange of vid1 is not waited
     call COMM_vars8( a3d, vid1, gid=gid1 )
     call COMM_vars_finalize( gid1 )
     call COMM_wait ( a3d, vid1, .false., gid=gid1 )
  case ( 'grp_wait_mem' )
     ! must abort at the wait of a3d: c3d is not posted
     call COMM_vars8( a3d, vid1, gid=gid1 )
     call COMM_wait ( a3d, vid1, .false., gid=gid1 )
  case ( 'grp_repost' )
     ! must abort at the second post of a3d: c3d is not waited
     call COMM_vars8( a3d, vid1, gid=gid1 )
     call COMM_vars8( c3d, vid3, gid=gid1 )
     call COMM_wait ( a3d, vid1, .false., gid=gid1 )
     call COMM_vars8( a3d, vid1, gid=gid1 )
  case ( 'grp_rel_round' )
     ! must abort at the release: c3d is not waited
     call COMM_vars8( a3d, vid1, gid=gid1 )
     call COMM_vars8( c3d, vid3, gid=gid1 )
     call COMM_wait ( a3d, vid1, .false., gid=gid1 )
     call COMM_vars_finalize( gid1 )
  case ( 'grp_dup' )
     ! must abort here: c3d is in the group
     call COMM_vars_group_init( 'check-group2', (/ vid4, vid3 /), gid=gid1 )
  end select

  call COMM_finalize
  !$acc exit data delete(a3d, b3d, c3d, d3d, e3d)
  call PRC_CARTESC_finalize
  call SCALE_finalize

end program test_comm_check
