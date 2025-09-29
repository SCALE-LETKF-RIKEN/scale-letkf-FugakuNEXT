program scaleles_pp_ens
  !-----------------------------------------------------------------------------

  use mpi
  use common_mpi, only: &
     nprocs, &
     myrank
  use common_nml
  use common_scale, only: &
     set_common_conf
  use common_mpi_scale, only: &
     myrank_to_mem, &
     myrank_to_pe, &
     nitmax, &
     set_mem_node_proc, &
     mpi_timer

  use scale_io, only: &
     H_LONG, &
     IO_L, &
     IO_FID_CONF, &
     IO_FID_STDOUT, &
     IO_set_universalrank 
  use scale_prc, only: &
     PRC_MPIstart, &
     PRC_UNIVERSAL_setup, &
     PRC_GLOBAL_setup, &
     PRC_MPIfinish, &
     PRC_MPIsplit_nest, &
     PRC_UNIVERSAL_myrank, &
     PRC_DOMAIN_nlim
  use mod_rm_prep

  implicit none

  integer :: it, im, idom, color, key, ierr

  integer :: universal_comm
  integer :: universal_nprocs
  logical :: universal_master
  integer :: universal_myrank
  integer :: global_comm
  integer :: local_comm
  integer :: intercomm_parent ! not used
  integer :: intercomm_child

  character(len=H_LONG) :: confname_domains(PRC_DOMAIN_nlim)
  character(len=H_LONG) :: confname

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  ! start MPI
  call PRC_MPIstart( universal_comm ) ! [OUT]

  call mpi_timer('', 1)

  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
                            universal_nprocs, & ! [OUT]
                            universal_myrank, & ! [OUT]
                            universal_master  ) ! [OUT]
  universal_myrank = PRC_UNIVERSAL_myrank
  nprocs = universal_nprocs
  myrank = universal_myrank

  call set_common_conf( myrank )

!-----------------------------------------------------------------------

  if (DET_RUN) then
    call set_mem_node_proc(MEMBER+2)
  else
    call set_mem_node_proc(MEMBER+1)
  end if

  call mpi_timer('INITIALIZE', 1, barrier=universal_comm)

!-----------------------------------------------------------------------
! Run SCALE-RM_pp
!-----------------------------------------------------------------------

  ! split MPI communicator for single members
  if (myrank_to_mem(1) >= 1) then
    color = myrank_to_mem(1) - 1
    key   = myrank_to_pe
  else
    color = MPI_UNDEFINED
    key   = MPI_UNDEFINED
  endif

  call MPI_COMM_SPLIT(universal_comm, color, key, global_comm, ierr)

  if (global_comm /= MPI_COMM_NULL) then

    call PRC_GLOBAL_setup( .false.,    & ! [IN]
                           global_comm ) ! [IN]

    do idom = 1, NUM_DOMAIN
      confname_domains(idom) = trim(CONF_FILES)
      call filename_replace_dom(confname_domains(idom), idom)
    end do

    !--- split for nesting
    ! communicator split for nesting domains
    call PRC_MPIsplit_nest( global_comm,      & ! [IN]
                            NUM_DOMAIN,       & ! [IN]
                            PRC_DOMAINS(:),   & ! [IN]
                            .false.,          & ! [IN]
                            COLOR_REORDER,    & ! [IN]
                            local_comm,       & ! [OUT]
                            idom)!,             & ! [OUT]
                            !intercomm_parent, & ! [OUT]           
                            !intercomm_child )   ! [OUT]

    do it = 1, nitmax
      im = myrank_to_mem(it)
      if (im >= 1 .and. im <= MEMBER_RUN) then
        confname = confname_domains(idom)
        if (CONF_FILES_SEQNUM) then
          call filename_replace_mem(confname, im)
        else
          if (im <= MEMBER) then
            call filename_replace_mem(confname, im)
          else if (im == MEMBER+1) then
            call filename_replace_mem(confname, memf_mean)
          else if (im == MEMBER+2) then
            call filename_replace_mem(confname, memf_mdet)
          end if
        end if
        if ( LOG_OUT ) WRITE(6,'(A,I6.6,2A)') 'MYRANK ',universal_myrank,' is running a model with configuration file: ', trim(confname)

        call IO_set_universalrank( universal_myrank, & ! [IN]
                                   im,               & ! [IN]
                                   idom              ) ! [IN]

        call rm_prep ( local_comm,       &
                       trim(confname),   &
                       "",               &
                       .false.  )
                       !intercomm_parent, &
                       !intercomm_child,  &
                       !trim(confname)  )
      end if
    end do ! [ it = 1, nitmax ]

  else ! [ global_comm /= MPI_COMM_NULL ]

    write (6, '(A,I6.6,A)') 'MYRANK=',universal_myrank,': This process is not used!'

  end if ! [ global_comm /= MPI_COMM_NULL ]

  call mpi_timer('SCALE_RM', 1, barrier=universal_comm)

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  if ( myrank == 0 ) then
    write(6,'(a)') 'scale_rm_pp_ens finished successfully'
  endif

!  call PRC_MPIfinish

  call MPI_Finalize(ierr)

  stop
end program scaleles_pp_ens
