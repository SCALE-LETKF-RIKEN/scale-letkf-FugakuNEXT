program efso
!=======================================================================
!
! [PURPOSE:] Main program of forecast sensitivity to observations using LETKF
!
! [HISTORY:]
!   09/29/2011 Yoichiro Ohta     created from main program of LETKF
!   07/01/2013 Daisuke Hotta     ported to GFS-LEKTF system
!   12/19/2013 Guo-Yuan Lien     merged to GFS-LETKF main development
!
!=======================================================================
!$use omp_lib
  use common
  use common_nml
  use common_scale
  use common_mpi_scale
  use efso_tools
  use letkf_obs
  use letkf_tools

  implicit none

  real(r_size), allocatable :: gues3d(:,:,:)
  real(r_size), allocatable :: gues2d(:,:)
  real(r_size), allocatable :: fcst3d(:,:,:,:)
  real(r_size), allocatable :: fcst2d(:,:,:)
  real(r_size), allocatable :: fcer3d(:,:,:)
  real(r_size), allocatable :: fcer2d(:,:)
  real(r_size), allocatable :: work3d(:,:,:)
  real(r_size), allocatable :: work2d(:,:)
  real(RP),     allocatable :: work3dg(:,:,:,:)
  real(RP),     allocatable :: work2dg(:,:,:)
!  real(r_size), allocatable :: uadf(:,:), vadf(:,:)
!  real(r_size), allocatable :: uada(:,:), vada(:,:)

  real(r_sngl) :: rtimer00,rtimer
  INTEGER :: n, ierr
  CHARACTER(8) :: stdoutf='NOUT-000'
  CHARACTER(4) :: fcstf='fc01'
  CHARACTER(9) :: analf='anal0.grd'
  CHARACTER(9) :: gmeanf='gmean.grd'
  CHARACTER(9) :: fmean0='fme00.grd'
  CHARACTER(9) :: fmean6='fme06.grd'
  CHARACTER(9) :: ameanf='amean.grd'

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  call mpi_timer('', 1)

  if ( LOG_OUT ) write(6,'(a)') 'Hello from EFSO'

  if (DET_RUN) then
    call set_mem_node_proc(MEMBER+2)
  else
    call set_mem_node_proc(MEMBER+1)
  end if
  call set_scalelib('LETKF')

  if (myrank_use) then

    call set_common_scale
    call set_common_mpi_scale
    call set_common_obs_scale

    call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_a)
  
    !-----------------------------------------------------------------------
    ! Read observations
    !-----------------------------------------------------------------------
  
    allocate (obs(OBS_IN_NUM))
    call read_obs_all_mpi(obs)

    call mpi_timer('READ_OBS', 1, barrier=MPI_COMM_a)

    !
    ! EFSO GRID setup
    !
    call set_common_mpi_grid
    call mpi_timer('SET_GRID', 1, barrier=MPI_COMM_a)

    allocate( gues3d(nij1,nlev,nv3d) )
    allocate( gues2d(nij1,nv2d) )
    allocate( fcst3d(nij1,nlev,nens,nv3d) )
    allocate( fcst2d(nij1,nens,nv2d) )
    allocate( fcer3d(nij1,nlev,nv3d) )
    allocate( fcer2d(nij1,nv2d) )
    allocate( work3d(nij1,nlev,nv3d) )
    allocate( work2d(nij1,nv2d) )
    allocate( work3dg(nlev,nlon,nlat,nv3d) )
    allocate( work2dg(nlon,nlat,nv2d) )

    !-----------------------------------------------------------------------
    ! Read observation diagnostics
    !-----------------------------------------------------------------------

    call set_efso_obs
    call mpi_timer('READ_OBSDIAG', 1, barrier=MPI_COMM_a)

    !-----------------------------------------------------------------------
    ! Read model data
    !-----------------------------------------------------------------------
    !
    ! Forecast ensemble
    !
    call read_ens_mpi(fcst3d, fcst2d, EFSO=.true.)
    !!! fcst3d,fcst2d: (xmean+X)^f_t [Eq.(6), Ota et al. 2013]

    !
    ! Forecast error at evaluation time
    !
    ! forecast from the ensemble mean of first guess
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_FCST_FROM_GUES_BASENAME), work3dg, work2dg)
      call state_trans(work3dg)
    endif
    call scatter_grd_mpi(mmean_rank_e,real(work3dg,RP),real(work2dg,RP),fcer3d,fcer2d)

    ! forecast from the analysis ensemble mean
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_FCST_FROM_ANAL_BASENAME), work3dg, work2dg)
      call state_trans(work3dg)
    endif
    call scatter_grd_mpi(mmean_rank_e,real(work3dg,RP),real(work2dg,RP),work3d,work2d)
    fcer3d(:,:,:) = 0.5_r_size * ( fcer3d(:,:,:) + work3d(:,:,:) )
    fcer2d(:,:) = 0.5_r_size * ( fcer2d(:,:) + work2d(:,:) )

    if ( myrank_a == 0 ) then
do n = 1, nv3d
  write(6,'(a,i7,2f10.1)')"DEBUG1 ", n, maxval( work3d(:,:,n) ), minval( work3d(:,:,n) )
  write(6,'(a,i7,2f10.1)')"DEBUG2 ", n, maxval( fcer3d(:,:,n) ), minval( fcer3d(:,:,n) )
enddo
    endif

    ! reference analysis ensemble mean
    if ( myrank_e == mmean_rank_e ) then  
      call read_restart( trim(EFSO_ANAL_IN_BASENAME), work3dg, work2dg)
      call state_trans(work3dg)
    endif
    call scatter_grd_mpi(mmean_rank_e,real(work3dg,RP),real(work2dg,RP),work3d,work2d)
    deallocate( work3dg, work2dg )

    if ( myrank_a == 0 ) then
do n = 1, nv3d
  write(6,'(a,i7,2f10.1)')"DEBUG3 ", n, maxval( work3d(:,:,n) ), minval( work3d(:,:,n) )
enddo
    endif

    !!! fcer3d,fcer2d: [1/2(K-1)](e^f_t+e^g_t) [Eq.(6), Ota et al. 2013]
    fcer3d(:,:,:) = ( fcer3d(:,:,:) - work3d(:,:,:) ) / real( MEMBER-1, r_size )
    fcer2d(:,:) = ( fcer2d(:,:) - work2d(:,:) ) / real( MEMBER-1, r_size )

    deallocate( work3d, work2d )

    if ( myrank_a == 0 ) then
do n = 1, nv3d
  write(6,'(a,i7,2f10.1)')"DEBUG4 ", n, maxval( fcer3d(:,:,n) ), minval( fcer3d(:,:,n) )
enddo
    endif
    !
    ! Norm
    !
    call lnorm( fcst3d, fcst2d, fcer3d, fcer2d )
    !! fcst3d,fcst2d: C^(1/2)*X^f_t
    !! fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)

    !-----------------------------------------------------------------------
    ! Winds for advection
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! EFSO computation
    !-----------------------------------------------------------------------
    call init_obsense
    call das_efso( gues3d, gues2d, fcst3d, fcst2d, fcer3d, fcer2d )

    deallocate( gues3d, gues2d )
    deallocate( fcst3d, fcst2d )
    deallocate( fcer3d, fcer2d )

    !-----------------------------------------------------------------------
    ! EFSO output
    !-----------------------------------------------------------------------
    !  IF(myrank == 0) CALL print_obsense()
    call destroy_obsense()
    !

  end if ! [ myrank_use ]

  call unset_scalelib

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  call finalize_mpi_scale

  stop

  !
  ! Forecast error at evaluation time
  !
  IF(myrank == 0) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',fmean0
    !CALL read_grd4(fmean0,work3dg,work2dg,1)
  END IF
  CALL scatter_grd_mpi(0,real(work3dg,RP),real(work2dg,RP),fcer3d,fcer2d)
  IF(myrank == 0) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',fmean6
    !CALL read_grd4(fmean6,work3dg,work2dg,1)
  END IF
  CALL scatter_grd_mpi(0,real(work3dg,RP),real(work2dg,RP),work3d,work2d)
  fcer3d(:,:,:) = 0.5_r_size * (fcer3d(:,:,:) + work3d(:,:,:))
  fcer2d(:,:) = 0.5_r_size * (fcer2d(:,:) + work2d(:,:))
  IF(myrank == 0) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',ameanf
    !CALL read_grd4(ameanf,work3dg,work2dg,1)
  END IF
  CALL scatter_grd_mpi(0,real(work3dg,RP),real(work2dg,RP),work3d,work2d)
  fcer3d(:,:,:) = (fcer3d(:,:,:) - work3d(:,:,:)) / real(MEMBER-1,r_size)
  fcer2d(:,:) = (fcer2d(:,:) - work2d(:,:)) / real(MEMBER-1,r_size)
  !!! fcer3d,fcer2d: [1/2(K-1)](e^f_t+e^g_t) [Eq.(6), Ota et al. 2013]
  !
  ! Norm
  !
  !CALL lnorm(fcst3d,fcst2d,fcer3d,fcer2d)
  !!! fcst3d,fcst2d: C^(1/2)*X^f_t
  !!! fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)
  !
  ! Guess mean for full-level pressure computation
  !
  IF(myrank == 0) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',gmeanf
   ! CALL read_grd4(gmeanf,work3dg,work2dg,0)
  END IF
  CALL scatter_grd_mpi(0,real(work3dg,RP),real(work2dg,RP),gues3d,gues2d)
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_FORECAST):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Winds for advection
!-----------------------------------------------------------------------
!  IF(ABS(locadv_rate) > TINY(locadv_rate)) THEN
!    ALLOCATE(uadf(nij1,nlev))
!    ALLOCATE(vadf(nij1,nlev))
!    ALLOCATE(uada(nij1,nlev))
!    ALLOCATE(vada(nij1,nlev))
!    uadf(:,:) = work3d(:,:,iv3d_u)
!    vadf(:,:) = work3d(:,:,iv3d_v)
!    IF(myrank == 0) THEN
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',analf
!      CALL read_grd4(analf,work3dg,work2dg,0)
!    END IF
!    CALL scatter_grd_mpi(0,real(work3dg,r_size),real(work2dg,r_size),work3d,work2d)
!    uada(:,:) = work3d(:,:,iv3d_u)
!    vada(:,:) = work3d(:,:,iv3d_v)
!    CALL loc_advection(uada,vada,uadf,vadf) ! ADVECTION for FSO
!    DEALLOCATE(uadf,vadf,uada,vada)
!!
!    CALL CPU_TIME(rtimer)
!    WRITE(6,'(A,2F10.2)') '### TIMER(WIND_ADVECTION):',rtimer,rtimer-rtimer00
!    rtimer00=rtimer
!  END IF
!  DEALLOCATE(work3d,work2d)
!-----------------------------------------------------------------------
! EFSO computation
!-----------------------------------------------------------------------
  !CALL init_obsense()
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !CALL das_efso(gues3d,gues2d,fcst3d,fcst2d,fcer3d,fcer2d)
  DEALLOCATE(gues3d,gues2d,fcst3d,fcst2d,fcer3d,fcer2d)
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(DAS_EFSO):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! EFSO output
!-----------------------------------------------------------------------
!  IF(myrank == 0) CALL print_obsense()
!  CALL destroy_obsense()
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(EFSO_OUTPUT):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  if ( myrank == 0 ) then
    write(6,'(a)') 'efso finished successfully'
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  stop
end program efso
