module radar_obs
!=======================================================================
!
! [PURPOSE:] Radar observations
!
!=======================================================================
!$USE OMP_LIB
  use common
  use common_nml
  use common_scale
  use common_obs_scale
  use common_mpi_scale, only: &
    mpi_timer
  use scale_precision, only: DP
  use radar_tools, only: & 
    MPI_COMM_o, nprocs_o, myrank_o, &
    radar_georeference, & 
    define_grid, &  
    radar_superobing 
  use dec_pawr_nml

  implicit none
  public
  integer, save :: utime_obs(6) = (/-1,-1,-1,-1,-1,-1/)

  ! time label
  character(len=24), save :: timelabel_obs  = ''
  
contains

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
subroutine read_obs_radar_toshiba(cfile, obs)
  use iso_c_binding
#ifdef MPW
  use read_toshiba_mpr_f
#else
  use read_toshiba_f
#endif

  use scale_atmos_grid_cartesC, only: &
      DX, DY
  use scale_time, only: &
      TIME_gettimelabel, &
      TIME_NOWDATE
  implicit none

  character(len=*), intent(in) :: cfile
  type(obs_info), intent(out) :: obs
  type(obs_info) :: obs_ref
!  REAL(r_sngl) :: wk(8)
!  INTEGER :: nrec
!  REAL(r_sngl) :: tmp
!  INTEGER :: n,iunit,ios

  character(len=1024) :: jitdt_place
#ifdef MPW
  integer, parameter :: n_type = 2
  character(len=4), parameter :: file_type_sfx(n_type) = &
   (/'.ze', '.vr'/)
  logical, parameter :: input_is_dbz = .true.
  integer, parameter :: opt_verbose = 1 !!! for MP-PAWR toshiba format    

  integer :: access !FILE INQUIRY
  integer :: ios
  integer(4),save :: shadow_na, shadow_ne
  integer(4), allocatable, save:: tmpshadow(:)
  integer(2), allocatable ,save :: shadow(:,:)
  real(8),save :: shadow_del_az
#else
  integer, parameter :: n_type = 3
  character(len=4), parameter :: file_type_sfx(n_type) = &
    (/'.ze ', '.vr ', '.qcf'/)
  logical, parameter :: input_is_dbz = .true.
#endif

#ifdef MPW
  type(c_mppawr_header) :: hd(n_type)
#else
  type(c_pawr_header) :: hd(n_type)
#endif
!  real(kind=c_float) :: az(AZDIM, ELDIM, n_type)
!  real(kind=c_float) :: el(AZDIM, ELDIM, n_type)
!  real(kind=c_float) :: rtdat(RDIM, AZDIM, ELDIM, n_type)
  real(kind=c_float), allocatable, save :: rtdat(:, :, :, :)
  real(kind=c_float), allocatable, save :: az(:, :, :)
  real(kind=c_float), allocatable, save :: el(:, :, :)
  integer :: j, ierr, ierr2
  character(len=3) :: fname
  integer, save::i=0

  real(r_size), allocatable :: ze(:, :, :), vr(:, :, :), qcflag(:, :, :), attenuation(:, :, :), rrange(:)
  real(r_size), allocatable :: radlon(:, :, :), radlat(:, :, :), radz(:, :, :)
  real(r_size), allocatable :: lon(:), lat(:), z(:)
  integer(8), allocatable :: grid_index(:), grid_count_ze(:), grid_count_vr(:)
  real(r_size), allocatable :: grid_ze(:), grid_vr(:)
  real(r_size), allocatable :: grid_lon_ze(:), grid_lat_ze(:), grid_z_ze(:)
  real(r_size), allocatable :: grid_lon_vr(:),  grid_lat_vr(:),  grid_z_vr(:)

  character(len=1024) :: input_fname(n_type)
  integer ia, ir, ie
  real(r_size) :: dlon, dlat
  integer :: nlon , nlat , nlev
  integer(8) nobs_sp

  integer,save :: na, nr, ne
  real(r_size),save :: lon0, lat0, z0
  real(r_size),save :: missing
  integer,save :: range_res

  real(r_size) :: max_obs_ze , min_obs_ze , max_obs_vr , min_obs_vr 
  integer :: nobs_ze, nobs_vr
  integer(8) :: idx, n, n_ref
  integer :: pos
  integer, parameter :: int1 = selected_int_kind(1) !1-BYTE INT
  integer(kind = int1) :: tmp_qcf, valid_qcf

!  integer,parameter :: qcf_mask(8)=(/ 1, 0, 0, 0, 0, 0, 0, 0 /) !! valid, shadow, clutter possible, clutter certain, interference, range sidelobe /
  integer,parameter :: qcf_mask(8)=(/ 0, 1, 1, 1, 1, 0, 0, 0 /) !! valid, shadow, clutter possible, clutter certain, interference, range sidelobe /
!!!  integer,parameter :: qcf_mask(8)=(/ 0, 0, 0, 0, 0, 0, 0, 0 /) !! valid, shadow, clutter possible, clutter certain, interference, range sidelobe /

  integer::qcf_count(0:255)

  character(len=90) :: plotname
  character(len=19) :: timelabel

  integer :: ii, jj, kk

  real(r_sngl), allocatable :: ref3d(:,:,:)
  real(r_sngl), allocatable :: vr3d(:,:,:)
  character(len=255) :: filename
  integer :: irec, iunit, iolen
  integer :: k

  call mpi_timer('', 3)

  RADAR_SO_SIZE_HORI = max( real( DX, kind=r_size ), RADAR_SO_SIZE_HORI )
  RADAR_SO_SIZE_HORI = max( real( DY, kind=r_size ), RADAR_SO_SIZE_HORI )

!!! MP-PAWR shadow masking 
#ifdef MPW
  if ( USE_PAWR_MASK .and. .not. (allocated(shadow)) ) then
    if (myrank_o == 0)then
        write(6, '("reading ", A)') trim(pawr_mask_file)
        open(99, file = trim(pawr_mask_file), status = "old", access = "stream", form = "unformatted", convert = "little_endian")
        read(99,iostat=ios) shadow_na, shadow_ne
        if ( ios == 0 )then
          allocate(shadow(shadow_na, shadow_ne))
          read(99,iostat=ios) shadow
          close(99)
          if( ios /= 0 ) shadow = 0
        else
          write(6,'(3A)') 'file ',trim(pawr_mask_file) ,' not found or unsupported format.'
          stop 1
        end if 
    end if
    if ( nprocs_o /= 1 )then
      call MPI_BCAST(shadow_na, 1, MPI_INTEGER4, 0, MPI_COMM_o, ierr)
      call MPI_BCAST(shadow_ne, 1, MPI_INTEGER4, 0, MPI_COMM_o, ierr)
      if (myrank_o /= 0) allocate(shadow(shadow_na,shadow_ne))
      call MPI_BCAST(shadow, shadow_na*shadow_ne, MPI_INTEGER2, 0, MPI_COMM_o, ierr)
    end if
  end if
#endif

    if ( .not. (allocated(rtdat)) ) allocate(rtdat(RDIM, AZDIM, ELDIM, n_type))
    if ( .not. (allocated(az)) ) allocate(az(AZDIM, ELDIM, n_type))
    if ( .not. (allocated(el)) ) allocate(el(AZDIM, ELDIM, n_type))

    if (LOG_LEVEL >= 3 .and. myrank_o == 0) then
      write(6, *) RDIM, AZDIM, ELDIM
      write(6, *) "dx = ", RADAR_SO_SIZE_HORI
      write(6, *) "dy = ", RADAR_SO_SIZE_HORI
      write(6, *) "dz = ", RADAR_SO_SIZE_VERT
    endif

        if (myrank_o == 0)then
          do j = 1, n_type
            input_fname(j) = trim(cfile)
            
            call str_replace(input_fname(j), '<type>', trim(file_type_sfx(j)), pos)

            if (pos == 0) then
              write (6, '(5A)') "[Error] Keyword '<type>' is not found in '", trim(cfile), "'."
              stop 1
            end if
          end do
        
          if (LOG_LEVEL >= 3) then
            write(6, *) "file1 = ", trim(input_fname(1))
            write(6, *) "file2 = ", trim(input_fname(2))
#ifndef MPW
            write(6, *) "file3 = ", trim(input_fname(3))
#endif
          endif

         endif

          do j = 1, n_type
            if (myrank_o == 0)then
#ifdef MPW
              ierr = read_toshiba_mpr(input_fname(j), opt_verbose, hd(j), az(:, :, j), el(:, :, j), rtdat(:, :, :, j))
#else
              ierr = read_toshiba(input_fname(j), hd(j), az(:, :, j), el(:, :, j), rtdat(:, :, :, j))
#endif
              if (LOG_LEVEL >= 3) then
                write(6, *) "return code = ", ierr
              endif
            end if
            call MPI_BCAST(ierr, 1, MPI_INTEGER, 0, MPI_COMM_o, ie)
            if (ierr /= 0) then
              obs%nobs = 0
              return
            endif
          end do
        
      call mpi_timer('read_obs_radar_toshiba:read_toshiba:', 2, barrier=MPI_COMM_o)
 
      ! Set obs information
      if (myrank_o == 0) then
        lon0 = hd(1)%longitude
        lat0 = hd(1)%latitude
        z0 = hd(1)%altitude
        missing = real(hd(1)%mesh_offset, r_size)
        range_res = hd(1)%range_res
      
        nr = hd(1)%range_num
#ifdef MPW
        na = hd(1)%ray_num
#else
        na = hd(1)%sector_num
#endif
        ne = hd(1)%el_num
    
        call jst2utc(hd(1)%s_yr, hd(1)%s_mn, hd(1)%s_dy, hd(1)%s_hr, hd(1)%s_mi, hd(1)%s_sc, 0.0_DP, utime_obs)
        if (myrank_o == 0 ) then
          write(6,'(a)') "get new data ..."
          write(6,'(a,i4.4,i2.2,i2.2,1x,i2.2,1a,i2.2,1a,i2.2)') "SCALE-LETKF:",&
                TIME_NOWDATE(1),TIME_NOWDATE(2),TIME_NOWDATE(3),&
                TIME_NOWDATE(4),":",TIME_NOWDATE(5),":",TIME_NOWDATE(6)
          write(6,'(a,i4.4,i2.2,i2.2,1x,i2.2,1a,i2.2,1a,i2.2)') "PAWR OBS:",&
                utime_obs(1),utime_obs(2),utime_obs(3),&
                utime_obs(4),":",utime_obs(5),":",utime_obs(6)
        endif
       
      endif ! [ myrank_o == 0 ]


      ! broadcast obs information
      call  pawr_toshiba_hd_mpi(lon0, lat0, z0, missing,&
                               range_res, na, nr, ne, &
                               AZDIM, ELDIM, n_type, RDIM, &
                               az, el, rtdat, utime_obs)

      
      call mpi_timer('read_obs_radar_toshiba:comm:', 2)
#ifdef MPW
   ierr=0
   if (myrank_o == 0 ) then
     if ((iand(hd(1)%status, 1) .ne. 0) .or. (iand(hd(1)%status, 2) .ne. 0) &
  & .or. (iand(hd(2)%status, 1) .ne. 0) .or. (iand(hd(2)%status, 2) .ne. 0)) then
       ierr=1
     end if 
   end if ! [ myrank_o == 0 ]
   call MPI_BCAST(ierr, 1, MPI_INTEGER, 0, MPI_COMM_o, ie)
   if (ierr /= 0) then
     write(6,'(a)') "MP-PAWR status error. Skip."
     obs%nobs = 0
     return
   end if
#endif

  if (LOG_LEVEL >= 2 .and. myrank_o == 0) then
!    write(*, '(I4.4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2, &
!         &     " -> ", I4.4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2)') &
!         & hd(1)%s_yr, hd(1)%s_mn, hd(1)%s_dy, hd(1)%s_hr, hd(1)%s_mi, hd(1)%s_sc, &
!         & hd(1)%e_yr, hd(1)%e_mn, hd(1)%e_dy, hd(1)%e_hr, hd(1)%e_mi, hd(1)%e_sc
    write(*, '(I4.4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2)') &
         & utime_obs(1), utime_obs(2), utime_obs(3), utime_obs(4), utime_obs(5), utime_obs(6)
     write(*, *) lon0, lat0, z0
    write(*, *)   na,   nr, ne
    write(*, *) "missing = ", missing
  endif


  allocate(ze(na, nr, ne), vr(na, nr, ne), qcflag(na, nr, ne), attenuation(na, nr, ne))

#ifdef MPW
  if ( USE_PAWR_MASK .and. .not. allocated(tmpshadow) ) then
     shadow_del_az = 360.0d0 / shadow_na
     allocate(tmpshadow(na))
  end if
#else
  valid_qcf = 0
  do j = 1, 8  
    if(qcf_mask(j) > 0) valid_qcf = ibset(valid_qcf, j - 1) 
  end do
#endif

!$omp parallel do &
#ifdef MPW
!$omp private(tmpshadow) &
#else
!$omp private(tmp_qcf) &
#endif
!$omp private(ia, ir, ie)

  do ie = 1, ne
     do ir = 1, nr
        do ia = 1, na
           ze(ia, ir, ie) = rtdat(ir, ia, ie, 1)
           vr(ia, ir, ie) = rtdat(ir, ia, ie, 2)                                                                 
     
           if(vr(ia, ir, ie) > RADAR_MAX_ABS_VR .or. vr(ia, ir, ie) < -RADAR_MAX_ABS_VR) vr(ia, ir, ie) = missing

#ifdef MPW
         end do
       end do

     do ir = 1, nr
        do ia = 1, na
           qcflag(ia,ir,ie) = 0.0d0 !valid
        end do
     end do
     if ( USE_PAWR_MASK .and. allocated(shadow) ) then
        if (ie <= shadow_ne) then
           do ia = 1, na
              tmpshadow(ia) = shadow(min(shadow_na,nint(az(ia, ie, 1) / shadow_del_az) + 1), ie)
           end do
           do ir = 1, nr
              do ia = 1, na
                 if ( tmpshadow(ia) /= 0 .and. ir >= tmpshadow(ia) ) then
                    qcflag(ia, ir, ie) = 1000.0d0  !invalid
                 end if
              end do
           end do
        end if
     end if

     do ir = 1, nr
        do ia = 1, na

#else
           tmp_qcf = int(rtdat(ir, ia, ie, 3), int1)
!          qcf_count(tmp_qcf)=qcf_count(tmp_qcf)+1
           if(iand(valid_qcf, tmp_qcf) == 0) then      
             qcflag(ia, ir, ie) = 0.0d0 !valid
           else
             qcflag(ia, ir, ie) = 1000.0d0 !invalid
           end if
#endif

           attenuation(ia, ir, ie) = 1.0d0 !not implemented yet
        end do
     end do
  end do
!$omp end parallel do
  deallocate(rtdat)


  allocate(rrange(nr))
!$omp parallel do private(ir)
  do ir = 1, nr
     rrange(ir) = (dble(ir) - 0.5d0) * range_res
  end do
!$omp end parallel do

  call mpi_timer('read_obs_radar_toshiba:qc:', 2)

  allocate(radlon(na, nr, ne), radlat(na, nr, ne), radz(na, nr, ne))

  call radar_georeference(lon0, lat0, z0, na, nr, ne, &                                   ! input
       &                  real(az(1:na, 1, 1), r_size), rrange, real(el(1, 1:ne, 1), r_size), & ! input (assume ordinary scan strategy)
       &                  radlon, radlat, radz, &                                     ! output  
       &                  MPI_COMM_o)

  call mpi_timer('read_obs_radar_toshiba:radar_georeference:', 2)


  call define_grid(lon0, lat0, nr, rrange, rrange(nr), RADAR_ZMAX, & ! input
                   RADAR_SO_SIZE_HORI, RADAR_SO_SIZE_HORI, RADAR_SO_SIZE_VERT, & ! input
       &           dlon, dlat, nlon, nlat, nlev, lon, lat, z)              ! output

  call mpi_timer('read_obs_radar_toshiba:define_grid:', 2)

  call radar_superobing(na, nr, ne, radlon, radlat, radz, ze, vr, & ! input spherical
        &                qcflag, attenuation, & ! input spherical
        &                nlon, nlat, nlev, lon, lat, z, dlon, dlat, RADAR_SO_SIZE_VERT, & ! input cartesian
        &                missing, input_is_dbz, & ! input param
        &                lon0, lat0, &
        &                nobs_sp, grid_index, & ! output array info
        &                grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze, & ! output ze
        &                grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr, & ! output vr
!        &                0) !!! dummy : not used
        &                MPI_COMM_o)

  if(allocated(ze)) deallocate(ze)
  if(allocated(vr)) deallocate(vr)
  if(allocated(qcflag)) deallocate(qcflag)
  if(allocated(attenuation)) deallocate(attenuation)
  if(allocated(rrange)) deallocate(rrange)

  call mpi_timer('read_obs_radar_toshiba:radar_superobing:', 2)


  call MPI_BCAST(grid_lon_ze, int(nobs_sp), MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(grid_lat_ze, int(nobs_sp), MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(grid_z_ze, int(nobs_sp), MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(grid_ze, int(nobs_sp), MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(grid_vr, int(nobs_sp), MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(grid_count_ze, int(nobs_sp), MPI_INTEGER8, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(grid_count_vr, int(nobs_sp), MPI_INTEGER8, 0, MPI_COMM_o, ierr)
  call mpi_timer('read_obs_radar_toshiba:plot_comm:', 2, barrier=MPI_COMM_o)


!!!!! check
!write(*,*) nlon,nlat,nlev,nobs_sp

  if ( OUT_PAWR_GRADS ) then
    if (.not. allocated(ref3d) ) allocate(ref3d(nlon,nlat,nlev))
    if (.not. allocated(vr3d) ) allocate(vr3d(nlon,nlat,nlev))
    ref3d = undef
    vr3d = undef
  endif

  obs%meta(1) = lon0
  obs%meta(2) = lat0
  obs%meta(3) = z0

  obs%nobs = 0
  obs_ref%nobs = 0
  do idx = 1, nobs_sp

    ! Thinning
    ii = nint( ( grid_lon_ze(idx) - lon(1) ) / dlon ) + 1
    jj = nint( ( grid_lat_ze(idx) - lat(1) ) / dlat ) + 1
    kk = nint( ( grid_z_ze(idx) - z(1) ) / RADAR_SO_SIZE_VERT ) + 1

    if ( mod(ii, RADAR_THIN_HORI) /= 0 .or. mod(jj, RADAR_THIN_HORI) /= 0 .or. &
         mod(kk, RADAR_THIN_VERT) /= 0 ) cycle

    if (grid_count_ze(idx) > 0) then
      obs%nobs = obs%nobs + 1

      ! Count refrectivity obs ( > MIN_RADAR_REF ) below RADAR_ZMAX
      if ( grid_ze(idx) > MIN_RADAR_REF .and. grid_z_ze(idx) < RADAR_ZMAX ) then
        obs_ref%nobs = obs_ref%nobs + 1
      end if
    end if
    if (grid_count_vr(idx) > 0) then
      obs%nobs = obs%nobs + 1
    end if
  end do
  call obs_info_allocate(obs, extended=.true.)
  call obs_info_allocate(obs_ref, extended=.true.)

  n = 0
  n_ref = 0
  nobs_ze = 0
  nobs_vr = 0
  min_obs_ze = huge(real(1.0,r_size))
  max_obs_ze = -huge(real(1.0,r_size))
  min_obs_vr = huge(real(1.0,r_size))
  max_obs_vr = -huge(real(1.0,r_size))
  do idx = 1, nobs_sp

    ii = nint( ( grid_lon_ze(idx) - lon(1) ) / dlon ) + 1
    jj = nint( ( grid_lat_ze(idx) - lat(1) ) / dlat ) + 1
    kk = nint( ( grid_z_ze(idx) - z(1) ) / RADAR_SO_SIZE_VERT ) + 1
    if ( OUT_PAWR_GRADS ) then
      if ( ii > 0 .and. ii <= nlon .and. &
           jj > 0 .and. jj <= nlat .and. &
           kk > 0 .and. kk <= nlev ) then
        if ( grid_count_ze(idx) > 0 .and. &
           grid_ze(idx) > 0.0_r_size ) then
          ref3d(ii,jj,kk) = 10.0*log10(grid_ze(idx))
        endif
        if ( grid_count_vr(idx) > 0 ) then
          vr3d(ii,jj,kk) = grid_vr(idx)
        endif
      endif
    endif

    ! Thinning
    if ( mod(ii, RADAR_THIN_HORI) /= 0 .or. mod(jj, RADAR_THIN_HORI) /= 0 .or. &
         mod(kk, RADAR_THIN_VERT) /= 0 ) cycle

    if (grid_count_ze(idx) > 0) then
      n = n + 1
      obs%elm(n) = id_radar_ref_obs
      obs%lon(n) = grid_lon_ze(idx)
      obs%lat(n) = grid_lat_ze(idx)
      obs%lev(n) = grid_z_ze(idx)
      obs%dat(n) = grid_ze(idx)
      ! Add RADAR_BIAS_CONST_DBZ in dBZ
      if ( RADAR_BIAS_COR_RAIN .and. grid_ze(idx) > MIN_RADAR_REF ) then 
        obs%dat(n) = grid_ze(idx) * RADAR_BIAS_RAIN_CONST
      elseif ( RADAR_BIAS_COR_CLR .and. grid_ze(idx) < MIN_RADAR_REF )  then
        obs%dat(n) = grid_ze(idx) * RADAR_BIAS_CLR_CONST
      endif
      obs%err(n) = OBSERR_RADAR_REF
      obs%typ(n) = 22
      obs%dif(n) = 0.0d0
      nobs_ze = nobs_ze + 1
      if (grid_ze(idx) > max_obs_ze) max_obs_ze = grid_ze(idx)
      if (grid_ze(idx) < min_obs_ze) min_obs_ze = grid_ze(idx)

      if ( grid_ze(idx) > MIN_RADAR_REF .and. grid_z_ze(idx) < RADAR_ZMAX ) then
        n_ref = n_ref + 1
        obs_ref%elm(n_ref) = id_radar_ref_obs
        obs_ref%lon(n_ref) = grid_lon_ze(idx)
        obs_ref%lat(n_ref) = grid_lat_ze(idx)
        obs_ref%lev(n_ref) = grid_z_ze(idx)
        if ( RADAR_BIAS_COR_RAIN ) then
          ! Add RADAR_BIAS_CONST_DBZ in dBZ
          obs_ref%dat(n_ref) = grid_ze(idx) * RADAR_BIAS_RAIN_CONST
        else
          obs_ref%dat(n_ref) = grid_ze(idx)
        end if
      end if
    end if

    if (grid_count_vr(idx) > 0) then
      n = n + 1
      obs%elm(n) = id_radar_vr_obs
      obs%lon(n) = grid_lon_ze(idx)
      obs%lat(n) = grid_lat_ze(idx)
      obs%lev(n) = grid_z_ze(idx)
      obs%dat(n) = grid_vr(idx)
      obs%err(n) = OBSERR_RADAR_VR
      obs%typ(n) = 22
      obs%dif(n) = 0.0d0
      nobs_vr = nobs_vr + 1
      if (grid_vr(idx) > max_obs_vr) max_obs_vr = grid_vr(idx)
      if (grid_vr(idx) < min_obs_vr) min_obs_vr = grid_vr(idx)
    end if


    write(6,'(2I8, 5F12.4)') n , obs%elm(n),  obs%lon(n), obs%lat(n), obs%lev(n), obs%dat(n), obs%err(n)

  end do

  if (LOG_LEVEL >= 2 .and. myrank_o == 0 ) then
    write (6, *) "Reflectivity obs. range = ", min_obs_ze, " to ", max_obs_ze
    write (6, *) "Radial vel. obs. range  = ", min_obs_vr, " to ", max_obs_vr
    write (6, *) "ze: ", nobs_ze, ", vr: ", nobs_vr, ", ze(rain): ", obs_ref%nobs
  endif

  call mpi_timer('read_obs_radar_toshiba:save_obs_info:', 2)

  if(allocated(radlon)) deallocate(radlon)
  if(allocated(radlat)) deallocate(radlat)
  if(allocated(radz)) deallocate(radz)
  deallocate(az, el)

  if ( OUT_PAWR_GRADS ) then
    if ( myrank_o == 0 ) then
      call TIME_gettimelabel(timelabel)
      filename = trim(OUT_PAWR_GRADS_PATH)//trim(timelabel(1:15))//".grd"
      iunit = 55
      inquire (iolength=iolen) iolen
      open(iunit, file=trim(filename), form='unformatted', access='direct', &
            status='unknown', convert='big_endian', recl=nlon*nlat*iolen)
      irec = 0

      do k = 1, nlev
        irec = irec + 1
        write(iunit, rec=irec) ref3d(:,:,k)
      enddo
      do k = 1, nlev
        irec = irec + 1
        write(iunit, rec=irec) vr3d(:,:,k)
      enddo
      write(6,'(a)') 'PAWR GrADS info'
      write(6,'(i5,2f13.8)') nlon, lon(1), dlon
      write(6,'(i5,2f13.8)') nlat, lat(1), dlat
      write(6,'(i5,2f13.8)') nlev, z(1), RADAR_SO_SIZE_VERT
      write(6,'(a)') ''

      close( iunit )
    endif
  endif


  if ( myrank_o == 0 ) then
    call TIME_gettimelabel(timelabel)
    filename=trim(OUT_PAWR_SUPEROB_PATH)//"_"//timelabel(1:8)//timelabel(10:15)//".dat"
    call write_obs_radar(trim(filename),obs,missing=.false.)
    write(6,*) 'processed PAWR obs data output : ', trim(filename) 
  end if
 
  call obs_info_deallocate( obs_ref )
  if (myrank_o /= 0) then
    call obs_info_deallocate( obs )
  endif


!  write(*, *) "writing LETKF data ..."
!  call output_letkf_obs(lon0, lat0, z0, 1, nobs_sp, &
!       &                grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze, &
!       &                grid_vr, grid_count_vr, &
!       &                OBSERR_RADAR_REF, OBSERR_RADAR_VR)
!  write(*, *) "done"

!  !!! OUTPUT CARTESIAN COORDINATE DATE FOR DEBUG !!!
!if (nobs_ze.gt.0) then 
!  write(*, *) "writing cartesian data ..."
!  i=i+1
!  write(fname, '(I03.3)') i
!  call output_grads_obs("super_" // fname // ".grd", nlon, nlat, nlev, nobs_sp, grid_index, &
!       &                grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze, &
!       &                grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr, missing)
!  write(*, *) "done"
! endif



! Assume radar geographical information does not change during DA cyclying
!  if(allocated(radlon)) deallocate(radlon)
!  if(allocated(radlat)) deallocate(radlat)
!  if(allocated(radz)) deallocate(radz)



  if(allocated(grid_index)) deallocate(grid_index)
  if(allocated(grid_ze)) deallocate(grid_ze)
  if(allocated(grid_lon_ze)) deallocate(grid_lon_ze)
  if(allocated(grid_lat_ze)) deallocate(grid_lat_ze)
  if(allocated(grid_z_ze)) deallocate(grid_z_ze)
  if(allocated(grid_count_ze)) deallocate(grid_count_ze)
  if(allocated(grid_vr)) deallocate(grid_vr)
  if(allocated(grid_lon_vr)) deallocate(grid_lon_vr)
  if(allocated(grid_lat_vr)) deallocate(grid_lat_vr)
  if(allocated(grid_z_vr)) deallocate(grid_z_vr)
  if(allocated(grid_count_vr)) deallocate(grid_count_vr)

  call mpi_timer('read_obs_radar_toshiba:deallocate_vars:', 2)

  return
end subroutine read_obs_radar_toshiba


!!! not used
subroutine read_obs_radar_jrc(cfile, obs)
  use netcdf
  use common_ncio
  implicit none

  character(*), intent(in) :: cfile
  type(obs_info), intent(inout) :: obs
  integer :: ncid, varid
  integer :: status
  integer :: xdim, ydim, zdim
  integer :: i, j, k
  integer :: n, nobs_vr, nobs_zh

  real(r_sngl), allocatable :: vr3d(:,:,:), zh3d(:,:,:)
  real(r_sngl), allocatable :: lon1d(:), lat1d(:), z1d(:)
  real(r_sngl) :: sf_vr, sf_zh ! scale factor for vr & Zh
  real(r_sngl) :: fill_vr, fill_zh ! fill values for vr & Zh
  real(r_size) :: radar_lon, radar_lat
#ifdef SINGLE_LETKF
  real(8) :: radar_lon_r8, radar_lat_r8
#endif

  real(r_sngl) :: max_obs_zh, min_obs_zh
  real(r_sngl) :: max_obs_vr, min_obs_vr

  call mpi_timer('', 3)

  status = nf90_open(trim(cfile), NF90_NOWRITE, ncid)
  if (status /= nf90_noerr) then
    obs%nobs = 0
    return
  endif

  call ncio_read_dim(ncid, "lon", xdim)
  call ncio_read_dim(ncid, "lat", ydim)
  call ncio_read_dim(ncid, "alt", zdim)

  allocate(vr3d(xdim,ydim,zdim), zh3d(xdim,ydim,zdim))
  allocate(lon1d(xdim), lat1d(ydim), z1d(zdim)) 

  ! get lon/lat/height information
  call ncio_check(nf90_inq_varid(ncid, "lon",  varid))
  call ncio_check(nf90_get_var(ncid, varid, lon1d(:), &
                               start = (/ 1/),    &
                               count = (/ xdim/)))

  call ncio_check(nf90_inq_varid(ncid, "lat",  varid))
  call ncio_check(nf90_get_var(ncid, varid, lat1d(:), &
                               start = (/ 1/),    &
                               count = (/ ydim/)))

  call ncio_check(nf90_inq_varid(ncid, "alt",  varid))
  call ncio_check(nf90_get_var(ncid, varid, z1d(:), &
                               start = (/ 1/),    &
                               count = (/ zdim/)))



  call ncio_check(nf90_inq_varid(ncid, "Vel",  varid))
  call ncio_check(nf90_get_var(ncid, varid, vr3d(:,:,:), &
                               start = (/ 1, 1, 1/),    &
                               count = (/ xdim, ydim, zdim /)))

  status = nf90_get_att(ncid,varid,"scale_factor",sf_vr)
  if (status /= nf90_noerr) then
    sf_vr = 1.0
  endif
  status = nf90_get_att(ncid,varid,"_FillValue",fill_vr)
  if (status /= nf90_noerr) then
    fill_vr = -32768.0
  endif

  call ncio_check(nf90_inq_varid(ncid, "Zh MTI",  varid))
  call ncio_check(nf90_get_var(ncid, varid, zh3d(:,:,:), &
                               start = (/ 1, 1, 1/),    &
                               count = (/ xdim, ydim, zdim /)))

  status = nf90_get_att(ncid,varid,"scale_factor",sf_zh)
  if (status /= nf90_noerr) then
    sf_zh = 1.0
  endif
  status = nf90_get_att(ncid,varid,"_FillValue",fill_zh)
  if (status /= nf90_noerr) then
    fill_zh = -32768.0
  endif

#ifdef SINGLE_LETKF
  call ncio_read_gattr_r8(ncid, "site_positions_center_latitude",  radar_lon_r8)
  call ncio_read_gattr_r8(ncid, "site_positions_center_longitude", radar_lat_r8)
 radar_lon=real(radar_lon_r8)
 radar_lat=real(radar_lat_r8)
#else
  call ncio_read_gattr_r8(ncid, "site_positions_center_latitude",  radar_lon)
  call ncio_read_gattr_r8(ncid, "site_positions_center_longitude", radar_lat)
#endif

  call ncio_close(ncid)


  obs%nobs = 0

  obs%meta(1) = real(radar_lon, kind=r_size)
  obs%meta(2) = real(radar_lat, kind=r_size)
  obs%meta(3) = 82.0_r_size

  ! count number of obs 
  do k = 1, zdim, RADAR_THIN_VERT
    do j = 1, ydim, RADAR_THIN_HORI
      do i = 1, xdim, RADAR_THIN_HORI
        if (zh3d(i,j,k) /= fill_zh .and. vr3d(i,j,k) /= fill_vr)then
          obs%nobs = obs%nobs + 2
        endif
      enddo ! i
    enddo ! j
  enddo ! k

  call obs_info_allocate(obs, extended=.true.)

  ! substitute obs info
  n = 0
  nobs_vr = 0
  nobs_zh = 0
  min_obs_zh = huge(real(1.0,r_size))
  max_obs_zh = -huge(real(1.0,r_size))
  min_obs_vr = huge(real(1.0,r_size))
  max_obs_vr = -huge(real(1.0,r_size))

  do k = 1, zdim, RADAR_THIN_VERT
    do j = 1, ydim, RADAR_THIN_HORI
      do i = 1, xdim, RADAR_THIN_HORI
        if (zh3d(i,j,k) /= fill_zh .and. vr3d(i,j,k) /= fill_vr)then
          ! zh
          n = n + 1
          obs%elm(n) = id_radar_ref_obs
          obs%lon(n) = real(lon1d(i), kind=r_size)
          obs%lat(n) = real(lat1d(j), kind=r_size)
          obs%lev(n) = real(z1d(k), kind=r_size)
          obs%dat(n) = real(zh3d(i,j,k) * sf_zh, kind=r_size)
          ! Add RADAR_BIAS_CONST_DBZ in dBZ
          if ( RADAR_BIAS_COR_RAIN .and. obs%dat(n) > MIN_RADAR_REF ) then
            obs%dat(n) = real(zh3d(i,j,k) * sf_zh, kind=r_size) * RADAR_BIAS_RAIN_CONST
          elseif ( RADAR_BIAS_COR_CLR .and. obs%dat(n) < MIN_RADAR_REF ) then
            obs%dat(n) = real(zh3d(i,j,k) * sf_zh, kind=r_size) * RADAR_BIAS_CLR_CONST
          endif
          obs%err(n) = OBSERR_RADAR_REF
          obs%typ(n) = 22
          obs%dif(n) = 0.0_r_size
          nobs_zh = nobs_zh + 1

          if (zh3d(i,j,k) > max_obs_zh) max_obs_zh = zh3d(i,j,k)
          if (zh3d(i,j,k) < min_obs_zh) min_obs_zh = zh3d(i,j,k)

          ! vr
          n = n + 1
          obs%elm(n) = id_radar_vr_obs
          obs%lon(n) = real(lon1d(i), kind=r_size)
          obs%lat(n) = real(lat1d(j), kind=r_size)
          obs%lev(n) = real(z1d(k), kind=r_size)
          obs%dat(n) = real(vr3d(i,j,k) * sf_vr, kind=r_size)
          obs%err(n) = OBSERR_RADAR_VR
          obs%typ(n) = 22
          obs%dif(n) = 0.0_r_size
          nobs_vr = nobs_vr + 1

          if (vr3d(i,j,k) > max_obs_vr) max_obs_vr = vr3d(i,j,k)
          if (vr3d(i,j,k) < min_obs_vr) min_obs_vr = vr3d(i,j,k)
        endif
      enddo ! i
    enddo ! j
  enddo ! k

  write (6, *) "Reflectivity obs. range = ", min_obs_zh * sf_zh, &
               " to ", max_obs_zh * sf_zh
  write (6, *) "Radial vel. obs. range  = ", min_obs_vr * sf_vr, &
               " to ", max_obs_vr * sf_vr
  write (6, *) "zh: ", nobs_zh, ", vr: ", nobs_vr

  deallocate(vr3d, zh3d)
  deallocate(lon1d, lat1d, z1d)

  call mpi_timer('read_obs_radar_jrc:read_obs_radar_jrc:', 2)

  return
end subroutine read_obs_radar_jrc

subroutine pawr_toshiba_hd_mpi(lon0, lat0, z0, missing, range_res, na, nr, ne, &
                            AZDIM, ELDIM, n_type, RDIM, az, el, rtdat, utime_obs)
  implicit none

  real(r_size), intent(inout) :: lon0, lat0, z0, missing
  integer, intent(inout) :: range_res, na, nr, ne
  integer, intent(in) :: AZDIM, ELDIM, n_type, RDIM
  real(r_sngl), intent(inout) :: az(AZDIM, ELDIM, n_type), el(AZDIM, ELDIM, n_type)
  real(r_sngl), intent(inout) :: rtdat(RDIM, AZDIM, ELDIM, n_type)
  integer, intent(inout) :: utime_obs(6)

  integer :: ierr

  if ( nprocs_o < 2 ) return

  call MPI_BCAST(lon0, 1, MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(lat0, 1, MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(z0, 1, MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(missing, 1, MPI_r_size, 0, MPI_COMM_o, ierr)

  call MPI_BCAST(range_res, 1, MPI_INTEGER, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(na, 1, MPI_INTEGER, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(nr, 1, MPI_INTEGER, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(ne, 1, MPI_INTEGER, 0, MPI_COMM_o, ierr)

  call MPI_BCAST(az,    AZDIM*ELDIM*n_type, MPI_REAL, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(el,    AZDIM*ELDIM*n_type, MPI_REAL, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(rtdat, RDIM*AZDIM*ELDIM*n_type, MPI_REAL, 0, MPI_COMM_o, ierr)

  call MPI_BCAST(utime_obs, 6, MPI_INTEGER, 0, MPI_COMM_o, ierr)

  return
end subroutine pawr_toshiba_hd_mpi

subroutine jst2utc(jyear, jmonth, jday, jhour, jminute, jsecond, jtime_ms, utime)
  use scale_calendar, only: &
      CALENDAR_date2daysec, &
      CALENDAR_daysec2date, &
      CALENDAR_adjust_daysec
  use scale_precision, only: DP
  implicit none

  integer, intent(in) :: jyear, jmonth, jday
  integer, intent(in) :: jhour, jminute, jsecond
  real(DP) :: jtime_ms
  integer, intent(out) :: utime(6)
  integer :: jtime(6)
  integer :: absday
  real(DP) :: abssec, utime_ms

  jtime(1) = jyear
  jtime(2) = jmonth
  jtime(3) = jday
  jtime(4) = jhour
  jtime(5) = jminute
  jtime(6) = jsecond

  call CALENDAR_date2daysec( absday,       & ! [OUT]
                             abssec,       & ! [OUT]
                             jtime,        & ! [IN]
                             jtime_ms,     & ! [IN]
                             0             ) ! [IN]

  abssec = abssec - real(3600*9, kind=DP)

  call CALENDAR_adjust_daysec( absday,   & ! [INOUT]
                               abssec )    ! [INOUT]

  call CALENDAR_daysec2date( utime,   & ! [OUT]
                             utime_ms, & ! [OUT]
                             absday,      & ! [IN]
                             abssec,      & ! [IN]
                             0            ) ! [IN]

  return
end subroutine jst2utc


subroutine get_timelabel_obs
  use scale_time, only: &
    TIME_gettimelabel
  implicit none

!  if (OBS_POSTFIX_TIMELABEL) then
    timelabel_obs = '_???????????????????.dat'
    call TIME_gettimelabel(timelabel_obs(2:20))
!  else
!    timelabel_obs = ''
!  end if
  if (LOG_LEVEL >= 4) then
    write (6, '(2A)') 'Timelabel for observation files: ', timelabel_obs
  end if

end subroutine get_timelabel_obs

!=======================================================================
end module radar_obs
