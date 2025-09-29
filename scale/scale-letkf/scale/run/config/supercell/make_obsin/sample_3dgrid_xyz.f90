program main
use common_ncio
implicit real(a-h,o-z)

real(4),allocatable::axlon(:,:),axlat(:,:),axz(:),pres(:,:,:)

integer,parameter::nelm=3
integer,parameter::elms(nelm)=(/2819,2820,3073/) !! U,V,T 
real(4),parameter::errs(nelm)=(/3.0,3.0,1.0/)    !! U,V,T 

integer,parameter::intv_x=4
integer,parameter::intv_y=4
integer,parameter::intv_z=2

real(4)::wk(8)
character(len=200)::cfile
character(len=200)::ncfile_in='history_merge.pe000000.nc'

integer::ncid, vidlon, vidlat,vidz

  call ncio_open( trim(ncfile_in), nf90_nowrite, ncid )
  call ncio_read_dim(ncid,"x",nlon)
  call ncio_read_dim(ncid,"y",nlat)
  call ncio_read_dim(ncid,"z",nlev)

  allocate(pres(nlon,nlat,nlev))
  allocate(axlon(nlon,nlat))
  allocate(axlat(nlon,nlat))
  allocate(axz(nlev))

  call ncio_read_const(ncid, "lon", nlon, nlat, axlon)
  call ncio_read_const(ncid, "lat", nlon, nlat, axlat)
  call ncio_read_const(ncid, "z", nlev, axz)
  call ncio_close( ncid ) 

cfile="test_obs_3d_xyz.dat"

open (21, file=trim(cfile), form='unformatted', access='sequential', convert='big_endian')

do ilon=1,nlon,intv_x
do ilat=1,nlat,intv_y
do ilev=1,nlev,intv_z
do ie=1,nelm
  wk(1)=real(elms(ie))  
  wk(2)=axlon(ilon,ilat)
  wk(3)=axlat(ilon,ilat)
  wk(4)=axz(iz)
  wk(5)=10.0  !!! dat
  wk(6)=errs(ie)   !!! err 
  wk(7)=1.0  !!! typ ADPUPA
  wk(8)=0.0   !!! dif
  write(21,iostat=ios) wk(1:8)
end do
!  write(*,'(F6.1,5F14.4)') wk(1:6)
end do
end do
end do

close(21)

stop
end program main
