program main
use common_ncio
implicit real(a-h,o-z)

real(4),allocatable::axlon(:,:),axlat(:,:),axz(:),pres(:,:,:)

integer,parameter::nelm=2
integer,parameter::elms(nelm)=(/4001,4004/) !! ze, ze(zero)
real(4),parameter::errs(nelm)=(/1.0,1.0/)     

real(4),parameter::er=6400.0e3
real(4),parameter::drad=3.141592/180.0

integer,parameter::intv_x=1
integer,parameter::intv_y=1
integer,parameter::intv_z=1

real(4)::wk(8)
character(len=200)::cfile
character(len=200)::ncfile_in='history_merge.pe000000.nc'

integer::ncid, vidlon, vidlat,vidz

real(4),parameter::radar_lon = 135.232
real(4),parameter::radar_lat = 34.662
real(4),parameter::radar_z   = 0.0  

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
  call ncio_read(ncid,"PRES",nlon,nlat,nlev,1,pres)
  call ncio_close( ncid ) 

cfile="test_obs_3d_radar.dat"

open (21, file=trim(cfile), form='unformatted', access='sequential', convert='big_endian')

! Radar header
write(21) radar_lon
write(21) radar_lat
write(21) radar_z

icount=0
do ilon=1,nlon,intv_x
do ilat=1,nlat,intv_y
do ilev=1,nlev,intv_z

!  vdist=  (er*drad*(axlat(ilon,ilat)-radar_lat))**2  &
!        + (er*cos(drad*radar_lat)*drad*(axlon(ilon,ilat)-radar_lon))**2 & 
!        + (axz(ilev)-radar_z)**2
!  vdist=sqrt(vdist)
!  if(vdist < vdist_limit)then
    do ie=1,nelm
      wk(1)=elms(ie)!!! elm radar ref
      wk(2)=axlon(ilon,ilat)
      wk(3)=axlat(ilon,ilat)
      wk(4)=axz(ilev)
      wk(5)=10.0  !!! dat (dummy)
      wk(6)=errs(ie)   !!! err
      wk(7)=22.0  !!! typ PHARAD
      wk(8)=0.0   !!! dif
      write(21,iostat=ios) wk(1:7)
    end do
!  else 
!  end if

end do
end do
end do


close(21)

stop
end program main
