program main
implicit real(a-h,o-z)

integer,parameter::nx=40,ny=40,nz=10

real(4),parameter::vlonl=134.4,vlonr=136.0
real(4),parameter::vlatl=34.0,vlatr=35.2
real(4),parameter::vzb=100.0,vzt=1000.0 !!! hPa

real(4)::vlons(nx)=(/( vlonl + (vlonr-vlonl) * (real(ix)-0.5)/real(nx) , ix=1,nx)/)
real(4)::vlats(ny)=(/( vlatl + (vlatr-vlatl) * (real(iy)-0.5)/real(ny) , iy=1,ny)/)
real(4)::vzs(nz)=(/( vzb + (vzt-vzb) * (real(iz)-0.5)/real(nz) , iz=1,nz)/)

integer,parameter::nelm=4
integer,parameter::elms(nelm)=(/2819,2820,3073,3331/) !! U,V,T,RH 
real(4),parameter::errs(nelm)=(/3.0,3.0,1.0,10.0/)    !! U,V,T,RH(%) 

real(4)::wk(8)
character(len=200)::cfile

cfile="test_obs_3d_sonde.dat"

open (21, file=trim(cfile), form='unformatted', access='sequential', convert='big_endian')

do iy=1,ny
do ix=1,nx
do iz=1,nz
do ie=1,nelm
  wk(1)=real(elms(ie))  
  wk(2)=vlons(ix)
  wk(3)=vlats(iy)
  wk(4)=vzs(iz)
  wk(5)=10.0  !!! dat
  wk(6)=errs(ie)   !!! err 
  wk(7)=1.0  !!! typ ADPUPA
  wk(8)=0.0   !!! dif
  write(21,iostat=ios) wk(1:8)
end do
end do
end do
end do

close(21)

stop
end program main
