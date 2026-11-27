! Regression checks for the vertical momentum fluxes of upstream e016ebc73 (#3754).
module test_atmos_dyn_fluxzj
  use scale_precision
  use scale_atmos_grid_cartesC_index
  use scale_atmos_grid_cartesC, only: CDZ => ATMOS_GRID_CARTESC_CDZ
  use scale_atmos_dyn_fvm_flux, only: ATMOS_DYN_FVM_flux_setup
  use dc_test, only: AssertLessThan
  implicit none
  private
  public :: test_atmos_dyn_fluxzj_run
contains
  ! The velocities have opposite signs, non-unit map factors and a prescribed total.
  ! A constant transported value checks every part analytically. An arbitrary value
  ! with zero total velocity must give cancelling parts, including the reduced-order
  ! levels and the full boundary levels. UD1 also has an explicit upwind-side oracle.
  subroutine test_atmos_dyn_fluxzj_run( scheme, tracer_scheme )
    character(len=*), intent(in) :: scheme, tracer_scheme
    character(len=12), parameter :: schemes(9) = &
         (/ 'UD1         ', 'CD2         ', 'UD3         ', 'UD3KOREN1993', 'CD4         ', &
            'UD5         ', 'CD6         ', 'UD7         ', 'CD8         ' /)
    real(RP) :: momz(KA,IA,JA), momx(KA,IA,JA), momy(KA,IA,JA), val(KA,IA,JA)
    real(RP) :: dens(KA,IA,JA), gsqrt(KA,IA,JA), f2hw(KA,IA,JA)
    real(RP) :: j13g(KA,IA,JA), j23g(KA,IA,JA), mapf(IA,JA,2), nd(KA,IA,JA)
    real(RP) :: fl(KA,IA,JA,3), ref(KA,IA,JA,3), dz(KA), fdz(KA-1)
    real(RP) :: c(3), target, fv, fac, ex(3), err, scl, sentinel
    integer :: n, comp, d, g, r, t, ii, iis, iie, k, i, j, p, kend
    logical :: twod, dispatch

    write(*,*) 'Test fused vertical momentum fluxes (#3754)'
    sentinel = -17.0_RP
    dens = 1.0_RP
    gsqrt = 1.0_RP
    f2hw = 0.37_RP
    momx = 3.0_RP
    momy = -2.0_RP
    j13g = 0.3_RP
    j23g = 0.2_RP
    mapf(:,:,1) = 0.8_RP
    mapf(:,:,2) = 1.4_RP
    nd = 0.0_RP
    do g = 1, 2
       do k = 1, KA
          dz(k) = CDZ(KS) * (1.0_RP + 0.03_RP * (g-1))**(k-KS)
       end do
       fdz = 0.5_RP * ( dz(1:KA-1) + dz(2:KA) )
       do comp = 1, 3
       do d = 1, 2
          twod = (d == 2)
          c(2) = merge(0.0_RP, 3.0_RP*0.3_RP*0.8_RP, twod)
          c(3) = -2.0_RP*0.2_RP*1.4_RP
          do n = 1, 9
             ! The dispatcher enforces the halo requirements; direct calls also
             ! exercise higher-order vertical stencils on this grid.
             dispatch = (n <= 5)
             if (dispatch) call ATMOS_DYN_FVM_flux_setup(schemes(n), tracer_scheme)
             do r = 1, 3
                target = real(r-2,RP) * 0.35_RP
                c(1) = target - c(2) - c(3)
                momz = c(1)
                do t = 1, 2
                   val = 1.0_RP
                   if (t == 2) then
                      do j = 1, JA
                      do i = 1, IA
                      do k = 1, KA
                         val(k,i,j) = sin(real(3*k+i+2*j,RP)) + 0.01_RP*k
                      end do
                      end do
                      end do
                   end if
                   fl = sentinel
                   ref = sentinel
                   do ii = 1, merge(IE-IS+1, 1, twod)
                      iis = IS
                      iie = IE
                      if (twod) then
                         iis = IS + ii - 1
                         iie = iis
                      end if
                      call fluxzj_new_run(n, comp, fl(:,:,:,1), fl(:,:,:,2), fl(:,:,:,3), &
                           momz, momx, momy, val, dens, gsqrt, f2hw, j13g, j23g, mapf, nd, &
                           1.0_RP, dz, fdz, twod, iis, iie)
                      if (dispatch) call fluxzj_wrap_run(comp, ref(:,:,:,1), ref(:,:,:,2), ref(:,:,:,3), &
                           momz, momx, momy, val, dens, gsqrt, f2hw, j13g, j23g, mapf, nd, &
                           dz, fdz, twod, iis, iie)
                   end do
                   if (dispatch) call AssertLessThan('fluxZJ dispatch', tiny(1.0_RP), maxval(abs(fl-ref)))
                   if (twod) call AssertLessThan('fluxZJ 2-D leaves J13 unchanged', tiny(1.0_RP), &
                                                 maxval(abs(fl(:,:,:,2)-sentinel)))
                   err = 0.0_RP
                   scl = maxval(abs(c)) * maxval(abs(val))
                   kend = KE
                   if (comp == 1) kend = KE-1
                   do j = JS, JE
                   do i = IS, IE
                   do k = KS-1, kend
                      if (r == 2) then
                         ex = fl(k,i,j,:)
                         if (twod) ex(2) = 0.0_RP
                         err = max(err, abs(sum(ex)))
                      end if
                      if (t /= 1 .and. n /= 1) cycle
                      fac = 1.0_RP
                      fv = 1.0_RP
                      if (t == 2) fv = merge(val(k,i,j), val(k+1,i,j), target >= 0.0_RP)
                      if (comp == 1) then
                         if (k == KS-1 .or. k == KE-1) then
                            fac = 0.5_RP
                            if (t == 2) fv = val(merge(KS, KE-1, k == KS-1),i,j)
                         end if
                      else
                         if (k == KS-1 .or. k == KE) fac = 0.0_RP
                      end if
                      do p = 1, 3
                         if (twod .and. p == 2) cycle
                         ! At exactly zero total velocity either side is valid;
                         ! only cancellation is asserted for a varying value.
                         if (r == 2 .and. t == 2) cycle
                         err = max(err, abs(fl(k,i,j,p) - fac*c(p)*fv))
                      end do
                   end do
                   end do
                   end do
                   call AssertLessThan('fluxZJ analytic parts and cancellation', &
                                        200.0_RP*epsilon(1.0_RP)*scl, err)
                end do
             end do
          end do
       end do
       end do
    end do
    call ATMOS_DYN_FVM_flux_setup(scheme, tracer_scheme)
  end subroutine test_atmos_dyn_fluxzj_run
subroutine fluxzj_new_run( n, comp, flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, &
                           J13G_l, J23G_l, MAPF_l, nd, J33G_l, CDZ_l, FDZ_l, TwoD, IIS, IIE )
  use scale_atmos_dyn_fvm_flux_ud1, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_ud1, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_ud1, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_ud1
  use scale_atmos_dyn_fvm_flux_cd2, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_cd2, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_cd2, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_cd2
  use scale_atmos_dyn_fvm_flux_ud3, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_ud3, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_ud3, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_ud3
  use scale_atmos_dyn_fvm_flux_ud3Koren1993, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_ud3Koren1993, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_ud3Koren1993, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_ud3Koren1993
  use scale_atmos_dyn_fvm_flux_cd4, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_cd4, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_cd4, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_cd4
  use scale_atmos_dyn_fvm_flux_ud5, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_ud5, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_ud5, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_ud5
  use scale_atmos_dyn_fvm_flux_cd6, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_cd6, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_cd6, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_cd6
  use scale_atmos_dyn_fvm_flux_ud7, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_ud7, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_ud7, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_ud7
  use scale_atmos_dyn_fvm_flux_cd8, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW_cd8, &
     ATMOS_DYN_FVM_fluxZJ_UYZ_cd8, &
     ATMOS_DYN_FVM_fluxZJ_XVZ_cd8
  integer,  intent(in)    :: n
  integer,  intent(in)    :: comp
  real(RP), intent(inout) :: flux(KA,IA,JA)
  real(RP), intent(inout) :: fluxJ13(KA,IA,JA)
  real(RP), intent(inout) :: fluxJ23(KA,IA,JA)
  real(RP), intent(in)    :: momz(KA,IA,JA)
  real(RP), intent(in)    :: momx(KA,IA,JA)
  real(RP), intent(in)    :: momy(KA,IA,JA)
  real(RP), intent(in)    :: val(KA,IA,JA)
  real(RP), intent(in)    :: DENS_l(KA,IA,JA)
  real(RP), intent(in)    :: GSQRT_l(KA,IA,JA)
  real(RP), intent(in)    :: F2HW_l(KA,IA,JA)
  real(RP), intent(in)    :: J13G_l(KA,IA,JA)
  real(RP), intent(in)    :: J23G_l(KA,IA,JA)
  real(RP), intent(in)    :: MAPF_l(IA,JA,2)
  real(RP), intent(in)    :: nd(KA,IA,JA)
  real(RP), intent(in)    :: J33G_l
  real(RP), intent(in)    :: CDZ_l(KA)
  real(RP), intent(in)    :: FDZ_l(KA-1)
  logical,  intent(in)    :: TwoD
  integer,  intent(in)    :: IIS, IIE

  integer :: sel

  sel = n * 3 + comp
  !$acc data copyin(momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J13G_l, J23G_l, MAPF_l, nd, CDZ_l, FDZ_l) &
  !$acc      copy(flux, fluxJ13, fluxJ23)
  select case( sel )
  case( 4 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_ud1( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 5 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_ud1( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 6 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_ud1( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 7 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_cd2( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 8 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_cd2( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 9 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_cd2( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 10 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_ud3( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 11 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_ud3( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 12 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_ud3( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 13 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_ud3Koren1993( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 14 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_ud3Koren1993( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 15 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_ud3Koren1993( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 16 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_cd4( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 17 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_cd4( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 18 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_cd4( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 19 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_ud5( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 20 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_ud5( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 21 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_ud5( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 22 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_cd6( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 23 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_cd6( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 24 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_cd6( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 25 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_ud7( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 26 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_ud7( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 27 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_ud7( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 28 )
     call ATMOS_DYN_FVM_fluxZJ_XYW_cd8( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 29 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ_cd8( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 30 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ_cd8( &
          flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J33G_l, J13G_l, J23G_l, MAPF_l, nd, &
          CDZ_l, TwoD, IIS, IIE, JS, JE )
  case default
     write(*,*) "fluxzj_new_run: invalid scheme: ", n
     call abort
  end select
  !$acc wait
  !$acc end data

end subroutine fluxzj_new_run

subroutine fluxzj_wrap_run( comp, flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, &
                            J13G_l, J23G_l, MAPF_l, nd, CDZ_l, FDZ_l, TwoD, IIS, IIE )
  use scale_atmos_dyn_fvm_flux, only: &
     ATMOS_DYN_FVM_fluxZJ_XYW, &
     ATMOS_DYN_FVM_fluxZJ_UYZ, &
     ATMOS_DYN_FVM_fluxZJ_XVZ
  integer,  intent(in)    :: comp
  real(RP), intent(inout) :: flux(KA,IA,JA)
  real(RP), intent(inout) :: fluxJ13(KA,IA,JA)
  real(RP), intent(inout) :: fluxJ23(KA,IA,JA)
  real(RP), intent(in)    :: momz(KA,IA,JA)
  real(RP), intent(in)    :: momx(KA,IA,JA)
  real(RP), intent(in)    :: momy(KA,IA,JA)
  real(RP), intent(in)    :: val(KA,IA,JA)
  real(RP), intent(in)    :: DENS_l(KA,IA,JA)
  real(RP), intent(in)    :: GSQRT_l(KA,IA,JA)
  real(RP), intent(in)    :: F2HW_l(KA,IA,JA)
  real(RP), intent(in)    :: J13G_l(KA,IA,JA)
  real(RP), intent(in)    :: J23G_l(KA,IA,JA)
  real(RP), intent(in)    :: MAPF_l(IA,JA,2)
  real(RP), intent(in)    :: nd(KA,IA,JA)
  real(RP), intent(in)    :: CDZ_l(KA)
  real(RP), intent(in)    :: FDZ_l(KA-1)
  logical,  intent(in)    :: TwoD
  integer,  intent(in)    :: IIS, IIE

  !$acc data copyin(momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, J13G_l, J23G_l, MAPF_l, nd, CDZ_l, FDZ_l) &
  !$acc      copy(flux, fluxJ13, fluxJ23)
  select case( comp )
  case( 1 )
     call ATMOS_DYN_FVM_fluxZJ_XYW( flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, &
                                    1.0_RP, J13G_l, J23G_l, MAPF_l, nd, CDZ_l, FDZ_l, 1.0_RP, TwoD, IIS, IIE, JS, JE )
  case( 2 )
     call ATMOS_DYN_FVM_fluxZJ_UYZ( flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, &
                                    1.0_RP, J13G_l, J23G_l, MAPF_l, nd, CDZ_l, TwoD, IIS, IIE, JS, JE )
  case( 3 )
     call ATMOS_DYN_FVM_fluxZJ_XVZ( flux, fluxJ13, fluxJ23, momz, momx, momy, val, DENS_l, GSQRT_l, F2HW_l, &
                                    1.0_RP, J13G_l, J23G_l, MAPF_l, nd, CDZ_l, TwoD, IIS, IIE, JS, JE )
  end select
  !$acc wait
  !$acc end data

end subroutine fluxzj_wrap_run
end module test_atmos_dyn_fluxzj
