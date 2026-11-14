module test_atmos_phy_mp_common
  use scale_precision, only: RP
  use dc_test, only: AssertEqual
  use scale_atmos_phy_mp_common, only: ATMOS_PHY_MP_precipitation_momentum
#ifdef _OPENACC
  use scale_atmos_phy_mp_common, only: ATMOS_PHY_MP_precipitation_momentum_3D
#endif
  implicit none
  private
  public :: test_atmos_phy_mp_common_run
contains
  !> Momentum flux divergence: column and grid entry points, top and halo cells.
  subroutine test_atmos_phy_mp_common_run
    implicit none
    integer, parameter :: KA=14, KS=3, KE=12, IA=7, JA=6
    real(RP), parameter :: sentinel=-99.0_RP, tol=64.0_RP*epsilon(1.0_RP)
    real(RP) :: dens(KA,IA,JA), momz(KA,IA,JA), u(KA,IA,JA), v(KA,IA,JA)
    real(RP) :: mflx(KA,IA,JA), rcdz(KA,IA,JA), rfdz(KA,IA,JA)
    real(RP) :: column(KA,IA,JA,3), grid(KA,IA,JA,3)
    real(RP) :: flux(KA,3), ref, scale, error
    integer :: k, i, j, n, mode, ie, je
    !-------------------------------------------------------------------------
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       dens(k,i,j) = 1.2_RP * exp(-0.02_RP*real(k,RP)) + 0.001_RP*real(i+j,RP)
       momz(k,i,j) = dens(k,i,j) * sin(0.1_RP*real(k+i+j,RP))
       u(k,i,j) = 10.0_RP * sin(0.1_RP*real(k+i,RP))
       v(k,i,j) = 10.0_RP * cos(0.1_RP*real(k+j,RP))
       mflx(k,i,j) = -1.0e-4_RP * (1.0_RP + 0.01_RP*real(mod(k+i+2*j,23),RP))
       rcdz(k,i,j) = 1.0_RP / (100.0_RP + real(k,RP))
       rfdz(k,i,j) = 1.0_RP / (105.0_RP + real(k,RP))
    end do
    end do
    end do
    ! No vertical momentum or precipitation flux through the model top.
    momz(KE,:,:) = 0.0_RP
    mflx(KE:KA,:,:) = 0.0_RP

    do mode = 1, 2
       ie = IA-1
       je = JA-1
       if ( mode == 2 ) then
          ie = 2
          je = 2
       end if
       column = sentinel
       grid = sentinel
       !$acc data copyin(dens, momz, u, v, mflx, rcdz, rfdz) copy(column, grid)
       !$acc parallel loop collapse(2) gang async present(dens, momz, u, v, mflx, rcdz, rfdz, column)
       do j = 2, je
       do i = 2, ie
          call ATMOS_PHY_MP_precipitation_momentum( &
               KA, KS, KE, dens(:,i,j), momz(:,i,j), u(:,i,j), v(:,i,j), &
               mflx(:,i,j), rcdz(:,i,j), rfdz(:,i,j), &
               column(:,i,j,1), column(:,i,j,2), column(:,i,j,3) )
       end do
       end do
       !$acc end parallel loop
#ifdef _OPENACC
       call ATMOS_PHY_MP_precipitation_momentum_3D( &
            KA, KS, KE, IA, 2, ie, JA, 2, je, dens, momz, u, v, mflx, rcdz, rfdz, &
            grid(:,:,:,1), grid(:,:,:,2), grid(:,:,:,3) )
#endif
       !$acc wait
       !$acc end data

       error = 0.0_RP
       do j = 2, je
       do i = 2, ie
          ! Independent reference: form face fluxes first, then their divergence.
          flux = 0.0_RP
          do k = KS, KE-1
             flux(k,1) = (mflx(k,i,j)+mflx(k-1,i,j))*momz(k,i,j) &
                       / (dens(k+1,i,j)+dens(k,i,j))
          end do
          do k = KS-1, KE-1
             flux(k,2) = mflx(k,i,j)*u(k+1,i,j)
             flux(k,3) = mflx(k,i,j)*v(k+1,i,j)
          end do
          do n = 1, 3
          do k = KS, KE
             if ( n == 1 ) then
                if ( k == KE ) then
                   call AssertEqual('momentum top: column', .true., column(k,i,j,n) == 0.0_RP)
#ifdef _OPENACC
                   call AssertEqual('momentum top: grid', .true., grid(k,i,j,n) == 0.0_RP)
#endif
                   cycle
                end if
                ref = -(flux(k+1,n)-flux(k,n))*rfdz(k,i,j)
                scale = (abs(flux(k+1,n))+abs(flux(k,n)))*rfdz(k,i,j)
             else
                ref = -(flux(k,n)-flux(k-1,n))*rcdz(k,i,j)
                scale = (abs(flux(k,n))+abs(flux(k-1,n)))*rcdz(k,i,j)
             end if
             ! Scale by both face fluxes, allowing FMA changes under cancellation.
             error = max(error, abs(column(k,i,j,n)-ref)/max(scale,tiny(1.0_RP)))
#ifdef _OPENACC
             error = max(error, abs(grid(k,i,j,n)-ref)/max(scale,tiny(1.0_RP)))
#endif
          end do
          end do
       end do
       end do
       call AssertEqual('momentum flux divergence within roundoff', .true., error <= tol)
#ifdef _OPENACC
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          if ( i >= 2 .and. i <= ie .and. j >= 2 .and. j <= je .and. k >= KS .and. k <= KE ) cycle
          call AssertEqual('momentum grid preserves halo', .true., all(grid(k,i,j,:) == sentinel))
       end do
       end do
       end do
#endif
    end do
  end subroutine test_atmos_phy_mp_common_run
end module test_atmos_phy_mp_common
