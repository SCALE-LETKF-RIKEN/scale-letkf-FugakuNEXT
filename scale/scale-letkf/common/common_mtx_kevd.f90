module common_mtx
  !-----------------------------------------------------------------------------
  !
  ! [purpose:] matrix functions
  !
  ! [created:] 07/20/2004 takemasa miyoshi
  ! [updated:] 10/16/2004 takemasa miyoshi
  !
  ! [public:]
  !   mtx_eigen  : eigenvalue decomposition
  !   mtx_inv    : real symmetric matrix inverse
  !   mtx_sqrt   : real symmetric matrix square root
  !
  ! [references:]
  !    core subroutines are adapted from netlib.org
  !
  ! [history:]
  !  07/20/2003 takemasa miyoshi  created at university of maryland, college park
  !  02/05/2019                   Kevd (special lib for small eigenvalue decomposition) ver.
  !  04/06/2021                   imported from NICAM-LETKF
  !
  !-----------------------------------------------------------------------------

  use common

  implicit none
  private
  !-----------------------------------------------------------------------------
  public :: mtx_setup
  public :: mtx_finalize
  public :: mtx_eigen

  integer, private :: lda
  integer, private :: lwork
  integer, private :: liwork
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine mtx_setup(n)
    implicit none
    integer, intent(in) :: n

    external :: dsyevdt_worksizes, ssyevdt_worksizes
    !---------------------------------------------------------------------------

!    write(6,*)
!    write(6,*) '+++ Module[matrix functions] / Category[common]'
!    write(6,*)
!    write(6,*) '--- use Kevd (optimized DSYEVD) for eigenvalue decomposition'

    lda = (n+15)/16*16
    if( mod(lda,512) == 0 ) lda = lda + 16

#ifdef SINGLE_EVP
    call ssyevdt_worksizes('V', 'U', n, lda, lwork, liwork);
#else
    call dsyevdt_worksizes('V', 'U', n, lda, lwork, liwork);
#endif

    return
  end subroutine mtx_setup

  !-----------------------------------------------------------------------------
  subroutine mtx_finalize
    implicit none
    !---------------------------------------------------------------------------

    ! nothing to do

    return
  end subroutine mtx_finalize

  !-----------------------------------------------------------------------------
  !  eigenvalue decomposition using subroutine dsyevd
  !    input
  !      integer   :: n          : dimension of matrix
  !      real(4/8) :: a(n,n)     : input matrix
  !    output
  !      real(4/8) :: eival(n)   : eiven values in decending order. i.e. eival(1) is the largest
  !      real(4/8) :: eivec(n,n) : eiven vectors
  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine mtx_eigen(n,a,eival,eivec)
    implicit none

    integer, intent(in)  :: n
    real(RP_EVP), intent(in)  :: a    (n,n)
    real(RP_EVP), intent(out) :: eival(n)
    real(RP_EVP), intent(out) :: eivec(n,n)

    real(RP_EVP) :: b    (lda,n)
    real(RP_EVP) :: w    (lda)
    real(RP_EVP) :: work (lwork)
    integer      :: iwork(liwork)

    real(RP_EVP) :: eival_inc

    integer, parameter :: simdlen = 64
    integer :: iblk, jblk
    integer :: imax, jmax
    integer :: ivec, jvec

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    do jblk = 1, n, simdlen
       jmax = min( n-jblk+1, simdlen )
       do iblk = 1, n, simdlen
          imax = min( n-iblk+1, simdlen )

          do jvec = 1, jmax
             j = jblk + jvec - 1
             do ivec = 1, imax
                i = iblk + ivec - 1

                b(i,j) = 0.5_RP_EVP * ( a(i,j) + a(j,i) )
             enddo
          enddo
       enddo
    enddo

#ifdef SINGLE_EVP
    call ssyevdt("V","U",n,b,lda,w,work,lwork,iwork,liwork,ierr)
#else
    call dsyevdt("V","U",n,b,lda,w,work,lwork,iwork,liwork,ierr)
#endif

    if ( ierr /= 0 ) then
       write(*,*) 'xxx [mtx_eigen/Kevd] Kevd/DSYEVDT error code is ', ierr, '! STOP.'
       write(*,*) 'input a'
       do j = 1, n
          write(*,*) j ,a(:,j)
       enddo
       write(*,*) 'output eival'
       write(*,*) w(:)
       write(*,*) 'output eivec'
       do j = 1, n
          write(*,*) j, b(:,j)
       enddo

       ! Temporary treatment because of the instability of SYEVD
       do jblk = 1, n, simdlen
          jmax = min( n-jblk+1, simdlen )
          do iblk = 1, n, simdlen
             imax = min( n-iblk+1, simdlen )
   
             do jvec = 1, jmax
                j = jblk + jvec - 1
                do ivec = 1, imax
                   i = iblk + ivec - 1
   
                   b(i,j) = 0.5_RP_EVP * ( a(i,j) + a(j,i) )
                enddo
             enddo
          enddo
       enddo

#ifdef SINGLE_EVP
       call ssyev ("V","L",n,b,lda,w,work,lwork,ierr)
#else
       call dsyev ("V","L",n,b,lda,w,work,lwork,ierr)
#endif

       if ( ierr /= 0 ) then
          write(*,*) 'xxx [mtx_eigen/LAPACK] LAPACK/SYEV error code is ', ierr, '! STOP.'
          stop
       endif

    endif

    if ( w(1) < 0.0_RP_EVP ) then
       eival_inc = w(n) / ( 1.E+5_RP_EVP - 1.0_RP_EVP )
       do i = 1, n
          w(i) = w(i) + eival_inc
       enddo
    elseif( w(n)/1.E+5_RP_EVP > w(1) ) then
       eival_inc = ( w(n) - w(1)*1.E+5_RP_EVP ) / ( 1.E+5_RP_EVP - 1.0_RP_EVP )
       do i = 1, n
          w(i) = w(i) + eival_inc
       enddo
    endif

    ! reverse order
    do i = 1, n
       eival(i) = w(i) 
    enddo

    do j = 1, n
    do i = 1, n
       eivec(i,j) = b(i,j)
    enddo
    enddo

    ! check zero
    if ( eival(n) <= 0.0_RP_EVP ) then
       write(*,*) 'xxx [mtx_eigen/Kevd] All eigenvalues are below 0! STOP.'
       stop
    endif

    return
  end subroutine mtx_eigen

end module common_mtx
