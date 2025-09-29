module common_mtx
!=======================================================================
!
! [PURPOSE:] Matrix Functions
!
! [CREATED:] 07/20/2004 Takemasa Miyoshi
! [UPDATED:] 10/16/2004 Takemasa Miyoshi
!
! [PUBLIC:]
!   mtx_eigen  : eigenvalue decomposition
!   mtx_inv    : real symmetric matrix inverse
!   mtx_sqrt   : real symmetric matrix square root
!
! [REFERENCES:]
!    Core subroutines are adapted from netlib.org
!
! [HISTORY:]
!  07/20/2003 Takemasa Miyoshi  Created at University of Maryland, College Park
!
!=======================================================================
  use common

  implicit none

  integer, private :: lda
  integer, private :: lwork
  integer, private :: liwork

  public :: mtx_setup
  public :: mtx_eigen !, mtx_inv, mtx_inv_rg

contains
  !-----------------------------------------------------------------------------
  subroutine mtx_setup(n)
    implicit none
    !---------------------------------------------------------------------------
    integer, intent(in) :: n

    lda = n

    lwork  = 2*n*n + 6*n + 1
    liwork = 5*n + 3
  
    return
  end subroutine mtx_setup
  !-----------------------------------------------------------------------------
  !  eigenvalue decomposition using subroutine dsyevd
  !    input
  !      integer   :: n          : dimension of matrix
  !      real(4/8) :: a(n,n)     : input matrix
  !    output
  !      real(4/8) :: eival(n)   : eiven values in decending order. i.e.
  !      eival(1) is the largest
  !      real(4/8) :: eivec(n,n) : eiven vectors
  !    work
  !      real(4/8) :: work(lwork)   : working array, size is lwork  = 2*n*n +
  !      6*n + 1
  !      integer   :: iwork(liwork) : working array, size is liwork = 5*n + 3
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

    integer :: nrank_eff

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

                b(i,j) = 0.5_r_size * ( a(i,j) + a(j,i) )
             enddo
          enddo
       enddo
    enddo

#ifdef SINGLE_EVP
    call ssyevd ("V","L",n,b,lda,w,work,lwork,iwork,liwork,ierr)
#else
    call dsyevd ("V","L",n,b,lda,w,work,lwork,iwork,liwork,ierr)
#endif

    if ( ierr /= 0 ) then
       write(*,*) 'xxx [mtx_eigen/LAPACK] LAPACK/DSYEVD error code is ', ierr,'! STOP.'
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

                   b(i,j) = 0.5_r_size * ( a(i,j) + a(j,i) )
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
          write(*,*) 'xxx [mtx_eigen/LAPACK] LAPACK/SYEV error code is ', ierr,'! STOP.'
          stop
       endif
    endif

    if ( w(1) < 0.0_r_size ) then
       eival_inc = w(n) / ( 1.E+5_r_size - 1.0_r_size )
       do i = 1, n
          w(i) = w(i) + eival_inc
       enddo
    elseif( w(n)/1.E+5_r_size > w(1) ) then
       eival_inc = ( w(n) - w(1)*1.E+5_r_size ) / ( 1.E+5_r_size - 1.0_r_size )
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
    if ( eival(n) <= 0.0_r_size ) then
       write(*,*) 'xxx [mtx_eigen/LAPACK] All eigenvalues are below 0! STOP.'
       stop
    endif

    nrank_eff = n
    if( eival(n) > 0.0_r_size ) then
      do i = 1, n
        if( eival(i) < abs(eival(n))*sqrt(epsilon(eival)) ) then
          nrank_eff = nrank_eff - 1
          eival(i) = 0.0_r_size
          eivec(:,i) = 0.0_r_size
        endif
      enddo
    else
      write(6,'(A)') '!!! ERROR (mtx_eigen): All Eigenvalues are below 0'
      stop 2
    endif

    return
  end subroutine mtx_eigen

!!=======================================================================
!!  Real symmetric matrix inversion using subroutine dspdi
!!    INPUT
!!      INTEGER :: n               : dimension of matrix
!!      REAL(r_size) :: a(n,n)     : input matrix (real symmetric)
!!    OUTPUT
!!      REAL(r_size) :: ainv(n,n)  : inverse of a
!!=======================================================================
!SUBROUTINE mtx_inv(n,a,ainv)
!  IMPLICIT NONE
!
!  INTEGER,INTENT(IN) :: n
!  REAL(r_size),INTENT(IN) :: a(1:n,1:n)
!  REAL(r_size),INTENT(OUT) :: ainv(1:n,1:n)
!
!  REAL(r_dble) :: acmp(n*(n+1)/2)
!  REAL(r_dble) :: det(2)
!  REAL(r_dble) :: work(n)
!  INTEGER :: kpvt(n)
!  INTEGER :: inert(3)
!  INTEGER :: info
!  INTEGER :: i,j,k
!
!  IF(n==1) THEN
!    ainv(1,1) = 1.0d0 / a(1,1)
!  ELSE
!
!!-----------------------------------------------------------------------
!!  Packed form of matrix
!!-----------------------------------------------------------------------
!  k=0
!  DO j=1,n
!    DO i=1,j
!      k = k+1
!      acmp(k) = a(i,j)
!    END DO
!  END DO
!!-----------------------------------------------------------------------
!!  dspfa
!!-----------------------------------------------------------------------
!  CALL dspfa(acmp,n,kpvt,info)
!  IF(info /= 0) THEN
!    WRITE(6,'(A,I4)') '!!! ERROR (mtx_inv): dspfa error code is ',info
!    STOP 3
!  END IF
!!-----------------------------------------------------------------------
!!  dspdi
!!-----------------------------------------------------------------------
!  CALL dspdi(acmp,n,kpvt,det,inert,work,001)
!!-----------------------------------------------------------------------
!!  unpack matrix
!!-----------------------------------------------------------------------
!  k=0
!  DO j=1,n
!    DO i=1,j
!      k = k+1
!      ainv(i,j) = acmp(k)
!    END DO
!  END DO
!
!  DO j=1,n
!    DO i=j+1,n
!      ainv(i,j) = ainv(j,i)
!    END DO
!  END DO
!
!  END IF
!
!  RETURN
!END SUBROUTINE mtx_inv
!!=======================================================================
!!  Compute square root of real symmetric matrix
!!    INPUT
!!      INTEGER :: n                : dimension of matrix
!!      REAL(r_size) :: a(n,n)      : input matrix (real symmetric)
!!    OUTPUT
!!      REAL(r_size) :: a_sqrt(n,n) : square root of a
!!=======================================================================
!SUBROUTINE mtx_sqrt(n,a,a_sqrt)
!  IMPLICIT NONE
!
!  INTEGER,INTENT(IN) :: n
!  REAL(r_size),INTENT(IN) :: a(1:n,1:n)
!  REAL(r_size),INTENT(OUT) :: a_sqrt(1:n,1:n)
!
!  REAL(r_size) :: eival(n)   ! holds eivenvalue of a
!  REAL(r_size) :: eivec(n,n) ! holds eivenvector of a
!  REAL(r_size) :: wk(n,n)
!  INTEGER :: i,j,k
!
!  CALL mtx_eigen(1,n,a,eival,eivec,i)
!
!  DO i=1,n
!    wk(:,i) = eivec(:,i) * SQRT( eival(i) )
!  END DO
!
!!  a_sqrt = matmul(wk,transpose(eivec))
!  DO j=1,n
!    DO i=1,n
!      a_sqrt(i,j) = wk(i,1)*eivec(j,1)
!      DO k=2,n
!        a_sqrt(i,j) = a_sqrt(i,j) + wk(i,k)*eivec(j,k)
!      END DO
!    END DO
!  END DO
!
!  RETURN
!END SUBROUTINE mtx_sqrt

end module common_mtx
