MODULE common_letkf
!=======================================================================
!
! [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
! [HISTORY:]
!  01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mtx

  IMPLICIT NONE

  PUBLIC

CONTAINS
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     ne               : ensemble size                                           !GYL
!     nobs             : array size, but only first nobsl elements are used
!     nobsl            : total number of observation assimilated at the point
!     hdxb(nobs,ne)    : obs operator times fcst ens perturbations
!     rdiag(nobs)      : observation error variance
!     rloc(nobs)       : localization weighting function
!     dep(nobs)        : observation departure (yo-Hxb)
!     parm_infl        : covariance inflation parameter
!     rdiag_wloc       : (optional) flag indicating that rdiag = rdiag / rloc    !GYL
!     infl_update      : (optional) flag to return updated inflation parameter   !GYL
!     depd(nobs)       : observation departure (yo-Hxb) for deterministic run    !GYL
!   OUTPUT
!     parm_infl        : updated covariance inflation parameter
!     trans(ne,ne)     : transformation matrix
!     transm(ne)       : (optional) transformation matrix mean                   !GYL
!     pao(ne,ne)       : (optional) analysis covariance matrix in ensemble space !GYL
!     transmd(ne)      : (optional) transformation matrix mean for deterministic run !GYL
!=======================================================================
!OCL SERIAL
SUBROUTINE letkf_core(ne,nobs,nobsl,hdxb,rdiag,rloc,dep,parm_infl,trans,transm,pao,rdiag_wloc,infl_update,depd,transmd)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ne                      !GYL
  INTEGER,INTENT(IN) :: nobs
  INTEGER,INTENT(IN) :: nobsl
  REAL(r_size),INTENT(IN) :: hdxb(1:nobs,1:ne)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobs)
  REAL(r_size),INTENT(IN) :: rloc(1:nobs)
  REAL(r_size),INTENT(IN) :: dep(1:nobs)
  REAL(r_size),INTENT(INOUT) :: parm_infl
  REAL(r_size),INTENT(OUT) :: trans(ne,ne)
  REAL(r_size),INTENT(OUT),OPTIONAL :: transm(ne)
  REAL(r_size),INTENT(OUT),OPTIONAL :: pao(ne,ne)
  LOGICAL,INTENT(IN),OPTIONAL :: rdiag_wloc     !GYL
  LOGICAL,INTENT(IN),OPTIONAL :: infl_update    !GYL
  REAL(r_size),INTENT(IN),OPTIONAL :: depd(1:nobs) !GYL
  REAL(r_size),INTENT(OUT),OPTIONAL :: transmd(ne)  !GYL

  REAL(r_size) :: hdxb_rinv(nobsl,ne)
  REAL(r_size) :: eivec(ne,ne)
  REAL(r_size) :: eival(ne)
  REAL(r_size) :: pa(ne,ne)
  REAL(r_size) :: work1(ne,ne)
  REAL(r_size) :: work2(ne)
  REAL(r_size) :: work2d(ne)
  REAL(r_size) :: work3(ne)
  REAL(r_size) :: rho
  REAL(r_size) :: parm(4),sigma_o,gain
  REAL(r_size),PARAMETER :: sigma_b = 0.04d0 !error stdev of parm_infl
  LOGICAL :: rdiag_wloc_
  LOGICAL :: infl_update_
  INTEGER :: i,j,k

  real(RP_EVP) :: work1_evp(ne,ne)
  real(RP_EVP) :: eivec_evp(ne,ne)
  real(RP_EVP) :: eival_evp(ne)

  rdiag_wloc_ = .FALSE.                               !GYL
  IF(present(rdiag_wloc)) rdiag_wloc_ = rdiag_wloc    !GYL
  infl_update_ = .FALSE.                              !GYL
  IF(present(infl_update)) infl_update_ = infl_update !GYL

  IF(nobsl == 0) THEN
    trans = 0.0d0
    DO i=1,ne
      trans(i,i) = SQRT(parm_infl)
    END DO
    IF (PRESENT(transm)) THEN   !GYL
      transm = 0.0d0            !GYL
    END IF                      !GYL
    IF (PRESENT(transmd)) THEN  !GYL
      transmd = 0.0d0           !GYL
    END IF                      !GYL
    IF (PRESENT(pao)) THEN                        !GYL
      pao = 0.0d0                                 !GYL
      DO i=1,ne                                   !GYL
        pao(i,i) = parm_infl / REAL(ne-1,r_size)  !GYL
      END DO                                      !GYL
    END IF                                        !GYL
    RETURN
  END IF
!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  IF(rdiag_wloc_) THEN                            !GYL
    DO j=1,ne                                     !GYL
      DO i=1,nobsl                                !GYL
        hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i)     !GYL
      END DO                                      !GYL
    END DO                                        !GYL
  ELSE                                            !GYL
    DO j=1,ne
      DO i=1,nobsl
        hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
      END DO
    END DO
  END IF                                          !GYL
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
#ifdef SINGLE_LETKF
  CALL sgemm('t','n',ne,ne,nobsl,1.0e0,hdxb_rinv,nobsl,hdxb(1:nobsl,:),&
    & nobsl,0.0e0,work1,ne)
#else
  CALL dgemm('t','n',ne,ne,nobsl,1.0d0,hdxb_rinv,nobsl,hdxb(1:nobsl,:),&
    & nobsl,0.0d0,work1,ne)
#endif
!  DO j=1,ne
!    DO i=1,ne
!      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
!      DO k=2,nobsl
!        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
  rho = 1.0d0 / parm_infl
  DO i=1,ne
    work1(i,i) = work1(i,i) + REAL(ne-1,r_size) * rho
  END DO
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
  work1_evp = real( work1, kind=RP_EVP )
#ifdef KEVD
  call mtx_eigen( ne, work1_evp, eival_evp, eivec_evp )
#else
  call mtx_eigen( ne, work1_evp, eival_evp, eivec_evp )
#endif
  eival = real( eival_evp, kind=r_size )
  eivec = real( eivec_evp, kind=r_size )
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
  DO j=1,ne
    DO i=1,ne
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO
#ifdef SINGLE_LETKF
  CALL sgemm('n','t',ne,ne,ne,1.0e0,work1,ne,eivec,&
    & ne,0.0e0,pa,ne)
#else
  CALL dgemm('n','t',ne,ne,ne,1.0d0,work1,ne,eivec,&
    & ne,0.0d0,pa,ne)
#endif
!  DO j=1,ne
!    DO i=1,ne
!      pa(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,ne
!        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  hdxb_rinv^T dep
!-----------------------------------------------------------------------
#ifdef SINGLE_LETKF
  call sgemv('t',nobsl,ne,1.0e0,hdxb_rinv,nobsl,dep,1,0.0e0,work2,1)
#else
  call dgemv('t',nobsl,ne,1.0d0,hdxb_rinv,nobsl,dep,1,0.0d0,work2,1)
#endif
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
#ifdef SINGLE_LETKF
  call sgemv('n',ne,ne,1.0e0,pa,ne,work2,1,0.0e0,work3,1)
#else
  call dgemv('n',ne,ne,1.0d0,pa,ne,work2,1,0.0d0,work3,1)
#endif

  IF (PRESENT(depd) .AND. PRESENT(transmd)) THEN 
#ifdef SINGLE_LETKF
    call sgemv('t',nobsl,ne,1.0e0,hdxb_rinv,nobsl,depd,1,0.0e0,work2d,1)
    call sgemv('n',ne,ne,1.0e0,pa,ne,work2d,1,0.0e0,transmd,1)
#else
    call dgemv('t',nobsl,ne,1.0d0,hdxb_rinv,nobsl,depd,1,0.0d0,work2d,1)
    call dgemv('n',ne,ne,1.0d0,pa,ne,work2d,1,0.0d0,transmd,1)
#endif

  END IF
!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------
  DO j=1,ne
    rho = SQRT( REAL(ne-1,r_size) / eival(j) )
    DO i=1,ne
      work1(i,j) = eivec(i,j) * rho
    END DO
  END DO
#ifdef SINGLE_LETKF
  CALL sgemm('n','t',ne,ne,ne,1.0e0,work1,ne,eivec,&
    & ne,0.0e0,trans,ne)
#else
  CALL dgemm('n','t',ne,ne,ne,1.0d0,work1,ne,eivec,&
    & ne,0.0d0,trans,ne)
#endif
!  DO j=1,ne
!    DO i=1,ne
!      trans(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,ne
!        trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  IF (PRESENT(transm)) THEN                !GYL - if transm is present,
    transm = work3                         !GYL - return both trans and transm without adding them
  ELSE                                     !GYL
    DO j=1,ne
      DO i=1,ne
        trans(i,j) = trans(i,j) + work3(i)
      END DO
    END DO
  END IF                                   !GYL
  IF (PRESENT(pao)) pao = pa               !GYL

  IF (.NOT. infl_update_) RETURN           !GYL - skip the following if no inflation update is required
!-----------------------------------------------------------------------
!  Inflation estimation
!-----------------------------------------------------------------------
  parm = 0.0d0
  IF(rdiag_wloc_) THEN                            !GYL
    DO i=1,nobsl                                  !GYL
      parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i)  !GYL
    END DO                                        !GYL
  ELSE                                            !GYL
    DO i=1,nobsl
      parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i) * rloc(i)
    END DO
  END IF                                          !GYL
  DO j=1,ne
    DO i=1,nobsl
      parm(2) = parm(2) + hdxb_rinv(i,j) * hdxb(i,j)
    END DO
  END DO
  parm(2) = parm(2) / REAL(ne-1,r_size)
  parm(3) = SUM(rloc(1:nobsl))
  parm(4) = (parm(1)-parm(3))/parm(2) - parm_infl
!  sigma_o = 1.0d0/REAL(nobsl,r_size)/MAXVAL(rloc(1:nobsl))
  sigma_o = 2.0d0/parm(3)*((parm_infl*parm(2)+parm(3))/parm(2))**2
  gain = sigma_b**2 / (sigma_o + sigma_b**2)
  parm_infl = parm_infl + gain * parm(4)

  RETURN
END SUBROUTINE letkf_core

END MODULE common_letkf
