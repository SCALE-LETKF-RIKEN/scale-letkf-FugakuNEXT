!-------------------------------------------------------------------------------
!> Solve tridiagonal systems in place with the partition method
!!
!! This file is included in a gang (one or more columns per gang) of an OpenACC
!! compute region, so that the solver is inlined; the work arrays of the caller
!! (e.g., in the shared memory) are accessed directly.
!! A subroutine version is slower: it is not inlined, the arrays are accessed
!! through generic and volatile loads, and the registers of the caller increase.
!!
!! Before including, define the following macros (they are undefined at the end):
!!   TDPM_UD(k,l), TDPM_MD(k,l), TDPM_LD(k,l) : upper, middle, and lower diagonals
!!   TDPM_X(k,l) : right hand side (in) and solution (out)
!!   TDPM_KS, TDPM_KE : range of the rows
!!   TDPM_NC          : number of the systems (columns, l = 1:TDPM_NC)
!!   TDPM_NPART       : upper limit of the number of the chunks per system
!! TDPM_UD and TDPM_LD are overwritten. TDPM_LD(TDPM_KS,:) and TDPM_UD(TDPM_KE,:)
!! are not read. The variables in matrix_solver_tridiagonal_pm_decl.h must be
!! declared in the caller.
!!
!! The rows of a system are divided into chunks by single separator rows, and a
!! vector lane handles a chunk.
!! 1. Each lane eliminates the rows in its chunk, giving the solution in the
!!    chunk as x = y + v * Xl + w * Xu, where Xl and Xu are the solutions at the
!!    separators below and above the chunk. y, v, and w overwrite X, LD, and UD.
!! 2. Substituting them into the separator rows gives a tridiagonal system for
!!    the separators, which is solved with one system on each lane.
!! 3. Each lane substitutes Xl and Xu back into its chunk.
!-------------------------------------------------------------------------------

    ! The first tdpm_nrem chunks have tdpm_nbase+1 rows and the others have tdpm_nbase rows.
    tdpm_npart = max( min( TDPM_NPART, ( TDPM_KE - TDPM_KS + 2 ) / 2 ), 1 )
    tdpm_nbase = ( TDPM_KE - TDPM_KS + 2 - tdpm_npart ) / tdpm_npart
    tdpm_nrem  = mod( TDPM_KE - TDPM_KS + 2 - tdpm_npart, tdpm_npart )

    ! 1. elimination in the chunks
    !$acc loop vector collapse(2) independent private(tdpm_k,tdpm_kb,tdpm_kt,tdpm_tmp)
    do tdpm_l = 1, TDPM_NC
    do tdpm_p = 0, tdpm_npart-1
       tdpm_kb = TDPM_KS + tdpm_p * ( tdpm_nbase + 1 ) + min( tdpm_p, tdpm_nrem )
       tdpm_kt = tdpm_kb + tdpm_nbase - 1 + merge( 1, 0, tdpm_p < tdpm_nrem )
       ! forward elimination: X, LD, and UD store y', v', and c'
       tdpm_tmp = 1.0_RP / TDPM_MD(tdpm_kb,tdpm_l)
       TDPM_X(tdpm_kb,tdpm_l) = TDPM_X(tdpm_kb,tdpm_l) * tdpm_tmp
       if ( tdpm_p > 0 ) then
          TDPM_LD(tdpm_kb,tdpm_l) = - TDPM_LD(tdpm_kb,tdpm_l) * tdpm_tmp
       else
          TDPM_LD(tdpm_kb,tdpm_l) = 0.0_RP
       end if
       if ( tdpm_kb < tdpm_kt ) TDPM_UD(tdpm_kb,tdpm_l) = TDPM_UD(tdpm_kb,tdpm_l) * tdpm_tmp
       !$acc loop seq
       do tdpm_k = tdpm_kb+1, tdpm_kt
          tdpm_tmp = 1.0_RP / ( TDPM_MD(tdpm_k,tdpm_l) - TDPM_LD(tdpm_k,tdpm_l) * TDPM_UD(tdpm_k-1,tdpm_l) )
          TDPM_X (tdpm_k,tdpm_l) = ( TDPM_X(tdpm_k,tdpm_l) - TDPM_LD(tdpm_k,tdpm_l) * TDPM_X(tdpm_k-1,tdpm_l) ) * tdpm_tmp
          TDPM_LD(tdpm_k,tdpm_l) = - TDPM_LD(tdpm_k,tdpm_l) * TDPM_LD(tdpm_k-1,tdpm_l) * tdpm_tmp
          if ( tdpm_k < tdpm_kt ) TDPM_UD(tdpm_k,tdpm_l) = TDPM_UD(tdpm_k,tdpm_l) * tdpm_tmp
       end do
       ! w' at the top row of the chunk
       if ( tdpm_p < tdpm_npart-1 ) then
          TDPM_UD(tdpm_kt,tdpm_l) = - TDPM_UD(tdpm_kt,tdpm_l) * tdpm_tmp
       else
          TDPM_UD(tdpm_kt,tdpm_l) = 0.0_RP
       end if
       ! back substitution: UD is overwritten from c' to w
       !$acc loop seq
       do tdpm_k = tdpm_kt-1, tdpm_kb, -1
          TDPM_X (tdpm_k,tdpm_l) = TDPM_X (tdpm_k,tdpm_l) - TDPM_UD(tdpm_k,tdpm_l) * TDPM_X (tdpm_k+1,tdpm_l)
          TDPM_LD(tdpm_k,tdpm_l) = TDPM_LD(tdpm_k,tdpm_l) - TDPM_UD(tdpm_k,tdpm_l) * TDPM_LD(tdpm_k+1,tdpm_l)
          TDPM_UD(tdpm_k,tdpm_l) =                        - TDPM_UD(tdpm_k,tdpm_l) * TDPM_UD(tdpm_k+1,tdpm_l)
       end do
    end do
    end do

    ! 2. tridiagonal system for the separators, by the Thomas algorithm.
    ! At the separator k, the sub- and super-diagonal elements are LD(k) * v(k-1)
    ! and UD(k) * w(k+1), and the diagonal element is MD(k) + LD(k) * w(k-1)
    ! + UD(k) * v(k+1). UD and X are overwritten by c' and the solution.
    !$acc loop vector independent private(tdpm_k,tdpm_p,tdpm_kb,tdpm_tmp,tdpm_Xl)
    do tdpm_l = 1, TDPM_NC
       !$acc loop seq
       do tdpm_p = 0, tdpm_npart-2
          tdpm_k = TDPM_KS + tdpm_p * ( tdpm_nbase + 1 ) + min( tdpm_p, tdpm_nrem ) &
                 + tdpm_nbase + merge( 1, 0, tdpm_p < tdpm_nrem )
          tdpm_tmp = TDPM_X(tdpm_k,tdpm_l) - TDPM_LD(tdpm_k,tdpm_l) * TDPM_X(tdpm_k-1,tdpm_l) &
                                           - TDPM_UD(tdpm_k,tdpm_l) * TDPM_X(tdpm_k+1,tdpm_l)
          tdpm_Xl = TDPM_MD(tdpm_k,tdpm_l) + TDPM_LD(tdpm_k,tdpm_l) * TDPM_UD(tdpm_k-1,tdpm_l) &
                                           + TDPM_UD(tdpm_k,tdpm_l) * TDPM_LD(tdpm_k+1,tdpm_l) ! diagonal element
          if ( tdpm_p > 0 ) then
             ! LD(k) * v(k-1) is the sub-diagonal element
             tdpm_Xl  = tdpm_Xl  - TDPM_LD(tdpm_k,tdpm_l) * TDPM_LD(tdpm_k-1,tdpm_l) * TDPM_UD(tdpm_kb,tdpm_l)
             tdpm_tmp = tdpm_tmp - TDPM_LD(tdpm_k,tdpm_l) * TDPM_LD(tdpm_k-1,tdpm_l) * TDPM_X (tdpm_kb,tdpm_l)
          end if
          tdpm_Xl = 1.0_RP / tdpm_Xl
          TDPM_UD(tdpm_k,tdpm_l) = TDPM_UD(tdpm_k,tdpm_l) * TDPM_UD(tdpm_k+1,tdpm_l) * tdpm_Xl
          TDPM_X (tdpm_k,tdpm_l) = tdpm_tmp * tdpm_Xl
          tdpm_kb = tdpm_k ! the previous separator
       end do
       !$acc loop seq
       do tdpm_p = tdpm_npart-3, 0, -1
          tdpm_k = TDPM_KS + tdpm_p * ( tdpm_nbase + 1 ) + min( tdpm_p, tdpm_nrem ) &
                 + tdpm_nbase + merge( 1, 0, tdpm_p < tdpm_nrem )
          TDPM_X(tdpm_k,tdpm_l) = TDPM_X(tdpm_k,tdpm_l) - TDPM_UD(tdpm_k,tdpm_l) * TDPM_X(tdpm_kb,tdpm_l)
          tdpm_kb = tdpm_k ! the next separator
       end do
    end do

    ! 3. back substitution of the separators into the chunks
    !$acc loop vector collapse(2) independent private(tdpm_k,tdpm_kb,tdpm_kt,tdpm_Xl,tdpm_Xu)
    do tdpm_l = 1, TDPM_NC
    do tdpm_p = 0, tdpm_npart-1
       tdpm_kb = TDPM_KS + tdpm_p * ( tdpm_nbase + 1 ) + min( tdpm_p, tdpm_nrem )
       tdpm_kt = tdpm_kb + tdpm_nbase - 1 + merge( 1, 0, tdpm_p < tdpm_nrem )
       tdpm_Xl = 0.0_RP
       tdpm_Xu = 0.0_RP
       if ( tdpm_p > 0 )            tdpm_Xl = TDPM_X(tdpm_kb-1,tdpm_l)
       if ( tdpm_p < tdpm_npart-1 ) tdpm_Xu = TDPM_X(tdpm_kt+1,tdpm_l)
       !$acc loop seq
       do tdpm_k = tdpm_kb, tdpm_kt
          TDPM_X(tdpm_k,tdpm_l) = TDPM_X(tdpm_k,tdpm_l) + TDPM_LD(tdpm_k,tdpm_l) * tdpm_Xl + TDPM_UD(tdpm_k,tdpm_l) * tdpm_Xu
       end do
    end do
    end do

#undef TDPM_UD
#undef TDPM_MD
#undef TDPM_LD
#undef TDPM_X
#undef TDPM_KS
#undef TDPM_KE
#undef TDPM_NC
#undef TDPM_NPART
