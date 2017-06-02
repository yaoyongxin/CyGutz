      module gprimme
      use gprec
      implicit none
!
      contains
!
      subroutine Zprimme(basismax, numemax, maxMatvecs, etol,  
     :                   printLevel,whichEvals, method, N, 
     :                   Nest, evals, evecs, AV)
      implicit none
      integer basismax, numemax, maxMatvecs, printLevel, whichEvals, 
     :        method, N, Nest
      real(8) etol, evals(numemax)
      COMPLEX(8) evecs(N*NUMEmax)
      external :: AV
! local
      integer(8) primme
      integer ierr, i
      real(8) rnorms(numemax)
      include 'primme_f77.h'
!
      ! Initialize PRIMME
      call primme_initialize_f77(primme)
      call primme_set_member_f77(primme, PRIMMEF77_n, N)
      call primme_set_member_f77(primme, PRIMMEF77_numEvals, NUMEmax)
      call primme_set_member_f77(primme, PRIMMEF77_maxBasisSize, 
     :                           BASISmax)
      call primme_set_member_f77(primme, PRIMMEF77_matrixMatvec, AV)
! Set the method to be used (after n, numEvals, and precondition have
! been set. Also after basisSize is set if desired.)
      call primme_set_method_f77(primme, method, ierr)
      call primme_set_member_f77(primme, PRIMMEF77_eps, ETOL)
      call primme_set_member_f77(primme, PRIMMEF77_target, whichEvals)
      call primme_set_member_f77(primme, PRIMMEF77_printLevel,
     :                           printLevel)
      call primme_set_member_f77(primme, PRIMMEF77_maxMatvecs, 
     :                           maxMatvecs)
      call primme_set_member_f77(primme, PRIMMEF77_initSize, Nest)
! Calling the PRIMME solver
      call Zprimme_f77(evals, evecs, rnorms, primme, ierr)
      if (ierr /= 0) then
        print *, ' ZPRIMME returned with error: ', ierr
      endif
      call primme_display_stats_f77(primme)
      call primme_free_f77(primme)
      write(*, '(" zprimme:")')
      do i = 1, numemax
        write (*, 9000) i, evals(i),rnorms(i)
      enddo
      return
9000   FORMAT (1x,'E(',i1,') = ',G24.16,4x,'residual norm =', E12.4)
      end subroutine Zprimme
      end module gprimme
