      SUBROUTINE DSDRV1(N,NEV,NCV,bmat,which,tol,sigma,rvec,WE,ZE,MAXITR,NCONV,AV)
      IMPLICIT NONE
!     Simple program to illustrate the idea of reverse communication
!     in regular mode for a standard symmetric eigenvalue problem.
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!     ... OP = A  and  B = I.
!     ... Assume "call av (n,x,y)" computes y = A*x.
!     ... Use mode 1 of DSAUPD.
! Routines called:
!     dsaupd  ARPACK reverse communication interface routine.
!     dseupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     av      Matrix vector multiplication routine that computes A*x.
      integer          N,NEV,NCV,MAXITR,NCONV
      character        bmat*1, which*2
      REAL(8)          tol, sigma
      logical          rvec
      REAL(8)          WE(NEV),ZE(N,NEV)
! LOCAL
      INTEGER maxn,maxnev,maxncv,ldv
      REAL(8),ALLOCATABLE::v(:,:),workl(:),workd(:),d(:,:),resid(:),ax(:)
      logical          select(ncv)
      integer          iparam(11), ipntr(11)
      integer          ido, lworkl, ierr, j, mode, ishfts, info
!     | BLAS & LAPACK routines used |
      REAL(8)          dnrm2
      external         dnrm2, daxpy, AV
!     %----------------------------------------------------%
!     | N is the 
!     | dimension of the matrix.  A standard eigenvalue    |
!     | problem is solved (BMAT = 'I'). NEV is the number  |
!     | of eigenvalues to be approximated.  The user can   |
!     | modify NEV, NCV, WHICH to solve problems of        |
!     | different sizes, and to get different parts of the |
!     | spectrum.  However, The following conditions must  |
!     | be satisfied:                                      |
!     |                   N <= MAXN,                       | 
!     |                 NEV <= MAXNEV,                     |
!     |             NEV + 1 <= NCV <= MAXNCV               | 
!     %----------------------------------------------------% 
!
      maxn=N; maxnev=NEV; maxncv=NCV; ldv=maxn
      ALLOCATE(v(ldv,maxncv),workl(maxncv*(maxncv+8)), &
              &workd(3*maxn),d(maxncv,2),resid(maxn),ax(maxn))
      v=0; workl=0; workd=0; d=0; resid=0; ax=0
!
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SDRV1: N is greater than MAXN '
         STOP
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
         STOP
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
         STOP
      end if
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
      lworkl = ncv*(ncv+8)
      info = 0
      ido = 0
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
      ishfts = 1
      mode   = 1
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
10    continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
     &                 ncv, v, ldv, iparam, ipntr, workd, workl, &
     &                 lworkl, info )
         if (ido .eq. -1 .or. ido .eq. 1) then
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
            call av (N,workd(ipntr(1)),workd(ipntr(2)))
            go to 10
         end if 
      if ( info .lt. 0 ) then
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
      else 
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        | Computed eigenvalues may be extracted.    |  
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma, &
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv, & 
     &        iparam, ipntr, workd, workl, lworkl, ierr )
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
         if ( ierr .ne. 0) then
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
         else
             nconv =  iparam(5)
             do 20 j=1, nconv
!               %---------------------------%
!               | Compute the residual norm |
!               |   ||  A*x - lambda*x ||   |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
                call av(n, v(1,j), ax)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
 20          continue
!             call dmout(7, nconv, 2, d, maxncv, -6, &
!     &            'Ritz values and relative residuals')
         end if
         WE=d(1:NEV,1); ZE=V(:,1:NEV)
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!         if ( info .eq. 1) then
!            print *, ' Maximum number of iterations reached.'
!         else if ( info .eq. 3) then
!            print *, ' No shifts could be applied during implicit', &
!     &               ' Arnoldi update, try increasing NCV.'
!         end if      
      end if
!
      end SUBROUTINE dsdrv1
!
!*********************************************************************************
      SUBROUTINE zhdrv1(N,NEV,NCV,bmat,which,tol,sigma,rvec,WE,ZE,MAXITR,NCONV,av)
      IMPLICIT NONE
!     for a standard complex Hermitian eigenvalue problem. 
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!     ... OP = A  and  B = I.
!     ... Assume "call av (n,x,y)" computes y = A*x
!     ... Use mode 1 of ZNAUPD .
!\Routines called
!     znaupd   ARPACK reverse communication interface routine.
!     zneupd   ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dznrm2   Level 1 BLAS that computes the norm of a complex vector.
!     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     av      Matrix vector multiplication routine that computes A*x.
!     %-----------------------------%
!     | Define maximum dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
      integer          N,NEV,NCV,MAXITR,NCONV
      character        bmat*1, which*2
      REAL(8)          tol
      Complex*16       sigma
      logical          rvec
      REAL(8)          WE(NEV)
      COMPLEX(8)       ZE(N,NEV)
! LOCAL
      integer           maxn, maxnev, maxncv, ldv
      integer           iparam(11), ipntr(14)
      logical           select(ncv)
      Complex*16        ax(n), d(ncv), &
     &                  v(n,ncv), workd(3*n), &
     &                  workev(3*ncv), resid(n), &
     &                  workl(3*ncv*ncv+5*ncv)
      Double precision  rwork(ncv), rd(ncv,3) 
      integer           ido, lworkl, info, j, &
     &                  ierr, ishfts, mode
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
      Double precision  dznrm2 , dlapy2 
      external          dznrm2 , zaxpy , dlapy2, av
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!     %--------------------------------------------------%
!     | The number N is the dimension of the             |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify N, NEV, NCV, WHICH to solve problems of   |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               | 
!     %--------------------------------------------------% 
!
      maxn=N; maxnev=NEV; maxncv=NCV; ldv=maxn
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD  as         | 
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         | 
!     %---------------------------------------------------%
      lworkl  = 3*ncv**2+5*ncv 
      ido    = 0
      info   = 0
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD .                                           |
!     %---------------------------------------------------%
      ishfts = 1
      mode   = 1
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode 
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
 10   continue
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD  and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
         call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv, &
     &        v, ldv, iparam, ipntr, workd, workl, lworkl, &
     &        rwork,info )
         if (ido .eq. -1 .or. ido .eq. 1) then
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
            CALL av(n,workd(ipntr(1)),workd(ipntr(2)))
            go to 10
         end if
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
      if ( info .lt. 0 ) then
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD   |
!        %--------------------------%
         WRITE(0,*) ' '
         WRITE(0,*) ' Error with _naupd, info = ', info
         WRITE(0,*) ' Check the documentation of znaupd'
         WRITE(0,*) ' '
      else 
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD .                |
!        | Computed eigenvalues may be extracted.    |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
         call zneupd  (rvec, 'A', select, d, v, ldv, sigma,  &
     &        workev, bmat, n, which, nev, tol, resid, ncv, &
     &        v, ldv, iparam, ipntr, workd, workl, lworkl, &
     &        rwork, ierr)
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  | 
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
         if ( ierr .ne. 0) then
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD . |
!           %------------------------------------%
             WRITE(0,*) ' '
             WRITE(0,*) ' Error with _neupd, info = ', ierr
             WRITE(0,*) ' Check the documentation of zneupd. '
             WRITE(0,*) ' '
         else
             nconv = iparam(5)
             do 20 j=1, nconv
!               %---------------------------%
!               | Compute the residual norm |
!               |   ||  A*x - lambda*x ||   |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
                CALL av(n,v(1,j),ax)
                call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble (d(j))
                rd(j,2) = dimag (d(j))
                rd(j,3) = dznrm2 (n, ax, 1)
                rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
                IF(ABS(rd(j,2))>1.D-6)THEN
                  WRITE(0,'(" ERROR IN zhdrv1: Eigen-value not real! imag(eval)=",E10.2)')rd(j,2)
                ENDIF
 20          continue
          end if
          WE=rd(1:NEV,1); ZE=V(:,1:NEV)
      end if
      RETURN
!
      end SUBROUTINE zhdrv1
