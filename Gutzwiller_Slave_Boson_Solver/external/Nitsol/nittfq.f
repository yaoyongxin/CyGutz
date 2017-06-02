      subroutine nittfq (n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, 
     &     ipar, ijacv, irpre, iksmax, ifdord, nfe, njve, nrpre, nli, r,
     &     rcgs, rtil, d, p, q, u, v, y, rwork1, rwork2, rsnrm, dinpr, 
     &                                                   dnorm, itrmks )

      implicit none

      integer ifdord, ijacv, iksmax, irpre, itrmks, n, nfe, njve,
     &        nrpre, nli

      integer ipar(*)

      double precision eta, fcnrm, rsnrm

      double precision d(n), fcur(n), p(n), q(n), rcgs(n), r(n),
     &                 rpar(*), rtil(n), rwork1(n), rwork2(n), step(n),
     &                 u(n), v(n), xcur(n), y(n)

      double precision dinpr, dnorm

      external f, jacv, dinpr, dnorm

c ------------------------------------------------------------------------
c
c This is nittfq v0.1, the TFQMR routine for determining (trial) inexact 
c Newton steps. The original reference is R. W. Freund, "A Transpose-Free
c Quasi-Minimal Residual Algorithm for Non-Hermitian Linear Systems", 
c SIAM J. Sci. Comput., 14 (1993), pp. 470-482.  The implementation here
c is based on the right preconditioned algorithm that appears in 
c J. N. Shadid and R. S. Tuminaro, "A Comparison of Preconditioned
c Nonsymmetric Krylov Methods on a Large-Scale MIMD Machine", SIAM J.
c Sci. Comput., 15 (1994), pp. 440-459.
c
c ------------------------------------------------------------------------
c
c Explanation:
c
c  n       = dimension of the problem.
c
c  xcur    = vector of length n, current approximate solution. 
c
c  fcur    = vector of length n, value of f at xcur.
c
c  fcnrm   = norm of fcur. 
c
c  step    = vector of length n, (trial) step. 
c
c  eta     = relative residual reduction factor. 
c
c  f      = name of user-supplied subroutine for evaluating the function 
c           the zero of which is sought; this routine has the form
c
c                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
c
c           where xcur is the array containing the current x value, fcur 
c           is f(xcur) on output, rpar and ipar are, respectively, real 
c           and integer parameter/work arrays for use by the subroutine,
c           and itrmf is an integer termination flag.  The meaning of
c           itrmf is as follows:
c             0 => normal termination; desired function value calculated.
c             1 => failure to produce f(xcur).
c 
c  jacv   = name of user-supplied subroutine for evaluating J*v or 
c           P(inverse)*v, where J is the Jacobian of f and P is a 
c           right preconditioning operator. If neither analytic J*v 
c           evaluations nor right preconditioning is used, this can 
c           be a dummy subroutine; if right preconditioning is used but 
c           not analytic J*v evaluations, this need only evaluate 
c           P(inverse)*v. The form is 
c
c           subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
c
c           where xcur and fcur are vectors of length n containing the 
c           current x and f values, ijob is an integer flag indicating 
c           which product is desired, v is a vector of length n to be 
c           multiplied, z is a vector of length n containing the desired 
c           product on output, rpar and ipar are, respectively, real 
c           and integer parameter/work arrays for use by the subroutine, 
c           and itrmjv is an integer termination 
c           flag. The meaning of ijob is as follows: 
c             0 => z = J*v
c             1 => z = P(inverse)*v 
c           The meaning of itrmjv is as follows:
c             0 => normal termination; desired product evaluated. 
c             1 => failure to produce J*v.
c             2 => failure to produce P(inverse)*v. 
c           This subroutine is called only from nitjv, and is always 
c           called with v .ne. z. 
c
c  rpar    = real parameter/work array passed to the f and jacv routines. 
c
c  ipar    = integer parameter/work array passed to the f and jacv routines. 
c
c  ijacv   = flag for determining method of J*v evaluation.
c              0 => finite-difference evaluation (default). 
c              1 => analytic evaluation. 
c
c  irpre   = flag for right preconditioning. 
c              0 => no right preconditioning
c              1 => right preconditioning
c
c  iksmax  = maximum allowable number of TFQMR iterations. 
c
c  ifdord  = order of the finite-difference formula used in BiCGSTAB 
c            when J*v products are evaluated using finite-differences. 
c            When ijacv = 0 on input to nitsol, ifdord is set to 1, 2, or 
c            4 in nitsol; otherwise, it is irrelevant. When ijacv = 0 on 
c            input to this subroutine, ifdord determines the order of the 
c            finite-difference formula used at each BiCGSTAB iteration 
c            (default 1). In this case, ijacv is set to -1 below to 
c            signal to nitjv that the order of the finite-difference 
c            formula is to be determined by ifdord. The original value 
c            ijacv = 0 is restored on return. 
c            
c  nfe     = number of function evaluations.
c
c  njve    = number of J*v evaluations. 
c
c  nrpre   = number of P(inverse)*v evaluations.
c
c  nli     = number of linear iterations.
c
c  r       = residual vector (for the QMR process)
c
c  rcgs    = residual vector (of the underlying CGS process)
c
c  rtil    = 'shadow' residual vector used in bi-orthogonalization
c
c  d       = vector used in TFQMR
c
c  p       = vector used in TFQMR
c
c  q       = vector used in TFQMR
c
c  u       = vector used in TFQMR
c
c  v       = vector used in TFQMR
c
c  y       = vector used in TFQMR
c
c  rwork1  = work vector, passed on to nitjv
c
c  rwork2  = work vector, passed on to nitjv
c
c  rsnrm   = TFQMR residual norm on return. 
c
c  dinpr   = inner-product routine, either user-supplied or blas ddot. 
c
c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
c
c  itrmks  = termination flag; values have the following meanings: 
c              0 => normal termination: acceptable step found. 
c              1 => J*v failure in nitjv. 
c              2 => P(inverse)*v failure in nitjv. 
c              3 => acceptable step not found in iksmax TFQMR iterations. 
c              4 => TFQMR breakdown.
c              5 => floating point error (the underlying CGS iteration
c                   has probably blown up)
c
c             Note: On return, nitsol terminates if itrmks is 1 or 2. 
c             If itrmks is 3 or 4, nitsol may terminate or continue. 
c             In this event, the step returned is a meaningful inexact 
c             Newton step only if the residual norm has been reduced. 
c             A decision on termination/continuation is made in nitdrv 
c             according to whether there is sufficient residual norm 
c             reduction, even though the desired inexact Newton condition 
c             may not hold.  
c
c -------------------------------------------------------------------------

c  Subroutines required by this and all called routines:

c    user supplied:  nitjv

c    nitsol routines: none

c    BLAS routines- dcopy, daxpy, dscal, dswap

c    LAPACK routines - dlamch

c    user supplied or BLAS:  dinpr, dnorm

c    explanation: In nitsol, dinpr and dnorm are set to either the BLAS 
c    ddot and dnrm2 routines or the user-supplied routines. 

c This subroutine called by: nitdrv

c Subroutines called by this subroutine: daxpy, dcopy, dscal, dswap, dinpr,
c    dlamch, dnorm, nitjv, dlamch

c Common block: 

      include 'nitprint.h'

c If diagnostic information is desired, include this common block in the 
c main program and set iplvl and ipunit according to the following: 
c
c     iplvl = 0 => no printout
c           = 1 => iteration numbers and F-norms
c           = 2 => ... + some stats, step norms, and linear model norms
c           = 3 => ... + some Krylov solver and backtrack information
c           = 4 => ... + more Krylov solver and backtrack information
c
c     ipunit = printout unit number.

c  Parameters-

      double precision zero,       one
      parameter      ( zero=0.0d0, one=1.0d0 )

c  Local variables-

      integer i
      integer itask
      integer itrmjv
      integer itfq
      integer k

      double precision alpha
      double precision abstol
      double precision beta
      double precision c
      double precision cgsnorm
      double precision qmreta
      double precision omega
      double precision rho
      double precision rho_old
      double precision sigma
      double precision t
      double precision tau
      double precision theta

      character*2 ab(0:1)
      data ab /'.a','.b'/

      double precision sfmin
      data sfmin /zero/

c  External subroutines-

      double precision dlamch

      external daxpy
      external dcopy
      external dlamch
      external dscal
      external dswap
      external nitjv

c  Intrinsics-

      intrinsic abs
      intrinsic dble
      intrinsic sqrt

c  Start of executable code-

c  Initialize sfmin only on first entry.

      if ( sfmin .eq. zero ) sfmin = dlamch( 's' )

c If finite-differences are used to evaluate J*v products (ijacv= 0), then 
c ijacv is set to -1 within this subroutine to signal to nitjv that the 
c order of the finite-difference formula is to be determined by ifdord. 
c The original value ijacv= 0 is restored on return. 

      if (ijacv .eq. 0) ijacv = -1 

c Set the stopping tolerance, initialize the step, etc. 

      rsnrm = fcnrm
      abstol = eta*rsnrm
      do 10 i = 1, n
         step(i) = zero
 10   continue
      itfq = 0

c For printing:

      if ( iplvl .ge. 3 ) then 
         write(ipunit,*) 
         write(ipunit,800) eta 
      endif
      if ( iplvl .ge. 4 ) then 
         write(ipunit,810) 
         write(ipunit,*) 
         write(ipunit,820) itfq, rsnrm 
      endif

c  Initialize residual and work vectors.

      call dcopy( n, fcur, 1, rcgs, 1 )
      call dscal( n, -one, rcgs, 1 )

c  Choice here is rtil = r.

      call dcopy( n, rcgs, 1, rtil, 1 )

      if ( irpre .eq. 0 ) then
         call dcopy( n, rcgs, 1, p, 1 )
      else
         itask = 2
         call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
     &               itask, nfe, njve, nrpre, rcgs, p, rwork1, rwork2,
     &                                                   dnorm, itrmjv )    
         if ( itrmjv .ne. 0 ) then
            itrmks = 2
            goto 900
         endif
      endif
      call dcopy( n, p, 1, u, 1 )
      rho = dinpr( n, rcgs, 1, rtil, 1 )
      do 20 i = 1, n
         d(i) = zero
 20   continue

      itask = 0
      call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
     &            itask, nfe, njve, nrpre, p, v, rwork1, rwork2,
     &                                                   dnorm, itrmjv )
      if ( itrmjv .ne. 0 ) then
         itrmks = 1
         goto 900
      end if

      alpha = zero
      omega = fcnrm
      tau = omega
      theta = zero
      qmreta = zero

c  Start iterations.

 100  continue
      itfq = itfq + 1
      nli = nli + 1
      sigma = dinpr( n, rtil, 1, v, 1 )

c  If sigma = 0 we have a serious breakdown.  We check this condition
c  by trying to detect whether division by sigma causes an overflow.

      if ( abs(sigma) .lt. sfmin*abs(rho) ) then
         itrmks = 4
         goto 900
      else
         alpha = rho/sigma
      endif

c  Need Pv for calculation of q.  First store result in q, then
c  swap some vectors to cast calculation of q as a SAXPY.

      if ( irpre .eq. 0 ) then
         call dcopy( n, v, 1, q, 1 )
      else
         itask = 2
         call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
     &               itask, nfe, njve, nrpre, v, q, rwork1, rwork2,
     &                                                   dnorm, itrmjv )    
         if ( itrmjv .ne. 0 ) then
            itrmks = 2
            goto 900
         endif
      endif
      call dcopy( n, u, 1, y, 1 )
      call dswap( n, q, 1, u, 1 )
      call daxpy( n, -alpha, u, 1, q, 1 )
      call dcopy( n, y, 1, u, 1 )

c  Update residual.

      do 30 i = 1, n
         y(i) = u(i) + q(i)
 30   continue
      itask = 0
      call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
     &            itask, nfe, njve, nrpre, y, v, rwork1, rwork2,
     &                                                   dnorm, itrmjv )
      if ( itrmjv .ne. 0 ) then
         itrmks = 1
         goto 900
      end if
      call daxpy( n, -alpha, v, 1, rcgs, 1 )
      cgsnorm = dnorm( n, rcgs, 1 )

c  Check for cgsnorm = NaN.

      if ( cgsnorm .ne. cgsnorm ) then
         itrmks = 5
         goto 900
      endif

c  QMR section.

      do 60 k = 0, 1

c  Use weighting strategy from (5.11) of Freund reference.

         t = qmreta*theta**2
         t = t/alpha

         if ( k .eq. 0 ) then
            do 40 i = 1, n
               d(i) = u(i) + t*d(i)
 40         continue
            omega = sqrt(omega*cgsnorm)
         else if ( k .eq. 1 ) then
            do 50 i = 1, n
               d(i) = q(i) + t*d(i)
 50         continue
            omega = cgsnorm
         endif

         theta = omega/tau
         c = one/sqrt(one + theta**2)
         tau = tau*theta*c
         qmreta = alpha*c**2

c For printing:

         if ( iplvl .ge. 4 .and. tau .gt. abstol ) then 
            write(ipunit,830) itfq, ab(k), tau, '     (estimated)'
         endif

         call daxpy( n, qmreta, d, 1, step, 1 )

c  Convergence check.  Do a cheap test to save on Jacobi-vector products.
c  In case residual history is requested by iplvl, we must calculate 
c  the QMR residual from scratch.  Note termination is always determined
c  by the smoothed residual, calculated from scratch.

         if ( tau .le. abstol ) then

            itask = 0
            call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv,
     &            ifdord, itask, nfe, njve, nrpre, step, r, rwork1,
     &                                           rwork2, dnorm, itrmjv )

c  This calculation of the QMR residual is off by a factor
c  of -1, but we don't care until we return from this routine.

            call daxpy( n, one, fcur, 1, r, 1 )
            rsnrm = dnorm( n, r, 1 )

c For printing:

            if ( iplvl .ge. 4 ) then 
               write(ipunit,830) itfq, ab(k), rsnrm, '  (from scratch)'
            endif

c  Check for rsnrm = NaN.

            if ( rsnrm .ne. rsnrm ) then
               itrmks = 5
               goto 900
            endif

c  If rsnrm is small enough, exit.

            if ( rsnrm .lt. abstol ) then
               itrmks = 0
               goto 900
            endif
         endif

 60   continue

      rho_old = rho
      rho = dinpr( n, rtil, 1, rcgs, 1 )

c  If rho_old = 0 we have a serious breakdown.  We check this condition
c  by trying to detect whether division by rho_old causes an overflow.

      if ( abs(rho_old) .lt. sfmin*abs(rho) ) then
         itrmks = 4
         goto 900
      else
         beta = rho/rho_old
      endif

      if ( irpre .eq. 0 ) then
         call dcopy( n, rcgs, 1, v, 1 )
      else
         itask = 2
         call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
     &               itask, nfe, njve, nrpre, rcgs, v, rwork1, rwork2,
     &                                                   dnorm, itrmjv )    
         if ( itrmjv .ne. 0 ) then
            itrmks = 2
            goto 900
         endif
      endif
      call daxpy( n, beta, q, 1, v, 1 )
      call dcopy( n, v, 1, u, 1 )
      call daxpy( n, beta, p, 1, q, 1 )
      call daxpy( n, beta, q, 1, v, 1 )
      call dcopy( n, v, 1, p, 1 )

      itask = 0
      call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
     &            itask, nfe, njve, nrpre, p, v, rwork1, rwork2,
     &                                                   dnorm, itrmjv )
      if ( itrmjv .ne. 0 ) then
         itrmks = 1
         goto 900
      end if
      if ( itfq .ge. iksmax ) then
         itrmks = 3
         goto 900
      end if

c  Do again

      goto 100

c  All returns made here.

  900 continue

c  If residual hasn't been updated, force
c  computation of residual from scratch.

      if ( rsnrm .eq. fcnrm ) then

         itask = 0
         call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv,
     &         ifdord, itask, nfe, njve, nrpre, step, r, rwork1,
     &                                        rwork2, dnorm, itrmjv )

         call daxpy( n, one, fcur, 1, r, 1 )
         rsnrm = dnorm( n, r, 1 )
         if ( rsnrm .le. abstol ) itrmks = 0

      end if

c  Correct residual before returning.

      call dscal( n, -one, r, 1 )

c If ijacv = -1, then restore it to the original value ijacv = 0. 

      if (ijacv .eq. -1) ijacv = 0 

c For printing:

      if ( iplvl .ge. 3 ) then 
         write(ipunit,*) 
         if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
            write(ipunit,840) itrmks, rsnrm 
         else
            write(ipunit,850) itrmks
         endif
      endif

      return

 800  format('nittfq:  eta =', 1pd10.3)
 810  format('nittfq:  TFQMR iteration no. (parts a and b),',
     &                                        ' linear residual norm')
 820  format(5x,i4,5x,1pd10.3)
 830  format(5x,i4,a2,3x,1pd10.3,a16)
 840  format('nittfq:  itrmks =', i2, '   final lin. res. norm =', 
     &                                                          1pd10.3)
 850  format('nittfq: itrmks:', i4) 

c  End of nittfq.

      end
