      subroutine nitstb (n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, 
     $     ipar, ijacv, irpre, iksmax, ifdord, nfe, njve, 
     $     nrpre, nli, r, rtil, p, phat, v, t, rwork1, rwork2, 
     $     rsnrm, dinpr, dnorm, itrmks)

      implicit none 

      integer n, ipar(*), ijacv, irpre, iksmax, ifdord, nfe, njve, 
     $     nrpre, nli, itrmks
      double precision xcur(n), fcur(n), fcnrm, step(n), eta, rpar(*), 
     $     r(n), rtil(n), p(n), phat(n), v(n), t(n), rwork1(n), 
     $     rwork2(n), rsnrm, dinpr, dnorm 
      external f, jacv, dinpr, dnorm 

c ------------------------------------------------------------------------
c
c This is nitstb v0.3, the BiCGSTAB routine for determining (trial) inexact 
c Newton steps. The original reference is H. van der Vorst, "Bi-CGSTAB: 
c A fast and smoothly converging variant of Bi-CG for the soluton of 
c nonsymmetric linear systems," SIAM J. Sci. Statist. Comput., 13 (1992), 
c pp. 631--644. 
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
c  iksmax  = maximum allowable number of BiCGSTAB iterations. 
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
c  r       = residual vector 
c
c  rtil    = "r-tilde" vector used in BiCGSTAB
c
c  p       = vector used in BiCGSTAB
c
c  phat    = vector used in BiCGSTAB
c
c  v       = vector used in BiCGSTAB
c
c  t       = vector used in BiCGSTAB
c
c  rwork1  = work vector, passed on to nitjv
c
c  rwork2  = work vector, passed on to nitjv
c
c  rsnrm   = BiCGSTAB residual norm on return. 
c
c  dinpr   = inner-product routine, either user-supplied or blas ddot. 
c
c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
c
c  itrmks  = termination flag; values have the following meanings: 
c              0 => normal termination: acceptable step found. 
c              1 => J*v failure in nitjv. 
c              2 => P(inverse)*v failure in nitjv. 
c              3 => acceptable step not found in iksmax BiCGSTAB iterations. 
c              4 => BiCGSTAB breakdown. 
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
c
c Subroutines required by this and all called routines: 
c
c    user supplied: f, jacv 
c
c    nitsol routines: nitjv, nitfd
c
c    blas routines: daxpy, dcopy, dscal
c
c    lapack routine:  dlamch
c
c    user supplied or blas: dinpr, dnorm 
c
c    explanation: In nitsol, dinpr and dnorm are set to either the blas 
c    ddot and dnrm2 routines or the user-supplied usrnpr and usrnrm 
c    routines. 
c
c This subroutine called by: nitdrv
c
c Subroutines called by this subroutine: daxpy, dcopy, dscal, dinpr, dlamch,
c    dnorm, nitjv
c
c Common block: 
c
      include 'nitprint.h'
c
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
c
c ------------------------------------------------------------------------
c
      double precision abstol, alpha, beta, omega, 
     $     rho, rhomns, tau, temp
      integer i, istb, itask, itrmjv

      double precision dlamch
      external dlamch
c
      double precision sfmin
      data sfmin /0.0d0/
c ------------------------------------------------------------------------
c
c  Initialize sfmin only on first entry.
c
      if ( sfmin .eq. 0.0d0 ) sfmin = dlamch( 's' )
c ------------------------------------------------------------------------
c If finite-differences are used to evaluate J*v products (ijacv= 0), then 
c ijacv is set to -1 within this subroutine to signal to nitjv that the 
c order of the finite-difference formula is to be determined by ifdord. 
c The original value ijacv= 0 is restored on return. 
c ------------------------------------------------------------------------
      if (ijacv .eq. 0) ijacv = -1 
c ------------------------------------------------------------------------
c Set the stopping tolerance, initialize the step, etc. 
c ------------------------------------------------------------------------
      rsnrm = fcnrm
      abstol = eta*rsnrm
      do 10 i = 1, n
         step(i) = 0.d0
 10   continue
      istb = 0
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 3) then 
         write(ipunit,*) 
         write(ipunit,800) eta 
      endif
 800  format('nitstb:  eta =', 1pd10.3)
      if (iplvl .ge. 4) then 
         write(ipunit,810) 
         write(ipunit,*) 
         write(ipunit,820) istb, rsnrm 
      endif
 810  format('nitstb:  BiCGSTAB iteration no. (parts a and b)',
     $     ' linear residual norm, ')
 820  format(5x,i4,5x,1pd10.3)
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Set up r and rtil. 
c ------------------------------------------------------------------------
      call dcopy(n,fcur,1,r,1)
      temp = -1.d0
      call dscal(n,temp,r,1)
      call dcopy(n,r,1,rtil,1)
c ------------------------------------------------------------------------
c Top of the iteration loop. 
c ------------------------------------------------------------------------
 100  continue
      istb = istb + 1
      nli = nli + 1
c ------------------------------------------------------------------------
c Perform the first "half-iteration". 
c ------------------------------------------------------------------------
      rho = dinpr(n,rtil,1,r,1)
      if (istb .eq. 1) then 
         call dcopy(n,r,1,p,1)
      else
         if ( abs(rhomns) .lt. sfmin*abs(rho) ) then
            itrmks = 4
            goto 900
         else
            beta = (rho/rhomns)*(alpha/omega)
            call daxpy(n,-omega,v,1,p,1)
            call dscal(n,beta,p,1)
            call daxpy(n,1.d0,r,1,p,1)
         endif
      endif
      if (irpre .eq. 0) then 
         call dcopy(n,p,1,phat,1)
      else
         itask = 2
         call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, 
     $        ifdord, itask, nfe, njve, nrpre, p, phat, 
     $        rwork1, rwork2, dnorm, itrmjv)
         if (itrmjv .gt. 0) then 
            itrmks = 2
            go to 900
         endif
      endif
      itask = 0
      call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, 
     $     ifdord, itask, nfe, njve, nrpre, phat, v, 
     $     rwork1, rwork2, dnorm, itrmjv)
      if (itrmjv .gt. 0) then 
         itrmks = 1
         go to 900
      endif
      tau = dinpr(n,rtil,1,v,1)
      if ( abs(tau) .lt. sfmin*abs(rho) ) then
         itrmks = 4
         goto 900
      else
         alpha = rho/tau
      endif
      call daxpy(n,-alpha,v,1,r,1)
      call daxpy(n,alpha,phat,1,step,1)
      rsnrm = dnorm(n,r,1)
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 4) then 
         write(ipunit,830) istb, rsnrm 
 830     format(5x,i4,'.a',3x,1pd10.3)
      endif
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Test for termination. 
c ------------------------------------------------------------------------
      if (rsnrm .le. abstol) then 
         itrmks = 0
         go to 900
      endif
c ------------------------------------------------------------------------
c Perform the second "half-iteration". 
c ------------------------------------------------------------------------
      if (irpre .eq. 0) then 
         call dcopy(n,r,1,phat,1)
      else
         itask = 2
         call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, 
     $        ifdord, itask, nfe, njve, nrpre, r, phat, 
     $        rwork1, rwork2, dnorm, itrmjv)
         if (itrmjv .gt. 0) then 
            itrmks = 2
            go to 900
         endif
      endif
      itask = 0
      call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, 
     $     ifdord, itask, nfe, njve, nrpre, phat, t, 
     $     rwork1, rwork2, dnorm, itrmjv)
      if (itrmjv .gt. 0) then 
         itrmks = 1
         go to 900
      endif
      tau = dnorm(n,t,1)
      tau = tau*tau
      temp = dinpr(n,t,1,r,1)
      if ( tau .le. sfmin*abs(temp) ) then
         itrmks = 4
         goto 900
      else
         omega = temp/tau
      endif
      if ( abs(omega) .lt. sfmin*abs(alpha) ) then 
         itrmks = 4
         go to 900
      endif
      call daxpy(n,-omega,t,1,r,1)
      call daxpy(n,omega,phat,1,step,1)
      rsnrm = dnorm(n,r,1)
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 4) then 
         write(ipunit,840) istb, rsnrm 
 840     format(5x,i4,'.b',3x,1pd10.3)
      endif
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Test for termination. 
c ------------------------------------------------------------------------
      if (rsnrm .le. abstol) then 
         itrmks = 0
         go to 900
      endif
      if (istb .ge. iksmax) then 
         itrmks = 3
         go to 900
      endif
c ------------------------------------------------------------------------
c If continuing, update and return to the top of the iteration loop. 
c ------------------------------------------------------------------------
      rhomns = rho
      go to 100
c ------------------------------------------------------------------------
c Bottom of the iteration loop. 
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 3) then 
         write(ipunit,*) 
         if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
            write(ipunit,850) itrmks, rsnrm 
 850        format('nitstb:  itrmks =', i2, '   final lin. res. norm =', 
     $           1pd10.3)
         else
            write(ipunit,860) itrmks
 860        format('nitstb: itrmks:', i4) 
         endif
      endif
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c If ijacv = -1, then restore it to the original value ijacv = 0. 
c ------------------------------------------------------------------------
      if (ijacv .eq. -1) ijacv = 0 
      return
      end
