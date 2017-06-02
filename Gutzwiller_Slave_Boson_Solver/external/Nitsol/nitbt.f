      subroutine nitbt(n, xcur, fcnrm, step, eta, xpls, fpls, fpnrm, 
     $     oftjs, redfac, nfe, ibt, ibtmax, f, rpar, ipar, dnorm, 
     $     itrmbt)

      implicit none 

      integer n, nfe, ibt, ibtmax, ipar(*), itrmbt
      double precision xcur(n), fcnrm, step(n), eta, 
     $     xpls(n), fpls(n), fpnrm, oftjs, redfac, rpar(*), dnorm
      external f, dnorm 

c ------------------------------------------------------------------------
c
c This is nitbt v0.3, the backtracking routine for the (inexact) Newton 
c iterations. 
c
c ------------------------------------------------------------------------
c 
c Explanation: 
c
c  n       = dimension of the problem
c
c  xcur    = vector of length n, current approximate solution (input). 
c
c  fcnrm   = norm of f(xcur) 
c
c  step    = vector of length n, initial (trial) step on input, final 
c            acceptable step on output. 
c
c  eta     = initial inexact Newton level on input, final inexact Newton 
c            level on output. 
c
c  xpls    = vector of length n, next approximate solution on output. 
c
c  fpls    = vector of length n, value of f at xpls. 
c
c  fpnrm   = norm of f(xpls) 
c
c  oftjs   = original value of f(transpose)*Js. 
c
c  redfac   = scalar factor by which the original step is reduced on output. 
c
c  nfe     = number of function evaluations.
c
c  ibt     = number of backtracks on this call. 
c
c  ibtmax   = maximum allowable number of backtracks (step reductions) 
c             per call to nitbt (default 10). 
c
c             USAGE NOTE: If ibtmax = -1, then backtracking is turned 
c             off. In this case, the only function of this subroutine 
c             is to update xpls, fpls, and fpnrm. 
c
c  f       = name of user-supplied subroutine for evaluating the function 
c            the zero of which is sought; this routine has the form 
c
c                  subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
c
c           where xcur is the array containing the current x value, fcur 
c           is f(xcur) on output, rpar and ipar are, respectively, real 
c           and integer parameter/work arrays for use by the subroutine,
c           and itrmf is an integer termination flag.  The meaning of
c           itrmf is as follows:
c             0 => normal termination; desired function value calculated.
c             1 => failure to produce f(xcur).
c 
c  rpar    = real parameter/work array passed to the f and jacv routines. 
c
c  ipar    = integer parameter/work array passed to the f and jacv routines. 
c
c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
c
c  itrmbt  = termination flag; values have the following meanings: 
c              0 => normal termination: acceptable step found. 
c              1 => acceptable step not found in ibtmax reductions.
c              2 => error in evaluation of f.
c 
c ------------------------------------------------------------------------
c
c Subroutines required by this and all called routines: 
c
c    user supplied: f
c
c    nitsol routines: none 
c
c    blas routine: dscal
c
c    user supplied or blas: dnorm 
c
c    explanation: In nitsol, dnorm is set to either the blas 
c    dnrm2 routine or the user-supplied usrnrm routine. 
c
c This subroutine called by: nitdrv
c
c Subroutines called by this subroutine: dnorm, dscal, f
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
      double precision t, theta
      integer i, itrmf

      include 'nitparam.h'

      external nitbd

c ------------------------------------------------------------------------ 
c
c ------------------------------------------------------------------------
c Initialize.
c ------------------------------------------------------------------------
      t = 1.d-4
      ibt = 0
      redfac = 1.d0
c ------------------------------------------------------------------------
c Backtracking loop. 
c ------------------------------------------------------------------------
 100  continue
      do 110 i = 1, n
         xpls(i) = xcur(i) + step(i) 
 110  continue
      call f(n, xpls, fpls, rpar, ipar, itrmf) 
      if (itrmf .ne. 0) then
         itrmbt = 2
         go to 900
      endif
      nfe = nfe + 1
      fpnrm = dnorm(n, fpls, 1) 
c ------------------------------------------------------------------------
c If t-condition is met or backtracking is turned off, return. 
c ------------------------------------------------------------------------
      if (fpnrm .le. (1.d0 - t*(1.d0-eta))*fcnrm .or. ibtmax .eq. -1) 
     $ then 
         itrmbt = 0
         go to 900 
      endif
c ------------------------------------------------------------------------
c Otherwise, ... 
c ------------------------------------------------------------------------
      ibt = ibt + 1
      if (ibt .gt. ibtmax) then 
         itrmbt = 1
         go to 900
      endif
c ------------------------------------------------------------------------
c ... choose theta ...
c ------------------------------------------------------------------------
      theta = -(oftjs*redfac)/(fpnrm**2 - fcnrm**2 - 2.d0*oftjs*redfac)
      if(theta .lt. thmin) theta = thmin
      if(theta .gt. thmax) theta = thmax
c ------------------------------------------------------------------------
c ... then reduce the step, increase eta, update redfac ... 
c ------------------------------------------------------------------------
      call dscal(n, theta, step, 1)
      eta = 1.d0 - theta*(1.d0 - eta)
      redfac = theta*redfac
c ------------------------------------------------------------------------
c ... and return to the top of the loop. 
c ------------------------------------------------------------------------
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 4) then 
         if (ibt .eq. 1) then 
            write(ipunit,*)
            write(ipunit,800)
            write(ipunit,*)
         endif
         write(ipunit,810) ibt, fpnrm, theta
      endif
 800  format('nitbt:  Step reduction no., trial F norm, current ',
     $     'reduction factor')
 810  format(5x,i4,2(5x,1pd10.3))
c ------------------------------------------------------------------------
      go to 100
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 3 .and. ibtmax .ne. -1) then 
         write(ipunit,*) 
         if (ibt .eq. 0) then 
            write(ipunit,820) 
         else
            write(ipunit,830) ibt, redfac
         endif
      endif
 820  format( 'nitbt:  no. of step reductions. = 0')
 830  format( 'nitbt:  no. of step reductions. =', i2, 
     $        '    total reduction factor =', 1pd10.3)
c ------------------------------------------------------------------------
      return
      end
