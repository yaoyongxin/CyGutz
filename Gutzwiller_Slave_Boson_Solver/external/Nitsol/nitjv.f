      subroutine nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, 
     $     ifdord, itask, nfe, njve, nrpre, v, z, 
     $     rwork1, rwork2, dnorm, itrmjv)

      implicit none  

      integer n, ipar(*), ijacv, ifdord, itask, nfe, njve, nrpre, 
     $     itrmjv
      double precision xcur(n), fcur(n), rpar(*), v(n), z(n), 
     $     rwork1(n), rwork2(n), dnorm 
      external f, jacv, dnorm

c ------------------------------------------------------------------------
c
c This is nitjv v0.3, the routine for controlling evaluation of products 
c J*v or J*P(inverse)*v or P(inverse)*v, where J is the Jacobian of f 
c and P is a right preconditioning operator. 
c
c ------------------------------------------------------------------------
c 
c Explanation: 
c
c  n       = dimension of the problem.
c
c  xcur    = vector of length n, initial guess on input and final 
c            approximate solution on output. 
c
c  fcur    = vector of length n, value of f at xcur. 
c
c  f       = name of user-supplied subroutine for evaluating the function 
c            the zero of which is sought; this routine has the form 
c
c                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
c
c            where xcur is the array containing the current x value, fcur 
c            is f(xcur) on output, rpar and ipar are, respectively, real 
c            and integer parameter/work arrays for use by the subroutine,
c            and itrmf is an integer termination flag.  The meaning of
c            itrmf is as follows:
c              0 => normal termination; desired function value calculated.
c              1 => failure to produce f(xcur).
c
c  jacv    = name of user-supplied subroutine for evaluating J*v or 
c            P(inverse)*v, where J is the Jacobian of f and P is a 
c            right preconditioning operator. If neither analytic J*v 
c            evaluations nor right preconditioning is used, this can 
c            be a dummy subroutine; if right preconditioning is used but 
c            not analytic J*v evaluations, this need only evaluate 
c            P(inverse)*v. The form is 
c
c            subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
c
c            where xcur and fcur are vectors of length n containing the 
c            current x and f values, ijob is an integer flag indicating 
c            which product is desired, v is a vector of length n to be 
c            multiplied, z is a vector of length n containing the desired 
c            product on output, rpar and ipar are, respectively, real 
c            and integer parameter/work arrays for use by the subroutine, 
c            and itrmjv is an integer termination 
c            flag. The meaning of ijob is as follows: 
c              0 => z = J*v
c              1 => z = P(inverse)*v 
c            The meaning of itrmjv is as follows:
c              0 => normal termination; desired product evaluated. 
c              1 => failure to produce J*v.
c              2 => failure to produce P(inverse)*v. 
c            This subroutine is called only from nitjv, and is always 
c            called with v .ne. z. 
c
c  rpar    = real parameter/work array passed to the f and jacv routines. 
c
c  ipar    = integer parameter/work array passed to the f and jacv routines. 
c
c  ijacv   = flag for determining method of J*v evaluation.
c             -1 => finite-difference evaluation of order ifdord. 
c              0 => first-order finite-difference evaluation. 
c              1 => analytic evaluation. 
c
c  ifdord  = order of the finite-difference formula (sometimes) used when 
c            J*v products are evaluated using finite-differences. When 
c            ijacv = 0 on input to nitsol, ifdord is set to 1, 2, or 4 in 
c            nitsol; subsequently, ijacv may be temporarily reset to -1 in 
c            the Krylov solver to force a finite-difference evaluation of 
c            order ifdord. If ijacv = 1 on input to nitsol, then ifdord is 
c            irrelevant. When ijacv = 0 on input to nitsol, the precise 
c            meaning of ifdord is as follows: 
c
c            If GMRES is used, then ifdord matters only if iresup = 1, in  
c            which case it determines the order of the finite-difference 
c            formula used in evaluating the initial residual at each GMRES 
c            restart (default 2). NOTE: This only affects initial residuals 
c            at restarts; first-order differences are always used within 
c            each GMRES cycle. Using higher-order differences at restarts 
c            only should give the same accuracy as if higher-order 
c            differences were used throughout; see K. Turner and H. F. 
c            Walker, "Efficient high accuracy solutions with GMRES(m)," 
c            SIAM J. Sci. Stat. Comput., 13 (1992), pp. 815--825. 
c               
c            If BiCGSTAB or TFQMR is used, then ifdord determines the 
c            order of the finite-difference formula used at each 
c            iteration (default 1). 
c
c  itask   = flag for determining which product is produced.
c              0 => z = J*v
c              1 => z = J*P(inverse)*v 
c              2 => z = P(inverse)*v 
c
c  nfe     = number of function evaluations.
c
c  njve    = number of J*v evaluations. 
c
c  nrpre   = number of P(inverse)*v evaluations.
c
c  v       = vector to be multiplied. 
c
c  z       = desired product. 
c
c  rwork1  = vector of length n, work array. 
c
c  rwork2  = vector of length n, work array. Note: rwork2 is only referenced 
c            when a J-product is evaluated with a finite-difference of 
c            order 2 or 4, i.e., when itask = 0 or 1, ijacv = 0, and 
c            ifdord = 2 or 4. 
c
c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
c
c  itrmjv  = termination flag; values have the following meanings: 
c              0 => normal termination; desired product evaluated. 
c              1 => failure to produce J*v.
c              2 => failure to produce P(inverse)*v. 
c
c ------------------------------------------------------------------------
c
c Subroutines required by this and all called routines: 
c
c    user supplied: f, jacv 
c
c    nitsol routines: nitfd 
c
c    blas routine: daxpy, dcopy, dscal 
c
c    user supplied or blas: dnorm 
c
c    explanation: In nitsol, dnorm is set to either the blas 
c    dnrm2 routine or the user-supplied usrnrm routine. 
c
c This subroutine called by: nitgm, nitstb, nittfq
c
c Subroutines called by this subroutine: dcopy, jacv, nitfd 
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
      integer ijob 
c
c ------------------------------------------------------------------------ 
c
c ------------------------------------------------------------------------
c If z = J*v is desired (itask = 0), then copy v into rwork1; if 
c z = J*P(inverse)*v or z = P(inverse)*v is desired (itask = 1,2), 
c then compute P(inverse)*v in rwork1. 
c ------------------------------------------------------------------------
      if (itask .eq. 0) then 
         call dcopy(n, v, 1, rwork1, 1)
      else
         ijob = 1
         call jacv(n, xcur, fcur, ijob, v, rwork1, rpar, ipar, itrmjv)
         nrpre = nrpre + 1
         if (itrmjv .ne. 0) go to 900
      endif
c ------------------------------------------------------------------------
c If only z = P(inverse)*v is desired (itask = 2), then copy rwork1 into 
c z and exit.
c ------------------------------------------------------------------------
      if (itask .eq. 2) then 
         call dcopy(n, rwork1, 1, z, 1)
         go to 900
      endif
c ------------------------------------------------------------------------
c If z = J*v or z = J*P(inverse)*v is desired (itask = 0, 1), then 
c compute J*rwork1 in z by either analytic evaluation (ijacv = 1) or 
c finite-differences (ijacv = 0, -1). 
c ------------------------------------------------------------------------
      if (ijacv .eq. 1) then 
         ijob = 0
         call jacv(n, xcur, fcur, ijob, rwork1, z, rpar, ipar, itrmjv)
      else	
         call nitfd(n, xcur, fcur, f, rpar, ipar, ijacv, ifdord, 
     $     nfe, rwork1, z, rwork2, dnorm, itrmjv)
      endif
      njve = njve + 1
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
      return
      end
      subroutine nitfd(n, xcur, fcur, f, rpar, ipar, ijacv, ifdord, 
     $     nfe, v, z, rwork, dnorm, itrmjv)

      implicit none  

      integer n, ipar(*), ijacv, ifdord, nfe, itrmjv
      double precision xcur(n), fcur(n), rpar(*), v(n), z(n), 
     $     rwork(n), dnorm 
      external f, dnorm

c ------------------------------------------------------------------------
c
c This is nitfd v0.3, the routine for finite-difference evaluation of 
c products z = J*v, where J is the Jacobian of f. 
c
c ------------------------------------------------------------------------
c 
c Explanation: 
c
c  n       = dimension of the problem.
c
c  xcur    = vector of length n, initial guess on input and final 
c            approximate solution on output. 
c
c  fcur    = vector of length n, value of f at xcur. 
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
c  rpar    = real parameter/work array passed to the f routine. 
c
c  ipar    = integer parameter/work array passed to the f routine. 
c
c  ijacv   = flag for determining method of J*v evaluation. In this 
c            subroutine, this should be 0 or -1 on input, as follows: 
c             -1 => finite-difference evaluation of order ifdord. 
c              0 => first-order finite-difference evaluation. 
c
c  ifdord  = order of the finite-difference formula used when ijacv = -1.  
c            When ijacv = 0 on input to nitsol, ifdord is set to 1, 2, or 
c            4 in nitsol; subsequently, ijacv may be temporarily reset to 
c            -1 in the Krylov solver to force a finite-difference evaluation 
c            of order ifdord. If ijacv = 1 on input to nitsol, then ifdord 
c            is irrelevant. When ijacv = 0 on input to nitsol, the precise 
c            meaning of ifdord is as follows: 
c
c            If GMRES is used, then ifdord matters only if iresup = 1, in  
c            which case it determines the order of the finite-difference 
c            formula used in evaluating the initial residual at each GMRES 
c            restart (default 2). NOTE: This only affects initial residuals 
c            at restarts; first-order differences are always used within 
c            each GMRES cycle. Using higher-order differences at restarts 
c            only should give the same accuracy as if higher-order 
c            differences were used throughout; see K. Turner and H. F. 
c            Walker, "Efficient high accuracy solutions with GMRES(m)," 
c            SIAM J. Sci. Stat. Comput., 13 (1992), pp. 815--825. 
c               
c            If BiCGSTAB or TFQMR is used, then ifdord determines the 
c            order of the finite-difference formula used at each 
c            iteration (default 1). 
c
c  nfe     = number of function evaluations.
c
c  v       = vector to be multiplied in the product J*v. 
c
c  z       = desired product J*v. 
c
c  rwork   = vector of length n, work array. 
c
c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
c
c  itrmjv  = termination flag; values have the following meanings: 
c              0 => normal termination; desired product evaluated. 
c              1 => failure to produce J*v.
c
c ------------------------------------------------------------------------
c
c Remark: In the selection of the difference step eps in the perturbation 
c x + eps*v for approximating J*v by a finite difference, we assume the 
c following: 
c     1. f and a few derivatives (up to five when ifdord = 4) all 
c        have about the same scale. 
c     2. The relative error in f-evaluations is about epsmach (machine 
c        epsilon). 
c     3. The computed value of the sum of two vectors y and z has error 
c        bounded by epsmach*(||y|| + ||z||). 
c
c The choice of eps is 
c
c             eps = {[(1 + ||x||)*epsmach]**(1/(ifdord+1)} /||v||, 
c
c which approximately minimizes a bound on the relative error in the 
c difference approximation. 
c
c ------------------------------------------------------------------------
c
c Subroutines required by this and all called routines: 
c
c    user supplied: f 
c
c    nitsol routines: none 
c
c    blas routine: daxpy, dscal 
c
c    lapack routine:  dlamch

c    user supplied or blas: dnorm 
c
c    explanation: In nitsol, dnorm is set to either the blas 
c    dnrm2 routine or the user-supplied usrnrm routine. 
c
c This subroutine called by: nitjv
c
c Subroutines called by this subroutine: daxpy, dlamch, dscal, dnorm, f 
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
      double precision eps, epsmach, temp
      integer i, itrmf, ncall
c
      double precision dlamch
      external dlamch
c
c ------------------------------------------------------------------------
      data ncall / 0 /
      save ncall, epsmach 
c ------------------------------------------------------------------------
c Set epsmach (machine epsilon) on first call. 
c ------------------------------------------------------------------------
      if (ncall .eq. 0) epsmach = 2.0d0*dlamch( 'e' )
      ncall = 1
c ------------------------------------------------------------------------
c Compute z = J*v by finite-differences: First, set eps = ||v||for later 
c use in computing the difference step; then evaluate the difference 
c formula according to ijacv and ifdord. 
c ------------------------------------------------------------------------
      eps = dnorm(n, v, 1)
      if (eps .eq. 0.d0) then 
         itrmjv = 1
         go to 900
      endif
c ------------------------------------------------------------------------
c Here ijacv = 0 or ifdord = 1 => first-order forward difference. 
c ------------------------------------------------------------------------
      if (ijacv .eq. 0 .or. ifdord .eq. 1) then 
         eps = dsqrt((1.d0 + dnorm(n,xcur,1))*epsmach)/eps
         do 100 i = 1, n
            v(i) = xcur(i) + eps*v(i)
 100     continue
         call f(n, v, z, rpar, ipar, itrmf)
         if (itrmf .ne. 0) then
            itrmjv = 1
            goto 900
         endif 
         nfe = nfe + 1
         do 110 i = 1, n
            z(i) = (z(i) - fcur(i))/eps
 110     continue
         itrmjv = 0 
         go to 900
      endif
c ------------------------------------------------------------------------
c Here ijacv = -1 and ifdord = 2 => second-order central difference. 
c ------------------------------------------------------------------------
      if (ifdord .eq. 2) then 
         eps = (((1.d0 + dnorm(n,xcur,1))*epsmach)**(1.d0/3.d0))/eps
         do 200 i = 1, n
            rwork(i) = xcur(i) + eps*v(i)
 200     continue
         call f(n, rwork, z, rpar, ipar, itrmf) 
         nfe = nfe + 1
         if (itrmf .ne. 0) then
            itrmjv = 1
            goto 900
         endif 
         do 210 i = 1, n
            rwork(i) = xcur(i) - eps*v(i)
 210     continue
         call f(n, rwork, v, rpar, ipar, itrmf) 
         nfe = nfe + 1
         if (itrmf .ne. 0) then
            itrmjv = 1
            goto 900
         endif 
         temp = 2.d0*eps
         do 220 i = 1, n
            z(i) = (z(i) - v(i))/temp
 220     continue
         itrmjv = 0 
         go to 900
      endif
c ------------------------------------------------------------------------
c Here ijacv = -1 and ifdord = 4 => fourth-order difference. 
c ------------------------------------------------------------------------
      if (ifdord .eq. 4) then 
         eps = (((1.d0 + dnorm(n,xcur,1))*epsmach)**(1.d0/5.d0))/eps
         do 300 i = 1, n
            rwork(i) = xcur(i) + eps*v(i)
 300     continue
         call f(n, rwork, z, rpar, ipar, itrmf) 
         nfe = nfe + 1
         if (itrmf .ne. 0) then
            itrmjv = 1
            goto 900
         endif 
         temp = -eps 
         call daxpy(n, temp, v, 1, xcur, 1)
         call f(n, xcur, rwork, rpar, ipar, itrmf) 
         nfe = nfe + 1
         if (itrmf .ne. 0) then
            itrmjv = 1
            goto 900
         endif 
         do 310 i = 1, n
            z(i) = rwork(i) - z(i)
 310     continue
         temp = eps/2.d0
         call daxpy(n, temp, v, 1, xcur, 1)
         call f(n, xcur, rwork, rpar, ipar, itrmf) 
         nfe = nfe + 1
         if (itrmf .ne. 0) then
            itrmjv = 1
            goto 900
         endif 
         temp = -8.d0
         call daxpy(n, temp, rwork, 1, z, 1)
         call daxpy(n, eps, v, 1, xcur, 1)
         call f(n, xcur, rwork, rpar, ipar, itrmf) 
         nfe = nfe + 1
         if (itrmf .ne. 0) then
            itrmjv = 1
            goto 900
         endif 
         temp = 8.d0
         call daxpy(n, temp, rwork, 1, z, 1)
         temp = 1.d0/(6.d0*eps)
         call dscal(n, temp, z, 1)
         temp = -eps/2.d0
         call daxpy(n, temp, v, 1, xcur, 1)
         itrmjv = 0 
      endif
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
      return
      end
