      subroutine nitgm (n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, 
     $     ipar, ijacv, irpre, iksmax, iresup, ifdord, nfe, njve,  
     $     nrpre, nli, kdmax, kdmaxp1, vv, rr, svbig, svsml, w, rwork, 
     $     rsnrm, dinpr, dnorm, itrmks)

      implicit none 

      integer n, ipar(*), ijacv, irpre, iksmax, iresup, ifdord, nfe, 
     $     njve, nrpre, nli, kdmax, kdmaxp1, itrmks
      double precision xcur(n), fcur(n), fcnrm, step(n), eta, rpar(*), 
     $     vv(n,kdmaxp1), rr(kdmax,kdmax), svbig(kdmax), svsml(kdmax), 
     $     w(kdmax), rwork(n), rsnrm, dinpr, dnorm 
      external f, jacv, dinpr, dnorm 

c ------------------------------------------------------------------------
c
c This is nitgm v0.3, the GMRES routine for determining (trial) inexact 
c Newton steps. This implementation is the "simpler" Gram-Schmidt GMRES 
c implementation from L. Zhou and H. F. Walker, "A simpler GMRES," 
c J. Numerical Lin. Alg. Appl., 1 (1994), pp. 571-581. 
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
c              0 => finite-difference evaluation (default). 
c              1 => analytic evaluation. 
c
c  irpre   = flag for right preconditioning. 
c              0 => no right preconditioning
c              1 => right preconditioning
c
c  iksmax  = maximum allowable number of GMRES iterations. 
c
c  iresup  = residual update flag; on GMRES restarts, the residual is 
c            updated as follows: 
c               0 => linear combination (default) 
c               1 => direct evaluation
c            The first is cheap (one n-vector saxpy) but may lose 
c            accuracy with extreme residual reduction; the second 
c            retains accuracy better but costs one J*v product 
c
c  ifdord  = order of the finite-difference formula (sometimes) used on 
c            GMRES restarts when J*v products are evaluated using finite-
c            differences. When ijacv = 0 on input to nitsol, ifdord is set 
c            to 1, 2, or 4 in nitsol; otherwise, it is irrelevant. When 
c            ijacv = 0 on input to this subroutine, the precise meaning is 
c            as follows: 
c
c            With GMRES, ifdord matters only if iresup = 1, in which case 
c            it determines the order of the finite-difference formula used 
c            in evaluating the initial residual at each GMRES restart 
c            (default 2). If iresup = 1 and ijacv = 0 on input to this 
c            subroutine, then ijacv is temporarily reset to -1 at each 
c            restart below to force a finite-difference evaluation of order 
c            ifdord. NOTE: This only affects initial residuals at restarts; 
c            first-order differences are always used within each GMRES 
c            cycle. Using higher-order differences at restarts only should 
c            give the same accuracy as if higher-order differences were 
c            used throughout; see K. Turner and H. F. Walker, "Efficient 
c            high accuracy solutions with GMRES(m)," SIAM J. Sci. Stat. 
c            Comput., 13 (1992), pp. 815--825. 
c               
c  nfe     = number of function evaluations.
c
c  njve    = number of J*v evaluations. 
c
c  nrpre   = number of P(inverse)*v evaluations.
c
c  nli     = number of linear iterations.
c
c  kdmax   = maximum Krylov subspace dimension; default 10. 
c
c  kdmaxp1 = kdmax + 1. 
c
c  vv      = n x (kdmax+1) matrix for storage of Krylov basis in GMRES;
c            on return, the residual vector is contained in the first 
c            column.
c
c  rr      = kdmax x kdmax matrix for storage of triangular matrix in GMRES. 
c
c  svbig   = vector of length kdmax for storage of estimate of singular 
c            vector of rr with largest singular value. 
c
c  svsml   = vector of length kdmax for storage of estimate of singular 
c            vector of rr with smallest singular value. 
c
c  w       = vector of length kdmax, contains right-hand side of 
c            triangular system and least-squares residual norm in GMRES. 
c
c  rwork   = vector of length n, work array. 
c
c  rsnrm   = GMRES residual norm on return. 
c
c  dinpr   = inner-product routine, either user-supplied or blas ddot. 
c
c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
c
c  itrmks  = termination flag; values have the following meanings: 
c              0 => normal termination: acceptable step found. 
c              1 => J*v failure in nitjv. 
c              2 => P(inverse)*v failure in nitjv. 
c              3 => acceptable step not found in iksmax GMRES iterations. 
c              4 => insufficient residual norm reduction over a cycle 
c                   of kdmax steps (stagnation) before an acceptable step 
c                   has been found. 
c              5 => dangerous ill-conditioning detected before an acceptable 
c                   step has been found. 
c
c             Note: On return, nitsol terminates if itrmks is 1 or 2. If  
c             itrmks is 3, 4, or 5, nitsol may terminate or continue. In 
c             this event, a meaningful inexact Newton step is returned, 
c             even though the desired inexact Newton condition may not 
c             hold, and a decision on termination/continuation is made 
c             in nitdrv according to whether there is sufficient residual 
c             norm reduction.
c
c -------------------------------------------------------------------------
c
c Subroutines required by this and all called routines: 
c
c    user supplied: f, jacv 
c
c    nitsol routines: nitjv, nitfd
c
c    lapack routine: dlaic1, dlamch
c
c    blas routines: daxpy, dcopy, dscal
c
c    user supplied or blas: dinpr, dnorm 
c
c    explanation: In nitsol, dinpr and dnorm are set to either the blas 
c    ddot and dnrm2 routines or the user-supplied usrnpr and usrnrm 
c    routines. 
c
c This subroutine called by: nitdrv
c
c Subroutines called by this subroutine: daxpy, dcopy, dscal, dlaic1, 
c    dlamch, dinpr, dnorm, nitjv
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
      double precision abstol, big, cs, cndmax, epsmach, rsnrm0, 
     $     sestpr, small, sn, temp 
      integer i, igm, ijob, itask, itrmjv, kd, kdp1, ncall

      double precision dlamch
      external dlamch

c ------------------------------------------------------------------------ 
c ------------------------------------------------------------------------
      data ncall / 0 /
      save ncall, cndmax, epsmach 
c ------------------------------------------------------------------------
c Initialize. 
c ------------------------------------------------------------------------
      if (ncall .eq. 0) then
         epsmach = 2.0d0*dlamch( 'e' )
         ncall = 1
         cndmax = 1.d0/(100.d0*epsmach)
      endif
      do 20 i = 1, n
         step(i) = 0.d0
 20   continue
      igm = 0
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 3) then 
         write(ipunit,*) 
         write(ipunit,800) eta 
      endif
 800  format('nitgm:  eta =', 1pd10.3)
      if (iplvl .ge. 4) then 
         write(ipunit,810) 
         write(ipunit,*) 
         write(ipunit,820) igm, fcnrm 
      endif
 810  format('nitgm:  GMRES iteration no., linear residual norm, ',
     $     'condition no. estimate')
 820  format(5x,i4,2(5x,1pd10.3))
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Set the stopping tolerance, etc. 
c ------------------------------------------------------------------------
      rsnrm0 = fcnrm
      abstol = eta*rsnrm0
c ------------------------------------------------------------------------
c Place the normalized initial residual in the first column of the vv array.
c ------------------------------------------------------------------------
      call dcopy(n, fcur, 1, vv(1,1), 1)
      temp = -1.d0/fcnrm
      call dscal(n, temp, vv(1,1), 1)
c ------------------------------------------------------------------------
c Top of the outer GMRES loop. 
c ------------------------------------------------------------------------
 100  continue
      kd = 0
      rsnrm = 1.d0
c ------------------------------------------------------------------------
c Top of the inner GMRES loop.
c ------------------------------------------------------------------------
 200  continue
      kd = kd + 1
      kdp1 = kd + 1
      nli = nli + 1
      igm = igm + 1
c ------------------------------------------------------------------------
c Evaluate J*(kd-th Krylov subspace basis vector) in vv(.,kdp1). 
c Note: rwork can be used for both work arrays in this call because 
c the second is not referenced within nitjv. 
c ------------------------------------------------------------------------
      if (irpre .eq. 0) then 
         itask = 0
      else
         itask = 1
      endif
      call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, 
     $     ifdord, itask, nfe, njve, nrpre, vv(1,kd), vv(1,kdp1), 
     $     rwork, rwork, dnorm, itrmjv)
      if (itrmjv .gt. 0) then 
         if (itrmjv .eq. 1) itrmks = 1
         if (itrmjv .eq. 2) itrmks = 2
         go to 900
      endif
c ------------------------------------------------------------------------
c Do modified Gram-Schmidt. 
c ------------------------------------------------------------------------
      do 210 i = 2, kd
         rr(i-1,kd) = dinpr(n, vv(1,i), 1, vv(1,kdp1), 1)
         call daxpy(n, -rr(i-1,kd), vv(1,i), 1, vv(1,kdp1), 1)
 210  continue
      rr(kd,kd) = dnorm(n, vv(1,kdp1), 1)
c ------------------------------------------------------------------------
c Update the estimates of the largest and smallest singular values. 
c ------------------------------------------------------------------------
      if (kd .eq. 1) then
         big = rr(1,1)
         small = big
         svbig(1) = 1.d0
         svsml(1) = 1.d0
      else
         ijob = 1
         call dlaic1(ijob, kd-1, svbig, big, rr(1,kd), rr(kd,kd),
     $        sestpr, sn, cs)
         big = sestpr
         call dscal(kd-1, sn, svbig, 1)
         svbig(kd) = cs
         ijob = 2
         call dlaic1(ijob, kd-1, svsml, small, rr(1,kd), rr(kd,kd),
     $        sestpr, sn, cs)
         small = sestpr
         call dscal(kd-1, sn, svsml, 1)
         svsml(kd) = cs
      endif
c ------------------------------------------------------------------------
c Terminate if the estimated condition number is too great. 
c ------------------------------------------------------------------------
      if (big .ge. small*cndmax) then 
         if (kd .eq. 1) then 
            itrmks = 5
            go to 900
         else 
            kdp1 = kd
            kd = kd - 1
            call daxpy(n, w(kd), vv(1,kdp1), 1, vv(1,1), 1)
            go to 300
         endif
      endif
c ------------------------------------------------------------------------
c Normalize vv(.,kdp1). 
c ------------------------------------------------------------------------
      temp = 1.d0/rr(kd,kd) 
      call dscal(n, temp, vv(1,kdp1), 1)
c ------------------------------------------------------------------------
c Update w and the residual norm by rsnrm <- rsnrm*dsin(dacos(w(kd)/rsnrm). 
c ------------------------------------------------------------------------
      w(kd) = dinpr(n, vv(1,1), 1, vv(1,kdp1), 1)
      temp = max(min(w(kd)/rsnrm,1.0D0),-1.0d0) 
      rsnrm = rsnrm*dsin(dacos(temp))
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 4) then 
         write(ipunit,820) igm, rsnrm*rsnrm0, big/small
c         if (kd .eq. kdmax) write(ipunit,*) 
      endif
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Test for termination of the inner loop.
c ------------------------------------------------------------------------
      if ( (rsnrm0*rsnrm .le. abstol) .or. (kd .eq. kdmax) .or.  
     $     (igm .ge. iksmax) )  go to 300
c ------------------------------------------------------------------------
c If not terminating the inner loop, update the residual vector 
c and go to the top of the inner loop. 
c ------------------------------------------------------------------------
      call daxpy(n, -w(kd), vv(1,kdp1), 1, vv(1,1), 1)
      go to 200
c ------------------------------------------------------------------------
c Bottom of inner loop.
c ------------------------------------------------------------------------
 300  continue
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 4) then 
         write(ipunit,*) 
      endif
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Compute the solution: 
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Use svbig for storage of the original components of w. 
c ------------------------------------------------------------------------
      call dcopy(kd, w, 1, svbig, 1)
c ------------------------------------------------------------------------
c Overwrite w with the solution of the upper triangular system.
c ------------------------------------------------------------------------
      do 310 i = kd, 1, -1
         w(i) = w(i)/rr(i,i)
         if (i .gt. 1) call daxpy(i-1, -w(i), rr(1,i), 1, w, 1)
 310  continue
c ------------------------------------------------------------------------
c Now form the linear combination to accumulate the correction in 
c the work vector.
c ------------------------------------------------------------------------
      call dcopy(n, vv(1,1), 1, rwork, 1)
      call dscal(n, w(1), rwork, 1)
      if (kd .gt. 1) then 
         call daxpy(kd-1, w(1), svbig, 1, w(2), 1)
         do 320 i = 2, kd
            call daxpy(n, w(i), vv(1,i), 1, rwork, 1)
 320     continue
      endif
c ------------------------------------------------------------------------
c If iresup .eq. 0, then update the residual vector by linear 
c combination. This frees vv(.,kdp1) for use as a work array. 
c ------------------------------------------------------------------------
      if (iresup .eq. 0) then 
         call daxpy(n, -svbig(kd), vv(1,kdp1), 1, vv(1,1), 1)
      endif
c ------------------------------------------------------------------------
c If right preconditioning is used, overwrite 
c correction <-- P(inverse)*correction, using vv(.,kdp1) as a work array. 
c Note: vv(.,kdp1) can be used for both work arrays in this call because 
c the second is not referenced within nitjv. 
c ------------------------------------------------------------------------
      if (irpre .gt. 0) then 
         itask = 2
         call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, 
     $        ifdord, itask, nfe, njve, nrpre, rwork, rwork, 
     $        vv(1,kdp1), vv(1,kdp1), dnorm, itrmjv)
         if (itrmjv .gt. 0) then 
            itrmks = 2
            go to 900
         endif
      endif
c ------------------------------------------------------------------------
c Update the step. This frees rwork for use as a work array.
c ------------------------------------------------------------------------
      call daxpy(n, rsnrm0, rwork, 1, step, 1)
c ------------------------------------------------------------------------
c If iresup .eq. 1, then update the residual vector by direct evaluation, 
c using rwork and vv(.,kdp1) as work arrays. Note: Two distinct work  
c arrays are needed in this call because both are referenced within nitjv 
c if the J*step product is evaluated with a finite-difference of order 
c two or higher. If finite-differences are used (ijacv= 0), then ijacv 
c is temporarily set to -1 to signal to nitjv that the order of the 
c finite-difference formula is to be determined by ifdord. 
c ------------------------------------------------------------------------
      if (iresup .eq. 1) then 
         itask = 0
         if (ijacv .eq. 0) ijacv = -1 
         call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, 
     $        ifdord, itask, nfe, njve, nrpre, step, vv(1,1), 
     $        rwork, vv(1,kdp1), dnorm, itrmjv)
         if (ijacv .eq. -1) ijacv = 0
         if (itrmjv .gt. 0) then 
            itrmks = 1
            go to 900
         endif
         do 330 i = 1, n
            vv(i,1) = -fcur(i) - vv(i,1)
 330     continue
      endif
c ------------------------------------------------------------------------
c Test for termination.  
c ------------------------------------------------------------------------
      if (rsnrm0*rsnrm .le. abstol) then 
         itrmks = 0
         go to 900
      endif
      if (igm .ge. iksmax) then 
         itrmks = 3
         go to 900
      endif
      if (big .ge. small*cndmax) then 
         itrmks = 5
         go to 900
      endif
      temp = dfloat(kd)*dlog(abstol/(rsnrm0*rsnrm))/
     $     dlog(rsnrm/(1.d0 + 10.d0*epsmach))
      if (temp .ge. 1000.d0*dfloat(iksmax - igm)) then 
         itrmks = 4
         go to 900
      endif
c ------------------------------------------------------------------------
c If not terminating, then normalize the initial residual, etc., and 
c return to the top of the outer loop. 
c ------------------------------------------------------------------------
      if (iresup .eq. 0) then 
         rsnrm0 = rsnrm0*rsnrm
         temp = 1.d0/rsnrm
      else
         rsnrm0 = dnorm(n, vv(1,1), 1)
         temp = 1.d0/rsnrm0
      endif
      call dscal(n, temp, vv(1,1), 1)
      go to 100
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
      if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
         if (iresup .eq. 0) then 
            call dscal(n, rsnrm0, vv(1,1), 1)
            rsnrm = rsnrm0*rsnrm
         else
            rsnrm = dnorm(n, vv(1,1), 1)
         endif
      endif
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 3) then 
         if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
            write(ipunit,830) itrmks, rsnrm 
 830        format('nitgm:  itrmks =', i2, '    final lin. res. norm =', 
     $           1pd10.3)
         else
            write(ipunit,840) itrmks
 840        format('nitgm: itrmks:', i4) 
         endif
      endif
c ------------------------------------------------------------------------
      return
      end
