      subroutine nitdrv(n, xcur, fcur, xpls, fpls, 
     $     step, f, jacv, rpar, ipar, ftol, stptol, nnimax, 
     $     ijacv, ikrysl, kdmax, irpre, iksmax, iresup, ifdord, 
     $     ibtmax, ieta, iterm, nfe, njve, nrpre, nli, nni, nbt, 
     $     rwork, dinpr, dnorm)

      implicit none  

      integer n, ipar(*), nnimax, ijacv, ikrysl, kdmax, irpre, 
     $     iksmax, iresup, ifdord, ibtmax, ieta, iterm, nfe, njve, 
     $     nrpre, nli, nni, nbt
      double precision xcur(n), fcur(n), xpls(n), fpls(n), 
     $     step(n), rpar(*), ftol, stptol, rwork(*), dinpr, dnorm
      external f, jacv, dinpr, dnorm

c ------------------------------------------------------------------------
c
c This is nitdrv v0.3, the driver routine for the Newton iterative 
c method.  
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
c  xpls    = vector of length n, next (trial) approximate solution. 
c            Also used as a work array where indicated.  
c
c  fpls    = vector of length n, value of f at xpls. Also used as a 
c            work array where indicated.  
c
c  step    = vector of length n, (trial) step. 
c
c  f       = name of user-supplied subroutine for evaluating the function 
c            the zero of which is sought; this routine has the form
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
c            and itrmjv is an integer termination flag. The meaning of 
c            ijob is as follows: 
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
c  ftol    = stopping tolerance on the f-norm.
c
c  stptol  = stopping tolerance on the steplength.
c
c  nnimax  = maximum allowable number of nonlinear iterations (default 200). 
c
c  ijacv   = flag for determining the method of J*v evaluation: 
c              0 => finite-difference evaluation (default) 
c              1 => analytic evaluation 
c
c  ikrysl  = flag for determining the Krylov solver: 
c              0 => GMRES (default) 
c              1 => BiCGSTAB
c              2 => TFQMR
c
c            For brief descriptions of the solvers plus references, 
c            see the subroutines nitgm, nitstb, and nittfq. 
c
c  kdmax   = maximum Krylov subspace dimension when GMRES is used 
c            (default 20).
c
c  irpre   = flag for right preconditioning: 
c              0 => no right preconditioning
c              1 => right preconditioning
c
c  iksmax  = maximum allowable number of iterations per call to the Krylov 
c            solver (default 1000). 
c
c  iresup  = residual update flag when GMRES is used; on GMRES restarts, 
c            the residual is updated as follows: 
c              0 => linear combination (default) 
c              1 => direct evaluation
c            The first is cheap (one n-vector saxpy) but may lose 
c            accuracy with extreme residual reduction; the second 
c            retains accuracy better but costs one J*v product per 
c            restart. 
c
c  ifdord  = order of the finite-difference formula (sometimes) used when 
c            ijacv = 0. When ijacv = 0, this is set to 1, 2, or 4 in nitsol; 
c            otherwise it is irrelevant. With ijacv = 0, the precise meaning 
c            is as follows: 
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
c  ibtmax  = maximum allowable number of backtracks (step 
c            reductions) per call to nitbt (default 10). 
c
c            USAGE NOTE: Backtracking can be turned off by setting 
c            ibtmax = -1. Other negative values of ibtmax are not 
c            valid. 
c
c
c  ieta    = flag determining the forcing term eta as follows: 
c               0 => (||fcur|| - ||fprev+Jprev*sprev||)/||fprev||
c                    (default)
c               1 => (||fcur||/||fprev||)**2
c               2 => gamma*(||fcur||/||fprev||)**alpha  
c                    for user-supplied gamma in (0,1] and alpha in (1,2] 
c               3 => user-supplied eta in [0,1). 
c                 3 => fixed (constant) eta in (0,1), either 0.1 (default) 
c		       or specified by the user (see USAGE NOTE below) 
c            Here, fcur = current f, fprev = previous f, etc. The Krylov 
c            iterations are terminated when an iterate s satisfies 
c            an inexact Newton condition ||F + J*s|| .le. eta*||F||.
c
c            USAGE NOTE: If ieta = 2, then alpha and gamma must be set 
c            in common block nitparam.h as described below. If 
c	     ieta = 3, then the desired constant eta may be similarly 
c	     set in nitparam.h if a value other than the default of 
c	     0.1 is desired. 
c               
c            The first three expressions above are from S. C. Eisenstat 
c            and H. F. Walker, "Choosing the forcing terms in an inexact 
c            Newton method", SIAM J. Scientific Computing, 17 (1996), 
c            pp. 16--32. (They may be modified according to certain 
c            safeguards below.) The first gives convergence that is 
c            q-superlinear and of r-order (1+sqrt(5))/2. The second gives 
c            convergence that is r-quadratic and of q-order p for every p 
c            in [1,2). The third gives convergence that is of q-order alpha 
c            when gamma < 1 and, when gamma = 1, of r-order alpha and 
c            q-order p for every p in [1,alpha). The fourth gives q-linear 
c            convergence with asymptotic rate constant alpha in a certain 
c            norm; see R. S. Dembo, S. C. Eisenstat, and T. Steihaug, 
c            "Inexact Newton methods", SIAM J. Numer. Anal., 18 (1982), 
c            pp. 400-408. 
c
c            Of these four choices, the 1st is usually satisfactory, 
c            the 2nd or 3rd is sometimes preferred, and the 4th may be 
c            useful in some situations, e.g., it may be desirable to 
c            choose a fairly large fixed eta in (0,1), such as eta = .1, 
c            when numerical inaccuracy prevents the Krylov solver 
c            from obtaining much residual reduction. 
c               
c  iterm   = termination flag; values have the following meanings: 
c              0 => normal termination: ||F||.le.ftol or ||step||.le.stptol.
c              1 => nnimax nonlinear iterations reached without success. 
c              2 => failure to evaluate F.
c              3 => in nitjv, J*v failure. 
c              4 => in nitjv, P(inverse)*v failure. 
c              5 => in nitdrv, insufficient initial model norm reduction 
c                   for adequate progress. NOTE: This can occur for several 
c                   reasons; examine itrmks on return from the Krylov 
c                   solver for further information. (This will be printed out 
c                   if iplvl .ge. 3, see the discussion of optional 
c                   common blocks below). 
c              6 => in nitbt, failure to reach an acceptable step through 
c                   backtracking. 
c 
c  nfe     = number of function evaluations.
c
c  njve    = number of J*v evaluations. 
c
c  nrpre   = number of P(inverse)*v evaluations.
c
c  nli     = number of linear iterations.
c
c  nni     = number of nonlinear iterations.
c
c  nbt     = number of backtracks. 
c
c  rwork   = real work vector for use by the Krylov solver. It is passed 
c            in as the tail of the rwork vector in nitsol. On input to 
c            nitsol, it should have length as follows: 
c
c             solver    rwork length
c             GMRES     n*(kdmax+5)+kdmax*(kdmax+3), where kdmax is the 
c                       maximum Krylov subspace dimension, either the 
c                       default value of 20 or another value specified 
c                       by the user. 
c             BiCGSTAB  11*n
c             TFQMR     14*n
c
c  dinpr   = inner-product routine, either user-supplied or BLAS ddot. 
c
c  dnorm   = norm routine, either user-supplied or BLAS dnrm2. 
c
c ------------------------------------------------------------------------
c
c Optional common blocks: 
c
c These can be used to control printing of diagnostic information by nitsol, 
c to pass information about the nonlinear iterations to jacv or other user 
c subroutines, or to control the default behavior of the nonlinear iterations. 
c
c For controlling printing of diagnostic information: 
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
c     ipunit = printout unit number, e.g., ipunit = 6 => standard output. 
c
c For passing information about the nonlinear iterations to user-supplied 
c subroutines: 
c
      include 'nitinfo.h'
c
c If information on the current state of the nonlinear iteration is
c desired in a user-supplied subroutine (for example, deciding 
c whether to update a preconditioner), include this common block
c in the subroutine. The variables are as follows: 
c
c     instep - inexact Newton step number. 
c
c    newstep - set to 0 at the beginning of an inexact Newton step.
c              This may be checked in a user-supplied jacv to decide
c              whether to update the preconditioner.  If you test on
c              newstep .eq. 0 to determine whether to take some 
c              special action at the beginning of a nonlinear iteration, 
c              you must also set newstep to some nonzero value to
c              subsequently avoid taking that action unnecessarily. 
c
c    krystat - status of the Krylov iteration; same as itrmks (see 
c              the nitsol documentation). 
c
c    avrate  - average rate of convergence of the Krylov solver during
c              the previous inexact Newton step.  This may be checked
c              in a user-supplied jacv to decide when to update the
c              preconditioner.
c
c    fcurnrm - ||f(xcur)||. 
c
c
c  For controlling the default behavior of the nonlinear iterations:
c
      include 'nitparam.h'

c nitparam contains some parameters that control the nonlinear
c iterations.  In some cases, the default values reflect prevailing  
c practice; in other cases, they are chosen to produce good 
c average-case behavior.  To change the default values, include this 
c common block in the main program and set the desired variables 
c according to the following:
c
c    choice1_exp -  parameter used in the update of the forcing term 
c                   eta when ieta = 0 (default).  This is the exponent
c                   for determining the etamin safeguard.  The default
c                   value is choice1_exp = (1+sqrt(5))/2.  A larger
c                   value will allow eta to decrease more rapidly,
c                   while a smaller value will result in a larger 
c                   value for the safeguard. 
c
c    choice2_exp  - parameter used in the update of the forcing term 
c                   eta when ieta = 2.  This is the exponent alpha 
c		    in the expression gamma*(||fcur||/||fprev||)**alpha; 
c		    it is also used to determine the etamin safeguard.  
c		    The default value is 2.0. Valid values are in the 
c		    range (1.0, 2.0].
c
c    choice2_coef - parameter used in the update of the forcing term eta 
c                   when ieta = 2.  This is the coefficient gamma used 
c		    in the expression gamma*(||fcur||/||fprev||)**alpha;
c                   it is also used to determine the etamin safeguard.
c                   The default value is 1.0. Valid values are in the 
c		    range (0.0, 1.0]. 
c
c    eta_cutoff   - parameter used to determine when to disable 
c                   safeguarding the update of the forcing term.  It
c                   only has meaning when ieta .ne. 3.  The default
c                   value is 0.1.  A value of 0.0 will enable 
c		    safeguarding always; a value of 1.0 will disable 
c		    safeguarding always. 
c
c    etamax       - parameter used to provide an upper bound on the 
c		    forcing terms when input(10) .ne. 3. This is 
c		    necessary to ensure convergence of the inexact Newton 
c		    iterates and is imposed whenever eta would otherwise 
c		    be too large. (An overly large eta can result from 
c		    the updating formulas when input(10) .ne. 3 or from 
c                   safeguarding when the previous forcing term has been 
c		    excessively increased during backtracking.) The 
c		    default value of etamax is 1.0 - 1.e-4.  When 
c		    backtracking occurs several times during a nonlinear 
c		    solve the forcing term can remain near etamax for several
c                   nonlinear steps and cause the nonlinear iterations
c                   to nearly stagnate.  In such cases a smaller value of 
c                   etamax may prevent this.  Valid values are in the 
c                   range (0.0, 1.0).
c
c    etafixed     - this is the user-supplied fixed eta when ieta = 3.
c                   The  default value is etafixed = 0.1.  Valid values
c                   are in the range (0.0,1.0).
c
c    thmin        - when backtracking occurs, this is the smallest
c                   reduction factor that will be applied to the current
c                   step in a single backtracking reduction.  The default
c                   value is 0.1.  Valid  values are in the range
c                   [0.0, thmax].
c
c    thmax        - when backtracking occurs, this is the largest
c                   reduction factor that will be applied to the current
c                   step in a single backtracking reduction.  The default
c                   value is 0.5.  Valid values are in the range
c                   [thmin, 1.0).
c
c  The values in this common block are not checked here.  We assume
c  that if you call nitdrv directly you know what you are doing.
c
c ------------------------------------------------------------------------
c
c Subroutines required by this and all called routines:
c
c    user supplied: f, jacv
c
c    nitsol routines: nitbd.f, nitbt, nitgm, nitjv, nitstb, nittfq, 
c
c    lapack routines: dlaic1, dlamch
c
c    blas routines: daxpy, dcopy, dscal, dswap 
c
c    user supplied or BLAS (see below): dinpr, dnorm 
c
c    Explanation: The call to nitsol specifies dinpr and dnorm as 
c    either user-supplied inner-product and norm routines or the 
c    BLAS ddot and dnrm2 routines. 
c
c This subroutine called by: nitsol
c
c Subroutines called by this subroutine: dcopy, dinpr, dlamch, dnorm, f,
c    nitbt, nitgm, nitstb, nittfq 
c
c ------------------------------------------------------------------------
c
c Remaining declarations: 
c
c NOTE: In nitinfo.h, instep, newstep, krystat are declared integer, 
c and avrate and fcurnrm are declared double precision. 
c 
      double precision alpha, epsmach, eta, etamin, fcnrm, 
     $     flmnrm, fpnrm, gamma, oftjs, oftlm, redfac, rsnrm, 
     $     stpnrm, temp 
      integer ibt, itrmbt, itrmf, itrmks, lrr, lsvbig, lsvsml, lvv, lw, 
     $     lr, ld, lrcgs, lrtil, lp, lphat, lq, lu, lv, lt, lrwork, ly, 
     $     kdmaxp1, ncall 

      double precision dlamch
      external dlamch

      external nitbd

c ------------------------------------------------------------------------
c For printing:
      integer infe, injve, inrpre, inli
c ------------------------------------------------------------------------
      data ncall / 0 /
      save ncall, alpha, gamma, epsmach 
c ------------------------------------------------------------------------
c Initialize.
c ------------------------------------------------------------------------
      if ( ncall .eq. 0 ) epsmach = 2.0d0*dlamch( 'e' )
      ncall = 1
      if (ieta .eq. 0) alpha = choice1_exp
      if (ieta .eq. 2) alpha = choice2_exp
      if (ieta .eq. 2) gamma = choice2_coef
      if (ieta .eq. 3) then 
         eta = etafixed
      else
         eta = .5d0
      endif
      nfe = 0
      njve = 0
      nrpre = 0
      nli = 0
      nni = 0
      nbt = 0
      avrate = 1.0d0
c ------------------------------------------------------------------------
c Evaluate f at initial x and initialize eta. 
c ------------------------------------------------------------------------
      call f(n, xcur, fcur, rpar, ipar, itrmf)
      if ( itrmf .ne. 0 ) then
         iterm = 2
         go to 900
      endif
      nfe = nfe + 1
      fcnrm = dnorm(n, fcur, 1) 
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 1) then 
         write(ipunit,*) 
         write(ipunit,*) 
         write(ipunit,*) 'nitdrv:  Beginning nonlinear iterations.'
      endif
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Nonlinear iteration loop.
c      tstrt = etime(dummy)
c ------------------------------------------------------------------------
 100  continue
c ------------------------------------------------------------------------ 
c For printing:
c      if (iplvl .eq. 0) then
c         if ( nni .eq. 0 ) then
c            write(ipunit,799) nni, njve, 0.0d0, fcnrm
c         else 
c            write(ipunit,799) nni, njve, etime(dummy)-tstrt, fcnrm
c         endif
c      endif
c 799  format( 2(2x,i4),2(2x,1pe12.3) )
      if (iplvl .ge. 1) then 
         write(ipunit,*)
         if (iplvl .eq. 1) write(ipunit,800) nni, fcnrm
         if (iplvl .ge. 2) write(ipunit,810) nni, fcnrm 
      endif
 800  format('    It. no.', i4, '      F norm =', 1pd12.3)
 810  format('--- It. no.', i4, '      F norm =', 1pd12.3,
     $     ' ---------------------------')
      if (iplvl .ge. 2) then 
         write(ipunit,*) 
         write(ipunit,820) nfe, njve, nrpre, nli
         infe = nfe
         injve = njve
         inrpre = nrpre 
         inli = nli
      endif
 820  format('  Initial totals:  nfe =', i4, '    njve =', i4, 
     $     '    nrpre =', i4, '    nli =', i4)
c ------------------------------------------------------------------------
c
c------------------------------------------------------------------------
c Test for stopping. 
c------------------------------------------------------------------------
      if (fcnrm .le. ftol) then 
         iterm = 0
         go to 900
      endif
      if (nni .gt. 0 .and. stpnrm .le. stptol .and. itrmks .eq. 0) then 
         iterm = 0
         go to 900
      endif
      if (nni .ge. nnimax) then 
         iterm = 1
         go to 900
      endif
c ------------------------------------------------------------------------
c Compute the (trial) inexact Newton step with the Krylov solver. 
c ------------------------------------------------------------------------
c Update data in nitinfo to mark the start of a new inexact Newton step.
c ------------------------------------------------------------------------
      newstep = 0
      instep = nni
      fcurnrm = fcnrm
c ------------------------------------------------------------------------
c If ikrysl = 0, apply GMRES, using fpls as a work array. 
c ------------------------------------------------------------------------
      if (ikrysl .eq. 0) then 
         kdmaxp1 = kdmax + 1 
         lvv = 1
         lrr = lvv + n*kdmaxp1
         lsvbig = lrr + kdmax*kdmax 
         lsvsml = lsvbig + kdmax
         lw = lsvsml + kdmax
         call nitgm(n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, 
     $        ipar, ijacv, irpre, iksmax, iresup, ifdord, nfe, njve,  
     $        nrpre, nli, kdmax, kdmaxp1, rwork(lvv), rwork(lrr), 
     $        rwork(lsvbig), rwork(lsvsml), rwork(lw), fpls, 
     $        rsnrm, dinpr, dnorm, itrmks)
      endif
c ------------------------------------------------------------------------
c If ikrysl = 1, apply BiCGSTAB, using fpls as a work array. 
c ------------------------------------------------------------------------
      if (ikrysl .eq. 1) then 
         lr = 1
         lrtil = lr + n
         lp = lrtil + n
         lphat = lp + n
         lv = lphat + n
         lt = lv + n
         lrwork = lt + n
         call nitstb (n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, 
     $        ipar, ijacv, irpre, iksmax, ifdord, nfe, njve, 
     $        nrpre, nli, rwork(lr), rwork(lrtil), rwork(lp), 
     $        rwork(lphat), rwork(lv), rwork(lt), rwork(lrwork), fpls, 
     $        rsnrm, dinpr, dnorm, itrmks)
      endif
c ------------------------------------------------------------------------
c If ikrysl = 2, apply TFQMR 
c ------------------------------------------------------------------------
      if (ikrysl .eq. 2) then 
         lr = 1
         lrcgs = lr + n
         lrtil = lrcgs + n
         ld = lrtil + n
         lp = ld + n
         lq = lp + n
         lu = lq + n
         lv = lu + n
         ly = lv + n
         lrwork = ly + n
         call nittfq (n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, 
     $        ipar, ijacv, irpre, iksmax, ifdord, nfe, njve, nrpre, nli,
     $        rwork(lr), rwork(lrcgs), rwork(lrtil), rwork(ld), 
     $        rwork(lp), rwork(lq), rwork(lu), rwork(lv), rwork(ly),
     $        rwork(lrwork), fpls, rsnrm, dinpr, dnorm, itrmks)
      endif
c ------------------------------------------------------------------------
c  Set values in nitinfo that reflect state of iterative solver.
c ------------------------------------------------------------------------
      krystat = itrmks
      avrate = (rsnrm/fcnrm)**(1.0d0/dble(nli))
c ------------------------------------------------------------------------
c Check itrmks and decide whether to terminate or continue: 
c        0 => continue, inexact Newton condition successfully met 
c   1 or 2 => terminate unconditionally, J*v or P(inverse)*v failure) 
c   .ge. 3 => terminate if the model norm has increased or if reduction at 
c             the current rate would at best require more than 1000 time 
c             the maximum remaining number of nonlinear iterations. 
c ------------------------------------------------------------------------
      if (itrmks .eq. 1 .or. itrmks .eq. 2) then
         iterm = itrmks + 2
         go to 900
      endif
      if (itrmks .ge. 3) then 
         if (rsnrm/fcnrm .ge. 1.d0) then 
            iterm = 5
            go to 900
         else
            temp = dlog(ftol/fcnrm)/
     $           dlog(rsnrm/((1.d0 + 10.d0*epsmach)*fcnrm))
            if (temp .gt. 1000.d0*dfloat(nnimax - nni)) then 
               iterm = 5
               go to 900
            endif
         endif
      endif
c ------------------------------------------------------------------------
c Compute the original value of f(transpose)*Js for backtracking; the 
c original value of f(transpose)*(linear model) is also computed for 
c later use. NOTE: The first n components of rwork contain the residual  
c vector for the Newton equation, which is -(linear model). 
c ------------------------------------------------------------------------
      oftlm = -dinpr(n, fcur, 1, rwork, 1)
      oftjs = oftlm - fcnrm**2
c ------------------------------------------------------------------------
c Determine an acceptable step via backtracking. 
c ------------------------------------------------------------------------
      call nitbt(n, xcur, fcnrm, step, eta, xpls, fpls, fpnrm, oftjs, 
     $     redfac, nfe, ibt, ibtmax, f, rpar, ipar, dnorm, itrmbt)
      if (itrmbt .eq. 1) then
         iterm = 6
         go to 900
      else if (itrmbt .eq. 2) then
         iterm = 2
         go to 900
      endif
      nbt = nbt + ibt
c ------------------------------------------------------------------------
c Set eta for next iteration. 
c ------------------------------------------------------------------------
      if (ieta .eq. 0) then 
         etamin = eta**alpha
         temp = 1.d0 - redfac
         flmnrm = dsqrt((temp*fcnrm)**2 + 2.d0*redfac*temp*oftlm + 
     $        (redfac*rsnrm)**2)
         eta = dabs(fpnrm - flmnrm)/fcnrm 
      endif
      if (ieta .eq. 1) then 
         etamin = eta**2
         eta = (fpnrm/fcnrm)**2
      endif
      if (ieta .eq. 2) then 
         etamin = gamma*eta**alpha
         eta = gamma*(fpnrm/fcnrm)**alpha
      endif
      if (ieta .ne. 3) then 
         if (etamin .le. eta_cutoff) etamin = 0.d0
chfw         if (etamin .le. 1.d-1) etamin = 1.d-1
         if (eta .lt. etamin) eta = etamin 
         if (eta .gt. etamax) eta = etamax 
         if (eta*fpnrm .le. 2.d0*ftol) eta = (.8d0*ftol)/fpnrm
      endif
      if (ieta .eq. 3) eta = etafixed
c ------------------------------------------------------------------------
c Update xcur, fcur, fcnrm, stpnrm, nni for next iteration.
c ------------------------------------------------------------------------
      call dcopy(n, xpls, 1, xcur, 1)
      call dcopy(n, fpls, 1, fcur, 1)
      fcnrm = fpnrm
      stpnrm = dnorm(n, step, 1)
      nni = nni + 1
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 2) then 
         infe = nfe - infe 
         injve = njve - injve 
         inrpre = nrpre  - inrpre 
         inli = nli - inli 
         if (ieta .gt. 0) then 
            temp = 1.d0 - redfac
            flmnrm = dsqrt((temp*fcnrm)**2 + 2.d0*redfac*temp*oftlm + 
     $           (redfac*rsnrm)**2)
         endif
         write(ipunit,*) 
         write(ipunit,830) infe, injve, inrpre, inli, stpnrm, flmnrm 
      endif
 830  format('  At this step:   nfe =', i4, '    njve =', i4, 
     $     '    nrpre =', i4, '    nli =', i4, /, 17x, 
     $     ' step norm =', 1pd10.3, 5x, 
     $     'final lin. model norm =', 1pd10.3)
c ------------------------------------------------------------------------
c
c ------------------------------------------------------------------------
c Return to top of loop for next iteration.
c ------------------------------------------------------------------------
      go to 100
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
c ------------------------------------------------------------------------ 
c For printing:
      if (iplvl .ge. 1) then 
         write(ipunit,*) 
         write(ipunit,*) 
         write(ipunit,*) 'nitdrv:  Terminating nonlinear iterations.'
      endif
c ------------------------------------------------------------------------
      return
      end
