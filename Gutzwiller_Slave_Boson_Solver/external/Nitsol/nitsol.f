      subroutine nitsol(n, x, f, jacv, ftol, stptol, 
     $     input, info, rpar, ipar, iterm)

      implicit none  

      integer n, input(10), info(6), ipar(*), iterm 
      double precision x(n), ftol, stptol, rpar(*) 
      real(8),allocatable :: rwork(:)
      real(8) ddot, dnrm2
      external f, jacv, ddot, dnrm2

c
c ------------------------------------------------------------------------
c
c This is nitsol v0.3. The basic algorithm is an inexact Newton method with 
c a backtracking globalization; the model is Algorithm INB of S. C. Eisenstat 
c and H. F. Walker, "Globally convergent inexact Newton methods", SIAM J. 
c Optimization, 4 (1994), pp. 393--422. Initial inexact Newton steps are 
c obtained by approximately solving the Newton equation with a transpose-free
c Krylov subspace method; the current choices are GMRES, BiCGSTAB, and 
c TFQMR. Jacobian-vector products are evaluated either by finite-difference 
c approximation or a user-supplied analytic-evaluation subroutine. An option 
c is provided for user-supplied right preconditioning. Left preconditioning 
c is not explicitly included as an option, but the user may provide this 
c in the subroutines for evaluating the function and Jacobian-vector 
c products. Various algorithmic options can be selected through the input 
c vector. Optional common blocks are also available for printing diagnostic 
c information, passing information about the nonlinear iterations to user 
c subroutines, and controlling the behavior of the nonlinear iterations. 
c Summary statistics are provided by the info vector on output. 
c
c This is the interface subroutine, which calls the driver subroutine nitdrv. 
c
c ------------------------------------------------------------------------
c
c Explanation: 
c
c  n      = dimension of the problem. 
c
c  x      = vector of length n, initial guess on input and final approximate 
c           solution on output. 
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
c  jacv   = name of user-supplied subroutine for optionally evaluating J*v 
c           or P(inverse)*v, where J is the Jacobian of f and P is a 
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
c           and itrmjv is an integer termination flag. The meaning of 
c           ijob is as follows: 
c             0 => z = J*v
c             1 => z = P(inverse)*v 
c           The meaning of itrmjv is as follows:
c             0 => normal termination; desired product evaluated. 
c             1 => failure to produce J*v.
c             2 => failure to produce P(inverse)*v. 
c           This subroutine is called only from nitjv, and is always 
c           called with v .ne. z. 
c
c  ftol   = stopping tolerance on the f-norm.
c
c  stptol = stopping tolerance on the steplength.
c
c  input  = integer vector of length 10 containing various user-specified 
c           inputs; see below. 
c
c  info   = integer vector of length 6 containing various outputs; 
c           see below. 
c
c  rwork  = real work vector with length as follows: 
c
c             solver    rwork length
c
c             GMRES     n*(kdmax+5)+kdmax*(kdmax+3), where kdmax is the 
c                       maximum Krylov subspace dimension, either the 
c                       default value of 20 or another value specified 
c                       by the user (see input(4) below). 
c
c             BiCGSTAB  11*n
c
c             TFQMR     14*n
c
c           The final f-value is contained in the first n components of 
c           rwork on return. 
c
c  rpar   = real parameter/work array passed to the f and jacv routines. 
c
c  ipar   = integer parameter/work array passed to the f and jacv routines. 
c               
c  iterm   = termination flag; values have the following meanings: 
c             -k => illegal value in input(k). 
c              0 => normal termination: ||F||.le.ftol or ||step||.le.stptol.
c              1 => nnimax nonlinear iterations reached without success. 
c              2 => failure to evaluate F.
c              3 => in nitjv, J*v failure. 
c              4 => in nitjv, P(inverse)*v failure. 
c              5 => in nitdrv, insufficient initial model norm reduction 
c                   for adequate progress. NOTE: This can occur for several 
c                   reasons; examine itrmks on return from the Krylov 
c                   solver for further information. (This will be printed out 
c                   if iplvl .ge. 3; see the discussion of optional 
c                   common blocks below). 
c              6 => in nitbt, failure to reach an acceptable step through 
c                   backtracking. 
c
c  dinpr  = name of user-supplied function for calculating vector inner
c           products.  This function must have the form
c
c              xdoty = dinpr( n, x, sx, y, sy )
c
c           where n is the length of the vectors, x and y are the
c           starting locations of the vectors, and sx (sy) is the stride
c           in memory between consecutive elements of x (y).  This is the
c           same signature as the BLAS routine ddot; if the Euclidean
c           inner product is desired the user can link to a local BLAS
c           library and provide the name ddot to nitsol.  dinpr must be
c           declared as an external function that returns a double
c           precision value in the calling program.
c
c  dnorm  = name of user-supplied function for calculating vector norms.
c           This function must have the form
c
c              xnorm = dnorm( n, x, sx )
c
c           where n is the length of the vector, x is the starting
c           location of the vector, and sx is the stride in memory
c           between consecutive elements of x.  This is the same
c           signature as the BLAS routine dnrm2; if the Euclidean
c           norm is desired the user can link to a local BLAS library
c           and provide the name dnrm2 to nitsol.  dnorm must be
c           declared as an external function that returns a double
c           precision value in the calling program.
c
c ------------------------------------------------------------------------
c 
c Further explanation of input: 
c
c This array allows the user to specify various options. It should be 
c declared an integer vector of length 11 in the calling program. To 
c specify an option, set the appropriate input component to the desired 
c value according to the specifications below. 
c
c USAGE NOTE: Setting a particular input component to zero gives the 
c default option for that component in all cases. 
c
c The first five input components are things that every user might wish 
c to modify; the remainder will usually be of interest only to more 
c experienced users. 
c
c Optional every-user input:
c
c    input(1) = nnimax = maximum number of nonlinear iterations (default 200).
c 
c    input(2) = ijacv = flag for determining the method of J*v evaluation:
c                 0 => finite-difference evaluation (default) 
c                 1 => analytic evaluation
c
c    input(3) = ikrysl = flag for determining the Krylov solver: 
c                 0 => GMRES (default)
c                 1 => BiCGSTAB
c                 2 => TFQMR
c
c               For brief descriptions of the solvers plus references, 
c               see the subroutines nitgm, nitstb, and nittfq. 
c
c    input(4) = kdmax = maximum Krylov subspace dimension when GMRES is used 
c               (default 20). 
c
c    input(5) = irpre = flag for right preconditioning: 
c                 0 => no right preconditioning
c                 1 => right preconditioning
c
c Optional experienced user input:
c
c    input(6) = iksmax = maximum allowable number of iterations per call 
c               to the Krylov solver routine (default 1000). 
c
c    input(7) = iresup = residual update flag when GMRES is used; on 
c               restarts, the residual is updated as follows: 
c                 0 => linear combination (default) 
c                 1 => direct evaluation
c               The first is cheap (one n-vector saxpy) but may lose 
c               accuracy with extreme residual reduction; the second 
c               retains accuracy better but costs one J*v product per 
c               restart. 
c
c    input(8) = ifdord = order of the finite-difference formula (sometimes) 
c               used when input(2) = ijacv = 0. When input(2) = ijacv = 0, 
c               this must be 0, 1, 2, or 4 on input; otherwise, it is 
c               irrelevant. With input(2) = ijacv = 0, the precise 
c               meaning is as follows: 
c
c               If GMRES is used, then ifdord matters only if input(7) = 
c               iresup = 1, in which case it determines the order of 
c               the finite-difference formula used in evaluating the 
c               initial residual at each GMRES restart (default 2); if 
c               ifdord = 0 on input, then it is set to 2 below. NOTE: This 
c               only affects initial residuals at restarts; first-order 
c               differences are always used within each GMRES cycle. Using 
c               higher-order differences at restarts only should give 
c               the same accuracy as if higher-order differences were 
c               used throughout; see K. Turner and H. F. Walker, "Efficient 
c               high accuracy solutions with GMRES(m)," SIAM J. Sci. 
c               Stat. Comput., 13 (1992), pp. 815--825. 
c               
c               If BiCGSTAB or TFQMR is used, then ifdord determines the 
c               order of the finite-difference formula used at each 
c               iteration (default 1); if ifdord = 0 on input, then it 
c               is set to 1 below. 
c
c    input(9) = ibtmax = maximum allowable number of backtracks (step 
c               reductions) per call to nitbt (default 10). 
c
c               USAGE NOTE: Backtracking can be turned off by setting 
c		ibtmax = -1. Other negative values of ibtmax are not 
c               valid. 
c
c    input(10) = ieta = flag determining the forcing term eta as follows: 
c                 0 => abs( ||fcur|| - ||fprev+Jprev*sprev|| )/||fprev||
c                      (default) 
c                 1 => (||fcur||/||fprev||)**2 
c                 2 => gamma*(||fcur||/||fprev||)**alpha 
c                      for user-supplied gamma in (0,1] and alpha in (1,2] 
c                 3 => fixed (constant) eta in (0,1), either 0.1 (default) 
c		       or specified by the user (see USAGE NOTE below) 
c               Here, fcur = current f, fprev = previous f, etc. The Krylov 
c               iterations are terminated when an iterate s satisfies 
c               an inexact Newton condition ||F + J*s|| .le. eta*||F||.
c
c               USAGE NOTE: If input(10) = ieta = 2, then alpha and gamma 
c               must be set in common block nitparam.h as described below. 
c		If input(10) = ieta = 3, then the desired constant eta may 
c		be similarly set in nitparam.h if a value other than the 
c		default of 0.1 is desired. 
c               
c               The first three expressions above are from S. C. Eisenstat 
c               and H. F. Walker, "Choosing the forcing terms in an inexact 
c               Newton method", SIAM J. Scientific Computing, 17 (1996), 
c               pp. 16--32. (They may be modified according to certain 
c               safeguards in subroutine nitdrv.) The first gives convergence 
c               that is q-superlinear and of r-order (1+sqrt(5))/2. The 
c               second gives convergence that is r-quadratic and of q-order 
c               p for every p in [1,2). The third gives convergence that is 
c               of q-order alpha when gamma < 1 and, when gamma = 1, of 
c               r-order alpha and q-order p for every p in [1,alpha). The 
c               fourth gives q-linear convergence with asymptotic rate 
c               constant eta in a certain norm; see R. S. Dembo, S. C. 
c		Eisenstat, and T. Steihaug, "Inexact Newton methods", 
c               SIAM J. Numer. Anal., 18 (1982), pp. 400-408. 
c
c               Of these four choices, the 1st is usually satisfactory, 
c               the 2nd or 3rd is sometimes preferred, and the 4th may be 
c               useful in some situations, e.g., it may be desirable to 
c               choose a fairly large fixed eta in (0,1), such as eta = .1, 
c               when numerical inaccuracy prevents the Krylov solver 
c               from obtaining much residual reduction. 
c
c ------------------------------------------------------------------------
c
c Further explanation of info: On output, the components of info are 
c as follows: 
c
c     info(1)   = nfe (number of function evaluations)
c     info(2)   = njve (number of J*v evaluations)
c     info(3)   = nrpre (number of P(inverse)*v evaluations)
c     info(4)   = nli (number of linear iterations)
c     info(5)   = nni (number of nonlinear iterations)
c     info(6)   = nbt (number of backtracks)
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
c              NOTE: If ipunit = 0 on input, then it is set to 6 below.
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
c        eta - forcing term. 
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
c  The values in this common block are checked once here in nitsol 
c  before the main solution driver is called.  If any parameter has
c  an invalid value, it is silently reset to the default value.
c
c ------------------------------------------------------------------------
c
c Subroutines required by this and all called routines: 
c
c     user supplied: f, jacv, dinpr, dnorm
c
c     nitsol routines: nitbd.f, nitbt, nitdrv, nitgm, nitjv, nitstb, nittfq
c
c     lapack routines: dlaic1, dlamch.f
c
c     blas routines: daxpy, dcopy, dscal, dswap 
c
c This subroutine called by: calling program 
c
c Subroutines called by this subroutine: nitdrv 
c
c ------------------------------------------------------------------------
c
c Remaining declarations:
c
      integer ibtmax, ieta, ifdord, ijacv, ikrysl, iksmax, iresup, 
     $     irpre, kdmax, lfcur, lfpls, lrwork, lstep, lxpls, 
     $     nbt, nfe, njve, nli, nni, nnimax, nrpre , dim_rwork

      include 'nitdflts.h'
 
      external nitbd

c ------------------------------------------------------------------------ 
c
c Begin executable code. 
c 
c ------------------------------------------------------------------------
c Check inputs and initialize parameters. 
c ------------------------------------------------------------------------
      if (ipunit .gt. 6) open( unit=ipunit, status='unknown' )
      if (ipunit .eq. 0) ipunit = 6
      if (input(1) .eq. 0) then
         nnimax = 200
      else
         if (input(1) .gt. 0) then 
            nnimax = input(1)
         else
            iterm = -1
            go to 900
         endif
      endif
      if (input(2) .eq. 0 .or. input(2) .eq. 1) then 
         ijacv = input(2) 
      else
         iterm = -2
         go to 900
      endif
      if (input(3) .ge. 0 .and. input(3) .le. 2) then 
         ikrysl = input(3)
      else 
         iterm = -3
         go to 900
      endif
      if (ikrysl .eq. 0) then 
         if (input(4) .eq. 0) then 
            kdmax = 20
         else
            if (input(4) .gt. 0) then 
               kdmax = input(4) 
            else
               iterm = -4
               go to 900
            endif
         endif
      endif
      if (input(5) .eq. 0 .or. input(5) .eq. 1) then 
         irpre = input(5) 
      else
         iterm = -5
         go to 900
      endif
      if (input(6) .eq. 0) then 
         iksmax = 1000
      else
         if (input(6) .gt. 0) then 
            iksmax = input(6)
         else
            iterm = -6
            go to 900
         endif
      endif
      if (ikrysl .eq. 0) then 
         if (input(7) .eq. 0 .or. input(7) .eq. 1) then 
            iresup = input(7)
         else
            iterm = -7
            go to 900
         endif
      endif
      if (ijacv .eq. 0) then 
         if (input(8) .eq. 0) then 
            if (ikrysl .eq. 0) then 
               ifdord = 2
            else
               ifdord = 1
            endif
         else
            if (input(8) .eq. 1 .or. input(8) .eq. 2 .or. 
     $           input(8) .eq. 4) then 
               ifdord = input(8)
            else
               iterm = -8
               go to 900
            endif
         endif
      endif
      if (input(9) .eq. 0) then 
         ibtmax = 10
      else
         if (input(9) .gt. 0 .or. input(9) .eq. -1) then 
            ibtmax = input(9)
         else
            iterm = -9 
            go to 900
         endif
      endif
      if (input(10) .ge. 0 .and. input(10) .le. 3) then 
         ieta = input(10)
      else
         iterm = -10
         go to 900
      endif
c ------------------------------------------------------------------------
c  Check possible invalid value for printout level.  In
c  case the value is invalid the default is restored.
c ------------------------------------------------------------------------
      if ( iplvl .lt. 0 .or. iplvl .gt. 4 ) iplvl = DFLT_PRLVL
c ------------------------------------------------------------------------
c  Check possible invalid values for various parameters.  In
c  case the values are invalid the defaults are restored.
c ------------------------------------------------------------------------
      if ( choice1_exp .le. 1.0d0 .or. choice1_exp .gt. 2.0d0 ) 
     &  choice1_exp = DFLT_CHOICE1_EXP
      if ( choice2_exp .le. 1.0d0 .or. choice2_exp .gt. 2.0d0 ) 
     &  choice2_exp = DFLT_CHOICE2_EXP
      if ( choice2_coef .lt. 0.0d0 .or. choice2_coef .gt. 1.0d0 ) 
     &  choice2_coef = DFLT_CHOICE2_COEF
      if ( etamax .le. 0.0d0 ) etamax = DFLT_ETA_MAX
      if ( thmin .lt. 0.0d0 .or. thmin .gt. thmax )
     &  thmin = DFLT_THMIN
      if ( thmax .gt. 1.0d0 .or. thmax .lt. thmin )
     &  thmax = DFLT_THMAX

      dim_rwork=max(n*(kdmax+5)+kdmax*(kdmax+3), 14*n)
      allocate(rwork(dim_rwork)); rwork=0
c ------------------------------------------------------------------------
c  Initialize some pointers into the rwork array.
c ------------------------------------------------------------------------
      lfcur = 1
      lxpls = lfcur + n
      lfpls = lxpls + n
      lstep = lfpls + n
      lrwork = lstep + n
c ------------------------------------------------------------------------
c Call nitdrv. 
c ------------------------------------------------------------------------
      call nitdrv(n, x, rwork(lfcur), rwork(lxpls), rwork(lfpls), 
     $        rwork(lstep), f, jacv, rpar, ipar, ftol, stptol, nnimax, 
     $        ijacv, ikrysl, kdmax, irpre, iksmax, iresup, ifdord, 
     $        ibtmax, ieta, iterm, nfe, njve, nrpre, nli, nni, nbt, 
     $        rwork(lrwork), ddot, dnrm2)
c ------------------------------------------------------------------------
c Set output for return. 
c ------------------------------------------------------------------------
      info(1) = nfe 
      info(2) = njve 
      info(3) = nrpre 
      info(4) = nli 
      info(5) = nni 
      info(6) = nbt 
      deallocate(rwork)
c ------------------------------------------------------------------------
c All returns made here.
c ------------------------------------------------------------------------
 900  continue
      if ( ipunit .gt. 6 ) close( unit=ipunit )
      return
      end
c -------------------- end of subroutine nitsol  --------------------------
