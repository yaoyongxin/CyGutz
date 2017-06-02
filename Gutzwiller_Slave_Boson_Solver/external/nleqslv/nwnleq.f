
      subroutine nwnleq(x0,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,cndtol,method,global,
     *                  xscalm,stepmx,delta,sigma,rjac,ldr,
     *                  rwork,lrwork,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,fjac,fvec,outopt,xp,
     *                  fp,gp,njcnt,nfcnt,iter,termcd)

      integer n,jacflg(*),maxit,njcnt,nfcnt,iter,termcd,method
      integer global,xscalm,ldr,lrwork,qrwsiz
      integer outopt(*)
      double precision  xtol,ftol,btol,cndtol,stepmx,delta,sigma
      double precision  xp(*),fp(*),gp(*),x0(*)
      double precision  rjac(ldr,*),rwork(*),rcdwrk(*),qrwork(*)
      double precision  scalex(*)
      integer           icdwrk(*)
      external fjac,fvec

c-------------------------------------------------------------------------
c
c     Solves systems of nonlinear equations using the Newton / Broyden
c     method with a global strategy either linesearch or double dogleg
c
c     In       x0      Real(*)         starting vector for x
c     In       n       Integer         dimension of problem
c     Inout    scalex  Real(*)         scaling factors x()
c     Inout    maxit   Integer         maximum number iterations
c     Inout    jacflg  Integer(*)      jacobian flag array
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting step allowed when
c                                                   singular or illconditioned
c     Inout    xtol    Real            x tolerance
c     Inout    ftol    Real            f tolerance
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    cndtol  Real            tolerance of test for ill conditioning
c     Inout    method  Integer         method to use
c                                        0 Newton
c                                        1 Broyden
c     In       global  Integer         global strategy to use
c                                        1 cubic linesearch
c                                        2 quadratic linesearch
c                                        3 geometric linesearch
c                                        4 double dogleg
c                                        5 powell dogleg
c                                        6 hookstep (More-Hebden Levenberg-Marquardt)
c     In       xscalm  Integer         scaling method
c                                        0 scale fixed and supplied by user
c                                        1 for scale from jac. columns a la Minpack
c     Inout    stepmx  Real            maximum stepsize
c     Inout    delta   Real            trust region radius
c                                        > 0.0 or special value for initial value
c                                        -1.0  ==> use min(Cauchy length, stepmx)
c                                        -2.0  ==> use min(Newton length, stepmx)
c     Inout    sigma   Real            reduction factor geometric linesearch
c     Inout    rjac    Real(ldr,*)     workspace jacobian
c                                         2*n*n for Broyden and n*n for Newton
c     In       ldr     Integer         leading dimension rjac
c     Out      rwork   Real(*)         real workspace (9n)
c     In       lrwork  Integer         size real workspace
c     In       rcdwrk  Real(*)         workspace for Dtrcon (3n)
c     In       icdwrk  Integer(*)      workspace for Dtrcon (n)
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c     In       fjac    Name            optional name of routine to calculate
c                                      user supplied jacobian
c     In       fvec    Name            name of routine to calculate f(x)
c     In       outopt  Integer(*)      output options
c                                       outopt(1)
c                                         0 no output
c                                         1 output an iteration report
c                                       outopt(2)
c                                         0 do not check any user supplied jacobian
c                                         1 check user supplied jacobian if supplied
c     Out      xp      Real(*)         final values for x()
c     Out      fp      Real(*)         final values for f(x)
c     Out      gp      Real(*)         gradient of f() at xp()
c     Out      njcnt   Integer         number of jacobian evaluations
c     Out      nfcnt   Integer         number of function evaluations
c     Out      iter    Integer         number of (outer) iterations
c     Out      termcd  Integer         termination code
c                                       > 0 process terminated
c                                             1  function criterion near zero
c                                             2  no better point found
c                                             3  x-values within tolerance
c                                             4  iteration limit exceeded
c                                             5  singular/ill-conditioned jacobian
c                                             6  totally singular jacobian
c                                                (when allowSingular=TRUE)
c
c                                       < 0 invalid input parameters
c                                            -1  n not positive
c                                            -2  insufficient workspace rwork
c                                            -3  cannot check user supplied jacobian (not supplied)
c
c    The subroutine fvec must be declared as
c
c!        subroutine fvec(x,f,n,flag)
c         double precision x(*), f(*)
c         integer  n, flag
c
c         x() are the x values for which to calculate the function values f(*)
c         The dimension of these vectors is n
c         The flag argument is set to
c            0  for calculation of function values
c           >0  indicating that jacobian column <flag> is being computed
c               so that fvec can abort.
c
c    The subroutine fjac must be declared as
c
c!        subroutine mkjac(rjac,ldr,x,n)
c         integer ldr
c         double precision rjac(ldr,*), x(*)
c         integer  n
c
c         The routine calculates the jacobian in point x(*) of the
c         function. If any illegal values are encountered during
c         calculation of the jacobian it is the responsibility of
c         the routine to quit.

c-------------------------------------------------------------------------

      double precision epsm

c     check input parameters

      call nwpchk(n,lrwork,xtol,ftol,btol,cndtol,maxit,
     *            jacflg,method,global,stepmx,delta,sigma,
     *            epsm,outopt,scalex,xscalm,termcd)
      if(termcd .lt. 0) then
         return
      endif

c     first argument of nwsolv/brsolv is leading dimension of rjac in those routines
c     should be at least n

      if( method .eq. 0 ) then

         call nwsolv(ldr,x0,n,scalex,maxit,jacflg,
     *               xtol,ftol,btol,cndtol,global,xscalm,
     *               stepmx,delta,sigma,
     *               rjac,
     *               rwork(1    ),rwork(1+  n),
     *               rwork(1+2*n),rwork(1+3*n),
     *               rwork(1+4*n),rwork(1+5*n),
     *               rwork(1+6*n),rwork(1+7*n),
     *               rwork(1+8*n),rcdwrk,icdwrk,qrwork,qrwsiz,epsm,
     *               fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)

      elseif( method .eq. 1 ) then

         call brsolv(ldr,x0,n,scalex,maxit,jacflg,
     *               xtol,ftol,btol,cndtol,global,xscalm,
     *               stepmx,delta,sigma,
     *               rjac,rjac(1,n+1),
     *               rwork(1    ),rwork(1+  n),
     *               rwork(1+2*n),rwork(1+3*n),
     *               rwork(1+4*n),rwork(1+5*n),
     *               rwork(1+6*n),rwork(1+7*n),
     *               rwork(1+8*n),rcdwrk,icdwrk,
     *               qrwork,qrwsiz,epsm,
     *               fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)

      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nwpchk(n,lrwk,
     *                  xtol,ftol,btol,cndtol,maxit,jacflg,method,
     *                  global,stepmx,delta,sigma,epsm,outopt,
     *                  scalex,xscalm,termcd)

      integer n,lrwk,jacflg(*)
      integer method,global,maxit,xscalm,termcd
      integer outopt(*)
      double precision  xtol,ftol,btol,cndtol,stepmx,delta,sigma,epsm
      double precision  scalex(*)

c-------------------------------------------------------------------------
c
c     Check input arguments for consistency and modify if needed/harmless
c
c     Arguments
c
c     In       n       Integer         dimension of problem
c     In       lrwk    Integer         size real workspace
c     Inout    xtol    Real            x tolerance
c     Inout    ftol    Real            f tolerance
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    cndtol  Real            tolerance of test for ill conditioning
c     Inout    maxit   Integer         maximum number iterations
c     Inout    jacflg  Integer(*)      jacobian flag
c     Inout    method  Integer         method to use (Newton/Broyden)
c     Inout    global  Integer         global strategy to use
c     Inout    stepmx  Real            maximum stepsize
c     Inout    delta     Real            trust region radius
c     Inout    sigma   Real            reduction factor geometric linesearch
c     Out      epsm                    machine precision
c     Inout    scalex  Real(*)         scaling factors x()
c     Inout    xscalm  integer         0 for fixed scaling, 1 for automatic scaling
c     Out      termcd  Integer         termination code (<0 on errors)
c
c-------------------------------------------------------------------------

      integer i,len
      double precision epsmch, dblhuge, Rhuge

      double precision Rzero, Rone, Rtwo, Rthree
      parameter(Rzero=0.0d0, Rone=1.0d0, Rtwo=2.0d0, Rthree=3.0d0)

      double precision Rhalf
      parameter(Rhalf = 0.5d0)

c     check that parameters only take on acceptable values
c     if not, set them to default values

c     initialize termcd to all ok

      termcd = 0

c     compute machine precision

      epsm = epsmch()

c     get largest double precision number
      Rhuge = dblhuge()

c     check dimensions of the problem

      if(n .le. 0) then
         termcd = -1
         return
      endif

c     check dimensions of workspace arrays

      len = 9*n
c      +2*n*n
      if(lrwk .lt. len) then
         termcd = -2
         return
      endif

c     check jacflg, method, and global

      if(jacflg(1) .gt. 3 .or. jacflg(1) .lt. 0) jacflg(1) = 0

      if(method .lt. 0 .or. method .gt. 1) method = 0

      if(global .lt. 0 .or. global .gt. 6) global = 4

c     set outopt to correct values

      if(outopt(1) .ne. 0 ) then
         outopt(1) = 1
      endif

      if(outopt(2) .ne. 0 ) then
         outopt(2) = 1
      endif

c     check scaling scale matrices

      if(xscalm .eq. 0) then
         do i = 1,n
            if(scalex(i) .lt. Rzero) scalex(i) = -scalex(i)
            if(scalex(i) .eq. Rzero) scalex(i) = Rone
         enddo
      else
         xscalm = 1
         do i = 1,n
            scalex(i) = Rone
         enddo
      endif
c     check step and function tolerances

      if(xtol .lt. Rzero) then
         xtol = epsm**(Rtwo/Rthree)
      endif

      if(ftol .lt. Rzero) then
         ftol = epsm**(Rtwo/Rthree)
      endif

      if( btol .lt. xtol ) btol = xtol

      cndtol = max(cndtol, epsm)

c     check reduction in geometric linesearch

      if( sigma .le. Rzero .or. sigma .ge. Rone ) then
         sigma = Rhalf
      endif

c     check iteration limit

      if(maxit .le. 0) then
         maxit = 150
      endif

c     set stepmx

      if(stepmx .le. Rzero) stepmx = Rhuge

c     check delta
      if(delta .le. Rzero) then
         if( delta .ne. -Rtwo ) then
            delta = -Rone
         endif
      elseif(delta .gt. stepmx) then
         delta = stepmx
      endif

      return
      end
