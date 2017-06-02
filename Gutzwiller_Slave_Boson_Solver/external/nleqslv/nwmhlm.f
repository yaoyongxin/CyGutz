
      subroutine nwmhlm(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol,
     *                  delta,qtf,scalex,fvec,d,xprev,
     *                  ssd,v,wa,fprev,xp,fp,fpnorm,retcd,gcnt,
     *                  priter,iter)

      integer ldr, n, retcd, gcnt, priter, iter
      double precision  fcnorm, stepmx, xtol, fpnorm, delta
      double precision  rjac(ldr,*), dn(*), g(*), xc(*), qtf(*)
      double precision  scalex(*), d(*)
      double precision  xprev(*), xp(*), fp(*)
      double precision  ssd(*), v(*), wa(*), fprev(*)
      external fvec

c-------------------------------------------------------------------------
c
c     Find a next iterate xp by the More-Hebden-Levenberg-Marquardt method
c
c     Arguments
c
c     In       n       Integer         size of problem: dimension x, f
c     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
c     In       ldr     Integer         leading dimension of Rjac
c     Inout    dn      Real(*)         newton direction
c     Inout    g       Real(*)         gradient at current point
c                                      trans(jac)*f()
c     In       xc      Real(*)         current iterate
c     In       fcnorm  Real            .5*||f(xc)||**2
c     In       stepmx  Real            maximum stepsize
c     In       xtol    Real            x-tolerance (stepsize)
c     Inout    delta   Real            on input: initial trust region radius
c                                                if -1 then set to something
c                                                reasonable
c                                      on output: final value
c                                      ! Do not modify between calls while
c                                        still iterating
c     In       qtf     Real(*)         trans(Q)*f(xc)
c     In       scalex  Real(*)         scaling factors for x()
c     In       fvec    Name            name of subroutine to evaluate f(x)
c                                      ! must be declared external in caller
c     Wk       d       Real(*)         work vector
c     Wk       xprev   Real(*)         work vector
c     Wk       ssd     Real(*)         work vector
c     Wk       v       Real(*)         work vector
c     Wk       wa      Real(*)         work vector
c     Wk       fprev   Real(*)         work vector
c     Out      xp      Real(*)         new x()
c     Out      fp      Real(*)         new f(xp)
c     Out      fpnorm  Real            new .5*||f(xp)||**2
c     Out      retcd   Integer         return code
c                                       0  new satisfactory x() found
c                                       1  no  satisfactory x() found
c     Out      gcnt    Integer         number of steps taken
c     In       priter  Integer         print flag
c                                       -1 no intermediate printing
c                                       >0 yes for print of intermediate results
c     In       iter    Integer         current iteration (only used for above)
c
c     All vectors at least size n
c
c-------------------------------------------------------------------------

      integer i
      double precision  dnlen,glen,ssdlen,alpha,beta,mu,fpred
      double precision  fpnsav,oarg(6)
      double precision  dnrm2
      logical nwtstep
      integer dtype

      integer idamax

      double precision Rone, Rtwo, Rhalf
      parameter(Rhalf=0.5d0)
      parameter(Rone=1.0d0, Rtwo=2.0d0)

c     length newton direction

      dnlen = dnrm2(n, dn, 1)

c     gradient length and steepest descent direction and length

      glen  = dnrm2(n,g,1)
      alpha = glen**2

      call dcopy(n, g, 1, d, 1)
      call dtrmv('U','N','N',n,rjac,ldr,d,1)
      beta = dnrm2(n,d,1)**2

      call dcopy(n, g, 1, ssd, 1)
      call dscal(n, -(alpha/beta), ssd, 1)

      ssdlen = alpha*glen/beta

c     set trust radius to ssdlen or dnlen if required

      if( delta .eq. -Rone ) then
         delta = min(ssdlen, stepmx)
      elseif( delta .eq. -Rtwo ) then
         delta = min(dnlen, stepmx)
      endif

      retcd = 4
      gcnt  = 0

      do while( retcd .gt. 1 )
c        find new step by More Hebden LM  algorithm
c        reuse ssd as sdiag

         call nwmhstep(Rjac,ldr,n,ssd,qtf,dn,dnlen,glen,delta,mu,
     *                 d, v, dtype)
         nwtstep = dtype .eq. 2
c        compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

         call dcopy(n,d,1,wa,1)
         call dtrmv('U','N','N',n,rjac,ldr,wa,1)
         call daxpy(n, Rone, qtf,1,wa,1)
         fpred = Rhalf * dnrm2(n,wa,1)**2

c        evaluate function at xp = xc + d

         do i=1,n
            xp(i) = xc(i) + d(i)
         enddo

         call nwfvec(xp,n,scalex,fvec,fp,fpnorm,wa)
         gcnt = gcnt + 1

c        check whether the global step is acceptable

         oarg(2) = delta
         call nwtrup(n,fcnorm,g,d,nwtstep,stepmx,xtol,delta,
     *               fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm)

         if( priter .gt. 0 ) then
            oarg(1) = mu
            oarg(3) = delta
            oarg(4) = dnrm2(n, d, 1)
            oarg(5) = fpnorm
            oarg(6) = abs(fp(idamax(n,fp,1)))
            call nwmhot(iter,dtype,retcd,oarg)
         endif

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine nwmhstep(R,ldr,n,sdiag,qtf,dn,dnlen,glen,delta,mu,
     *                  d, work, dtype)
      integer ldr, n
      double precision R(ldr,*)
      double precision sdiag(*), qtf(*), dn(*), d(*), work(*)
      double precision  dnlen, glen, delta, mu
      integer dtype

c-------------------------------------------------------------------------
c
c     Find a new step by the More Hebden Levemberg Marquardt algorithm
c     Internal routine for nwmhlm
c
c     Arguments
c
c     In       R       Real(ldr,*)     R of QR-factored jacobian
c     In       ldr     Integer         leading dimension of R
c     In       n       Integer         size of problem
c     Out      sdiag   Real(*)         diagonal of LM lower triangular modified R
c     In       qtf     Real(*)         trans(Q)*f(xc)
c     In       dn      Real(*)         current newton step
c     Out      dnlen   Real            length dn()
c     In       glen    Real            length gradient
c     In       delta   Real            current trust region radius
c     Inout    mu      Real            Levenberg-Marquardt parameter
c     Out      d       Real(*)         new step for x()
c     Work     work    Real(*)         work vector for limhpar
c     Out      dtype   Integer         steptype
c                                       1 LM step
c                                       2 full newton direction
c
c-----------------------------------------------------------------------

      double precision Rone
      parameter(Rone=1.0D0)

      if(dnlen .le. delta) then

c        Newton step smaller than trust radius ==> take it

         call dcopy(n, dn, 1, d, 1)
         delta = dnlen
         dtype = 2

      else

c        calculate LM step
         call limhpar(R, ldr, n, sdiag, qtf, dn, dnlen, glen, delta,
     *                mu, d, work)
c        change sign of step d (limhpar solves for trans(R)*R+mu*I)=qtf instead of -qtf)
         call dscal(n,-Rone,d,1)
         dtype = 1
      endif

      return
      end
