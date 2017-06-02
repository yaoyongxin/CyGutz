
      subroutine nwddlg(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol,
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
c     Find a next iterate xp by the double dogleg method
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
      double precision  dnlen,ssdlen,alpha,beta,lambda,fpred
      double precision  sqalpha,eta,gamma,fpnsav,oarg(7)
      double precision  dnrm2,ddot
      logical nwtstep
      integer dtype

      integer idamax

      double precision Rone, Rtwo, Rten, Rhalf, Rp2, Rp8
      parameter(Rhalf=0.5d0)
      parameter(Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)
      parameter(Rp2 = Rtwo/Rten, Rp8 = Rone - Rp2)
      double precision Rzero
      parameter(Rzero=0.0d0)

c     length newton direction

      dnlen = dnrm2(n, dn, 1)

c     steepest descent direction and length

      sqalpha = dnrm2(n,g,1)
      alpha   = sqalpha**2

      call dcopy(n, g, 1, d, 1)
      call dtrmv('U','N','N',n,rjac,ldr,d,1)
      beta = dnrm2(n,d,1)**2

      call dcopy(n, g, 1, ssd, 1)
      call dscal(n, -(alpha/beta), ssd, 1)

      ssdlen = alpha*sqalpha/beta

c     set trust radius to ssdlen or dnlen if required

      if( delta .eq. -Rone ) then
         delta = min(ssdlen, stepmx)
      elseif( delta .eq. -Rtwo ) then
         delta = min(dnlen, stepmx)
      endif

c     calculate double dogleg parameter

      gamma = alpha*alpha/(-beta*ddot(n,g,1,dn,1))
c      call dgdbg(gamma, alpha*alpha, -beta*ddot(n,g,1,dn,1))
c     precautionary (just in case)
      eta = max(Rzero, min(Rone,Rp2 + Rp8*gamma))

      retcd = 4
      gcnt  = 0

      do while( retcd .gt. 1 )
c        find new step by double dogleg algorithm

         call ddlgstp(n,dn,dnlen,delta,v,
     *               ssd,ssdlen,eta,d,dtype,lambda)
         nwtstep = dtype .eq. 4
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
            oarg(1) = lambda
            oarg(3) = delta
            oarg(4) = eta
            oarg(5) = fpnorm
            oarg(6) = abs(fp(idamax(n,fp,1)))
            call nwdgot(iter,dtype,retcd,oarg)
         endif

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine ddlgstp(n,dn,dnlen,delta,v,
     *                  ssd,ssdlen,eta,d,dtype,lambda)
      integer n
      double precision  dn(*), ssd(*), v(*), d(*)
      double precision  dnlen, delta, ssdlen, eta, lambda
      integer dtype

c-------------------------------------------------------------------------
c
c     Find a new step by the double dogleg algorithm
c     Internal routine for nwddlg
c
c     Arguments
c
c     In       n       Integer         size of problem
c     In       dn      Real(*)         current newton step
c     Out      dnlen   Real            length dn()
c     In       delta   Real            current trust region radius
c     Out      v       Real(*)         (internal) eta * dn() - ssd()
c     In       ssd     Real(*)         (internal) steepest descent direction
c     In       ssdlen  Real            (internal) length ssd
c     In       eta     Real            (internal) double dogleg parameter
c     Out      d       Real(*)         new step for x()
c     Out      dtype   Integer         steptype
c                                       1 steepest descent
c                                       2 combination of dn and ssd
c                                       3 partial newton step
c                                       4 full newton direction
c     Out      lambda  Real            weight of eta*dn() in d()
c                                      closer to 1 ==> more of eta*dn()
c
c-----------------------------------------------------------------------

      integer i
      double precision vssd, vlen
      double precision dnrm2, ddot

      if(dnlen .le. delta) then

c        Newton step smaller than trust radius ==> take it

         call dcopy(n, dn, 1, d, 1)
         delta = dnlen
         dtype = 4

      elseif(eta*dnlen .le. delta) then

c        take partial step in newton direction

         call dcopy(n, dn, 1, d, 1)
         call dscal(n, delta / dnlen, d, 1)
         dtype = 3

      elseif(ssdlen .ge. delta) then

c        take step in steepest descent direction

         call dcopy(n, ssd, 1, d, 1)
         call dscal(n, delta / ssdlen, d, 1)
         dtype = 1

      else

c        calculate convex combination of ssd and eta*dn with length delta

         do i=1,n
            v(i) = eta*dn(i) - ssd(i)
         enddo

         vssd = ddot(n,v,1,ssd,1)
         vlen = dnrm2(n,v,1)**2

         lambda =(-vssd+sqrt(vssd**2-vlen*(ssdlen**2-delta**2)))/vlen
         call dcopy(n, ssd, 1, d, 1)
         call daxpy(n, lambda, v, 1, d, 1)
         dtype = 2

      endif

      return
      end
