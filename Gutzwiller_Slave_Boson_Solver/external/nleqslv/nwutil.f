
      subroutine nwtcvg(xplus,fplus,xc,xtol,retcd,ftol,iter,
     *                  maxit,n,ierr,termcd)

      integer n,iter,maxit,ierr,termcd,retcd
      double precision  xtol,ftol
      double precision  xplus(*),fplus(*),xc(*)

c-------------------------------------------------------------------------
c
c     Decide whether to terminate the nonlinear algorithm
c
c     Arguments
c
c     In       xplus   Real(*)         new x values
c     In       fplus   Real(*)         new f values
c     In       xc      Real(*)         current x values
c     In       xtol    Real            stepsize tolerance
c     In       retcd   Integer         return code from global search routines
c     In       ftol    Real            function tolerance
c     In       iter    Integer         iteration number
c     In       maxit   Integer         maximum number of iterations allowed
c     In       n       Integer         size of x and f
c     In       ierr    Integer         return code of cndjac (condition estimation)
c
c     Out      termcd  Integer         termination code
c                                        0 no termination criterion satisfied
c                                          ==> continue iterating
c                                        1 norm of scaled function value is
c                                          less than ftol
c                                        2 scaled distance between last
c                                          two steps less than xtol
c                                        3 unsuccessful global strategy
c                                          ==> cannot find a better point
c                                        4 iteration limit exceeded
c                                        5 Jacobian too ill-conditioned
c                                        6 Jacobian singular
c
c-------------------------------------------------------------------------

      double precision  fmax,rsx, nuxnrm
      integer idamax

c     check whether function values are within tolerance

      termcd = 0

      if( ierr .ne. 0 ) then
         termcd = 4 + ierr
         return
      endif

      fmax = abs(fplus(idamax(n,fplus,1)))
      if( fmax .le. ftol) then
         termcd = 1
         return
      endif

c     initial check at start so there is no xplus
c     so only a check of function values is useful
      if(iter .eq. 0) return

      if(retcd .eq. 1) then
         termcd = 3
         return
      endif

c     check whether relative step length is within tolerance
c     Dennis Schnabel Algorithm A7.2.3

      rsx = nuxnrm(n, xplus, xc)
      if(rsx .le. xtol) then
        termcd = 2
        return
      endif

c     check iteration limit

      if(iter .ge. maxit) then
         termcd = 4
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,iter)
      double precision xc(*),fc(*),fcnorm,xp(*),fp(*),fpnorm
      integer n, gcnt, priter, iter

c-------------------------------------------------------------------------
c
c     calling routine got an error in decomposition/update of Jacobian/Broyden
c     jacobian an singular or too ill-conditioned
c     prepare return arguments
c
c     Arguments
c
c     In       n       Integer         size of x
c     In       xc      Real(*)         current (starting) x values
c     In       fc      Real(*)         function values f(xc)
c     In       fcnorm  Real            norm fc
c     Out      xp      Real(*)         final x values
c     Out      fp      Real(*)         function values f(xp)
c     Out      fpnorm  Real            final norm fp
c     Out      gcnt    Integer         # of backtracking steps (here set to 0)
c     In       priter  Integer         flag for type of output
c     In       iter    Integer         iteration counter
c
c-------------------------------------------------------------------------

      call dcopy(n,xc,1,xp,1)
      call dcopy(n,fc,1,fp,1)
      fpnorm = fcnorm
      gcnt   = 0
      if( priter .gt. 0 ) then
         call nwjerr(iter)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine chkjac1(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd)

      integer lda,n,termcd
      double precision  A(lda,*),xc(*),fc(*)
      double precision  epsm,scalex(*)
      double precision  fz(*),wa(*),xw(*)
      external fvec

c-------------------------------------------------------------------------
c
c     Check the user supplied jacobian against its finite difference approximation
c
c     Arguments
c
c     In       A       Real(lda,*)     user supplied jacobian
c     In       lda     Integer         leading dimension of ajanal
c     In       xc      Real(*)         vector of x values
c     In       fc      Real(*)         function values f(xc)
c     In       n       Integer         size of x
c     In       epsm    Real            machine precision
c     In       scalex  Real(*)         scaling vector for x()
c     Wk       fz      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c     Wk       xw      Real(*)         workspace
c     In       fvec    Name            name of routine to evaluate f(x)
c     Out      termcd  Integer         return code
c                                        0  user supplied jacobian ok
c                                      -10  user supplied jacobian NOT ok
c
c-------------------------------------------------------------------------

      integer i,j,errcnt
      double precision  ndigit,p,h,xcj,dinf
      double precision  tol
      double precision  rnudif
      integer idamax

      integer MAXERR
      parameter(MAXERR=10)

      double precision Rquart, Rten
      parameter(Rquart=0.25d0, Rten=10.0d0)

      termcd = 0

c     compute the finite difference jacobian and check it against
c     the analytic one

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))
      tol    = epsm**Rquart

      errcnt = 0
      call dcopy(n,xc,1,xw,1)
      call vunsc(n,xw,scalex)

      do j=1,n
         h = p + p * abs(xw(j))
         xcj   = xw(j)
         xw(j) = xcj + h

c        avoid (small) rounding errors
c        h = xc(j) - xcj but not here to avoid clever optimizers

         h = rnudif(xw(j), xcj)

         call fvec(xw,fz,n,j)
         xw(j) = xcj

         do i=1,n
            wa(i) = (fz(i)-fc(i))/h
         enddo

         dinf = abs(wa(idamax(n,wa,1)))

         do i=1,n
            if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
               errcnt = errcnt + 1
               if( errcnt .gt. MAXERR ) then
                  termcd = -10
                  return
               endif
               call nwckot(i,j,A(i,j),wa(i))
            endif
         enddo
      enddo

c      call vscal(n,xc,scalex)

      if( errcnt .gt. 0 ) then
         termcd = -10
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine chkjac2(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd,
     *                   dsub,dsuper)

      integer lda,n,termcd,dsub,dsuper
      double precision  A(lda,*),xc(*),fc(*)
      double precision  epsm,scalex(*)
      double precision  fz(*),wa(*),xw(*)
      external fvec

c-------------------------------------------------------------------------
c
c     Check the user supplied jacobian against its finite difference approximation
c
c     Arguments
c
c     In       A       Real(lda,*)     user supplied jacobian
c     In       lda     Integer         leading dimension of ajanal
c     In       xc      Real(*)         vector of x values
c     In       fc      Real(*)         function values f(xc)
c     In       n       Integer         size of x
c     In       epsm    Real            machine precision
c     In       scalex  Real(*)         scaling vector for x()
c     Wk       fz      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c     Wk       xw      Real(*)         workspace
c     In       fvec    Name            name of routine to evaluate f(x)
c     Out      termcd  Integer         return code
c                                        0  user supplied jacobian ok
c                                      -10  user supplied jacobian NOT ok
c
c-------------------------------------------------------------------------

      integer i,j,k,dsum,errcnt
      double precision  ndigit,p,h,dinf
      double precision  tol
      double precision w(n),xstep(n)

      integer MAXERR
      parameter(MAXERR=10)

      double precision Rquart, Rten, Rzero
      parameter(Rquart=0.25d0, Rten=10.0d0, Rzero=0.0d0)

      dsum = dsub + dsuper + 1

      termcd = 0

c     compute the finite difference jacobian and check it against
c     the user supplied one

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))
      tol    = epsm**Rquart

      errcnt = 0
      call dcopy(n,xc,1,xw,1)
      call vunsc(n,xw,scalex)

      do j=1,n
          xstep(j) = p + p * abs(xw(j))
          w(j) = xw(j)
      enddo

      do k=1,dsum
         do j=k,n,dsum
            xw(j) = xw(j) + xstep(j)
         enddo

c        for non finite values error message will be wrong
         call fvec(xw,fz,n,n+k)

         do j=k,n,dsum
             h = xstep(j)
             xw(j) = w(j)
             dinf = Rzero
             do i=max(j-dsuper,1),min(j+dsub,n)
                wa(i) = (fz(i)-fc(i)) / h
                if(abs(wa(i)).gt.dinf) dinf = abs(wa(i))
             enddo

             do i=max(j-dsuper,1),min(j+dsub,n)
                if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
                   errcnt = errcnt + 1
                   if( errcnt .gt. MAXERR ) then
                      termcd = -10
                      return
                   endif
                   call nwckot(i,j,A(i,j),wa(i))
                endif
             enddo
         enddo
      enddo

c      call vscal(n,xc,scalex)

      if( errcnt .gt. 0 ) then
         termcd = -10
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine chkjac(A,lda,xc,fc,n,epsm,jacflg,scalex,fz,wa,xw,fvec,
     *                  termcd)

      integer lda,n,termcd,jacflg(*)
      double precision  A(lda,*),xc(*),fc(*)
      double precision  epsm,scalex(*)
      double precision  fz(*),wa(*),xw(*)
      external fvec

c-------------------------------------------------------------------------
c
c     Check the user supplied jacobian against its finite difference approximation
c
c     Arguments
c
c     In       A       Real(lda,*)     user supplied jacobian
c     In       lda     Integer         leading dimension of ajanal
c     In       xc      Real(*)         vector of x values
c     In       fc      Real(*)         function values f(xc)
c     In       n       Integer         size of x
c     In       epsm    Real            machine precision
c     In       jacflg  Integer(*)      indicates how to compute jacobian
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting jacobian allowed when
c                                                   singular or illconditioned
c     In       scalex  Real(*)         scaling vector for x()
c     Wk       fz      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c     Wk       xw      Real(*)         workspace
c     In       fvec    Name            name of routine to evaluate f(x)
c     Out      termcd  Integer         return code
c                                        0  user supplied jacobian ok
c                                      -10  user supplied jacobian NOT ok
c
c-------------------------------------------------------------------------

      if(jacflg(1) .eq. 3) then
c        user supplied and banded
         call chkjac2(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd,
     *                jacflg(2),jacflg(3))
      else
         call chkjac1(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine fdjac0(xc,fc,n,epsm,fvec,fz,rjac,ldr)

      integer ldr,n
      double precision  epsm
      double precision  rjac(ldr,*),fz(*),xc(*),fc(*)
      external fvec

c-------------------------------------------------------------------------
c
c     Compute the finite difference jacobian at the current point xc
c
c     Arguments
c
c     In       xc      Real(*)         current point
c     In       fc      Real(*)         function values at current point
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       fvec    Name            name of routine to evaluate f(x)
c     Wk       fz      Real(*)         workspace
c     Out      rjac    Real(ldr,*)     jacobian matrix at x
c                                        entry [i,j] is derivative of
c                                        f(i) wrt to x(j)
c     In       ldr     Integer         leading dimension of rjac
c
c-------------------------------------------------------------------------

      integer i,j
      double precision  ndigit,p,h,xcj
      double precision  rnudif

      double precision Rten
      parameter(Rten=10d0)

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))

      do j=1,n
         h = p + p * abs(xc(j))

c        or as alternative h  = p * max(Rone, abs(xc(j)))

         xcj   = xc(j)
         xc(j) = xcj + h

c        avoid (small) rounding errors
c        h = xc(j) - xcj  but not here to avoid clever optimizers

         h = rnudif(xc(j), xcj)
         call fvec(xc,fz,n,j)
         xc(j) = xcj
         do i=1,n
            rjac(i,j) = (fz(i)-fc(i)) / h
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine fdjac2(xc,fc,n,epsm,fvec,fz,rjac,ldr,dsub,dsuper,
     *                  w,xstep)

      integer ldr,n,dsub,dsuper
      double precision  epsm
      double precision  rjac(ldr,*),fz(*),xc(*),fc(*)
      double precision  w(*), xstep(*)
      external fvec

c-------------------------------------------------------------------------
c
c     Compute a banded finite difference jacobian at the current point xc
c
c     Arguments
c
c     In       xc      Real(*)         current point
c     In       fc      Real(*)         function values at current point
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       fvec    Name            name of routine to evaluate f(x)
c     Wk       fz      Real(*)         workspace
c     Out      rjac    Real(ldr,*)     jacobian matrix at x
c                                        entry [i,j] is derivative of
c                                        f(i) wrt to x(j)
c     In       ldr     Integer         leading dimension of rjac
c     In       dsub    Integer         number of subdiagonals
c     In       dsuper  Integer         number of superdiagonals
c     Internal w       Real(*)         for temporary saving of xc
c     Internal xstep   Real(*)         stepsizes
c
c-------------------------------------------------------------------------

      integer i,j,k
      double precision  ndigit,p,h
      double precision  rnudif

      double precision Rten
      parameter(Rten=10d0)

      integer dsum

      dsum = dsub + dsuper + 1

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))

      do k=1,n
         xstep(k) = p + p * abs(xc(k))
      enddo

      do k=1,dsum
         do j=k,n,dsum
            w(j) = xc(j)
            xc(j) = xc(j) + xstep(j)
         enddo

         call fvec(xc,fz,n,n+k)
         do j=k,n,dsum
             call nuzero(n,rjac(1,j))
c            fdjac0 for why
c            doing this ensures that results for fdjac2 and fdjac0 will be identical
             h = rnudif(xc(j),w(j))
             xc(j) = w(j)
             do i=max(j-dsuper,1),min(j+dsub,n)
                rjac(i,j) = (fz(i)-fc(i)) / h
             enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      function nudnrm(n, d, x)
      integer n
      double precision  d(*), x(*)
      double precision nudnrm

c-------------------------------------------------------------------------
c
c     calculate  max( abs(d[*]) / max(x[*],1) )
c
c     Arguments
c
c     In   n        Integer       number of elements in d() and x()
c     In   d        Real(*)       vector d
c     In   x        Real(*)       vector x
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do i=1,n
         t = max(t, abs(d(i)) / max(abs(x(i)),Rone))
      enddo
      nudnrm = t

      return
      end

c-----------------------------------------------------------------------

      function nuxnrm(n, xn, xc)
      integer n
      double precision  xn(*), xc(*)
      double precision nuxnrm

c-------------------------------------------------------------------------
c
c     calculate  max( abs(xn[*]-xc[*]) / max(xn[*],1) )
c
c     Arguments
c
c     In   n        Integer       number of elements in xn() and xc()
c     In   xn       Real(*)       vector xn
c     In   xc       Real(*)       vector xc
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do i=1,n
         t = max(t, abs(xn(i)-xc(i)) / max(abs(xn(i)),Rone))
      enddo
      nuxnrm = t

      return
      end

c-----------------------------------------------------------------------

      function rnudif(x, y)
      double precision x, y
      double precision rnudif

c-------------------------------------------------------------------------
c
c     Return difference of x and y (x - y)
c
c     Arguments
c
c     In   x  Real      argument 1
c     In   y  Real      argument 2
c
c-------------------------------------------------------------------------

      rnudif = x - y
      return
      end

c-----------------------------------------------------------------------

      subroutine compmu(r,ldr,n,mu,y,ierr)

      integer ldr,n,ierr
      double precision r(ldr,*),mu,y(*)

c-------------------------------------------------------------------------
c
c     Compute a small perturbation mu for the (almost) singular matrix R.
c     mu is used in the computation of the Levenberg-Marquardt step.
c
c     Arguments
c
c     In       R       Real(ldr,*)     upper triangular matrix from QR
c     In       ldr     Integer         leading dimension of R
c     In       n       Integer         column dimension of R
c     Out      mu      Real            sqrt(l1 norm of R * infinity norm of R
c                                      * n * epsm * 100) designed to make
c                                        trans(R)*R + mu * I not singular
c     Wk       y       Real(*)         workspace for dlange
c     Out      ierr    Integer         0 indicating mu ok
c                                      3 indicating mu much too small
c
c-------------------------------------------------------------------------

      double precision aifnrm,al1nrm,epsm
      double precision dlantr
      double precision epsmch

      double precision Rhund
      parameter(Rhund=100d0)

c     get the infinity norm of R
c     get the l1 norm of R
      ierr = 0
      aifnrm = dlantr('I','U','N',n,n,r,ldr,y)
      al1nrm = dlantr('1','U','N',n,n,r,ldr,y)
      epsm = epsmch()
      mu = sqrt(n*epsm*Rhund)*aifnrm*al1nrm
c     matrix consists of zero's or near zero's
c     LM correction in liqrev will not work
      if( mu .le. Rhund*epsm ) then
         ierr = 3
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine cndjac(n,r,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)
      integer n,ldr,icdwrk(*),ierr
      double precision cndtol,rcond,r(ldr,*),rcdwrk(*)

c---------------------------------------------------------------------
c
c     Check r for singularity and/or ill conditioning
c
c     Arguments
c
c     In       n       Integer         dimension of problem
c     In       r       Real(ldr,*)     upper triangular R from QR decomposition
c     In       ldr     Integer         leading dimension of rjac
c     In       cndtol  Real            tolerance of test for ill conditioning
c                                       when rcond <= cndtol then ierr is set to 1
c                                       cndtol should be >= machine precision
c     Out      rcond   Real            inverse condition  of r
c     Wk       rcdwrk  Real(*)         workspace (for dtrcon)
c     Wk       icdwrk  Integer(*)      workspace (fordtrcon)
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian too ill-conditioned
c                                      2 indicating Jacobian completely singular
c
c---------------------------------------------------------------------

      integer i,info
      logical rsing
      double precision Rzero
      parameter(Rzero=0.0d0)

      ierr = 0

      rsing = .false.
      do i=1,n
         if( r(i,i) .eq. Rzero ) then
             rsing = .true.
         endif
      enddo

      if( rsing ) then
         ierr = 2
         rcond = Rzero
      else
         call dtrcon('1','U','N',n,r,ldr,rcond,rcdwrk,icdwrk,info)
         if( rcond .eq. Rzero ) then
             ierr = 2
         elseif( rcond .le. cndtol ) then
             ierr = 1
         endif
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nwfjac(x,scalex,f,fq,n,epsm,jacflg,fvec,mkjac,rjac,
     *                  ldr,xw,w,xstep)

      integer ldr,n,jacflg(*)
      double precision  epsm
      double precision  x(*),f(*),scalex(*),xw(*),w(*),xstep(*)
      double precision  rjac(ldr,*),fq(*)
      external fvec,mkjac

c-------------------------------------------------------------------------
c
c     Calculate the jacobian  matrix
c
c     Arguments
c
c     In       x       Real(*)         (scaled) current x values
c     In       scalex  Real(*)         scaling factors x
c     In       f       Real(*)         function values f(x)
c     Wk       fq      Real(*)         (internal) workspace
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       jacflg  Integer(*)      indicates how to compute jacobian
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting jacobian allowed when
c                                                   singular or illconditioned
c     In       fvec    Name            name of routine to evaluate f()
c     In       mkjac   Name            name of routine to evaluate jacobian
c     Out      rjac    Real(ldr,*)     jacobian matrix (unscaled)
c     In       ldr     Integer         leading dimension of rjac
c     Internal xw      Real(*)         used for storing unscaled x
c     Internal w       Real(*)         workspace for banded jacobian
c     Internal xstep   Real(*)         workspace for banded jacobian
c
c-------------------------------------------------------------------------

c     compute the finite difference or analytic jacobian at x

      call dcopy(n,x,1,xw,1)
      call vunsc(n,xw,scalex)
      if(jacflg(1) .eq. 0) then
         call fdjac0(xw,f,n,epsm,fvec,fq,rjac,ldr)
      elseif(jacflg(1) .eq. 2) then
         call fdjac2(xw,f,n,epsm,fvec,fq,rjac,ldr,jacflg(2),jacflg(3),
     *               w,xstep)
      else
         call mkjac(rjac,ldr,xw,n)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nwscjac(n,rjac,ldr,scalex)
      integer n, ldr
      double precision rjac(ldr,*), scalex(*)

c-------------------------------------------------------------------------
c
c     Scale jacobian
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     Inout    rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     In       scalex  Real(*)         scaling factors for x
c
c-------------------------------------------------------------------------

      integer j
      double precision t, Rone
      parameter(Rone=1.0d0)

      do j = 1,n
         t = Rone/scalex(j)
         call dscal(n,t,rjac(1,j),1)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine nwunscjac(n,rjac,ldr,scalex)
      integer n, ldr
      double precision rjac(ldr,*), scalex(*)

c-------------------------------------------------------------------------
c
c     Unscale jacobian
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     Inout    rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     In       scalex  Real(*)         scaling factors for x
c
c-------------------------------------------------------------------------

      integer j
      double precision t

      do j = 1,n
         t = scalex(j)
         call dscal(n,t,rjac(1,j),1)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine nwcpsx(n,rjac,ldr,scalex,epsm, mode)

      integer ldr,n,mode
      double precision  epsm
      double precision  scalex(*)
      double precision  rjac(ldr,*)

c-------------------------------------------------------------------------
c
c     Calculate scaling factors from the jacobian  matrix
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     In       rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     Out      scalex  Real(*)         scaling factors for x
c     In       epsm    Real            machine precision
c     In       mode    Integer         1: initialise, >1: adjust
c-------------------------------------------------------------------------

      integer k
      double precision  dnrm2

      if( mode .eq. 1 ) then
         do k=1,n
            scalex(k) = dnrm2(n,rjac(1,k),1)
            if( scalex(k) .le. epsm ) scalex(k) = 1
         enddo
      else if( mode .gt. 1 ) then
         do k=1,n
            scalex(k) = max(scalex(k),dnrm2(n,rjac(1,k),1))
         enddo
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine nwcpmt(n, x, scalex, factor, wrk, stepsiz)
      double precision x(*), scalex(*), wrk(*)
      double precision factor, stepsiz
      integer n

c-------------------------------------------------------------------------
c
c     Calculate maximum stepsize
c
c     Arguments
c
c     In       n       Integer     size of x
c     In       x       Real(*)     x-values
c     In       scalex  Real(*)     scaling factors for x
c     In       factor  Real        multiplier
c     Inout    wrk     Real(*)     workspace
c     Out      stepsiz Real        stepsize
c
c     Currently not used
c     Minpack uses this to calculate initial trust region size
c     Not (yet) used in this code because it doesn't seem to help
c     Manually setting an initial trust region size works better
c
c-------------------------------------------------------------------------

      double precision Rzero
      parameter(Rzero=0.0d0)

      double precision  dnrm2

      call dcopy(n,x,1,wrk,1)
      call vscal(n,wrk,scalex)
      stepsiz = factor * dnrm2(n,wrk,1)
      if( stepsiz .eq. Rzero ) stepsiz = factor
      return
      end

c-----------------------------------------------------------------------

      subroutine vscal(n,x,sx)

      integer n
      double precision  x(*),sx(*)

c-------------------------------------------------------------------------
c
c     Scale a vector x
c
c     Arguments
c
c     In       n       Integer         size of x
c     Inout    x       Real(*)         vector to scale
c     In       sx      Real(*)         scaling vector
c
c-------------------------------------------------------------------------

      integer i

      do i = 1,n
         x(i) = sx(i) * x(i)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine vunsc(n,x,sx)

      integer n
      double precision  x(*),sx(*)

c-------------------------------------------------------------------------
c
c     Unscale a vector x
c
c     Arguments
c
c     In       n       Integer         size of x
c     Inout    x       Real(*)         vector to unscale
c     In       sx      Real(*)         scaling vector
c
c-------------------------------------------------------------------------

      integer i

      do i = 1,n
         x(i) = x(i) / sx(i)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine nwfvec(x,n,scalex,fvec,f,fnorm,xw)

      integer n
      double precision  x(*),xw(*),scalex(*),f(*),fnorm
      external fvec

c-------------------------------------------------------------------------
c
c     Evaluate the function at current iterate x and scale its value
c
c     Arguments
c
c     In       x       Real(*)         x
c     In       n       Integer         size of x
c     In       scalex  Real(*)         scaling vector for x
c     In       fvec    Name            name of routine to calculate f(x)
c     Out      f       Real(*)         f(x)
c     Out      fnorm   Real            .5*||f(x)||**2
c     Internal xw      Real(*)         used for storing unscaled xc
c
c-------------------------------------------------------------------------

      double precision dnrm2

      double precision Rhalf
      parameter(Rhalf=0.5d0)

      call dcopy(n,x,1,xw,1)
      call vunsc(n,xw,scalex)
      call fvec(xw,f,n,0)

      fnorm = Rhalf * dnrm2(n,f,1)**2

      return
      end

c-----------------------------------------------------------------------

      function epsmch()

c     Return machine precision
c     Use Lapack routine

      double precision epsmch
      double precision dlamch
      external dlamch

c     dlamch('e') returns negeps (1-eps)
c     dlamch('p') returns 1+eps

      epsmch = dlamch('p')

      return
      end

c-----------------------------------------------------------------------

      function dblhuge()

c     Return largest double precision number
c     Use Lapack routine

      double precision dblhuge
      double precision dlamch
      external dlamch

c     dlamch('o') returns max double precision

      dblhuge = dlamch('o')

      return
      end
