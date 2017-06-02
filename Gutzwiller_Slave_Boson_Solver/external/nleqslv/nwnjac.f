      subroutine nwnjac(rjac,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg,
     *                  wrk1,wrk2,wrk3,
     *                  xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn,
     *                  qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

c-----------------------------------------------------------------------
c
c     Compute Jacobian matrix in xc, fc
c     scale it, compute gradient in xc and generate QR decomposition
c     calculate Newton step
c
c     Arguments
c
c     Out      rjac    Real(ldr,*)     jacobian (n columns)
c     In       ldr     Integer         leading dimension of rjac
c     In       n       Integer         dimensions of problem
c     In       xc      Real(*)         initial estimate of solution
c     Inout    fc      Real(*)         function values f(xc)
c     Wk       fq      Real(*)         workspace
c     In       fjac    Name            name of routine to calculate jacobian
c                                      (optional)
c     In       fvec    Name            name of routine to calculate f()
c     In       epsm    Real            machine precision
c     In       jacflg  Integer(*)      jacobian flag array
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting step allowed when
c                                                   singular or illconditioned
c     Wk       wrk1    Real(*)         workspace
c     Wk       wrk2    Real(*)         workspace
c     Wk       wrk3    Real(*)         workspace
c     In       xscalm  Integer         x scaling method
c                                        1 from column norms of first jacobian
c                                          increased if needed after first iteration
c                                        0 scaling user supplied
c     Inout    scalex  Real(*)         scaling factors x(*)
c     Out      gp      Real(*)         gradient at xp()
c     In       cndtol  Real            tolerance of test for ill conditioning
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     Out      dn      Real(*)         Newton step
c     Out      qtf     Real(*)         workspace for nwnstp
c     Out      rcond   Real            estimated inverse condition of R from QR
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c     Out      njcnt   Integer         number of jacobian evaluations
c     In       iter    Integer         iteration counter (used in scaling)
c     Inout    fstjac  logical         .true. if initial jacobian is available
c                                      on exit set to .false.
c     Out      ierr    Integer         error code
c                                        0 no error
c                                       >0 error in nwnstp (singular ...)
c
c-----------------------------------------------------------------------

      integer ldr,n,iter, njcnt, ierr
      integer jacflg(*),xscalm,qrwsiz
      logical fstjac
      double precision  epsm, cndtol, rcond
      double precision  rjac(ldr,*)
      double precision  xc(*),fc(*),dn(*)
      double precision  wrk1(*),wrk2(*),wrk3(*)
      double precision  qtf(*),gp(*),fq(*)
      double precision  scalex(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      external fjac,fvec

      logical stepadj
      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

c     evaluate the jacobian at the current iterate xc

      if( .not. fstjac ) then
         call nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,fvec,fjac,rjac,
     *               ldr,wrk1,wrk2,wrk3)
         njcnt = njcnt + 1
      else
         fstjac = .false.
      endif

c     if requested calculate x scale from jacobian column norms a la Minpack

      if( xscalm .eq. 1 ) then
         call vunsc(n,xc,scalex)
         call nwcpsx(n,rjac,ldr,scalex,epsm,iter)
         call vscal(n,xc,scalex)
      endif

      call nwscjac(n,rjac,ldr,scalex)

c     evaluate the gradient at the current iterate xc
c     gp = trans(Rjac) * fc
      call dgemv('T',n,n,Rone,rjac,ldr,fc,1,Rzero,gp,1)

c     get newton step
      stepadj = jacflg(4) .eq. 1
      call dcopy(n,fc,1,fq,1)
      call nwnstp(rjac,ldr,fq,n,cndtol, stepadj,
     *            wrk1,dn,qtf,ierr,rcond,
     *            rcdwrk,icdwrk,qrwork,qrwsiz)

c     save some data about jacobian for later output
      call nwsnot(0,ierr,rcond)

      return
      end

c-----------------------------------------------------------------------

      subroutine nwnstp(rjac,ldr,fn,n,cndtol, stepadj,
     *                  qraux,dn,qtf,ierr,rcond,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz)

      integer ldr,n,ierr,qrwsiz
      double precision  cndtol,rjac(ldr,*),qraux(*),fn(*)
      double precision  dn(*),qtf(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      double precision  rcond
      logical           stepadj

c-----------------------------------------------------------------------
c
c     Calculate the newton step
c
c     Arguments
c
c     Inout    rjac    Real(ldr,*)     jacobian matrix at current iterate
c                                      overwritten with QR decomposition
c     In       ldr     Integer         leading dimension of rjac
c     In       fn      Real(*)         function values at current iterate
c     In       n       Integer         dimension of problem
c     In       cndtol  Real            tolerance of test for ill conditioning
c     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
c     Inout    qraux   Real(*)         QR info from liqrfa (calling Lapack dgeqrf)
c     Out      dn      Real(*)         Newton direction
c     Out      qtf     Real(*)         trans(Q)*f()
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian ill-conditioned
c                                      2 indicating Jacobian completely singular
c                                      3 indicating almost zero LM correction
c     Out      rcond   Real            inverse condition of upper triangular R of QR
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c
c-----------------------------------------------------------------------

      integer info,k

      double precision Rone
      parameter(Rone=1.0d0)
      double precision mu

c     perform a QR factorization of rjac (simple Lapack routine)
c     check for singularity or ill conditioning
c     form qtf = trans(Q) * fn

      call liqrfa(rjac,ldr,n,qraux,qrwork,qrwsiz,ierr)

c     check for singularity or ill conditioning

      call cndjac(n,rjac,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)

c     compute qtf = trans(Q)*fn

      call dcopy(n,fn,1,qtf,1)
      call liqrqt(rjac, ldr, n, qraux, qtf, qrwork, qrwsiz, info)

      if( ierr .eq. 0 ) then
c         Normal Newton step
c         solve rjac*dn  =  -fn
c         ==> R*dn = - qtf

          call dcopy(n,qtf,1,dn,1)
          call dtrsv('U','N','N',n,rjac,ldr,dn,1)
          call dscal(n, -Rone, dn, 1)

      elseif( stepadj ) then
c         Adjusted Newton step
c         approximately from pseudoinverse(Jac+)
c         use mu to solve (trans(R)*R + mu*I*mu*I) * x = - trans(R) * fn
c         directly from the QR decomposition of R stacked with mu*I
c         a la Levenberg-Marquardt
          call compmu(rjac,ldr,n,mu,rcdwrk,ierr)
          if( ierr .eq. 0 ) then
             call liqrev(n,rjac,ldr,mu,qtf,dn,
     *                   rcdwrk(1+n),rcdwrk(2*n+1))
             call dscal(n, -Rone, dn, 1)

c            copy lower triangular Rjac to upper triangular
             do k=1,n
                call dcopy (n-k+1,rjac(k,k),1,rjac(k,k),ldr)
                rjac(k,k) = rcdwrk(1+n+k-1)
             enddo
          endif
      endif

      return
      end
