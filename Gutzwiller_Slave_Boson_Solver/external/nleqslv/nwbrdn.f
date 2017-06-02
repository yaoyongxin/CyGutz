
      subroutine brsolv(ldr,xc,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,cndtol,global,xscalm,
     *                  stepmx,delta,sigma,
     *                  rjac,r,wrk1,wrk2,wrk3,wrk4,fc,fq,dn,d,qtf,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,epsm,
     *                  fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter,
     *                  termcd)

      integer ldr,n,termcd,njcnt,nfcnt,iter
      integer maxit,jacflg(*),global,xscalm,qrwsiz
      integer outopt(*)
      double precision  xtol,ftol,btol,cndtol
      double precision  stepmx,delta,sigma,fpnorm,epsm
      double precision  rjac(ldr,*),r(ldr,*)
      double precision  xc(*),fc(*),xp(*),fp(*),dn(*),d(*)
      double precision  wrk1(*),wrk2(*),wrk3(*),wrk4(*)
      double precision  qtf(*),gp(*),fq(*)
      double precision  scalex(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      external fjac,fvec

c-----------------------------------------------------------------------
c
c     Solve system of nonlinear equations with Broyden and global strategy
c
c
c     Arguments
c
c     In       ldr     Integer         leading dimension of rjac
c     In       xc      Real(*)         initial estimate of solution
c     In       n       Integer         dimensions of problem
c     Inout    scalex  Real(*)         scaling factors x(*)
c     In       maxit   Integer         maximum number of allowable iterations
c     In       jacflg  Integer(*)      jacobian flag array
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting step allowed when
c                                                   singular or illconditioned
c     In       xtol    Real            tolerance at which successive iterates x()
c                                      are considered close enough to
c                                      terminate algorithm
c     In       ftol    Real            tolerance at which function values f()
c                                      are considered close enough to zero
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    cndtol  Real            tolerance of test for ill conditioning
c     In       global  Integer         global strategy to use
c                                        1 cubic linesearch
c                                        2 quadratic linesearch
c                                        3 geometric linesearch
c                                        4 double dogleg
c                                        5 powell dogleg
c                                        6 hookstep (More-Hebden Levenberg-Marquardt)
c     In       xscalm  Integer         x scaling method
c                                        1 from column norms of first jacobian
c                                          increased if needed after first iteration
c                                        0 scaling user supplied
c     In       stepmx  Real            maximum allowable step size
c     In       delta   Real            trust region radius
c     In       sigma   Real            reduction factor geometric linesearch
c     Inout    rjac    Real(ldr,*)     jacobian (n columns)(compact QR decomposition/Q matrix)
c     Inout    r       Real(ldr,*)     stored R from QR decomposition
c     Wk       wrk1    Real(*)         workspace
c     Wk       wrk2    Real(*)         workspace
c     Wk       wrk3    Real(*)         workspace
c     Wk       wrk4    Real(*)         workspace
c     Inout    fc      Real(*)         function values f(xc)
c     Wk       fq      Real(*)         workspace
c     Wk       dn      Real(*)         workspace
c     Wk       d       Real(*)         workspace
c     Wk       qtf     Real(*)         workspace
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c     In       epsm    Real            machine precision
c     In       fjac    Name            name of routine to calculate jacobian
c                                      (optional)
c     In       fvec    Name            name of routine to calculate f()
c     In       outopt  Integer(*)      output options
c     Out      xp      Real(*)         final x()
c     Out      fp      Real(*)         final f(xp)
c     Out      gp      Real(*)         gradient at xp()
c     Out      njcnt   Integer         number of jacobian evaluations
c     Out      nfcnt   Integer         number of function evaluations
c     Out      iter    Integer         number of (outer) iterations
c     Out      termcd  Integer         termination code
c
c-----------------------------------------------------------------------

      integer gcnt,retcd,ierr
      double precision  dum(2),dlt0,fcnorm,rcond
      logical fstjac
      logical jacevl,jacupd
      logical stepadj
      integer priter

      integer idamax

      double precision Rone
      parameter(Rone=1.0d0)

c     initialization

      retcd = 0
      iter  = 0
      njcnt = 0
      nfcnt = 0
      ierr  = 0

      dum(1) = 0
      dlt0 = delta

      if( outopt(1) .eq. 1 ) then
         priter = 1
      else
         priter = -1
      endif

c     evaluate function

      call vscal(n,xc,scalex)
      call nwfvec(xc,n,scalex,fvec,fc,fcnorm,wrk1)

c     evaluate user supplied or finite difference jacobian and check user supplied
c     jacobian, if requested

      fstjac = .false.
      if(mod(jacflg(1),2) .eq. 1) then

        if( outopt(2) .eq. 1 ) then
           fstjac = .true.
           njcnt = njcnt + 1
           call nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,fvec,fjac,rjac,
     *                 ldr,wrk1,wrk2,wrk3)
           call chkjac(rjac,ldr,xc,fc,n,epsm,jacflg,scalex,
     *                 fq,wrk1,wrk2,fvec,termcd)
           if(termcd .lt. 0) then
c              copy initial values
               call dcopy(n,xc,1,xp,1)
               call dcopy(n,fc,1,fp,1)
               call vunsc(n,xp,scalex)
               fpnorm = fcnorm
               return
           endif
        endif

      endif

c     check stopping criteria for input xc

      call nwtcvg(xc,fc,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

      if(termcd .gt. 0) then
          call dcopy(n,xc,1,xp,1)
          call dcopy(n,fc,1,fp,1)
          fpnorm = fcnorm
          if( outopt(3) .eq. 1 .and. .not. fstjac ) then
             njcnt = njcnt + 1
             call nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,fvec,fjac,rjac,
     *                   ldr,wrk1,wrk2,wrk3)
          endif
          return
      endif

      if( priter .gt. 0 ) then

         dum(1) = fcnorm
         dum(2) = abs(fc(idamax(n,fc,1)))

         if( global .eq. 0 ) then
            call nwprot(iter, -1, dum)
         elseif( global .le. 3 ) then
            call nwlsot(iter,-1,dum)
         elseif( global .eq. 4 ) then
            call nwdgot(iter,-1,0,dum)
         elseif( global .eq. 5 ) then
            call nwpwot(iter,-1,0,dum)
         elseif( global .eq. 6 ) then
            call nwmhot(iter,-1,0,dum)
         endif

      endif

      jacevl  = .true.
      stepadj = jacflg(4) .eq. 1

      do while( termcd .eq. 0 )
         iter = iter+1

         if( jacevl ) then

            call nwbjac(rjac,r,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg,
     *                  wrk1,wrk2,wrk3,
     *                  xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn,
     *                  qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

         else

c          - get broyden step
c          - calculate approximate gradient

            call dcopy(n,fc,1,fq,1)
            call brodir(rjac,ldr,r,fq,n,cndtol, stepadj,
     *                  dn,qtf,ierr,rcond,rcdwrk,icdwrk)

            if( ierr .eq. 0 ) then
               call dcopy(n,qtf,1,gp,1)
               call dtrmv('U','T','N',n,r,ldr,gp,1)
            endif
         endif
c      - choose the next iterate xp by a global strategy

         if( ierr .gt. 0 ) then
c           jacobian singular or too ill-conditioned
            call nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,iter)
         elseif(global .eq. 0) then
            call nwpure(n,xc,dn,stepmx,scalex,
     *                  fvec,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 1) then
            call nwclsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex,
     *                  fvec,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 2) then
            call nwqlsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex,
     *                  fvec,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 3) then
            call nwglsh(n,xc,fcnorm,dn,gp,sigma,stepmx,btol,scalex,
     *                  fvec,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 4) then
            call nwddlg(n,r,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  fvec,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         elseif(global .eq. 5) then
            call nwpdlg(n,r,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  fvec,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         elseif(global .eq. 6) then
            call nwmhlm(n,r,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  fvec,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         endif

         nfcnt = nfcnt + gcnt

c      - check stopping criteria for the new iterate xp

         call nwtcvg(xp,fp,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

         if( termcd .eq. 3 .and. .not. jacevl ) then
c           global strategy failed but jacobian is out of date
c           try again with proper jacobian
c           reset trust region radius

            jacevl = .true.
            jacupd = .false.
            delta = dlt0
            termcd = 0

         elseif(termcd .gt. 0) then
            jacupd = .false.
         else
            jacupd = .true.
            jacevl = .false.
         endif

         if( jacupd ) then
c           perform Broyden update of current jacobian
c           update xc, fc, and fcnorm
            call brupdt(n,rjac,r,ldr,xc,xp,fc,fp,epsm,
     *                  wrk1,wrk2,wrk3)
            call dcopy(n,xp,1,xc,1)
            call dcopy(n,fp,1,fc,1)
            fcnorm = fpnorm
         endif

      enddo

      if( outopt(3) .eq. 1 ) then
c        final update of jacobian
         call brupdt(n,rjac,r,ldr,xc,xp,fc,fp,epsm,
     *               wrk1,wrk2,wrk3)
c        reconstruct Broyden matrix
c        calculate Q * R where Q is overwritten by result
c        Q is in rjac and R is in r
         call dtrmm('R','U','N','N',n,n,Rone,r,n,rjac,n)
c        unscale
         call nwunscjac(n,rjac,ldr,scalex)
      endif

      call vunsc(n,xp,scalex)

      return
      end

c-----------------------------------------------------------------------

      subroutine brupdt(n,q,r,ldr,xc,xp,fc,fp,epsm,dx,df,wa)
      integer n,ldr
      double precision  q(ldr,*),r(ldr,*)
      double precision  xc(*),xp(*),fc(*),fp(*),dx(*),df(*),wa(*)
      double precision  epsm

c-----------------------------------------------------------------------
c
c     Calculate new Q and R from rank-1 update with xp-xc and fp-fc
c     using Broyden method
c
c     Arguments
c
c     In       n       Integer         size of xc() etc.
c     Inout    Q       Real(ldr,n)     orthogonal matrix Q from QR
c                                       On output updated Q
c     Inout    R       Real(ldr,n)     upper triangular R from QR
c                                       On output updated R
c     In       ldr     Integer         leading dimension of Q and R
c     In       xc      Real(*)         current x() values
c     In       xp      Real(*)         new     x() values
c     In       fc      Real(*)         current f(xc)
c     In       fp      Real(*)         new     f(xp)
c     In       epsm    Real            machine precision
c     Wk       dx      Real(*)         workspace
c     Wk       df      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c
c-----------------------------------------------------------------------

      integer i
      double precision  eta,sts
      double precision  dnrm2
      logical doupdt

      double precision Rzero, Rone, Rtwo, Rhund
      parameter(Rzero=0.0d0, Rone=1.0d0, Rtwo=2.0d0, Rhund=100d0)

      eta    = Rhund * Rtwo * epsm
      doupdt = .false.

      do i=1,n
         dx(i) = xp(i) - xc(i)
         df(i) = fp(i) - fc(i)
      enddo

c     clear lower triangle

      do i=1,n-1
         call nuzero(n-i,r(i+1,i))
      enddo

c     calculate df - B*dx = df - Q*R*dx
c     wa = R*dx
c     df = df - Q*(R*dx) (!not really needed if qrupdt were to be changed)
c     do not update with noise

      call dcopy(n,dx,1,wa,1)
      call dtrmv('U','N','N',n,r,ldr,wa,1)
      call dgemv('N',n,n,-Rone,q,ldr,wa,1,Rone,df,1)

      do i=1,n
         if( abs(df(i)) .gt. eta*( abs(fp(i)) + abs(fc(i)) ) ) then
            doupdt = .true.
         else
            df(i)  = Rzero
         endif
      enddo

      if( doupdt ) then
c        equation 8.3.1 from Dennis and Schnabel (page 187)(Siam edition)
         sts = dnrm2(n,dx,1)
         call dscal(n,Rone/sts,dx,1)
         call dscal(n,Rone/sts,df,1)
         call liqrup(q,ldr,n,r,ldr,df,dx,wa)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine brodir(q,ldr,r,fn,n,cndtol,stepadj,dn,qtf,
     *                  ierr,rcond,rcdwrk,icdwrk)

      integer ldr,n,ierr
      double precision  cndtol,q(ldr,*),r(ldr,*),fn(*)
      double precision  dn(*),qtf(*)
      double precision  rcdwrk(*)
      integer           icdwrk(*)
      double precision  rcond
      logical           stepadj

c-----------------------------------------------------------------------
c
c     Calculate the approximate newton direction
c
c     Arguments
c
c     Inout    Q       Real(ldr,*)     Q part from QR at current iterate
c     In       ldr     Integer         leading dimension of Q and R
c     In       R       Real(ldr,*)     upper triangular R from QR decomposition
c     In       fn      Real(*)         function values at current iterate
c     In       n       Integer         dimension of problem
c     In       cndtol  Real            tolerance of test for ill conditioning
c     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
c     Out      dn      Real(*)         Newton direction
c     Out      qtf     Real(*)         trans(Q)*f()
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian ill-conditioned
c                                      2 indicating Jacobian completely singular
c                                      3 indicating almost zero LM correction
c     Out      rcond   Real            inverse condition of matrix
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c
c     QR decomposition with no pivoting.
c
c-----------------------------------------------------------------------

      integer k
      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)
      double precision mu

c     check for singularity or ill conditioning

      call cndjac(n,r,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)

      if( ierr .eq. 0 ) then
c         form qtf = trans(Q) * fn

          call dgemv('T',n,n,Rone,q,ldr,fn,1,Rzero,qtf,1)

c         solve rjac*dn  =  -fn
c         ==> R*dn = - qtf

          call dcopy(n,qtf,1,dn,1)
          call dtrsv('U','N','N',n,r,ldr,dn,1)
          call dscal(n, -Rone, dn, 1)

      elseif( stepadj ) then
c         call intpr('ierr brodir', 12,ierr,1)
c         Adjusted Newton step
c         approximately from pseudoinverse(Jac+)
c         compute qtf = trans(Q)*fn

c         form qtf = trans(Q) * fn

          call dgemv('T',n,n,Rone,q,ldr,fn,1,Rzero,qtf,1)

c         use mu to solve (trans(R)*R + mu*I*mu*I) * x = - trans(R) * fn
c         directly from the QR decomposition of R stacked with mu*I
c         a la Levenberg-Marquardt
          call compmu(r,ldr,n,mu,rcdwrk,ierr)
          if( ierr .eq. 0 ) then
             call liqrev(n,r,ldr,mu,qtf,dn,
     *                   rcdwrk(1+n),rcdwrk(2*n+1))
             call dscal(n, -Rone, dn, 1)

c            copy lower triangular Rjac to upper triangular
             do k=1,n
                call dcopy (n-k+1,r(k,k),1,r(k,k),ldr)
                r(k,k) = rcdwrk(1+n+k-1)
             enddo
          endif
      endif
      call nwsnot(1,ierr,rcond)

      return
      end
