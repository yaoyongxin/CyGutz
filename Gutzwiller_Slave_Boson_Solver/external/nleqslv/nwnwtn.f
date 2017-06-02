
      subroutine nwsolv(ldr,xc,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,cndtol,global,xscalm,
     *                  stepmx,delta,sigma,
     *                  rjac,wrk1,wrk2,wrk3,wrk4,fc,fq,dn,d,qtf,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,epsm,
     *                  fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter,
     *                  termcd)

      integer ldr,n,termcd,njcnt,nfcnt,iter
      integer maxit,jacflg(*),global,xscalm,qrwsiz
      integer outopt(*)
      double precision  xtol,ftol,btol,cndtol
      double precision  stepmx,delta,sigma,fpnorm,epsm
      double precision  rjac(ldr,*)
      double precision  xc(*),fc(*),xp(*),fp(*),dn(*),d(*)
      double precision  wrk1(*),wrk2(*),wrk3(*),wrk4(*)
      double precision  qtf(*),gp(*),fq(*)
      double precision  scalex(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      external fjac,fvec

c-----------------------------------------------------------------------
c
c     Solve system of nonlinear equations with Newton and global strategy
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
c     Inout    rjac    Real(ldr,*)     jacobian (n columns)
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
      double precision  dum(2),fcnorm,rcond
      logical fstjac
      integer priter

      integer idamax

c     initialization

      retcd = 0
      iter  = 0
      njcnt = 0
      nfcnt = 0
      ierr  = 0

      dum(1) = 0

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

      do while( termcd .eq. 0 )
         iter = iter + 1

         call nwnjac(rjac,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg,wrk1,
     *               wrk2,wrk3,
     *               xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn,
     *               qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)
c        - choose the next iterate xp by a global strategy

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
            call nwddlg(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  fvec,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         elseif(global .eq. 5) then
            call nwpdlg(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  fvec,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         elseif(global .eq. 6) then
            call nwmhlm(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  fvec,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         endif

         nfcnt = nfcnt + gcnt

c        - check stopping criteria for the new iterate xp

         call nwtcvg(xp,fp,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

         if(termcd .eq. 0) then
c           update xc, fc, and fcnorm
            call dcopy(n,xp,1,xc,1)
            call dcopy(n,fp,1,fc,1)
            fcnorm = fpnorm
         endif

      enddo

      if( outopt(3) .eq. 1 ) then
         call nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,fvec,fjac,rjac,
     *               ldr,wrk1,wrk2,wrk3)
      endif

      call vunsc(n,xp,scalex)

      return
      end
