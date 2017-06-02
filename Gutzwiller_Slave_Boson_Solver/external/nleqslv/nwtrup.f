
      subroutine nwtrup(n,fcnorm,g,sc,nwtstep,stepmx,xtol,delta,
     *                  fpred,retcd,xprev,fpnsav,fprev,xp,fp,
     *                  fpnorm)

      integer n,retcd
      double precision  fcnorm,stepmx,xtol,delta,fpred,fpnsav,fpnorm
      double precision  xp(*),g(*)
      double precision  sc(*),xprev(*),fprev(*),fp(*)
      logical nwtstep

c-------------------------------------------------------------------------
c
c     Decide whether to accept xp=xc+sc as the next iterate
c     and updates the trust region delta
c
c     Arguments
c
c     In       n       Integer         size of xc()
c     In       fcnorm  Real            .5*||f(xc)||**2
c     In       g       Real(*)         gradient at xc()
c     In       sc      Real(*)         current step
c     In       nwtstep Logical         true if sc is newton direction
c     In       stepmx  Real            maximum step size
c     In       xtol    Real            minimum step tolerance
c     Inout    delta   Real            trust region radius
c     In       fpred   Real            predicted value of .5*||f()||**2
c
c     Inout    retcd   Integer         return code
c                                       0 xp accepted as next iterate;
c                                         delta trust region for next iteration.
c
c                                       1 xp unsatisfactory but
c                                         accepted as next iterate because
c                                         xp-xc .lt. smallest allowable
c                                         step length.
c
c                                       2 f(xp) too large.
c                                         continue current iteration with
c                                         new reduced delta.
c
c                                       3 f(xp) sufficiently small, but
c                                         quadratic model predicts f(xp)
c                                         sufficiently well to continue current
c                                         iteration with new doubled delta.
c
c                                      On first entry, retcd must be 4
c
c     Wk       xprev   Real(*)         (internal) work
c     Wk       fpnsav  Real            (internal) work
c     Wk       fprev   Real(*)         (internal) work
c     Inout    xp      Real(*)         new iterate x()
c     Inout    fp      Real(*)         new f(xp)
c     Inout    fpnorm  Real            new .5*||f(xp)||**2
c
c-------------------------------------------------------------------------

      double precision  ared,pred,slope,sclen,rln,dltmp
      double precision  dnrm2,ddot,nudnrm
      logical ret3ok

      double precision Rone, Rtwo, Rthree, Rfour, Rten
      double precision Rhalf, Rpten
      parameter(Rpten = 0.1d0)
      parameter(Rhalf=0.5d0)
      parameter(Rone=1.0d0, Rtwo=2.0d0, Rthree=3.0d0, Rfour=4.0d0)
      parameter(Rten=10.0d0)

      double precision Rp99,Rp4, Rp75
      parameter(Rp99=Rone-Rten**(-2), Rp4=Rten**(-4), Rp75=Rthree/Rfour)

      double precision alpha
      parameter(alpha = Rp4)

c     pred measures the predicted reduction in the function value

      ared  = fpnorm - fcnorm
      pred  = fpred  - fcnorm
      slope = ddot(n,g,1,sc,1)

      if(retcd .ne. 3) then
         ret3ok = .false.
      else
         ret3ok = fpnorm .ge. fpnsav .or. ared .gt. alpha * slope
      endif

      if(retcd .eq. 3 .and. ret3ok) then

c        reset xp,fp,fpnorm to saved values and terminate global step

         retcd = 0
         call dcopy(n,xprev,1,xp,1)
         call dcopy(n,fprev,1,fp,1)
         fpnorm = fpnsav
c        reset delta to initial value
c        but beware
c           if the trial step was a Newton step then delta is reset to
c           .5 * length(Newton step) which will be smaller
c           because delta is set to length(Newton step) elsewhere
c           see nwddlg.f and nwpdlg.f
         delta  = Rhalf*delta

      elseif(ared .gt. alpha * slope) then

c        fpnorm too large (decrease not sufficient)

         rln = nudnrm(n,sc,xp)
         if(rln .lt. xtol) then

c           cannot find satisfactory xp sufficiently distinct from xc

            retcd = 1

         else

c           reduce trust region and continue current global step

            retcd = 2
            sclen = dnrm2(n,sc,1)
            dltmp = -slope*sclen/(Rtwo*(ared-slope))

            if(dltmp .lt. Rpten*delta) then
               delta = Rpten*delta
            else
               delta = min(Rhalf*delta, dltmp)
            endif

         endif

      elseif(retcd .ne. 2 .and. (abs(pred-ared) .le. Rpten*abs(ared))
     *      .and. (.not. nwtstep) .and. (delta .le. Rp99*stepmx)) then

c        pred predicts ared very well
c        attempt a doubling of the trust region and continue global step
c        when not taking a newton step and trust region not at maximum

         call dcopy(n,xp,1,xprev,1)
         call dcopy(n,fp,1,fprev,1)
         fpnsav = fpnorm
         delta  = min(Rtwo*delta,stepmx)
         retcd  = 3

      else

c        fpnorm sufficiently small to accept xp as next iterate.
c        Choose new trust region.

         retcd = 0
         if(ared .ge. Rpten*pred) then

c           Not good enough. Decrease trust region for next iteration

            delta = Rhalf*delta
         elseif( ared .le. Rp75*pred ) then

c           Wonderful. Increase trust region for next iteration

            delta = min(Rtwo*delta,stepmx)
         endif

      endif

      return
      end
