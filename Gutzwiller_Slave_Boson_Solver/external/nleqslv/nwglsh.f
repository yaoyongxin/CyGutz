
      subroutine nwglsh(n,xc,fcnorm,d,g,sigma,stepmx,xtol,scalex,fvec,
     *                  xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)

      integer n,retcd,gcnt
      double precision  sigma,stepmx,xtol,fcnorm,fpnorm
      double precision  xc(*)
      double precision  d(*),g(*),xp(*),fp(*),xw(*)
      double precision  scalex(*)
      external fvec

      integer priter,iter

c-------------------------------------------------------------------------
c
c     Find a next acceptable iterate using geometric line search
c     along the newton direction
c
c     Arguments
c
c     In       n       Integer          dimension of problem
c     In       xc      Real(*)          current iterate
c     In       fcnorm  Real             0.5 * || f(xc) ||**2
c     In       d       Real(*)          newton direction
c     In       g       Real(*)          gradient at current iterate
c     In       sigma   Real             reduction factor for lambda
c     In       stepmx  Real             maximum stepsize
c     In       xtol    Real             relative step size at which
c                                       successive iterates are considered
c                                       close enough to terminate algorithm
c     In       scalex  Real(*)          diagonal scaling matrix for x()
c     In       fvec    Name             name of routine to calculate f()
c     In       xp      Real(*)          new x()
c     In       fp      Real(*)          new f(x)
c     In       fpnorm  Real             .5*||fp||**2
c     Out      xw      Real(*)          workspace for unscaling x
c
c     Out      retcd   Integer          return code
c                                         0 new satisfactory x() found
c                                         1 no  satisfactory x() found
c                                           sufficiently distinct from xc()
c
c     Out      gcnt    Integer          number of steps taken
c     In       priter  Integer           >0 if intermediate steps to be printed
c                                        -1 if no printing
c
c-------------------------------------------------------------------------

      integer i
      double precision  alpha,slope,rsclen,oarg(4)
      double precision  lambda,lamhi,lamlo
      double precision  ddot,dnrm2, nudnrm, ftarg
      double precision  dlen

      integer idamax

      parameter (alpha = 1.0d-4)

      double precision Rone
      parameter(Rone=1.0d0)

c     safeguard initial step size

      dlen = dnrm2(n,d,1)
      if( dlen .gt. stepmx ) then
          lamhi  = stepmx / dlen
      else
          lamhi  = Rone
      endif

c     compute slope  =  g-trans * d

      slope = ddot(n,g,1,d,1)

c     compute the smallest value allowable for the damping
c     parameter lambda ==> lamlo

      rsclen = nudnrm(n,d,xc)
      lamlo  = xtol / rsclen

c     initialization of retcd and lambda (linesearch length)

      retcd  = 2
      lambda = lamhi
      gcnt   = 0

      do while( retcd .eq. 2 )

c        compute next x

         do i=1,n
            xp(i) = xc(i) + lambda*d(i)
         enddo

c        evaluate functions and the objective function at xp

         call nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)
         gcnt = gcnt + 1
         ftarg = fcnorm + alpha * lambda * slope

         if( priter .gt. 0) then
            oarg(1) = lambda
            oarg(2) = ftarg
            oarg(3) = fpnorm
            oarg(4) = abs(fp(idamax(n,fp,1)))
            call nwlsot(iter,1,oarg)
         endif

c        test if the standard step produces enough decrease
c        of the objective function.
c        If not update lambda and compute a new next iterate

         if( fpnorm .le. ftarg ) then
            retcd = 0
         else
            lambda  = sigma * lambda
            if(lambda .lt. lamlo) then
               retcd = 1
            endif
         endif

      enddo

      return
      end
