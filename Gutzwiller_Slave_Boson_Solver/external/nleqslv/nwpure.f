
      subroutine nwpure(n,xc,d,stepmx,scalex,fvec,
     *                  xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)

      integer n,retcd,gcnt
      double precision  stepmx,fpnorm
      double precision  xc(*)
      double precision  d(*),xp(*),fp(*),xw(*)
      double precision  scalex(*)
      external fvec

      integer priter,iter

c-------------------------------------------------------------------------
c
c     Find a next iterate using geometric line search
c     along the newton direction
c
c     Arguments
c
c     In       n       Integer          dimension of problem
c     In       xc      Real(*)          current iterate
c     In       d       Real(*)          newton direction
c     In       stepmx  Real             maximum stepsize
c     In       scalex  Real(*)          diagonal scaling matrix for x()
c     In       fvec    Name             name of routine to calculate f()
c     In       xp      Real(*)          new x()
c     In       fp      Real(*)          new f(x)
c     In       fpnorm  Real             .5*||fp||**2
c     Out      xw      Real(*)          workspace for unscaling x
c
c     Out      retcd   Integer          return code
c                                         0 new satisfactory x() found (!always)
c
c     Out      gcnt    Integer          number of steps taken
c     In       priter  Integer           >0 if intermediate steps to be printed
c                                        -1 if no printing
c
c-------------------------------------------------------------------------

      integer i
      double precision  oarg(3)
      double precision  lambda
      double precision  dnrm2
      double precision  dlen

      integer idamax

      double precision Rone
      parameter(Rone=1.0d0)

c     safeguard initial step size

      dlen = dnrm2(n,d,1)
      if( dlen .gt. stepmx ) then
          lambda = stepmx / dlen
      else
          lambda = Rone
      endif

      retcd  = 0
      gcnt   = 1

c     compute the next iterate xp

      do i=1,n
         xp(i) = xc(i) + lambda*d(i)
      enddo

c     evaluate functions and the objective function at xp

      call nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)

      if( priter .gt. 0) then
         oarg(1) = lambda
         oarg(2) = fpnorm
         oarg(3) = abs(fp(idamax(n,fp,1)))
         call nwprot(iter,1,oarg)
      endif

      return
      end
