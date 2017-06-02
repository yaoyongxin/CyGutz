      subroutine limhpar(R, ldr, n, sdiag, qtf, dn, dnlen,
     *                   glen, delta, mu, d, work)
      integer ldr, n
      double precision R(ldr,*), sdiag(*)
      double precision qtf(*),dn(*), dnlen, glen, d(*),work(*)
      double precision delta, mu

c-----------------------------------------------------------------------------
c
c     Arguments
c
c       Inout  R      Real(ldr,*)     upper triangular matrix R from QR (unaltered)
c                                     strict lower triangle contains
c                                        transposed strict upper triangle of the upper
c                                        triangular matrix S.
c
c       In     n      Integer         order of R.
c
c       In     ldr    Integer         leading dimension of the array R.
c
c       Out    sdiag  Real(*)         vector of size n, containing the
c                                     diagonal elements of the upper
c                                     triangular matrix S.
c
c       In     qtr    Real(*)         trans(Q)*f vector of size n
c       In     dn     Real(*)         Newton step
c       In     dnlen  Real            length Newton step
c       In     glen   Real            length gradient vector
c
c       Inout  mu     Real            Levenberg-Marquardt parameter
c       In     delta  Real            size of trust region (euclidian norm)
c
c       Out    d      Real(*)         vector with step with norm very close to delta
c       Out    work   Real(*)         workspace of size n.
c
c     Description
c
c     determine Levenberg-Marquardt parameter mu such that
c     norm[(R**T R + mu*I)**(-1) * qtf] - delta approximately 0
c     See description in liqrev.f for further details
c
c     Algorithm comes from More: The Levenberg-Marquardt algorithm, Implementation and Theory
c     Lecture Notes in Mathematics, 1978, no. 630.
c     uses liqrev (in file liqrev.f) which is based on Nocedal's method (see comments in file)
c-----------------------------------------------------------------------------

      double precision phi, pnorm, qnorm, mulo, muhi,dmu, sqmu
      integer iter
      logical done
      double precision dnrm2

      double precision Rone
      parameter(Rone=1.0D0)

      phi = dnlen - delta
      muhi = glen/delta

      call dcopy(n,dn,1,d,1)
      call dscal(n, Rone/dnlen, d, 1)

c     solve R**T * x = dn
      call dtrsv("U","T","N",n,R,ldr,d,1)
      qnorm = dnrm2(n,d,1)
      mulo = (phi/dnlen)/qnorm**2
      mu = mulo

      iter = 0
      done = .false.
      do while( .not. done )
          iter = iter + 1
          sqmu = sqrt(mu)
          call liqrev(n, R, ldr, sqmu, qtf, d, sdiag, work)
          pnorm = dnrm2(n,d,1)
          call dcopy(n,d,1,work,1)
          call dtrstt(R, ldr, n, sdiag, work)
          done = abs(pnorm-delta) .le. .1d0*delta .or. iter .gt. 5
          if( .not. done ) then
              qnorm = dnrm2(n,work,1)
              if( pnorm .gt. delta ) then
                  mulo = max(mulo,mu)
              else if( pnorm .lt. delta ) then
                  muhi = min(muhi,mu)
              endif
              dmu = (pnorm-delta)/delta * (pnorm/qnorm)**2
              mu = max(mulo, mu + dmu)
          endif
      enddo
      return
      end
