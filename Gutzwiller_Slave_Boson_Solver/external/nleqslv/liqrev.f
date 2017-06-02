
c-----------------------------------------------------------------------------

      subroutine liqrev(n,r,ldr,diag,b,x,sdiag,wrk)
      integer n,ldr
      double precision  r(ldr,*),b(*),x(*),sdiag(*),wrk(*)
      double precision diag

c-----------------------------------------------------------------------------
c
c     Arguments
c
c       In     n      Integer         order of R.
c       Inout  R      Real(ldr,*)     upper triangular matrix R from QR
c                                     unaltered
c                                     strict lower triangle contains
c                                        transposed strict upper triangle of the upper
c                                        triangular matrix S.
c
c       In     diag   Real            scalar for matrix D
c
c       In     ldr    Integer         leading dimension of the array R.
c       In     b      Real(*)         vector of size n
c
c       Out    x      Real(*)         vector of size n
c                                     on output contains least squares solution
c                                     of the system R*x = b, D*x = 0.
c
c       Out    sdiag  Real(*)         vector of size n, containing the
c                                     diagonal elements of the upper
c                                     triangular matrix S.
c
c       Out    wrk    Real(*)         workspace of size n.
c
c     Description
c
c     Given an n by n upper triangular matrix R, a diagonal matrix D with positive entries
c     and an n-vector b, determine an x which solves the system
c
c         |R*x| = |b|
c         |D*x| = |0|
c
c     in the least squares sense where D=diag*I.
c     This routine can be used for two different purposes.
c     The first is to provide a method of slightly modifying a singular or ill-conditioned matrix.
c     The second is for calculating a least squares solution to the above problem within
c     the context of e.g. a Levenberg-Marquardt algorithm combined with a More-Hebden algorithm
c     to determine a value of D (diagonal mu) such that x has a predetermined 2-norm.
c
c     The routine could also be used when the matrix R from the QR decomposition of a Jacobian
c     is ill-conditioned (or singular). Then it is difficult to calculate a Newton step
c     accurately (Dennis and Schnabel). D&S advise perturbing trans(J)*J with a positive
c     diagonal matrix.
c
c     The idea is  to solve (J^T * J + mu*I)x=b where mu is a small positive number.
c     Calculation of mu must be done in the calling routine.
c     Using a QR decomposition of J solving this system
c     is equivalent solving (R^T*R + mu*I)x=b, where R comes from the QR decomposition.
c     Solving this system is equivalent to solving the above least squares problem with the
c     elements of the matrix D set to sqrt(mu) which should be done in the calling routine.
c
c     On output the routine also provides an upper triangular matrix S such that
c     (see description of arguments above for the details)
c
c         (trans(R)*R + D*D) = trans(S)*S .
c
c     Method used here is described in
c     Nocedal and Wright, 2006, Numerical Optimization, Springer, ISBN 978-0-387-30303-1
c     page 258--261 (second edition)
c-----------------------------------------------------------------------------

      integer j,k
      double precision  bj,c,s,sum,temp
      double precision  ddot
      double precision Rzero
      parameter(Rzero=0.0d0)

c     copy R and b to preserve input and initialise S.
c     Save the diagonal elements of R in wrk.
c     Beware: the algorithm operates on an upper triangular matrix,
c     which is stored in lower triangle of R.
c
      do j=1,n
         call dcopy(n-j+1,r(j,j),ldr,r(j,j),1)
         wrk(j) = r(j,j)
      enddo
      call dcopy(n,b,1,x,1)

c     eliminate the diagonal matrix D using givens rotations.
c     Nocedal method: start at the bottom right
c     at end of loop R contains diagonal of S
c     save in sdiag and restore original diagonal of R

      do j=n,1,-1

c        initialise the row of D to be eliminated

         call nuzero(n-j+1,sdiag(j))
         sdiag(j) = diag

c        the transformations to eliminate the row of D

         bj = Rzero
         do k=j,n

c           determine a givens rotation which eliminates the
c           appropriate element in the current row of D.
c           accumulate the transformation in the row of S.

c           eliminate the diagonal element in row j of D
c           this generates fill-in in columns [j+1 .. n] of row j of D
c           successively eliminate the fill-in with givens rotations
c           for R[j+1,j+1] and D[j,j+1].
c           rows of R have been copied into the columns of R initially (see above)
c           perform all operations on those columns to preserve the original R

            if (sdiag(k) .ne. Rzero) then

               call nuvgiv(r(k,k),sdiag(k),c,s)
               if( k .lt. n ) then
                   call drot(n-k,r(k+1,k),1,sdiag(k+1),1,c,s)
               endif

c              compute the modified element of (b,0).

               temp =  c*x(k) + s*bj
               bj   = -s*x(k) + c*bj
               x(k) = temp

            endif

         enddo

      enddo

c     retrieve diagonal of S from diagonal of R
c     restore original diagonal of R

      do k=1,n
         sdiag(k) = r(k,k)
         r(k,k) = wrk(k)
      enddo

c     x now contains modified b
c     solve trans(S)*x = x
c     still to be done: guard against division by 0 to be absolutely safe
c     call dblepr('liqrev sdiag', 12, sdiag, n)
      x(n) = x(n) / sdiag(n)
      do j=n-1,1,-1
         sum  = ddot(n-j,r(j+1,j),1,x(j+1),1)
         x(j) = (x(j) - sum)/sdiag(j)
      enddo

      return
      end

c ----------------------------------------------------------------------

      subroutine dtrstt(S,ldr,n,sdiag,x)
      integer ldr, n
      double precision S(ldr,*), sdiag(*), x(*)
      integer j
      double precision sum, ddot

c     solve S*x = x where x is the result from subroutine liqrev
c     S is a lower triangular matrix with diagonal entries in sdiag()
c     and here it is in the lower triangular part of R as returned by liqrev

      x(1) = x(1) / sdiag(1)
      do j=2,n
         sum  = ddot(j-1,S(j,1),n,x,1)
         x(j) = (x(j) - sum)/sdiag(j)
      enddo

      return
      end

c ----------------------------------------------------------------------

      subroutine nuvgiv(x,y,c,s)
      double precision x,y,c,s

c     Parameters
c
c     Inout   x     Real       x input / c*x+s*y on output
c     Inout   y     Real       y input / 0       on output
c     Out     c     Real       c of tranformation (cosine)
c     Out     s     Real       s of tranformation (  sine)
c
c     Description
c
c     Nuvgiv calculates the givens rotator
c
c             |  c   s |
c         G = |        |
c             | -s   c |
c
c     with  c*c+s*s=1
c
c     for which G * | x | = | z |
c                   | y |   | 0 |
c
c     resulting in
c
c            c * x + s * y = z
c           -s * x + c * y = 0   ==>  s/c = y/x or c/s = x/y
c
c     Use Lapack dlartg routine
c     return c and s and the modified x and y

      double precision t

      double precision Rzero
      parameter(Rzero=0.0d0)

      call dlartg(x,y,c,s,t)
      x = t
      y = Rzero
      return
      end
