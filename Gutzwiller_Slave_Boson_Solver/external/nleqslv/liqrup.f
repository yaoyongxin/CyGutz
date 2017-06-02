      subroutine liqrup(Q,ldq,n,R,ldr,u,v,wk)
      integer ldq,n,ldr
      double precision Q(ldq,*),R(ldr,*),u(*),v(*),wk(*)

c-----------------------------------------------------------------------------
c
c     Arguments
c
c     Inout  Q       Real(ldq,n)      orthogonal matrix from QR
c     In     ldq     Integer          leading dimension of Q
c     In     n       Integer          order of Q and R
c     Inout  R       Real(ldr,n)      upper triangular matrix R from QR
c     In     ldr     Integer          leading dimension of R
c     In     u       Real(*)          vector u of size n
c     In     v       Real(*)          vector v of size n
c     Out    wk      Real(*)          workspace of size n
c
c     on return
c
c        Q       Q is the matrix with orthonormal columns in a QR
c                decomposition of the matrix B = A + u*v'
c
c        R       R is the upper triangular matrix in a QR
c                decomposition of the matrix B = A + u*v'
c
c     Description
c
c     The matrices Q and R are a QR decomposition of a square matrix
c     A = Q*R.
c     Given Q and R, qrupdt computes a QR decomposition of the rank one
c     modification B = A + u*trans(v) of A. Here u and v are vectors.
c
c     Source : procedure outlined in Dennis & Schnabel (Appendix A)
c              Algorithm 3.1.4 and 3.4.1a
c              modified (to use Lapack routines and more)
c
c-----------------------------------------------------------------------------

c     Local variables and functions

      integer k,i
      double precision  ddot

c     calculate wk = trans(Q)*u

      do i=1,n
         wk(i) = ddot(n,Q(1,i),1,u,1)
      enddo

c     zero components wk(n),wk(n-1)...wk(2)
c     and apply rotators to R and Q.

      do k=n-1,1,-1
         call jacrot(wk(k),wk(k+1),k,n,Q,ldq,R,ldr,k)
      enddo

c     r(1,1:n) += wk(1)*v(1:n)
      call daxpy(n,wk(1),v,1,R(1,1),ldr)

c     R is of upper hessenberg form. Triangularize R.
c      kr argument == k+1 to start applying rotation at column k+1
c      otherwise R(k,k) will be rotated twice and this way it also
c      avoids tiny roundoff errors.

      do k=1,n-1
         call jacrot(R(k,k),R(k+1,k),k,n,Q,ldq,R,ldr,k+1)
      enddo

      return
      end

c-----------------------------------------------------------------------------

      subroutine jacrot(a,b,k,n,Q,ldq,R,ldr,kr)

      double precision a, b
      integer k,n,ldr,ldq,kr
      double precision Q(ldq,*), R(ldr,*)

c-----------------------------------------------------------------------------
c
c     Arguments
c
c     Inout  a       Real             rotate argument
c     Inout  b       Real             rotate argument to rotate to zero
c     In     k       Integer          row/column number for rotation
c     In     n       Integer          order of Q and R
c     Inout  Q       Real(ldq,n)      orthogonal matrix from QR
c     In     ldq     Integer          leading dimension of Q
c     Inout  R       Real(ldr,n)      upper triangular matrix R from QR
c     In     ldr     Integer          leading dimension of R
c     In     u       Real(*)          vector u of size n
c     In     v       Real(*)          vector v of size n
c     In     kr      Integer          start R rotation in column kr
c                                     (should be k or k+1)
c
c-----------------------------------------------------------------------------

      double precision t
      double precision c,s
      double precision Rzero
      parameter(Rzero=0.0d0)

      call dlartg(a,b,c,s,t)
      a = t
      b = Rzero
      call drot(n-kr+1,R(k,kr),ldr,R(k+1,kr),ldr,c,s)
      call drot(n     ,Q(1,k) ,1  ,Q(1,k+1) ,1  ,c,s)

      return
      end
