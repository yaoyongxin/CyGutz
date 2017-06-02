c
c-----------------------------------------------------------------------
c
c     Count the number of elements equal to a specified integer value.
c
      integer function icnteq (n, array, value)
c
      integer    n, value
      integer    array(*)
c
      k = 0
      do 10 i = 1, n
         if (array(i) .eq. value) k = k + 1
   10 continue
      icnteq = k
c
      return
      end
*--------------------------------------------------------------------
*\Documentation
*
*\Name: ICOPY
*
*\Description:
*     ICOPY copies an integer vector lx to an integer vector ly.
*
*\Usage:
*     call icopy ( n, lx, inc, ly, incy )
*
*\Arguments:
*    n        integer (input)
*             On entry, n is the number of elements of lx to be
c             copied to ly.
*
*    lx       integer array (input)
*             On entry, lx is the integer vector to be copied.
*
*    incx     integer (input)
*             On entry, incx is the increment between elements of lx.
*
*    ly       integer array (input)
*             On exit, ly is the integer vector that contains the
*             copy of lx.
*
*    incy     integer (input)
*             On entry, incy is the increment between elements of ly.
*
*\Enddoc
*
*--------------------------------------------------------------------
*
      subroutine icopy( n, lx, incx, ly, incy )
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      integer    incx, incy, n
      integer    lx( 1 ), ly( 1 )
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer           i, ix, iy
*
*     --------------------------
*     First executable statement
*     --------------------------
      if( n.le.0 )
     $   return
      if( incx.eq.1 .and. incy.eq.1 )
     $   go to 20
c
c.....code for unequal increments or equal increments
c     not equal to 1
      ix = 1
      iy = 1
      if( incx.lt.0 )
     $   ix = ( -n+1 )*incx + 1
      if( incy.lt.0 )
     $   iy = ( -n+1 )*incy + 1
      do 10 i = 1, n
         ly( iy ) = lx( ix )
         ix = ix + incx
         iy = iy + incy
   10 continue
      return
c
c.....code for both increments equal to 1
c
   20 continue
      do 30 i = 1, n
         ly( i ) = lx( i )
   30 continue
      return
      end
c
c-----------------------------------------------------------------------
c
c     Only work with increment equal to 1 right now.
c
      subroutine iset (n, value, array, inc)
c
      integer    n, value, inc
      integer    array(*)
c
      do 10 i = 1, n
         array(i) = value
   10 continue
c
      return
      end
      subroutine iswap (n,sx,incx,sy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      integer sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
   50 continue
      return
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dgetv0
c
c\Description: 
c  Generate a random initial residual vector for the Arnoldi process.
c  Force the residual vector to be in the range of the operator OP.  
c
c\Usage:
c  call dgetv0
c     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, 
c       IPNTR, WORKD, IERR )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to dgetv0.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B in the (generalized)
c          eigenvalue problem A*x = lambda*B*x.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  ITRY    Integer.  (INPUT)
c          ITRY counts the number of times that dgetv0 is called.  
c          It should be set to 1 on the initial call to dgetv0.
c
c  INITV   Logical variable.  (INPUT)
c          .TRUE.  => the initial residual vector is given in RESID.
c          .FALSE. => generate a random initial residual vector.
c
c  N       Integer.  (INPUT)
c          Dimension of the problem.
c
c  J       Integer.  (INPUT)
c          Index of the residual vector to be generated, with respect to
c          the Arnoldi process.  J > 1 in case of a "restart".
c
c  V       Double precision N by J array.  (INPUT)
c          The first J-1 columns of V contain the current Arnoldi basis
c          if this is a "restart".
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          Initial residual vector to be generated.  If RESID is 
c          provided, force RESID into the range of the operator OP.
c
c  RNORM   Double precision scalar.  (OUTPUT)
c          B-norm of the generated residual.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c
c  WORKD   Double precision work array of length 2*N.  (REVERSE COMMUNICATION).
c          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
c
c  IERR    Integer.  (OUTPUT)
c          =  0: Normal exit.
c          = -1: Cannot generate a nontrivial restarted residual vector
c                in the range of the operator OP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dlarnv  LAPACK routine for generating a random vector.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\SCCS Information: @(#) 
c FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dgetv0 
     &   ( ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, 
     &     ipntr, workd, ierr )
c 
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      logical    initv
      integer    ido, ierr, itry, j, ldv, n
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           resid(n), v(ldv,j), workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    first, inits, orth
      integer    idist, iseed(4), iter, jj
      Double precision
     &           rnorm0
      save       first, iseed, inits, iter, orth, rnorm0
!$omp threadprivate(first, iseed, inits, iter, orth, rnorm0)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dlarnv, dcopy, dgemv
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2
      external   ddot, dnrm2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, sqrt
c
c     %-----------------%
c     | Data Statements |
c     %-----------------%
c
      data       inits /.true./
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
      NOPX=0; NBX=0
c     %-----------------------------------%
c     | Initialize the seed of the LAPACK |
c     | random number generator           |
c     %-----------------------------------%
c
      if (inits) then
          iseed(1) = 1
          iseed(2) = 3
          iseed(3) = 5
          iseed(4) = 7
          inits = .false.
      end if
c
      if (ido .eq.  0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        %-------------------------------%
c
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
c
c        %-----------------------------------------------------%
c        | Possibly generate a random starting vector in RESID |
c        | Use a LAPACK random number generator used by the    |
c        | matrix generation routines.                         |
c        |    idist = 1: uniform (0,1)  distribution;          |
c        |    idist = 2: uniform (-1,1) distribution;          |
c        |    idist = 3: normal  (0,1)  distribution;          |
c        %-----------------------------------------------------%
c
         if (.not.initv) then
            idist = 2
            call dlarnv (idist, iseed, n, resid)
         end if
c 
c        %----------------------------------------------------------%
c        | Force the starting vector into the range of OP to handle |
c        | the generalized problem when B is possibly (singular).   |
c        %----------------------------------------------------------%
c
         if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call dcopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         end if
      end if
c 
c     %-----------------------------------------%
c     | Back from computing OP*(initial-vector) |
c     %-----------------------------------------%
c
      if (first) go to 20
c
c     %-----------------------------------------------%
c     | Back from computing B*(orthogonalized-vector) |
c     %-----------------------------------------------%
c
      if (orth)  go to 40
c 
c     %------------------------------------------------------%
c     | Starting vector is now in the range of OP; r = OP*r; |
c     | Compute B-norm of starting vector.                   |
c     %------------------------------------------------------%
c
      first = .TRUE.
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call dcopy (n, workd(n+1), 1, resid, 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call dcopy (n, resid, 1, workd, 1)
      end if
c 
   20 continue
c
      first = .FALSE.
      if (bmat .eq. 'G') then
          rnorm0 = ddot (n, resid, 1, workd, 1)
          rnorm0 = sqrt(abs(rnorm0))
      else if (bmat .eq. 'I') then
           rnorm0 = dnrm2(n, resid, 1)
      end if
      rnorm  = rnorm0
c
c     %---------------------------------------------%
c     | Exit if this is the very first Arnoldi step |
c     %---------------------------------------------%
c
      if (j .eq. 1) go to 50
c 
c     %----------------------------------------------------------------
c     | Otherwise need to B-orthogonalize the starting vector against |
c     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
c     | This is the case where an invariant subspace is encountered   |
c     | in the middle of the Arnoldi factorization.                   |
c     |                                                               |
c     |       s = V^{T}*B*r;   r = r - V*s;                           |
c     |                                                               |
c     | Stopping criteria used for iter. ref. is discussed in         |
c     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
c     %---------------------------------------------------------------%
c
      orth = .TRUE.
   30 continue
c
      call dgemv ('T', n, j-1, one, v, ldv, workd, 1, 
     &            zero, workd(n+1), 1)
      call dgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1, 
     &            one, resid, 1)
c 
c     %----------------------------------------------------------%
c     | Compute the B-norm of the orthogonalized starting vector |
c     %----------------------------------------------------------%
c
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call dcopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call dcopy (n, resid, 1, workd, 1)
      end if
c 
   40 continue
c
      if (bmat .eq. 'G') then
         rnorm = ddot (n, resid, 1, workd, 1)
         rnorm = sqrt(abs(rnorm))
      else if (bmat .eq. 'I') then
         rnorm = dnrm2(n, resid, 1)
      end if
c
c     %--------------------------------------%
c     | Check for further orthogonalization. |
c     %--------------------------------------%
c
      if (rnorm .gt. 0.717*rnorm0) go to 50
c 
      iter = iter + 1
      if (iter .le. 5) then
c
c        %-----------------------------------%
c        | Perform iterative refinement step |
c        %-----------------------------------%
c
         rnorm0 = rnorm
         go to 30
      else
c
c        %------------------------------------%
c        | Iterative refinement step "failed" |
c        %------------------------------------%
c
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = zero
         ierr = -1
      end if
c 
   50 continue
c
      ido = 99
c 
 9000 continue
      return
c
c     %---------------%
c     | End of dgetv0 |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dlaqrb
c
c\Description:
c  Compute the eigenvalues and the Schur decomposition of an upper 
c  Hessenberg submatrix in rows and columns ILO to IHI.  Only the
c  last component of the Schur vectors are computed.
c
c  This is mostly a modification of the LAPACK routine dlahqr.
c  
c\Usage:
c  call dlaqrb
c     ( WANTT, N, ILO, IHI, H, LDH, WR, WI,  Z, INFO )
c
c\Arguments
c  WANTT   Logical variable.  (INPUT)
c          = .TRUE. : the full Schur form T is required;
c          = .FALSE.: only eigenvalues are required.
c
c  N       Integer.  (INPUT)
c          The order of the matrix H.  N >= 0.
c
c  ILO     Integer.  (INPUT)
c  IHI     Integer.  (INPUT)
c          It is assumed that H is already upper quasi-triangular in
c          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
c          ILO = 1). SLAQRB works primarily with the Hessenberg
c          submatrix in rows and columns ILO to IHI, but applies
c          transformations to all of H if WANTT is .TRUE..
c          1 <= ILO <= max(1,IHI); IHI <= N.
c
c  H       Double precision array, dimension (LDH,N).  (INPUT/OUTPUT)
c          On entry, the upper Hessenberg matrix H.
c          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
c          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
c          standard form. If WANTT is .FALSE., the contents of H are
c          unspecified on exit.
c
c  LDH     Integer.  (INPUT)
c          The leading dimension of the array H. LDH >= max(1,N).
c
c  WR      Double precision array, dimension (N).  (OUTPUT)
c  WI      Double precision array, dimension (N).  (OUTPUT)
c          The real and imaginary parts, respectively, of the computed
c          eigenvalues ILO to IHI are stored in the corresponding
c          elements of WR and WI. If two eigenvalues are computed as a
c          complex conjugate pair, they are stored in consecutive
c          elements of WR and WI, say the i-th and (i+1)th, with
c          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
c          eigenvalues are stored in the same order as on the diagonal
c          of the Schur form returned in H, with WR(i) = H(i,i), and, if
c          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
c          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
c
c  Z       Double precision array, dimension (N).  (OUTPUT)
c          On exit Z contains the last components of the Schur vectors.
c
c  INFO    Integer.  (OUPUT)
c          = 0: successful exit
c          > 0: SLAQRB failed to compute all the eigenvalues ILO to IHI
c               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
c               elements i+1:ihi of WR and WI contain those eigenvalues
c               which have been successfully computed.
c
c\Remarks
c  1. None.
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dlabad  LAPACK routine that computes machine constants.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dlanv2  LAPACK routine that computes the Schur factorization of
c             2 by 2 nonsymmetric matrix in standard form.
c     dlarfg  LAPACK Householder reflection construction routine.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     drot    Level 1 BLAS that applies a rotation to a 2 by 2 matrix.
!
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c               Modified from the LAPACK routine dlahqr so that only the
c               last component of the Schur vectors are computed.
c
c\SCCS Information: @(#) 
c FILE: laqrb.F   SID: 2.2   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dlaqrb ( wantt, n, ilo, ihi, h, ldh, wr, wi,
     &                    z, info )
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      logical    wantt
      integer    ihi, ilo, info, ldh, n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h( ldh, * ), wi( * ), wr( * ), z( * )
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           zero, one, dat1, dat2
      parameter (zero = 0.0D+0, one = 1.0D+0, dat1 = 7.5D-1, 
     &           dat2 = -4.375D-1)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, i1, i2, itn, its, j, k, l, m, nh, nr
      Double precision
     &           cs, h00, h10, h11, h12, h21, h22, h33, h33s,
     &           h43h34, h44, h44s, ovfl, s, smlnum, sn, sum,
     &           t1, t2, t3, tst1, ulp, unfl, v1, v2, v3
      Double precision
     &           v( 3 ), work( 1 )
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch, dlanhs
      external   dlamch, dlanhs
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dlabad, dlanv2, dlarfg, drot
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      info = 0
c
c     %--------------------------%
c     | Quick return if possible |
c     %--------------------------%
c
      if( n.eq.0 )
     &   return
      if( ilo.eq.ihi ) then
         wr( ilo ) = h( ilo, ilo )
         wi( ilo ) = zero
         return
      end if
c 
c     %---------------------------------------------%
c     | Initialize the vector of last components of |
c     | the Schur vectors for accumulation.         |
c     %---------------------------------------------%
c
      do 5 j = 1, n-1
         z(j) = zero
  5   continue 
      z(n) = one
c 
      nh = ihi - ilo + 1
c
c     %-------------------------------------------------------------%
c     | Set machine-dependent constants for the stopping criterion. |
c     | If norm(H) <= sqrt(OVFL), overflow should not occur.        |
c     %-------------------------------------------------------------%
c
      unfl = dlamch( 'safe minimum' )
      ovfl = one / unfl
      call dlabad( unfl, ovfl )
      ulp = dlamch( 'precision' )
      smlnum = unfl*( nh / ulp )
c
c     %---------------------------------------------------------------%
c     | I1 and I2 are the indices of the first row and last column    |
c     | of H to which transformations must be applied. If eigenvalues |
c     | only are computed, I1 and I2 are set inside the main loop.    |
c     | Zero out H(J+2,J) = ZERO for J=1:N if WANTT = .TRUE.          |
c     | else H(J+2,J) for J=ILO:IHI-ILO-1 if WANTT = .FALSE.          |
c     %---------------------------------------------------------------%
c
      if( wantt ) then
         i1 = 1
         i2 = n
         do 8 i=1,i2-2
            h(i1+i+1,i) = zero
 8       continue
      else
         do 9 i=1, ihi-ilo-1
            h(ilo+i+1,ilo+i-1) = zero
 9       continue
      end if
c 
c     %---------------------------------------------------%
c     | ITN is the total number of QR iterations allowed. |
c     %---------------------------------------------------%
c
      itn = 30*nh
c 
c     ------------------------------------------------------------------
c     The main loop begins here. I is the loop index and decreases from
c     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
c     with the active submatrix in rows and columns L to I.
c     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
c     H(L,L-1) is negligible so that the matrix splits.
c     ------------------------------------------------------------------
c 
      i = ihi
   10 continue
      l = ilo
      if( i.lt.ilo )
     &   go to 150
!
c     %--------------------------------------------------------------%
c     | Perform QR iterations on rows and columns ILO to I until a   |
c     | submatrix of order 1 or 2 splits off at the bottom because a |
c     | subdiagonal element has become negligible.                   |
c     %--------------------------------------------------------------%
!
      do 130 its = 0, itn
c
c        %----------------------------------------------%
c        | Look for a single small subdiagonal element. |
c        %----------------------------------------------%
c
         do 20 k = i, l + 1, -1
            tst1 = abs( h( k-1, k-1 ) ) + abs( h( k, k ) )
            if( tst1.eq.zero )
     &         tst1 = dlanhs( '1', i-l+1, h( l, l ), ldh, work )
            if( abs( h( k, k-1 ) ).le.max( ulp*tst1, smlnum ) )
     &         go to 30
   20    continue
   30    continue
         l = k
         if( l.gt.ilo ) then
c
c           %------------------------%
c           | H(L,L-1) is negligible |
c           %------------------------%
c
            h( l, l-1 ) = zero
         end if
c
c        %-------------------------------------------------------------%
c        | Exit from loop if a submatrix of order 1 or 2 has split off |
c        %-------------------------------------------------------------%
c
         if( l.ge.i-1 )
     &      go to 140
c
c        %---------------------------------------------------------%
c        | Now the active submatrix is in rows and columns L to I. |
c        | If eigenvalues only are being computed, only the active |
c        | submatrix need be transformed.                          |
c        %---------------------------------------------------------%
c
         if( .not.wantt ) then
            i1 = l
            i2 = i
         end if
c 
         if( its.eq.10 .or. its.eq.20 ) then
c
c           %-------------------%
c           | Exceptional shift |
c           %-------------------%
c
            s = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
            h44 = dat1*s
            h33 = h44
            h43h34 = dat2*s*s
c
         else
c
c           %-----------------------------------------%
c           | Prepare to use Wilkinson's double shift |
c           %-----------------------------------------%
c
            h44 = h( i, i )
            h33 = h( i-1, i-1 )
            h43h34 = h( i, i-1 )*h( i-1, i )
         end if
c
c        %-----------------------------------------------------%
c        | Look for two consecutive small subdiagonal elements |
c        %-----------------------------------------------------%
c
         do 40 m = i - 2, l, -1
c
c           %---------------------------------------------------------%
c           | Determine the effect of starting the double-shift QR    |
c           | iteration at row M, and see if this would make H(M,M-1) |
c           | negligible.                                             |
c           %---------------------------------------------------------%
c
            h11 = h( m, m )
            h22 = h( m+1, m+1 )
            h21 = h( m+1, m )
            h12 = h( m, m+1 )
            h44s = h44 - h11
            h33s = h33 - h11
            v1 = ( h33s*h44s-h43h34 ) / h21 + h12
            v2 = h22 - h11 - h33s - h44s
            v3 = h( m+2, m+1 )
            s = abs( v1 ) + abs( v2 ) + abs( v3 )
            v1 = v1 / s
            v2 = v2 / s
            v3 = v3 / s
            v( 1 ) = v1
            v( 2 ) = v2
            v( 3 ) = v3
            if( m.eq.l )
     &         go to 50
            h00 = h( m-1, m-1 )
            h10 = h( m, m-1 )
            tst1 = abs( v1 )*( abs( h00 )+abs( h11 )+abs( h22 ) )
            if( abs( h10 )*( abs( v2 )+abs( v3 ) ).le.ulp*tst1 )
     &         go to 50
   40    continue
   50    continue
c
c        %----------------------%
c        | Double-shift QR step |
c        %----------------------%
c
         do 120 k = m, i - 1
c 
c           ------------------------------------------------------------
c           The first iteration of this loop determines a reflection G
c           from the vector V and applies it from left and right to H,
c           thus creating a nonzero bulge below the subdiagonal.
c
c           Each subsequent iteration determines a reflection G to
c           restore the Hessenberg form in the (K-1)th column, and thus
c           chases the bulge one step toward the bottom of the active
c           submatrix. NR is the order of G.
c           ------------------------------------------------------------
c 
            nr = min( 3, i-k+1 )
            if( k.gt.m )
     &         call dcopy( nr, h( k, k-1 ), 1, v, 1 )
            call dlarfg( nr, v( 1 ), v( 2 ), 1, t1 )
            if( k.gt.m ) then
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = zero
               if( k.lt.i-1 )
     &            h( k+2, k-1 ) = zero
            else if( m.gt.l ) then
               h( k, k-1 ) = -h( k, k-1 )
            end if
            v2 = v( 2 )
            t2 = t1*v2
            if( nr.eq.3 ) then
               v3 = v( 3 )
               t3 = t1*v3
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
               do 60 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j ) + v3*h( k+2, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
                  h( k+2, j ) = h( k+2, j ) - sum*t3
   60          continue
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
               do 70 j = i1, min( k+3, i )
                  sum = h( j, k ) + v2*h( j, k+1 ) + v3*h( j, k+2 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
                  h( j, k+2 ) = h( j, k+2 ) - sum*t3
   70          continue
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
               sum      = z( k ) + v2*z( k+1 ) + v3*z( k+2 )
               z( k )   = z( k ) - sum*t1
               z( k+1 ) = z( k+1 ) - sum*t2
               z( k+2 ) = z( k+2 ) - sum*t3
!
            else if( nr.eq.2 ) then
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
               do 90 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
   90          continue
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
               do 100 j = i1, i
                  sum = h( j, k ) + v2*h( j, k+1 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
  100          continue
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
               sum      = z( k ) + v2*z( k+1 )
               z( k )   = z( k ) - sum*t1
               z( k+1 ) = z( k+1 ) - sum*t2
            end if
  120    continue
!
  130 continue
c
c     %-------------------------------------------------------%
c     | Failure to converge in remaining number of iterations |
c     %-------------------------------------------------------%
c
      info = i
      return
!
  140 continue
!
      if( l.eq.i ) then
c
c        %------------------------------------------------------%
c        | H(I,I-1) is negligible: one eigenvalue has converged |
c        %------------------------------------------------------%
c
         wr( i ) = h( i, i )
         wi( i ) = zero
!
      else if( l.eq.i-1 ) then
c
c        %--------------------------------------------------------%
c        | H(I-1,I-2) is negligible;                              |
c        | a pair of eigenvalues have converged.                  |
c        |                                                        |
c        | Transform the 2-by-2 submatrix to standard Schur form, |
c        | and compute and store the eigenvalues.                 |
c        %--------------------------------------------------------%
c
         call dlanv2( h( i-1, i-1 ), h( i-1, i ), h( i, i-1 ),
     &                h( i, i ), wr( i-1 ), wi( i-1 ), wr( i ), wi( i ),
     &                cs, sn )
!
         if( wantt ) then
c
c           %-----------------------------------------------------%
c           | Apply the transformation to the rest of H and to Z, |
c           | as required.                                        |
c           %-----------------------------------------------------%
c
            if( i2.gt.i )
     &         call drot( i2-i, h( i-1, i+1 ), ldh, h( i, i+1 ), ldh,
     &                    cs, sn )
            call drot( i-i1-1, h( i1, i-1 ), 1, h( i1, i ), 1, cs, sn )
            sum      = cs*z( i-1 ) + sn*z( i )
            z( i )   = cs*z( i )   - sn*z( i-1 )
            z( i-1 ) = sum
         end if
      end if
c
c     %---------------------------------------------------------%
c     | Decrement number of remaining iterations, and return to |
c     | start of the main loop with new value of I.             |
c     %---------------------------------------------------------%
c
      itn = itn - its
      i = l - 1
      go to 10
!
  150 continue
      return
c
c     %---------------%
c     | End of dlaqrb |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnaitr
c
c\Description: 
c  Reverse communication interface for applying NP additional steps to 
c  a K step nonsymmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in dnaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call dnaitr
c     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and do not need to be
c          recompute in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.  See dnaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current size of V and H.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  NB      Integer.  (INPUT)
c          Blocksize to be used in the recurrence.          
c          Only work for NB = 1 right now.  The goal is to have a 
c          program that implement both the block and non-block method.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Double precision scalar.  (INPUT/OUTPUT)
c          B-norm of the starting residual on input.
c          B-norm of the updated residual r_{k+p} on output.
c
c  V       Double precision N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K 
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
c          H is used to store the generated upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not 
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On input, WORKD(1:N) = B*RESID and is used to save some 
c          computation at the first step.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of the spanning invariant subspace of OP found.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dgetv0  ARPACK routine to generate the initial vector.
c     dlabad  LAPACK routine that computes machine constants.
c     dlamch  LAPACK routine that determines machine constants.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c 
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: naitr.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c  
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        ( At present tol is zero )
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in dnaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        H(:,j) = w_{j};
c        H(j,j-1) = rnorm
c        rnorm = || r_(j) ||
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};   
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf 
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaitr
     &   (ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh, 
     &    ipntr, workd, info)
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, nb, np
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           h(ldh,k+np), resid(n), v(ldv,k+np), workd(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    first, orth1, orth2, rstart, step3, step4
      integer    ierr, i, infol, ipj, irj, ivj, iter, itry, j, 
     &           jj
      Double precision
     &           betaj, ovfl, temp1, rnorm1, smlnum, tst1, ulp, unfl, 
     &           wnorm
      save       first, orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, ovfl,
     &           betaj, rnorm1, smlnum, ulp, unfl, wnorm
!$omp threadprivate(first, orth1, orth2, rstart, step3, step4)
!$omp threadprivate(ierr, ipj, irj, ivj, iter, itry, j, ovfl)
!$omp threadprivate(betaj, rnorm1, smlnum, ulp, unfl, wnorm)
c
c     %-----------------------%
c     | Local Array Arguments | 
c     %-----------------------%
c
      Double precision
     &           xtemp(2)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dgemv, dgetv0, dlabad
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlanhs, dlamch
      external   ddot, dnrm2, dlanhs, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, sqrt
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
      data      first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      NOPX=0; NBX=0; NRORTH=0
      if (first) then
c
c        %-----------------------------------------%
c        | Set machine-dependent constants for the |
c        | the splitting and deflation criterion.  |
c        | If norm(H) <= sqrt(OVFL),               |
c        | overflow should not occur.              |
c        | REFERENCE: LAPACK subroutine dlahqr     |
c        %-----------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = one / unfl
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        %-------------------------------%
c
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
         j      = k + 1
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
c 
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true. when ....                        |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         dgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %-----------------------------%
c     | Else this is the first step |
c     %-----------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
!
 1000 continue
c
c        %---------------------------------------------------%
c        | STEP 1: Check if the B norm of j-th residual      |
c        | vector is zero. Equivalent to determing whether   |
c        | an exact j-step Arnoldi factorization is present. |
c        %---------------------------------------------------%
c
         betaj = rnorm
         if (rnorm .gt. zero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted.       |
c           %---------------------------------------------%
c 
            betaj  = zero
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
            call dgetv0 (ido, bmat, itry, .false., n, j, v, ldv, 
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               ido = 99
               go to 9000
            end if
c 
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call dcopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. unfl) then
             temp1 = one / rnorm
             call dscal (n, temp1, v(1,j), 1)
             call dscal (n, temp1, workd(ipj), 1)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
             call dlascl ('General', i, i, rnorm, one, n, 1, 
     &                    v(1,j), n, infol)
             call dlascl ('General', i, i, rnorm, one, n, 1, 
     &                    workd(ipj), n, infol)
         end if
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         nopx  = nopx + 1
         call dcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c 
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c 
         go to 9000 
   50    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
c        | if step3 = .true.                |
c        %----------------------------------%
c
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call dcopy (n, workd(irj), 1, resid, 1)
c 
c        %---------------------------------------%
c        | STEP 4:  Finish extending the Arnoldi |
c        |          factorization to length j.   |
c        %---------------------------------------%
c
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   60    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
c        | if step4 = .true.                |
c        %----------------------------------%
c
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
         if (bmat .eq. 'G') then  
             wnorm = ddot (n, resid, 1, workd(ipj), 1)
             wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c 
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, h(1,j), 1)
c
c        %--------------------------------------%
c        | Orthogonalize r_{j} against V_{j}.   |
c        | RESID contains OP*v_{j}. See STEP 3. | 
c        %--------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, h(1,j), 1,
     &               one, resid, 1)
c
         if (j .gt. 1) h(j,j-1) = betaj
c
         orth1 = .true.
c
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if 
   70    continue
c 
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
c 
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        | The following test determines whether the sine of the     |
c        | angle between  OP*x and the computed residual is less     |
c        | than or equal to 0.717.                                   |
c        %-----------------------------------------------------------%
c
         if (rnorm .gt. 0.717*wnorm) go to 100
         iter  = 0
         nrorth = nrorth + 1
c 
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c 
   80    continue
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1, 
     &               zero, workd(irj), 1)
c
c        %---------------------------------------------%
c        | Compute the correction to the residual:     |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
c        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
c        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
c        %---------------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, workd(irj), 1, 
     &               one, resid, 1)
         call daxpy (j, one, workd(irj), 1, h(1,j), 1)
c 
         orth2 = .true.
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if 
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c 
         if (bmat .eq. 'G') then         
             rnorm1 = ddot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if (rnorm1 .gt. 0.717*rnorm) then
c
c           %---------------------------------------%
c           | No need for further refinement.       |
c           | The cosine of the angle between the   |
c           | corrected residual vector and the old |
c           | residual vector is greater than 0.717 |
c           | In other words the corrected residual |
c           | and the old residual vector share an  |
c           | angle of less than arcCOS(0.717)      |
c           %---------------------------------------%
c
            rnorm = rnorm1
c 
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required.    |
c           %-------------------------------------------%
c
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
c 
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c 
  100    continue
c 
         rstart = .false.
         orth2  = .false.
c 
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            ido = 99
            do 110 i = max(1,k), k+np-1
c     
c              %--------------------------------------------%
c              | Check for splitting and deflation.         |
c              | Use a standard test as in the QR algorithm |
c              | REFERENCE: LAPACK subroutine dlahqr        |
c              %--------------------------------------------%
c     
               tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
               if( tst1.eq.zero )
     &              tst1 = dlanhs( '1', k+np, h, ldh, workd(n+1) )
               if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) 
     &              h(i+1,i) = zero
 110        continue
c     
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c 
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of dnaitr |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP implicit shifts resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix which is the product of rotations
c  and reflections resulting from the NP bulge chage sweeps.
c  The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call dnapps
c     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ, 
c       WORKL, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. size of matrix A.
c
c  KEV     Integer.  (INPUT/OUTPUT)
c          KEV+NP is the size of the input matrix H.
c          KEV is the size of the updated matrix HNEW.  KEV is only 
c          updated on ouput when fewer than NP shifts are applied in
c          order to keep the conjugate pair together.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFTR, Double precision array of length NP.  (INPUT)
c  SHIFTI  Real and imaginary part of the shifts to be applied.
c          Upon, entry to dnapps, the shifts must be sorted so that the 
c          conjugate pairs are in consecutive locations.
c
c  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, V contains the current KEV+NP Arnoldi vectors.
c          On OUTPUT, V contains the updated KEV Arnoldi vectors
c          in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, H contains the current KEV+NP by KEV+NP upper 
c          Hessenber matrix of the Arnoldi factorization.
c          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
c          matrix in the KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT, RESID contains the the residual vector r_{k+p}.
c          On OUTPUT, RESID is the update residual vector rnew_{k} 
c          in the first KEV locations.
c
c  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations and reflections
c          during the bulge chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length (KEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     dlabad  LAPACK routine that computes machine constants.
c     dlacpy  LAPACK matrix copy routine.
c     dlamch  LAPACK routine that determines machine constants. 
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlarf   LAPACK routine that applies Householder reflection to
c             a matrix.
c     dlarfg  LAPACK Householder reflection construction routine.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dscal   Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: napps.F   SID: 2.4   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the sublocks of
c     the Hessenberg matrix H and not just to the submatrix that it
c     comes from. Deflation as in LAPACK routine dlahqr (QR algorithm
c     for upper Hessenberg matrices ) is used.
c     The subdiagonals of H are enforced to be non-negative.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnapps
     &   ( n, kev, np, shiftr, shifti, v, ldv, h, ldh, resid, q, ldq, 
     &     workl, workd )
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h(ldh,kev+np), resid(n), shifti(np), shiftr(np), 
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, iend, ir, istart, j, jj, kplusp, nr
      logical    cconj, first
      Double precision
     &           c, f, g, h11, h12, h21, h22, h32, ovfl, r, s, sigmai, 
     &           sigmar, smlnum, ulp, unfl, u(3), t, tau, tst1
      save       first, ovfl, smlnum, ulp, unfl 
!$omp threadprivate(first, ovfl, smlnum, ulp, unfl)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dlacpy, dlarfg, dlarf,
     &           dlaset, dlabad, dlartg
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch, dlanhs, dlapy2
      external   dlamch, dlanhs, dlapy2
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs, max, min
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine dlahqr           |
c        %-----------------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = one / unfl
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------%
c
      kplusp = kev + np 
c 
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
      call dlaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
      cconj = .false.
      do 110 jj = 1, np
         sigmar = shiftr(jj)
         sigmai = shifti(jj)
c
c        %-------------------------------------------------%
c        | The following set of conditionals is necessary  |
c        | in order that complex conjugate pairs of shifts |
c        | are applied together or not at all.             |
c        %-------------------------------------------------%
c
         if ( cconj ) then
c
c           %-----------------------------------------%
c           | cconj = .true. means the previous shift |
c           | had non-zero imaginary part.            |
c           %-----------------------------------------%
c
            cconj = .false.
            go to 110
         else if ( jj .lt. np .and. abs( sigmai ) .gt. zero ) then
c
c           %------------------------------------%
c           | Start of a complex conjugate pair. |
c           %------------------------------------%
c
            cconj = .true.
         else if ( jj .eq. np .and. abs( sigmai ) .gt. zero ) then
c
c           %----------------------------------------------%
c           | The last shift has a nonzero imaginary part. |
c           | Don't apply it; thus the order of the        |
c           | compressed H is order KEV+1 since only np-1  |
c           | were applied.                                |
c           %----------------------------------------------%
c
            kev = kev + 1
            go to 110
         end if
         istart = 1
   20    continue
c
c        %--------------------------------------------------%
c        | if sigmai = 0 then                               |
c        |    Apply the jj-th shift ...                     |
c        | else                                             |
c        |    Apply the jj-th and (jj+1)-th together ...    |
c        |    (Note that jj < np at this point in the code) |
c        | end                                              |
c        | to the current block of H. The next do loop      |
c        | determines the current block ;                   |
c        %--------------------------------------------------%
c
         do 30 i = istart, kplusp-1
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine dlahqr    |
c           %----------------------------------------%
c
            tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
            if( tst1.eq.zero )
     &         tst1 = dlanhs( '1', kplusp-jj+1, h, ldh, workl )
            if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) then
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        %------------------------------------------------%
c
         if ( istart .eq. iend ) go to 100
c
c        %------------------------------------------------------%
c        | If istart + 1 = iend then no reason to apply a       |
c        | complex conjugate pair of shifts on a 2 by 2 matrix. |
c        %------------------------------------------------------%
c
         if ( istart + 1 .eq. iend .and. abs( sigmai ) .gt. zero ) 
     &      go to 100
c
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         if ( abs( sigmai ) .le. zero ) then
c
c           %---------------------------------------------%
c           | Real-valued shift ==> apply single shift QR |
c           %---------------------------------------------%
c
            f = h11 - sigmar
            g = h21
c 
            do 80 i = istart, iend-1
c
c              %-----------------------------------------------------%
c              | Contruct the plane rotation G to zero out the bulge |
c              %-----------------------------------------------------%
c
               call dlartg (f, g, c, s, r)
               if (i .gt. istart) then
c
c                 %-------------------------------------------%
c                 | The following ensures that h(1:iend-1,1), |
c                 | the first iend-2 off diagonal of elements |
c                 | H, remain non negative.                   |
c                 %-------------------------------------------%
c
                  if (r .lt. zero) then
                     r = -r
                     c = -c
                     s = -s
                  end if
                  h(i,i-1) = r
                  h(i+1,i-1) = zero
               end if
c
c              %---------------------------------------------%
c              | Apply rotation to the left of H;  H <- G'*H |
c              %---------------------------------------------%
c
               do 50 j = i, kplusp
                  t        =  c*h(i,j) + s*h(i+1,j)
                  h(i+1,j) = -s*h(i,j) + c*h(i+1,j)
                  h(i,j)   = t   
   50          continue
c
c              %---------------------------------------------%
c              | Apply rotation to the right of H;  H <- H*G |
c              %---------------------------------------------%
c
               do 60 j = 1, min(i+2,iend)
                  t        =  c*h(j,i) + s*h(j,i+1)
                  h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
                  h(j,i)   = t   
   60          continue
c
c              %----------------------------------------------------%
c              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c              %----------------------------------------------------%
c
               do 70 j = 1, min( i+jj, kplusp ) 
                  t        =   c*q(j,i) + s*q(j,i+1)
                  q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                  q(j,i)   = t   
   70          continue
c
c              %---------------------------%
c              | Prepare for next rotation |
c              %---------------------------%
c
               if (i .lt. iend-1) then
                  f = h(i+1,i)
                  g = h(i+2,i)
               end if
   80       continue
c
c           %-----------------------------------%
c           | Finished applying the real shift. |
c           %-----------------------------------%
c 
         else
c
c           %----------------------------------------------------%
c           | Complex conjugate shifts ==> apply double shift QR |
c           %----------------------------------------------------%
c
            h12 = h(istart,istart+1)
            h22 = h(istart+1,istart+1)
            h32 = h(istart+2,istart+1)
c
c           %---------------------------------------------------------%
c           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
c           %---------------------------------------------------------%
c
            s    = 2.0*sigmar
            t = dlapy2 ( sigmar, sigmai ) 
            u(1) = ( h11 * (h11 - s) + t * t ) / h21 + h12
            u(2) = h11 + h22 - s 
            u(3) = h32
c
            do 90 i = istart, iend-1
c
               nr = min ( 3, iend-i+1 )
c
c              %-----------------------------------------------------%
c              | Construct Householder reflector G to zero out u(1). |
c              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
c              %-----------------------------------------------------%
c
               call dlarfg ( nr, u(1), u(2), 1, tau )
c
               if (i .gt. istart) then
                  h(i,i-1)   = u(1)
                  h(i+1,i-1) = zero
                  if (i .lt. iend-1) h(i+2,i-1) = zero
               end if
               u(1) = one
c
c              %--------------------------------------%
c              | Apply the reflector to the left of H |
c              %--------------------------------------%
c
               call dlarf ('Left', nr, kplusp-i+1, u, 1, tau,
     &                     h(i,i), ldh, workl)
c
c              %---------------------------------------%
c              | Apply the reflector to the right of H |
c              %---------------------------------------%
c
               ir = min ( i+3, iend )
               call dlarf ('Right', ir, nr, u, 1, tau,
     &                     h(1,i), ldh, workl)
c
c              %-----------------------------------------------------%
c              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
c              %-----------------------------------------------------%
c
               call dlarf ('Right', kplusp, nr, u, 1, tau, 
     &                     q(1,i), ldq, workl)
c
c              %----------------------------%
c              | Prepare for next reflector |
c              %----------------------------%
c
               if (i .lt. iend-1) then
                  u(1) = h(i+1,i)
                  u(2) = h(i+2,i)
                  if (i .lt. iend-2) u(3) = h(i+3,i)
               end if
c
   90       continue
c
c           %--------------------------------------------%
c           | Finished applying a complex pair of shifts |
c           | to the current block                       |
c           %--------------------------------------------%
c 
         end if
c
  100    continue
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
         istart = iend + 1
         if (iend .lt. kplusp) go to 20
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
  110 continue
c
c     %--------------------------------------------------%
c     | Perform a similarity transformation that makes   |
c     | sure that H will have non negative sub diagonals |
c     %--------------------------------------------------%
c
      do 120 j=1,kev
         if ( h(j+1,j) .lt. zero ) then
              call dscal( kplusp-j+1, -one, h(j+1,j), ldh )
              call dscal( min(j+2, kplusp), -one, h(1,j+1), 1 )
              call dscal( min(j+np+1,kplusp), -one, q(1,j+1), 1 )
         end if
 120  continue
c
      do 130 i = 1, kev
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine dlahqr        |
c        %--------------------------------------------%
c
         tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
         if( tst1.eq.zero )
     &       tst1 = dlanhs( '1', kev, h, ldh, workl )
         if( h( i+1,i ) .le. max( ulp*tst1, smlnum ) ) 
     &       h(i+1,i) = zero
 130  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
      if (h(kev+1,kev) .gt. zero)
     &    call dgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero, 
     &                workd(n+1), 1)
c 
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
      do 140 i = 1, kev
         call dgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call dcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call dlacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
c 
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
      if (h(kev+1,kev) .gt. zero)
     &   call dcopy (n, workd(n+1), 1, v(1,kev+1), 1)
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call dscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,kev) .gt. zero)
     &   call daxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
c
 9000 continue
      return
c
c     %---------------%
c     | End of dnapps |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: dnaup2
c
c\Description: 
c  Intermediate level interface called by dnaupd.
c
c\Usage:
c  call dnaup2
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS, 
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dnaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in dnaupd.
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during 
c          each Arnoldi iteration.  
c          If ISHIFT=1, NP is adjusted dynamically at each iteration 
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector 
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV for two reasons. The first, is
c          to keep complex conjugate pairs of "wanted" Ritz values 
c          together. The 2nd, is that a leading block of the current
c          upper Hessenberg matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number 
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV 
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT)
c  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
c          imaginary) part of the computed Ritz values of OP.
c
c  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to 
c          the computed Ritz values.
c          
c  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length at least 
c          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c          On exit, the last 3*(NEV+NP) locations of WORKL contain
c          the Ritz values (real,imaginary) and associated Ritz
c          estimates of the current Hessenberg matrix.  They are
c          listed in the same order as returned from dneigh.
c
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
c          of WORKL are used in reverse communication to hold the user 
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in DNAUPD.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.  
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dgetv0  ARPACK initial vector generation routine. 
c     dnaitr  ARPACK Arnoldi factorization routine.
c     dnapps  ARPACK application of implicit shifts routine.
c     dnconv  ARPACK convergence of Ritz values routine.
c     dneigh  ARPACK compute Ritz values and error bounds routine.
c     dngets  ARPACK reorder Ritz values and error bounds routine.
c     dsortc  ARPACK sorting routine.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dswap   Level 1 BLAS that swaps two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c 
c\SCCS Information: @(#) 
c FILE: naup2.F   SID: 2.8   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, 
     &     ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds, 
     &     q, ldq, workl, ipntr, workd, info )
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Double precision
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n),
     &           ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np), 
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  wprime*2
      logical    cnorm , getv0, initv, update, ushift
      integer    ierr  , iter , j    , kplusp, nconv, 
     &           nevbef, nev0 , np0  , nptemp, numcnv
      Double precision
     &           rnorm , temp , eps23
      save       cnorm , getv0, initv, update, ushift,
     &           rnorm , iter , eps23, kplusp, nconv , 
     &           nevbef, nev0 , np0  , numcnv
!$omp threadprivate(cnorm , getv0, initv, update, ushift)
!$omp threadprivate(rnorm , iter , eps23, kplusp, nconv)
!$omp threadprivate(nevbef, nev0 , np0  , numcnv)
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(4)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy , dgetv0, dnaitr, dnconv, dneigh, 
     &           dngets, dnapps
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlapy2, dlamch
      external   ddot, dnrm2, dlapy2, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min, max, abs, sqrt
c
      NBX=0
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c 
c        %-------------------------------------%
c        | Get the machine dependent constant. |
c        %-------------------------------------%
c
         eps23 = dlamch('Epsilon-Machine')
         eps23 = eps23**(2.0D+0 / 3.0D+0)
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 0
c 
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c 
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call dgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. zero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c 
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call dnaitr (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, 
     &             h, ldh, ipntr, workd, info)
c 
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c 
 1000 continue
c
         iter = iter + 1
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine dnapps.                  |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call dnaitr (ido  , bmat, n  , nev, np , mode , resid, 
     &                rnorm, v   , ldv, h  , ldh, ipntr, workd,
     &                info)
c 
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         call dneigh (rnorm, kplusp, h, ldh, ritzr, ritzi, bounds,
     &                q, ldq, workl, ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained from dneigh.                       |
c        %----------------------------------------------------%
c
         call dcopy(kplusp, ritzr, 1, workl(kplusp**2+1), 1)
         call dcopy(kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZR,   |
c        | RITZI and BOUNDS respectively. The variables NEV  |
c        | and NP may be updated if the NEV-th wanted Ritz   |
c        | value has a non zero imaginary part. In this case |
c        | NEV is increased by one and NP decreased by one.  |
c        | NOTE: The last two arguments of dngets are no     |
c        | longer used as of version 2.1.                    |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         numcnv = nev
         call dngets (ishift, which, nev, np, ritzr, ritzi, 
     &                bounds, workl, workl(np+1))
         if (nev .eq. nev0+1) numcnv = nev0+1
c 
c        %-------------------%
c        | Convergence test. | 
c        %-------------------%
c
         call dcopy (nev, bounds(np+1), 1, workl(2*np+1), 1)
         call dnconv (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1), 
     &        tol, nconv)
c 
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c     
         if ( (nconv .ge. NUMCNV_FIX) .or.
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to _neupd if needed               |
c           %------------------------------------------%
!
            h(3,1) = rnorm
c
c           %----------------------------------------------%
c           | To be consistent with dngets, we first do a  |
c           | pre-processing sort in order to keep complex |
c           | conjugate pairs together.  This is similar   |
c           | to the pre-processing sort used in dngets    |
c           | except that the sort is done in the opposite |
c           | order.                                       |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SR'
            if (which .eq. 'SM') wprime = 'LR'
            if (which .eq. 'LR') wprime = 'SM'
            if (which .eq. 'SR') wprime = 'LM'
            if (which .eq. 'LI') wprime = 'SM'
            if (which .eq. 'SI') wprime = 'LM'
c
            call dsortc (wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %----------------------------------------------%
c           | Now sort Ritz values so that converged Ritz  |
c           | values appear within the first NEV locations |
c           | of ritzr, ritzi and bounds, and the most     |
c           | desired one appears at the front.            |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call dsortc(wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, numcnv
                temp = max(eps23,dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | esitmates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LR'
            call dsortc(wprime, .true., numcnv, bounds, ritzr, ritzi)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, numcnv
                temp = max(eps23, dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %------------------------------------------------%
c           | Sort the converged Ritz values again so that   |
c           | the "threshold" value appears at the front of  |
c           | ritzr, ritzi and bound.                        |
c           %------------------------------------------------%
c
            call dsortc(which, .true., nconv, ritzr, ritzi, bounds)
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. numcnv) info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. numcnv) .and. (ishift .eq. 1) ) then
c     
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call dngets (ishift, which, nev, np, ritzr, ritzi, 
     &              bounds, workl, workl(np+1))
c
         end if              
c     
         if (ishift .eq. 0) then
c
c           %-------------------------------------------------------%
c           | User specified shifts: reverse comminucation to       |
c           | compute the shifts. They are returned in the first    |
c           | 2*NP locations of WORKL.                              |
c           %-------------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
c 
   50    continue
c
c        %------------------------------------%
c        | Back from reverse communication;   |
c        | User specified shifts are returned |
c        | in WORKL(1:2*NP)                   |
c        %------------------------------------%
c
         ushift = .false.
c
         if ( ishift .eq. 0 ) then
c 
c            %----------------------------------%
c            | Move the NP shifts from WORKL to |
c            | RITZR, RITZI to free up WORKL    |
c            | for non-exact shift case.        |
c            %----------------------------------%
c
             call dcopy (np, workl,       1, ritzr, 1)
             call dcopy (np, workl(np+1), 1, ritzi, 1)
         end if
c
c
c        %---------------------------------------------------------%
c        | Apply the NP implicit shifts by QR bulge chasing.       |
c        | Each shift is applied to the whole upper Hessenberg     |
c        | matrix H.                                               |
c        | The first 2*N locations of WORKD are used as workspace. |
c        %---------------------------------------------------------%
c
         call dnapps (n, nev, np, ritzr, ritzi, v, ldv, 
     &                h, ldh, resid, q, ldq, workl, workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dnaitr.  |
c        %---------------------------------------------%
c
         cnorm = .true.
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd, 1)
         end if
c 
  100    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
         cnorm = .false.
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = numcnv
c     
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
 9000 continue
c
c     %---------------%
c     | End of dnaup2 |
c     %---------------%
c
      return
      end
c\BeginDoc
c
c\Name: dnaupd
c
c\Description: 
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This subroutine computes approximations to a few eigenpairs 
c  of a linear operator "OP" with respect to a semi-inner product defined by 
c  a symmetric positive semi-definite real matrix B. B may be the identity 
c  matrix. NOTE: If the linear operator "OP" is real and symmetric 
c  with respect to the real positive semi-definite symmetric matrix B, 
c  i.e. B*OP = (OP`)*B, then subroutine dsaupd should be used instead.
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  dnaupd is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*x = amu*x, then 
c           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
c           Note: If sigma is real, i.e. imaginary part of sigma is zero;
c                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M 
c                 amu == 1/(lambda-sigma). 
c  
c  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*x = amu*x, then 
c           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
c
c  Both mode 3 and 4 give the same enhancement to eigenvalues close to
c  the (complex) shift sigma.  However, as lambda goes to infinity,
c  the operator OP in mode 4 dampens the eigenvalues more strongly than
c  does OP defined in mode 3.
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call dnaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to dnaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          dnaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3 and 4, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) real and imaginary parts 
c                    of the shifts where INPTR(14) is the pointer
c                    into WORKL for placing the shifts. See Remark
c                    5 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT/OUTPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 2 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Arnoldi vectors are generated, the algorithm generates 
c          approximately NCV-NEV Arnoldi vectors at each subsequent update 
c          iteration. Most of the cost in generating each Arnoldi vector is 
c          in the matrix-vector operation OP*x. 
c          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz 
c          values are kept together. (See remark 4 below)
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The real and imaginary
c                      parts of the NCV eigenvalues of the Hessenberg
c                      matrix H are returned in the part of the WORKL 
c                      array corresponding to RITZR and RITZI. See remark 
c                      5 below.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to 
c                      restarting the iteration with a starting vector
c                      that is a linear combination of approximate Schur
c                      vectors associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced.
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.  
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4; See under \Description of dnaupd for the 
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), dnaupd returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          5 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.        
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
c                    H in WORKL.
c          IPNTR(6): pointer to the real part of the ritz value array 
c                    RITZR in WORKL.
c          IPNTR(7): pointer to the imaginary part of the ritz value array
c                    RITZI in WORKL.
c          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZR and RITZI in WORKL.
c
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by dneupd. See Remark 2 below.
c
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the 
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of 
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
c          associated with the converged Ritz values is desired, see remark
c          2 below, subroutine dneupd uses this output.
c          See Data Distribution Note below.  
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 6*NCV.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when
c     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
c     original problem may be obtained with the ARPACK subroutine dneupd.
c
c  2. If a basis for the invariant subspace corresponding to the converged Ritz 
c     values is needed, the user must call dneupd immediately following 
c     completion of dnaupd. This is new starting with release 2 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular 
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
c     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.  The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically. 
c     See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
c     NP = IPARAM(8) real and imaginary parts of the shifts in locations 
c         real part                  imaginary part
c         -----------------------    --------------
c     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
c     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
c                        .                          .
c                        .                          .
c                        .                          .
c     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
c
c     Only complex conjugate pairs of shifts may be applied and the pairs 
c     must be placed in consecutive locations. The real part of the 
c     eigenvalues of the current upper Hessenberg matrix are located in 
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part 
c     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
c     according to the order defined by WHICH. The complex conjugate
c     pairs are kept together and the associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note: 
c
c  Fortran-D syntax:
c  ================
c  Double precision resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Double precision  resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c  
c  CM2/CM5 syntax:
c  ==============
c  
c-----------------------------------------------------------------------
c
c     include   'ex-nonsym.doc'
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     dnaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/16/93: Version '1.1'
c
c\SCCS Information: @(#) 
c FILE: naupd.F   SID: 2.10   DATE OF SID: 08/23/02   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &     ipntr, workd, workl, lworkl, info )
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw, 
     &           ldh, ldq, levec, mode,  mxiter, nb,
     &           nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, iupd, iw, ldh, ldq,
     &           levec, mode, mxiter, nb, nev0, next,
     &           np, ritzi, ritzr
!$omp threadprivate(bounds, ih, iq, ishift, iupd, iw, ldh, ldq)
!$omp threadprivate(levec, mode, mxiter, nb, nev0, next)
!$omp threadprivate(np, ritzi, ritzr)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dnaup2
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch
      external   dlamch
c
      NOPX=0; NBX=0; NRORTH=0
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        %-------------------------------%
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         ierr   = 0
         ishift = iparam(1)
c         levec  = iparam(2)
         mxiter = iparam(3)
c         nb     = iparam(4)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev+1 .or.  ncv .gt. n) then
            ierr = -3
         else if (mxiter .le.          0) then
            ierr = 4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 4) then
            ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
            ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
            ierr = -12
         end if
c 
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c 
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)				nb = 1
         if (tol .le. zero)			tol = dlamch('EpsMach')
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev 
c 
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
  10     continue
c 
c        %-------------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
c        | etc... and the remaining workspace.                         |
c        | Also update pointer to be used on output.                   |
c        | Memory is laid out as follows:                              |
c        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
c        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
c        |                                   parts of ritz values      |
c        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
c        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
c        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
c        | The final workspace is needed by subroutine dneigh called   |
c        | by dnaup2. Subroutine dneigh calls LAPACK routines for      |
c        | calculating eigenvalues and the last row of the eigenvector |
c        | matrix.                                                     |
c        %-------------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw 
c
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
      call dnaup2 
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr), 
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw), 
     &     ipntr, workd, info )
c 
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c 
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within dnaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of dnaupd |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnconv
c
c\Description: 
c  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
c
c\Usage:
c  call dnconv
c     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.
c
c  RITZR,  Double precision arrays of length N.  (INPUT)
c  RITZI   Real and imaginary parts of the Ritz values to be checked
c          for convergence.
!
c  BOUNDS  Double precision array of length N.  (INPUT)
c          Ritz estimates for the Ritz values in RITZR and RITZI.
c
c  TOL     Double precision scalar.  (INPUT)
c          Desired backward error for a Ritz value to be considered
c          "converged".
c
c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnconv (n, ritzr, ritzi, bounds, tol, nconv)
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    n, nconv
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
!
      Double precision
     &           ritzr(n), ritzi(n), bounds(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i
      Double precision
     &           temp, eps23
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2, dlamch
      external   dlapy2, dlamch
!
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %-------------------------------------------------------------%
c     | Convergence test: unlike in the symmetric code, I am not    |
c     | using things like refined error bounds and gap condition    |
c     | because I don't know the exact equivalent concept.          |
c     |                                                             |
c     | Instead the i-th Ritz value is considered "converged" when: |
c     |                                                             |
c     |     bounds(i) .le. ( TOL * | ritz | )                       |
c     |                                                             |
c     | for some appropriate choice of norm.                        |
c     %-------------------------------------------------------------%
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
      nconv  = 0
      do 20 i = 1, n
         temp = max( eps23, dlapy2( ritzr(i), ritzi(i) ) )
         if (bounds(i) .le. tol*temp)   nconv = nconv + 1
   20 continue
c 
      return
c
c     %---------------%
c     | End of dnconv |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dneigh
c
c\Description:
c  Compute the eigenvalues of the current upper Hessenberg matrix
c  and the corresponding Ritz estimates given the current residual norm.
c
c\Usage:
c  call dneigh
c     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          Residual norm corresponding to the current upper Hessenberg 
c          matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the matrix H.
c
c  H       Double precision N by N array.  (INPUT)
c          H contains the current upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZR,  Double precision arrays of length N.  (OUTPUT)
c  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real 
c          (respectively imaginary) parts of the eigenvalues of H.
c
c  BOUNDS  Double precision array of length N.  (OUTPUT)
c          On output, BOUNDS contains the Ritz estimates associated with
c          the eigenvalues RITZR and RITZI.  This is equal to RNORM 
c          times the last components of the eigenvectors corresponding 
c          to the eigenvalues in RITZR and RITZI.
c
c  Q       Double precision N by N array.  (WORKSPACE)
c          Workspace needed to store the eigenvectors of H.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  This is needed to keep the full Schur form
c          of H and also in the calculation of the eigenvectors of H.
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from dlaqrb or dtrevc.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dlaqrb  ARPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix and last row of the Schur vectors.
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dneigh (rnorm, n, h, ldh, ritzr, ritzi, bounds, 
     &                   q, ldq, workl, ierr)
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    ierr, n, ldh, ldq
      Double precision     
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision     
     &           bounds(n), h(ldh,n), q(ldq,n), ritzi(n), ritzr(n),
     &           workl(n*(n+3))
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision     
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c 
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    select(1)
      integer    i, iconj
      Double precision     
     &           temp, vl(1)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dlacpy, dlaqrb, dtrevc
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2, dnrm2
      external   dlapy2, dnrm2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic  abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------%
c
c     %-----------------------------------------------------------%
c     | 1. Compute the eigenvalues, the last components of the    |
c     |    corresponding Schur vectors and the full Schur form T  |
c     |    of the current upper Hessenberg matrix H.              |
c     | dlaqrb returns the full Schur form of H in WORKL(1:N**2)  |
c     | and the last components of the Schur vectors in BOUNDS.   |
c     %-----------------------------------------------------------%
c
      call dlacpy ('All', n, n, h, ldh, workl, n)
      call dlaqrb (.true., n, 1, n, workl, n, ritzr, ritzi, bounds,
     &             ierr)
      if (ierr .ne. 0) go to 9000
c
c     %-----------------------------------------------------------%
c     | 2. Compute the eigenvectors of the full Schur form T and  |
c     |    apply the last components of the Schur vectors to get  |
c     |    the last components of the corresponding eigenvectors. |
c     | Remember that if the i-th and (i+1)-st eigenvalues are    |
c     | complex conjugate pairs, then the real & imaginary part   |
c     | of the eigenvector components are split across adjacent   |
c     | columns of Q.                                             |
c     %-----------------------------------------------------------%
c
      call dtrevc ('R', 'A', select, n, workl, n, vl, n, q, ldq,
     &             n, n, workl(n*n+1), ierr)
c
      if (ierr .ne. 0) go to 9000
c
c     %------------------------------------------------%
c     | Scale the returning eigenvectors so that their |
c     | euclidean norms are all one. LAPACK subroutine |
c     | dtrevc returns each eigenvector normalized so  |
c     | that the element of largest magnitude has      |
c     | magnitude 1; here the magnitude of a complex   |
c     | number (x,y) is taken to be |x| + |y|.         |
c     %------------------------------------------------%
c
      iconj = 0
      do 10 i=1, n
         if ( abs( ritzi(i) ) .le. zero ) then
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c    
            temp = dnrm2( n, q(1,i), 1 )
            call dscal ( n, one / temp, q(1,i), 1 )
         else
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we further normalize by the      |
c           | square root of two.                       |
c           %-------------------------------------------%
c
            if (iconj .eq. 0) then
               temp = dlapy2( dnrm2( n, q(1,i), 1 ), 
     &                        dnrm2( n, q(1,i+1), 1 ) )
               call dscal ( n, one / temp, q(1,i), 1 )
               call dscal ( n, one / temp, q(1,i+1), 1 )
               iconj = 1
            else
               iconj = 0
            end if
         end if         
   10 continue
c
      call dgemv ('T', n, n, one, q, ldq, bounds, 1, zero, workl, 1)
c
c     %----------------------------%
c     | Compute the Ritz estimates |
c     %----------------------------%
c
      iconj = 0
      do 20 i = 1, n
         if ( abs( ritzi(i) ) .le. zero ) then
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c    
            bounds(i) = rnorm * abs( workl(i) )
         else
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we need to take the magnitude    |
c           | of the last components of the two vectors |
c           %-------------------------------------------%
c
            if (iconj .eq. 0) then
               bounds(i) = rnorm * dlapy2( workl(i), workl(i+1) )
               bounds(i+1) = bounds(i)
               iconj = 1
            else
               iconj = 0
            end if
         end if
   20 continue
c
 9000 continue
      return
c
c     %---------------%
c     | End of dneigh |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: dneupd 
c
c\Description: 
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to DNAUPD .  DNAUPD  must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  See documentation in the header of the subroutine DNAUPD  for 
c  definition of OP as well as other terms and the relation of computed
c  Ritz values and Ritz vectors of OP with respect to the given problem
c  A*z = lambda*B*z.  For a brief description, see definitions of 
c  IPARAM(7), MODE and WHICH in the documentation of DNAUPD .
c
c\Usage:
c  call dneupd  
c     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT, 
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, 
c       LWORKL, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether a basis for the invariant subspace corresponding 
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
c                                See Remarks below. 
c 
c  HOWMNY  Character*1  (INPUT) 
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors; 
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
c
c  DR      Double precision  array of dimension NEV+1.  (OUTPUT)
c          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains 
c          the real part of the Ritz  approximations to the eigenvalues of 
c          A*z = lambda*B*z. 
c          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
c          DR contains the real part of the Ritz values of OP computed by 
c          DNAUPD . A further computation must be performed by the user
c          to transform the Ritz values computed for OP by DNAUPD  to those
c          of the original system A*z = lambda*B*z. See remark 3 below.
c
c  DI      Double precision  array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value 
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR.
c
c          NOTE: When Ritz values are complex, they will come in complex 
c                conjugate pairs.  If eigenvectors are requested, the 
c                corresponding Ritz vectors will also come in conjugate 
c                pairs and the real and imaginary parts of these are 
c                represented in two consecutive columns of the array Z 
c                (see below).
c
c  Z       Double precision  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
c          Z represent approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem 
c          A*z = lambda*B*z. 
c 
c          The complex Ritz vector associated with the Ritz value 
c          with positive imaginary part is stored in two consecutive 
c          columns.  The first column holds the real part of the Ritz 
c          vector and the 2nd column holds the imaginary part.  The 
c          Ritz vector associated with the Ritz value with negative 
c          imaginary part is simply the complex conjugate of the Ritz vector 
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by DNAUPD .  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
c
c  SIGMAR  Double precision   (INPUT)
c          If IPARAM(7) = 3 or 4, represents the real part of the shift. 
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  SIGMAI  Double precision   (INPUT)
c          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift. 
c          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
c
c  WORKEV  Double precision  work array of dimension 3*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DNAUPD  that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DNEUPD  following the last call
c         to DNAUPD .  These arguments MUST NOT BE MODIFIED between
c         the the last call to DNAUPD  and the call to DNEUPD .
c
c  Three of these parameters (V, WORKL, INFO) are also output parameters:
c
c  V       Double precision  N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by DNAUPD  .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.  See Remark 2 below.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
c          dnaupd .  They are not changed by dneupd .
c          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
c          real and imaginary part of the untransformed Ritz values,
c          the upper quasi-triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by dneupd .
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd  if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine dlahqr 
c                could not be reordered by LAPACK routine dtrsen .
c                Re-enter subroutine dneupd  with IPARAM(5)=NCV and 
c                increase the size of the arrays DR and DI to have 
c                dimension at least dimension NCV and allocate at least NCV 
c                columns for Z. NOTE: Not necessary if Z and V share 
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine dlahqr .
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine dtrevc .
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: DNAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: DNEUPD got a different count of the number of converged
c                 Ritz values than DNAUPD got.  This indicates the user
c                 probably made an error in passing data from DNAUPD to
c                 DNEUPD or that the data was modified before entering
c                 DNEUPD
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     dgeqr2   LAPACK routine that computes the QR factorization of 
c             a matrix.
c     dlacpy   LAPACK matrix copy routine.
c     dlahqr   LAPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlaset   LAPACK matrix initialization routine.
c     dorm2r   LAPACK routine that applies an orthogonal matrix in 
c             factored form.
c     dtrevc   LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form.
c     dtrsen   LAPACK routine that re-orders the Schur form.
c     dtrmm    Level 3 BLAS matrix times an upper triangular matrix.
c     dger     Level 2 BLAS rank one update to a matrix.
c     dcopy    Level 1 BLAS that copies one vector to another .
c     ddot     Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dscal    Level 1 BLAS that scales a vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c     Let trans(X) denote the transpose of X.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .TRUE. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c     trans(V(:,1:IPARAM(5))) * V(:,1:IPARAM(5)) = I are approximately 
c     satisfied. Here T is the leading submatrix of order IPARAM(5) of the 
c     real upper quasi-triangular matrix stored workl(ipntr(12)). That is,
c     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; 
c     each 2-by-2 diagonal block has its diagonal elements equal and its
c     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
c     diagonal block is a complex conjugate pair of Ritz values. The real
c     Ritz values are stored on the diagonal of T.
c
c  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
c     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
c     values computed by DNAUPD  for OP to those of A*z = lambda*B*z. 
c     Set RVEC = .true. and HOWMNY = 'A', and
c     compute 
c           trans(Z(:,I)) * A * Z(:,I) if DI(I) = 0.
c     If DI(I) is not equal to zero and DI(I+1) = - D(I), 
c     then the desired real and imaginary parts of the Ritz value are
c           trans(Z(:,I)) * A * Z(:,I) +  trans(Z(:,I+1)) * A * Z(:,I+1),
c           trans(Z(:,I)) * A * Z(:,I+1) -  trans(Z(:,I+1)) * A * Z(:,I), 
c     respectively.
c     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
c     compute trans(V(:,1:IPARAM(5))) * A * V(:,1:IPARAM(5)) and then an upper
c     quasi-triangular matrix of order IPARAM(5) is computed. See remark
c     2 above.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics          
c     Rice University           
c     Houston, Texas            
c 
c\SCCS Information: @(#) 
c FILE: neupd.F   SID: 2.7   DATE OF SID: 09/20/00   RELEASE: 2 
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine dneupd (rvec , howmny, select, dr    , di,    
     &                   z    , ldz   , sigmar, sigmai, workev,
     &                   bmat , n     , which , nev   , tol,
     &                   resid, ncv   , v     , ldv   , iparam,
     &                   ipntr, workd , workl , lworkl, info)
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision      
     &           sigmar, sigmai, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision 
     &           dr(nev+1)    , di(nev+1), resid(n)  , 
     &           v(ldv,ncv)   , z(ldz,*) , workd(3*n), 
     &           workl(lworkl), workev(3*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision 
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds, ierr  , ih    , ihbds   , 
     &           iheigr, iheigi, iconj , nconv   , 
     &           invsub, iuptri, iwev  , iwork(1),
     &           j     , k     , ldh   , ldq     ,
     &           mode  , outncv, ritzr   ,
     &           ritzi , wri   , wrr   , irr     ,
     &           iri   , ibd   , ishift, numcnv  ,
     &           np    , jj 
      logical    reord
      Double precision 
     &           conds  , rnorm, sep  , temp,
     &           vl(1,1), temp1, eps23
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dger   , dgeqr2 , dlacpy , 
     &           dlahqr , dlaset , dorm2r , 
     &           dtrevc , dtrmm  , dtrsen , dscal
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dlapy2 , dnrm2 , dlamch , ddot 
      external   dlapy2 , dnrm2 , dlamch , ddot 
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, min, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch ('Epsilon-Machine')
      eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c     %--------------%
c     | Quick return |
c     %--------------%
c
      ierr = 0
c
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev+1 .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 .and. sigmai .eq. zero) then
         type = 'SHIFTI'
      else if (mode .eq. 3 ) then
         type = 'REALPT'
      else if (mode .eq. 4 ) then
         type = 'IMAGPT'
      else 
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')    ierr = -11
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c 
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
c     |                                   parts of ritz values |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by DNEUPD .                  |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                             real part of the Ritz values. |
c     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
c     |                        imaginary part of the Ritz values. |
c     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
c     |                           error bounds of the Ritz values |
c     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
c     |                             quasi-triangular matrix for H |
c     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
c     |       associated matrix representation of the invariant   |
c     |       subspace for H.                                     |
c     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
c     %-----------------------------------------------------------%
c     
      ih     = ipntr(5)
      ritzr  = ipntr(6)
      ritzi  = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheigr = bounds + ldh
      iheigi = iheigr + ldh
      ihbds  = iheigi + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheigr
      ipntr(10) = iheigi
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wrr = 1
      wri = ncv + 1
      iwev = wri + ncv
c
c     %-----------------------------------------%
c     | irr points to the REAL part of the Ritz |
c     |     values computed by _neigh before    |
c     |     exiting _naup2.                     |
c     | iri points to the IMAGINARY part of the |
c     |     Ritz values computed by _neigh      |
c     |     before exiting _naup2.              |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irr = ipntr(14)+ncv*ncv
      iri = irr+ncv
      ibd = iri+ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c
      if (rvec) then
c     
         reord = .false.
c
c        %---------------------------------------------------%
c        | Use the temporary bounds array to store indices   |
c        | These will be used to mark the select array later |
c        %---------------------------------------------------%
c
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
c
c        %-------------------------------------%
c        | Select the wanted Ritz values.      |
c        | Sort the Ritz values so that the    |
c        | wanted ones appear at the tailing   |
c        | NEV positions of workl(irr) and     |
c        | workl(iri).  Move the corresponding |
c        | error estimates in workl(bound)     |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call dngets (ishift       , which     , nev       , 
     &                np           , workl(irr), workl(iri),
     &                workl(bounds), workl     , workl(np+1))
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23,
     &                 dlapy2 ( workl(irr+ncv-j), workl(iri+ncv-j) ))
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          workl(ibd+jj-1) .le. tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by dnaupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the dnaupd data.           |
c        %-----------------------------------------------------------%
c
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine dlahqr  to compute the real Schur form |
c        | of the upper Hessenberg matrix returned by DNAUPD .        |
c        | Make a copy of the upper Hessenberg matrix.               |
c        | Initialize the Schur vector matrix Q to the identity.     |
c        %-----------------------------------------------------------%
c     
         call dcopy (ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call dlaset ('All', ncv, ncv, 
     &                zero , one, workl(invsub),
     &                ldq)
         call dlahqr (.true., .true.       , ncv, 
     &                1     , ncv          , workl(iuptri), 
     &                ldh   , workl(iheigr), workl(iheigi),
     &                1     , ncv          , workl(invsub), 
     &                ldq   , ierr)
         call dcopy (ncv         , workl(invsub+ncv-1), ldq, 
     &               workl(ihbds), 1)
c     
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c     
         if (reord) then
c     
c           %-----------------------------------------------------%
c           | Reorder the computed upper quasi-triangular matrix. | 
c           %-----------------------------------------------------%
c     
            call dtrsen ('None'       , 'V'          , 
     &                   select       , ncv          ,
     &                   workl(iuptri), ldh          , 
     &                   workl(invsub), ldq          , 
     &                   workl(iheigr), workl(iheigi), 
     &                   nconv        , conds        ,
     &                   sep          , workl(ihbds) , 
     &                   ncv          , iwork        ,
     &                   1            , ierr)
c
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
c
         end if
c
c        %---------------------------------------%
c        | Copy the last row of the Schur vector |
c        | into workl(ihbds).  This will be used |
c        | to compute the Ritz estimates of      |
c        | converged Ritz values.                |
c        %---------------------------------------%
c
         call dcopy (ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
c
c        %----------------------------------------------------%
c        | Place the computed eigenvalues of H into DR and DI |
c        | if a spectral transformation was not used.         |
c        %----------------------------------------------------%
c
         if (type .eq. 'REGULR') then 
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
         end if
c     
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(invsub,ldq).                            |
c        %----------------------------------------------------------%
c     
         call dgeqr2 (ncv, nconv , workl(invsub), 
     &               ldq, workev, workev(ncv+1),
     &               ierr)
c
c        %---------------------------------------------------------%
c        | * Postmultiply V by Q using dorm2r .                     |   
c        | * Copy the first NCONV columns of VQ into Z.            |
c        | * Postmultiply Z by R.                                  |
c        | The N by NCONV matrix Z is now a matrix representation  |
c        | of the approximate invariant subspace associated with   |
c        | the Ritz values in workl(iheigr) and workl(iheigi)      |
c        | The first NCONV columns of V are now approximate Schur  |
c        | vectors associated with the real upper quasi-triangular |
c        | matrix of order NCONV in workl(iuptri)                  |
c        %---------------------------------------------------------%
c     
         call dorm2r ('Right', 'Notranspose', n            , 
     &                ncv   , nconv        , workl(invsub),
     &                ldq   , workev       , v            , 
     &                ldv   , workd(n+1)   , ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
c
         do 20 j=1, nconv
c     
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | quasi-triangular form of workl(iuptri,ldq)        |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones           |
c           %---------------------------------------------------%
c     
            if (workl(invsub+(j-1)*ldq+j-1) .lt. zero) then
               call dscal (nconv, -one, workl(iuptri+j-1), ldq)
               call dscal (nconv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
c     
 20      continue
c     
         if (howmny .eq. 'A') then
c     
c           %--------------------------------------------%
c           | Compute the NCONV wanted eigenvectors of T | 
c           | located in workl(iuptri,ldq).              |
c           %--------------------------------------------%
c     
            do 30 j=1, ncv
               if (j .le. nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
c
            call dtrevc ('Right', 'Select'     , select       , 
     &                   ncv    , workl(iuptri), ldq          , 
     &                   vl     , 1            , workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev , ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c     
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | dtrevc  returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1;                                   |
c           %------------------------------------------------%
c     
            iconj = 0
            do 40 j=1, nconv
c
               if ( workl(iheigi+j-1) .eq. zero ) then
c     
c                 %----------------------%
c                 | real eigenvalue case |
c                 %----------------------%
c     
                  temp = dnrm2 ( ncv, workl(invsub+(j-1)*ldq), 1 )
                  call dscal ( ncv, one / temp, 
     &                 workl(invsub+(j-1)*ldq), 1 )
c
               else
c     
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 | columns, we further normalize by the      |
c                 | square root of two.                       |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     temp = dlapy2 (dnrm2 (ncv, 
     &                                   workl(invsub+(j-1)*ldq), 
     &                                   1),
     &                             dnrm2 (ncv, 
     &                                   workl(invsub+j*ldq),
     &                                   1))  
                     call dscal (ncv, one/temp, 
     &                           workl(invsub+(j-1)*ldq), 1 )
                     call dscal (ncv, one/temp, 
     &                           workl(invsub+j*ldq), 1 )
                     iconj = 1
                  else
                     iconj = 0
                  end if
c
               end if
c
 40         continue
c
            call dgemv ('T', ncv, nconv, one, workl(invsub),
     &                 ldq, workl(ihbds), 1, zero,  workev, 1)
c
            iconj = 0
            do 45 j=1, nconv
               if (workl(iheigi+j-1) .ne. zero) then
c
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     workev(j) = dlapy2 (workev(j), workev(j+1))
                     workev(j+1) = workev(j)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               end if
 45         continue
c
c           %---------------------------------------%
c           | Copy Ritz estimates into workl(ihbds) |
c           %---------------------------------------%
c
            call dcopy (nconv, workev, 1, workl(ihbds), 1)
c
c           %---------------------------------------------------------%
c           | Compute the QR factorization of the eigenvector matrix  |
c           | associated with leading portion of T in the first NCONV |
c           | columns of workl(invsub,ldq).                           |
c           %---------------------------------------------------------%
c     
            call dgeqr2 (ncv, nconv , workl(invsub), 
     &                   ldq, workev, workev(ncv+1),
     &                   ierr)
c     
c           %----------------------------------------------%
c           | * Postmultiply Z by Q.                       |   
c           | * Postmultiply Z by R.                       |
c           | The N by NCONV matrix Z is now contains the  | 
c           | Ritz vectors associated with the Ritz values |
c           | in workl(iheigr) and workl(iheigi).          |
c           %----------------------------------------------%
c     
            call dorm2r ('Right', 'Notranspose', n            ,
     &                   ncv  , nconv        , workl(invsub),
     &                   ldq  , workev       , z            ,
     &                   ldz  , workd(n+1)   , ierr)
c     
            call dtrmm ('Right'   , 'Upper'       , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
c     
         end if
c     
      else 
c
c        %------------------------------------------------------%
c        | An approximate invariant subspace is not needed.     |
c        | Place the Ritz values computed DNAUPD  into DR and DI |
c        %------------------------------------------------------%
c
         call dcopy (nconv, workl(ritzr), 1, dr, 1)
         call dcopy (nconv, workl(ritzi), 1, di, 1)
         call dcopy (nconv, workl(ritzr), 1, workl(iheigr), 1)
         call dcopy (nconv, workl(ritzi), 1, workl(iheigi), 1)
         call dcopy (nconv, workl(bounds), 1, workl(ihbds), 1)
      end if
c 
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
         if (rvec) 
     &      call dscal (ncv, rnorm, workl(ihbds), 1)     
c     
      else 
c     
c        %---------------------------------------%
c        |   A spectral transformation was used. |
c        | * Determine the Ritz estimates of the |
c        |   Ritz values in the original system. |
c        %---------------------------------------%
c     
         if (type .eq. 'SHIFTI') then
c
            if (rvec) 
     &         call dscal (ncv, rnorm, workl(ihbds), 1)
c
            do 50 k=1, ncv
               temp = dlapy2 ( workl(iheigr+k-1), 
     &                        workl(iheigi+k-1) )
               workl(ihbds+k-1) = abs( workl(ihbds+k-1) ) 
     &                          / temp / temp
 50         continue
c
         else if (type .eq. 'REALPT') then
c
            do 60 k=1, ncv
 60         continue
c
         else if (type .eq. 'IMAGPT') then
c
            do 70 k=1, ncv
 70         continue
c
         end if
c     
c        %-----------------------------------------------------------%
c        | *  Transform the Ritz values back to the original system. |
c        |    For TYPE = 'SHIFTI' the transformation is              |
c        |             lambda = 1/theta + sigma                      |
c        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
c        |    Rayleigh quotients or a projection. See remark 3 above.| 
c        | NOTES:                                                    |
c        | *The Ritz vectors are not affected by the transformation. |
c        %-----------------------------------------------------------%
c     
         if (type .eq. 'SHIFTI') then 
c
            do 80 k=1, ncv
               temp = dlapy2 ( workl(iheigr+k-1), 
     &                        workl(iheigi+k-1) )
               workl(iheigr+k-1) = workl(iheigr+k-1)/temp/temp 
     &                           + sigmar   
               workl(iheigi+k-1) = -workl(iheigi+k-1)/temp/temp
     &                           + sigmai   
 80         continue
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         else if (type .eq. 'REALPT' .or. type .eq. 'IMAGPT') then
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         end if
c
      end if
c
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 2.                                   |
c     %-------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------%
c
         iconj = 0
         do 110 j=1, nconv
            if (workl(iheigi+j-1) .eq. zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) /
     &                      workl(iheigr+j-1)
            else if (iconj .eq. 0) then
               temp = dlapy2 ( workl(iheigr+j-1), workl(iheigi+j-1) )
               workev(j) = ( workl(invsub+(j-1)*ldq+ncv-1) * 
     &                       workl(iheigr+j-1) +
     &                       workl(invsub+j*ldq+ncv-1) * 
     &                       workl(iheigi+j-1) ) / temp / temp
               workev(j+1) = ( workl(invsub+j*ldq+ncv-1) * 
     &                         workl(iheigr+j-1) -
     &                         workl(invsub+(j-1)*ldq+ncv-1) * 
     &                         workl(iheigi+j-1) ) / temp / temp
               iconj = 1
            else
               iconj = 0
            end if
 110     continue
c
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
      return
c     
c     %---------------%
c     | End of DNEUPD  |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dngets
c
c\Description: 
c  Given the eigenvalues of the upper Hessenberg matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of 
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: call this even in the case of user specified shifts in order
c  to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call dngets
c     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> want the KEV eigenvalues of largest magnitude.
c          'SM' -> want the KEV eigenvalues of smallest magnitude.
c          'LR' -> want the KEV eigenvalues of largest real part.
c          'SR' -> want the KEV eigenvalues of smallest real part.
c          'LI' -> want the KEV eigenvalues of largest imaginary part.
c          'SI' -> want the KEV eigenvalues of smallest imaginary part.
c
c  KEV      Integer.  (INPUT/OUTPUT)
c           INPUT: KEV+NP is the size of the matrix H.
c           OUTPUT: Possibly increases KEV by one to keep complex conjugate
c           pairs together.
c
c  NP       Integer.  (INPUT/OUTPUT)
c           Number of implicit shifts to be computed.
c           OUTPUT: Possibly decreases NP by one to keep complex conjugate
c           pairs together.
c
c  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary 
c          parts of the eigenvalues of H.
c          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
c          eigenvalues are in the first NP locations and the wanted
c          portion is in the last KEV locations.  When exact shifts are 
c          selected, the unwanted part corresponds to the shifts to 
c          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
c          are further sorted so that the ones with largest Ritz values
c          are first.
c
c  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
c  
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dsortc  ARPACK sorting routine.
c     dcopy   Level 1 BLAS that copies one vector to another .
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dngets ( ishift, which, kev, np, ritzr, ritzi, bounds,
     &                    shiftr, shifti )
c
c     %----------------------------------------------------%
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      integer    ishift, kev, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           bounds(kev+np), ritzr(kev+np), ritzi(kev+np), 
     &           shiftr(1), shifti(1)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dsortc
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------%
c 
c     %----------------------------------------------------%
c     | LM, SM, LR, SR, LI, SI case.                       |
c     | Sort the eigenvalues of H into the desired order   |
c     | and apply the resulting order to BOUNDS.           |
c     | The eigenvalues are sorted so that the wanted part |
c     | are always in the last KEV locations.              |
c     | We first do a pre-processing sort in order to keep |
c     | complex conjugate pairs together                   |
c     %----------------------------------------------------%
c
      if (which .eq. 'LM') then
         call dsortc ('LR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SM') then
         call dsortc ('SR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'LR') then
         call dsortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SR') then
         call dsortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'LI') then
         call dsortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SI') then
         call dsortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      end if
c      
      call dsortc (which, .true., kev+np, ritzr, ritzi, bounds)
c     
c     %-------------------------------------------------------%
c     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
c     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero |
c     | Accordingly decrease NP by one. In other words keep   |
c     | complex conjugate pairs together.                     |
c     %-------------------------------------------------------%
c     
      if (       ( ritzr(np+1) - ritzr(np) ) .eq. zero
     &     .and. ( ritzi(np+1) + ritzi(np) ) .eq. zero ) then
         np = np - 1
         kev = kev + 1
      end if
c
      if ( ishift .eq. 1 ) then
c     
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when they shifts |
c        | are applied in subroutine dnapps.                     |
c        | Be careful and use 'SR' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c     
         call dsortc ( 'SR', .true., np, bounds, ritzr, ritzi )
      end if
c     
      return
c     
c     %---------------%
c     | End of dngets |
c     %---------------%
c     
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaitr
c
c\Description: 
c  Reverse communication interface for applying NP additional steps to 
c  a K step symmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in dsaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call dsaitr
c     ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and does not need to be
c          recomputed in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of matrix B that defines the
c          semi-inner product for the operator OP.  See dsaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current order of H and the number of columns of V.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  MODE    Integer.  (INPUT)
c          Signifies which form for "OP". If MODE=2 then
c          a reduction in the number of B matrix vector multiplies
c          is possible since the B-norm of OP*x is equivalent to
c          the inv(B)-norm of A*x.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Double precision scalar.  (INPUT/OUTPUT)
c          On INPUT the B-norm of r_{k}.
c          On OUTPUT the B-norm of the updated residual r_{k+p}.
c
c  V       Double precision N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K 
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (K+NP) by 2 array.  (INPUT/OUTPUT)
c          H is used to store the generated symmetric tridiagonal matrix
c          with the subdiagonal in the first column starting at H(2,1)
c          and the main diagonal in the 2nd column.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not 
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On INPUT, WORKD(1:N) = B*RESID where RESID is associated
c          with the K step Arnoldi factorization. Used to save some 
c          computation at the first step. 
c          On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
c          with the K+NP step Arnoldi factorization.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of an invariant subspace of OP is found that is
c               less than K + NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dgetv0  ARPACK routine to generate the initial vector.
c     dlamch  LAPACK routine that determines machine constants.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     xx/xx/93: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: saitr.F   SID: 2.6   DATE OF SID: 8/28/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c  
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in dsaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        alphaj <- j-th component of w_{j}
c        rnorm = || r_{j} ||
c        betaj+1 = rnorm
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};   
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf 
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsaitr
     &   (ido, bmat, n, k, np, mode, resid, rnorm, v, ldv, h, ldh, 
     &    ipntr, workd, info)
c
c     %----------------------------------------------------%
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, mode, np
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           h(ldh,2), resid(n), v(ldv,k+np), workd(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    first, orth1, orth2, rstart, step3, step4
      integer    i, ierr, ipj, irj, ivj, iter, itry, j, 
     &           infol, jj
      Double precision
     &           rnorm1, wnorm, safmin, temp1
      save       orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, 
     &           rnorm1, safmin, wnorm
!$omp threadprivate(orth1, orth2, rstart, step3, step4)
!$omp threadprivate(ierr, ipj, irj, ivj, iter, itry, j )
!$omp threadprivate(rnorm1, safmin, wnorm)
c
c     %-----------------------%
c     | Local Array Arguments | 
c     %-----------------------%
c
      Double precision
     &           xtemp(2)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dgemv, dgetv0, 
     &           dlascl
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlamch
      external   ddot, dnrm2, dlamch
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
c      data      first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      NOPX=0; NRORTH=0; NBX=0
      first=.true.
      if (first) then
         first = .false.
c
c        %--------------------------------%
c        | safmin = safe minimum is such  |
c        | that 1/sfmin does not overflow |
c        %--------------------------------%
c
         safmin = dlamch('safmin')
      end if
c
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        %-------------------------------%
c
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
c 
c        %--------------------------------%
c        | Pointer to the current step of |
c        | the factorization to build     |
c        %--------------------------------%
c
         j      = k + 1
c 
c        %------------------------------------------%
c        | Pointers used for reverse communication  |
c        | when using WORKD.                        |
c        %------------------------------------------%
c
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
c 
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true.                                  |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         dgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %------------------------------%
c     | Else this is the first step. |
c     %------------------------------%
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
c
 1000 continue
c
c        %---------------------------------------------------------%
c        | Check for exact zero. Equivalent to determing whether a |
c        | j-step Arnoldi factorization is present.                |
c        %---------------------------------------------------------%
c
         if (rnorm .gt. zero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted.         |
c           %---------------------------------------------%
c
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
            call dgetv0 (ido, bmat, itry, .false., n, j, v, ldv, 
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               ido = 99
               go to 9000
            end if
c 
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call dcopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. safmin) then
             temp1 = one / rnorm
             call dscal (n, temp1, v(1,j), 1)
             call dscal (n, temp1, workd(ipj), 1)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
             call dlascl ('General', i, i, rnorm, one, n, 1, 
     &                    v(1,j), n, infol)
             call dlascl ('General', i, i, rnorm, one, n, 1, 
     &                    workd(ipj), n, infol)
         end if
c 
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         nopx  = nopx + 1
         call dcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c 
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c 
         go to 9000
   50    continue
c 
c        %-----------------------------------%
c        | Back from reverse communication;  |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
c        %-----------------------------------%
c
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call dcopy (n, workd(irj), 1, resid, 1)
c 
c        %-------------------------------------------%
c        | STEP 4:  Finish extending the symmetric   |
c        |          Arnoldi to length j. If MODE = 2 |
c        |          then B*OP = B*inv(B)*A = A and   |
c        |          we don't need to compute B*OP.   |
c        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
c        | assumed to have A*v_{j}.                  |
c        %-------------------------------------------%
c
         if (mode .eq. 2) go to 65
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
              call dcopy(n, resid, 1 , workd(ipj), 1)
         end if
   60    continue
c 
c        %-----------------------------------%
c        | Back from reverse communication;  |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
c        %-----------------------------------%
c
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
   65    continue
         if (mode .eq. 2) then
c
c           %----------------------------------%
c           | Note that the B-norm of OP*v_{j} |
c           | is the inv(B)-norm of A*v_{j}.   |
c           %----------------------------------%
c
            wnorm = ddot (n, resid, 1, workd(ivj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'G') then         
            wnorm = ddot (n, resid, 1, workd(ipj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c
         if (mode .ne. 2 ) then
            call dgemv('T', n, j, one, v, ldv, workd(ipj), 1, zero, 
     &                  workd(irj), 1)
         else if (mode .eq. 2) then
            call dgemv('T', n, j, one, v, ldv, workd(ivj), 1, zero, 
     &                  workd(irj), 1)
         end if
c
c        %--------------------------------------%
c        | Orthgonalize r_{j} against V_{j}.    |
c        | RESID contains OP*v_{j}. See STEP 3. | 
c        %--------------------------------------%
c
         call dgemv('N', n, j, -one, v, ldv, workd(irj), 1, one, 
     &               resid, 1)
c
c        %--------------------------------------%
c        | Extend H to have j rows and columns. |
c        %--------------------------------------%
c
         h(j,2) = workd(irj + j - 1)
         if (j .eq. 1  .or.  rstart) then
            h(j,1) = zero
         else
            h(j,1) = rnorm
         end if
         orth1 = .true.
         iter  = 0
c 
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
c 
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        %-----------------------------------------------------------%
c
         if (rnorm .gt. 0.717*wnorm) go to 100
         nrorth = nrorth + 1
c 
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c
   80    continue
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1, 
     &               zero, workd(irj), 1)
c
c        %----------------------------------------------%
c        | Compute the correction to the residual:      |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
c        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
c        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
c        | H(j,j) is updated.                           |
c        %----------------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, workd(irj), 1, 
     &               one, resid, 1)
c
         if (j .eq. 1  .or.  rstart) h(j,1) = zero
         h(j,2) = h(j,2) + workd(irj + j - 1)
c 
         orth2 = .true.
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c 
         if (bmat .eq. 'G') then         
             rnorm1 = ddot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if (rnorm1 .gt. 0.717*rnorm) then
c
c           %--------------------------------%
c           | No need for further refinement |
c           %--------------------------------%
c
            rnorm = rnorm1
c 
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required.    |
c           %-------------------------------------------%
c
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
c 
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c
  100    continue
c 
         rstart = .false.
         orth2  = .false.
c 
c        %----------------------------------------------------------%
c        | Make sure the last off-diagonal element is non negative  |
c        | If not perform a similarity transformation on H(1:j,1:j) |
c        | and scale v(:,j) by -1.                                  |
c        %----------------------------------------------------------%
c
         if (h(j,1) .lt. zero) then
            h(j,1) = -h(j,1)
            if ( j .lt. k+np) then 
               call dscal(n, -one, v(1,j+1), 1)
            else
               call dscal(n, -one, resid, 1)
            end if
         end if
c 
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            ido = 99
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c 
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of dsaitr |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP shifts implicitly resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix of order KEV+NP. Q is the product of 
c  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi 
c  factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call dsapps
c     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.
c
c  KEV     Integer.  (INPUT)
c          INPUT: KEV+NP is the size of the input matrix H.
c          OUTPUT: KEV is the size of the updated matrix HNEW.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFT   Double precision array of length NP.  (INPUT)
c          The shifts to be applied.
c
c  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
c          INPUT: V contains the current KEV+NP Arnoldi vectors.
c          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
c          are in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (KEV+NP) by 2 array.  (INPUT/OUTPUT)
c          INPUT: H contains the symmetric tridiagonal matrix of the
c          Arnoldi factorization with the subdiagonal in the 1st column
c          starting at H(2,1) and the main diagonal in the 2nd column.
c          OUTPUT: H contains the updated tridiagonal matrix in the 
c          KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Double precision array of length (N).  (INPUT/OUTPUT)
c          INPUT: RESID contains the the residual vector r_{k+p}.
c          OUTPUT: RESID is the updated residual vector rnew_{k}.
c
c  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations during the bulge
c          chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dlamch  LAPACK routine that determines machine constants.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlacpy  LAPACK matrix copy routine.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dscal   Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     12/16/93: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: sapps.F   SID: 2.6   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the subblocks of
c     the tridiagonal matrix H and not just to the submatrix that it 
c     comes from. This routine assumes that the subdiagonal elements 
c     of H that are stored in h(1:kev+np,1) are nonegative upon input
c     and enforce this condition upon output. This version incorporates
c     deflation. See code for documentation.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsapps
     &   ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workd )
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h(ldh,2), q(ldq,kev+np), resid(n), shift(np), 
     &           v(ldv,kev+np), workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, iend, istart, itop, j, jj, kplusp
      logical    first
      Double precision
     &           a1, a2, a3, a4, big, c, epsmch, f, g, r, s
      save       epsmch, first
!$omp threadprivate(epsmch, first)
c
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dlacpy, dlartg, dlaset,  
     &           dgemv
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch
      external   dlamch
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
         epsmch = dlamch('Epsilon-Machine')
         first = .false.
      end if
      itop = 1
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------%
c
      kplusp = kev + np 
c 
c     %----------------------------------------------%
c     | Initialize Q to the identity matrix of order |
c     | kplusp used to accumulate the rotations.     |
c     %----------------------------------------------%
c
      call dlaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c 
c     %----------------------------------------------------------%
c     | Apply the np shifts implicitly. Apply each shift to the  |
c     | whole matrix and not just to the submatrix from which it |
c     | comes.                                                   |
c     %----------------------------------------------------------%
c
      do 90 jj = 1, np
c 
         istart = itop
c
c        %----------------------------------------------------------%
c        | Check for splitting and deflation. Currently we consider |
c        | an off-diagonal element h(i+1,1) negligible if           |
c        |         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   |
c        | for i=1:KEV+NP-1.                                        |
c        | If above condition tests true then we set h(i+1,1) = 0.  |
c        | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
c        %----------------------------------------------------------%
c
   20    continue
c
c        %------------------------------------------------%
c        | The following loop exits early if we encounter |
c        | a negligible off diagonal element.             |
c        %------------------------------------------------%
c
         do 30 i = istart, kplusp-1
            big   = abs(h(i,2)) + abs(h(i+1,2))
            if (h(i+1,1) .le. epsmch*big) then
               h(i+1,1) = zero
               iend = i
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
         if (istart .lt. iend) then
c 
c           %--------------------------------------------------------%
c           | Construct the plane rotation G'(istart,istart+1,theta) |
c           | that attempts to drive h(istart+1,1) to zero.          |
c           %--------------------------------------------------------%
c
             f = h(istart,2) - shift(jj)
             g = h(istart+1,1)
             call dlartg (f, g, c, s, r)
c 
c            %-------------------------------------------------------%
c            | Apply rotation to the left and right of H;            |
c            | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
c            | This will create a "bulge".                           |
c            %-------------------------------------------------------%
c
             a1 = c*h(istart,2)   + s*h(istart+1,1)
             a2 = c*h(istart+1,1) + s*h(istart+1,2)
             a4 = c*h(istart+1,2) - s*h(istart+1,1)
             a3 = c*h(istart+1,1) - s*h(istart,2) 
             h(istart,2)   = c*a1 + s*a2
             h(istart+1,2) = c*a4 - s*a3
             h(istart+1,1) = c*a3 + s*a4
c 
c            %----------------------------------------------------%
c            | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c            %----------------------------------------------------%
c
             do 60 j = 1, min(istart+jj,kplusp)
                a1            =   c*q(j,istart) + s*q(j,istart+1)
                q(j,istart+1) = - s*q(j,istart) + c*q(j,istart+1)
                q(j,istart)   = a1
   60        continue
c
c
c            %----------------------------------------------%
c            | The following loop chases the bulge created. |
c            | Note that the previous rotation may also be  |
c            | done within the following loop. But it is    |
c            | kept separate to make the distinction among  |
c            | the bulge chasing sweeps and the first plane |
c            | rotation designed to drive h(istart+1,1) to  |
c            | zero.                                        |
c            %----------------------------------------------%
c
             do 70 i = istart+1, iend-1
c 
c               %----------------------------------------------%
c               | Construct the plane rotation G'(i,i+1,theta) |
c               | that zeros the i-th bulge that was created   |
c               | by G(i-1,i,theta). g represents the bulge.   |
c               %----------------------------------------------%
c
                f = h(i,1)
                g = s*h(i+1,1)
c
c               %----------------------------------%
c               | Final update with G(i-1,i,theta) |
c               %----------------------------------%
c
                h(i+1,1) = c*h(i+1,1)
                call dlartg (f, g, c, s, r)
c
c               %-------------------------------------------%
c               | The following ensures that h(1:iend-1,1), |
c               | the first iend-2 off diagonal of elements |
c               | H, remain non negative.                   |
c               %-------------------------------------------%
c
                if (r .lt. zero) then
                   r = -r
                   c = -c
                   s = -s
                end if
c 
c               %--------------------------------------------%
c               | Apply rotation to the left and right of H; |
c               | H <- G * H * G',  where G = G(i,i+1,theta) |
c               %--------------------------------------------%
c
                h(i,1) = r
c 
                a1 = c*h(i,2)   + s*h(i+1,1)
                a2 = c*h(i+1,1) + s*h(i+1,2)
                a3 = c*h(i+1,1) - s*h(i,2)
                a4 = c*h(i+1,2) - s*h(i+1,1)
c 
                h(i,2)   = c*a1 + s*a2
                h(i+1,2) = c*a4 - s*a3
                h(i+1,1) = c*a3 + s*a4
c 
c               %----------------------------------------------------%
c               | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c               %----------------------------------------------------%
c
                do 50 j = 1, min( i+jj, kplusp )
                   a1       =   c*q(j,i) + s*q(j,i+1)
                   q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                   q(j,i)   = a1
   50           continue
c
   70        continue
c
         end if
c
c        %--------------------------%
c        | Update the block pointer |
c        %--------------------------%
c
         istart = iend + 1
c
c        %------------------------------------------%
c        | Make sure that h(iend,1) is non-negative |
c        | If not then set h(iend,1) <-- -h(iend,1) |
c        | and negate the last column of Q.         |
c        | We have effectively carried out a        |
c        | similarity on transformation H           |
c        %------------------------------------------%
c
         if (h(iend,1) .lt. zero) then
             h(iend,1) = -h(iend,1)
             call dscal(kplusp, -one, q(1,iend), 1)
         end if
c
c        %--------------------------------------------------------%
c        | Apply the same shift to the next block if there is any |
c        %--------------------------------------------------------%
c
         if (iend .lt. kplusp) go to 20
c
c        %-----------------------------------------------------%
c        | Check if we can increase the the start of the block |
c        %-----------------------------------------------------%
c
         do 80 i = itop, kplusp-1
            if (h(i+1,1) .gt. zero) go to 90
            itop  = itop + 1
   80    continue
c
c        %-----------------------------------%
c        | Finished applying the jj-th shift |
c        %-----------------------------------%
c
   90 continue
c
c     %------------------------------------------%
c     | All shifts have been applied. Check for  |
c     | more possible deflation that might occur |
c     | after the last shift is applied.         |                               
c     %------------------------------------------%
c
      do 100 i = itop, kplusp-1
         big   = abs(h(i,2)) + abs(h(i+1,2))
         if (h(i+1,1) .le. epsmch*big) then
            h(i+1,1) = zero
         end if
 100  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is not necessary if h(kev+1,1) = 0.         |
c     %-------------------------------------------------%
c
      if ( h(kev+1,1) .gt. zero ) 
     &   call dgemv ('N', n, kplusp, one, v, ldv,
     &                q(1,kev+1), 1, zero, workd(n+1), 1)
c 
c     %-------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order    |
c     | taking advantage that Q is an upper triangular matrix |    
c     | with lower bandwidth np.                              |
c     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
c     %-------------------------------------------------------%
c
      do 130 i = 1, kev
         call dgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call dcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  130 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call dlacpy ('All', n, kev, v(1,np+1), ldv, v, ldv)
c 
c     %--------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the |
c     | appropriate place if h(kev+1,1) .ne. zero. |
c     %--------------------------------------------%
c
      if ( h(kev+1,1) .gt. zero ) 
     &     call dcopy (n, workd(n+1), 1, v(1,kev+1), 1)
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call dscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,1) .gt. zero) 
     &   call daxpy (n, h(kev+1,1), v(1,kev+1), 1, resid, 1)
c
 9000 continue 
      return
c
c     %---------------%
c     | End of dsapps |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaup2
c
c\Description: 
c  Intermediate level interface called by dsaupd.
c
c\Usage:
c  call dsaup2 
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dsaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in dsaupd.
c  
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during 
c          each Arnoldi/Lanczos iteration.  
c          If ISHIFT=1, NP is adjusted dynamically at each iteration 
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector 
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV since a leading block of the current
c          upper Tridiagonal matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number 
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Lanczos basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (NEV+NP) by 2 array.  (OUTPUT)
c          H is used to store the generated symmetric tridiagonal matrix
c          The subdiagonal is stored in the first column of H starting 
c          at H(2,1).  The main diagonal is stored in the 2 column
c          of H starting at H(1,2). If dsaup2 converges store the 
c          B-norm of the final residual vector in H(1,1).
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZ    Double precision array of length NEV+NP.  (OUTPUT)
c          RITZ(1:NEV) contains the computed Ritz values of OP.
c
c  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.
c
c  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the 
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c          
c  WORKL   Double precision array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in the computation of the 
c          tridiagonal eigenvalue problem, the calculation and
c          application of the shifts and convergence checking.
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations
c          of WORKL are used in reverse communication to hold the user 
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in one of  
c                    the spectral transformation modes.  X is the current
c                    operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Lanczos iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in dsaupd.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: All possible eigenvalues of OP has been found.  
c                   NP returns the size of the invariant subspace
c                   spanning the operator OP. 
c          =     2: No shifts could be applied.
c          =    -8: Error return from trid. eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Lanczos factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c
c\Routines called:
c     dgetv0  ARPACK initial vector generation routine. 
c     dsaitr  ARPACK Lanczos factorization routine.
c     dsapps  ARPACK application of implicit shifts routine.
c     dsconv  ARPACK convergence of Ritz values routine.
c     dseigt  ARPACK compute Ritz values and error bounds routine.
c     dsgets  ARPACK reorder Ritz values and error bounds routine.
c     dsortr  ARPACK sorting routine.
c     dlamch  LAPACK routine that determines machine constants.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     dswap   Level 1 BLAS that swaps two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/15/93: Version ' 2.4'
c     xx/xx/95: Version ' 2.4'.  (R.B. Lehoucq)
c
c\SCCS Information: @(#) 
c FILE: saup2.F   SID: 2.7   DATE OF SID: 5/19/98   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, 
     &     ishift, mxiter, v, ldv, h, ldh, ritz, bounds, 
     &     q, ldq, workl, ipntr, workd, info )
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, ldh, ldq, ldv, mxiter,
     &           n, mode, nev, np
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           bounds(nev+np), h(ldh,2), q(ldq,nev+np), resid(n), 
     &           ritz(nev+np), v(ldv,nev+np), workd(3*n), 
     &           workl(3*(nev+np))
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  wprime*2
      logical    cnorm, getv0, initv, update, ushift
      integer    ierr, iter, j, kplusp, nconv, nevbef, nev0, 
     &           np0, nptemp, nevd2, nevm2, kp(3) 
      Double precision
     &           rnorm, temp, eps23
      save       cnorm, getv0, initv, update, ushift,
     &           iter, kplusp, nconv, nev0, np0,
     &           rnorm, eps23
!$omp threadprivate(cnorm, getv0, initv, update, ushift)
!$omp threadprivate(iter, kplusp, nconv, nev0, np0)
!$omp threadprivate(rnorm, eps23)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dgetv0, dsaitr, dscal, dsconv, dseigt, dsgets, 
     &           dsapps, dsortr, dswap
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlamch
      external   ddot, dnrm2, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min
c
      NBX=0
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        %-------------------------------%
c
c        %---------------------------------%
c        | Set machine dependent constant. |
c        %---------------------------------%
c
         eps23 = dlamch('Epsilon-Machine')
         eps23 = eps23**(2.0D+0/3.0D+0)
c
c        %-------------------------------------%
c        | nev0 and np0 are integer variables  |
c        | hold the initial values of NEV & NP |
c        %-------------------------------------%
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev0 + np0
         nconv  = 0
         iter   = 0
c 
c        %--------------------------------------------%
c        | Set flags for computing the first NEV steps |
c        | of the Lanczos factorization.              |
c        %--------------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c        %--------------------------------------------%
c        | User provides the initial residual vector. |
c        %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c 
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call dgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. zero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1200
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %------------------------------------------------------------%
c     | Back from reverse communication: continue with update step |
c     %------------------------------------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c 
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Lanczos factorization |
c     %----------------------------------------------------------%
c
      call dsaitr (ido, bmat, n, 0, nev0, mode, resid, rnorm, v, ldv, 
     &             h, ldh, ipntr, workd, info)
c 
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
c
c        %-----------------------------------------------------%
c        | dsaitr was unable to build an Lanczos factorization |
c        | of length NEV0. INFO is returned with the size of   |
c        | the factorization built. Exit main loop.            |
c        %-----------------------------------------------------%
c
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  LANCZOS  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Lanczos     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c 
 1000 continue
c
         iter = iter + 1
c
c        %------------------------------------------------------------%
c        | Compute NP additional steps of the Lanczos factorization. |
c        %------------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call dsaitr (ido, bmat, n, nev, np, mode, resid, rnorm, v, 
     &                ldv, h, ldh, ipntr, workd, info)
c 
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
c
c           %-----------------------------------------------------%
c           | dsaitr was unable to build an Lanczos factorization |
c           | of length NEV0+NP0. INFO is returned with the size  |  
c           | of the factorization built. Exit main loop.         |
c           %-----------------------------------------------------%
c
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current symmetric tridiagonal matrix.           |
c        %--------------------------------------------------------%
c
         call dseigt (rnorm, kplusp, h, ldh, ritz, bounds, workl, ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained from _seigt.                       |
c        %----------------------------------------------------%
c
         call dcopy(kplusp, ritz, 1, workl(kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The selection is based on the requested number of |
c        | eigenvalues instead of the current NEV and NP to  |
c        | prevent possible misconvergence.                  |
c        | * Wanted Ritz values := RITZ(NP+1:NEV+NP)         |
c        | * Shifts := RITZ(1:NP) := WORKL(1:NP)             |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         call dsgets (ishift, which, nev, np, ritz, bounds, workl)
c 
c        %-------------------%
c        | Convergence test. |
c        %-------------------%
c
         call dcopy (nev, bounds(np+1), 1, workl(np+1), 1)
         call dsconv (nev, ritz(np+1), workl(np+1), tol, nconv)
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c 
         if ( (nconv .ge. NEV0_FIX) .or.
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c     
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP since we don't want to |
c           | swap overlapping locations.                    |
c           %------------------------------------------------%
c
            if (which .eq. 'BE') then
c
c              %-----------------------------------------------------%
c              | Both ends of the spectrum are requested.            |
c              | Sort the eigenvalues into algebraically decreasing  |
c              | order first then swap low end of the spectrum next  |
c              | to high end in appropriate locations.               |
c              | NOTE: when np < floor(nev/2) be careful not to swap |
c              | overlapping locations.                              |
c              %-----------------------------------------------------%
c
               wprime = 'SA'
               call dsortr (wprime, .true., kplusp, ritz, bounds)
               nevd2 = nev0 / 2
               nevm2 = nev0 - nevd2 
               if ( nev .gt. 1 ) then
                  call dswap ( min(nevd2,np), ritz(nevm2+1), 1,
     &                 ritz( max(kplusp-nevd2+1,kplusp-np+1) ), 1)
                  call dswap ( min(nevd2,np), bounds(nevm2+1), 1,
     &                 bounds( max(kplusp-nevd2+1,kplusp-np+1)), 1)
               end if
c
            else
c
c              %--------------------------------------------------%
c              | LM, SM, LA, SA case.                             |
c              | Sort the eigenvalues of H into the an order that |
c              | is opposite to WHICH, and apply the resulting    |
c              | order to BOUNDS.  The eigenvalues are sorted so  |
c              | that the wanted part are always within the first |
c              | NEV locations.                                   |
c              %--------------------------------------------------%
c
               if (which .eq. 'LM') wprime = 'SM'
               if (which .eq. 'SM') wprime = 'LM'
               if (which .eq. 'LA') wprime = 'SA'
               if (which .eq. 'SA') wprime = 'LA'
c
               call dsortr (wprime, .true., kplusp, ritz, bounds)
c
            end if
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, nev0
               temp = max( eps23, abs(ritz(j)) )
               bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | esitmates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LA'
            call dsortr(wprime, .true., nev0, bounds, ritz)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, nev0
                temp = max( eps23, abs(ritz(j)) )
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %--------------------------------------------------%
c           | Sort the "converged" Ritz values again so that   |
c           | the "threshold" values and their associated Ritz |
c           | estimates appear at the appropriate position in  |
c           | ritz and bound.                                  |
c           %--------------------------------------------------%
c
            if (which .eq. 'BE') then
c
c              %------------------------------------------------%
c              | Sort the "converged" Ritz values in increasing |
c              | order.  The "threshold" values are in the      |
c              | middle.                                        |
c              %------------------------------------------------%
c
               wprime = 'LA'
               call dsortr(wprime, .true., nconv, ritz, bounds)
c
            else
c
c              %----------------------------------------------%
c              | In LM, SM, LA, SA case, sort the "converged" |
c              | Ritz values according to WHICH so that the   |
c              | "threshold" value appears at the front of    |
c              | ritz.                                        |
c              %----------------------------------------------%
!
               call dsortr(which, .true., nconv, ritz, bounds)
c
            end if
c
c           %------------------------------------------%
c           |  Use h( 1,1 ) as storage to communicate  |
c           |  rnorm to _seupd if needed               |
c           %------------------------------------------%
c
            h(1,1) = rnorm
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. nev) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. nev0) info = 2
c
            np = nconv
            go to 1100
c
         else if (nconv .lt. nev .and. ishift .eq. 1) then
c
c           %---------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.    |
c           | To prevent possible stagnation, adjust the number |
c           | of Ritz values and the shifts.                    |
c           %---------------------------------------------------%
c
            nevbef = nev
            nev = nev + min (nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 2) then
               nev = 2
            end if
            np  = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call dsgets (ishift, which, nev, np, ritz, bounds,
     &              workl)
c
         end if
c
         if (ishift .eq. 0) then
c
c           %-----------------------------------------------------%
c           | User specified shifts: reverse communication to     |
c           | compute the shifts. They are returned in the first  |
c           | NP locations of WORKL.                              |
c           %-----------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
c
   50    continue
c
c        %------------------------------------%
c        | Back from reverse communication;   |
c        | User specified shifts are returned |
c        | in WORKL(1:*NP)                   |
c        %------------------------------------%
c
         ushift = .false.
c 
c 
c        %---------------------------------------------------------%
c        | Move the NP shifts to the first NP locations of RITZ to |
c        | free up WORKL.  This is for the non-exact shift case;   |
c        | in the exact shift case, dsgets already handles this.   |
c        %---------------------------------------------------------%
c
         if (ishift .eq. 0) call dcopy (np, workl, 1, ritz, 1)
c
c        %---------------------------------------------------------%
c        | Apply the NP0 implicit shifts by QR bulge chasing.      |
c        | Each shift is applied to the entire tridiagonal matrix. |
c        | The first 2*N locations of WORKD are used as workspace. |
c        | After dsapps is done, we have a Lanczos                 |
c        | factorization of length NEV.                            |
c        %---------------------------------------------------------%
c
         call dsapps (n, nev, np, ritz, v, ldv, h, ldh, resid, q, ldq,
     &        workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dsaitr.  |
c        %---------------------------------------------%
c
         cnorm = .true.
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd, 1)
         end if
c 
  100    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
         cnorm = .false.
  130    continue
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c 
 1100 continue
c
      mxiter = iter
      nev = nconv
c 
 1200 continue
      ido = 99
c
c     %------------%
c     | Error exit |
c     %------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of dsaup2 |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaupd 
c
c\Description: 
c
c  Reverse communication interface for the Implicitly Restarted Arnoldi 
c  Iteration.  For symmetric problems this reduces to a variant of the Lanczos 
c  method.  This method has been designed to compute approximations to a 
c  few eigenpairs of a linear operator OP that is real and symmetric 
c  with respect to a real positive semi-definite symmetric matrix B, 
c  i.e.
c                   
c       B*OP = (OP`)*B.  
c
c  Another way to express this condition is 
c
c       < x,OPy > = < OPx,y >  where < z,w > = z`Bw  .
c  
c  In the standard eigenproblem B is the identity matrix.  
c  ( A` denotes transpose of A)
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  dsaupd  is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x, A symmetric 
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
c           ===> OP = (inv[K - sigma*M])*M  and  B = M. 
c           ===> Shift-and-Invert mode
c
c  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite, 
c           KG symmetric indefinite
c           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
c           ===> Buckling mode
c
c  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
c           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
c           ===> Cayley transformed mode
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call dsaupd  
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to dsaupd .  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          dsaupd  with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          (If Mode = 2 see remark 5 below)
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3,4 and 5, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) shifts where
c                    IPNTR(11) is the pointer into WORKL for
c                    placing the shifts. See remark 6 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          Specify which of the Ritz values of OP to compute.
c
c          'LA' - compute the NEV largest (algebraic) eigenvalues.
c          'SA' - compute the NEV smallest (algebraic) eigenvalues.
c          'LM' - compute the NEV largest (in magnitude) eigenvalues.
c          'SM' - compute the NEV smallest (in magnitude) eigenvalues. 
c          'BE' - compute NEV eigenvalues, half from each end of the
c                 spectrum.  When NEV is odd, compute one more from the
c                 high end than from the low end.
c           (see remark 1 below)
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N.
c
c  TOL     Double precision  scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH ('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH ).
c
c  RESID   Double precision  array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector. 
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          This will indicate how many Lanczos vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Lanczos vectors are generated, the algorithm generates 
c          NCV-NEV Lanczos vectors at each subsequent update iteration.
c          Most of the cost in generating each Lanczos vector is in the 
c          matrix-vector product OP*x. (See remark 4 below).
c
c  V       Double precision  N by NCV array.  (OUTPUT)
c          The NCV columns of V contain the Lanczos basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The NCV eigenvalues of
c                      the current tridiagonal matrix T are returned in
c                      the part of WORKL array corresponding to RITZ.
c                      See remark 6 below.
c          ISHIFT = 1: exact shifts with respect to the reduced 
c                      tridiagonal matrix T.  This is equivalent to 
c                      restarting the iteration with a starting vector 
c                      that is a linear combination of Ritz vectors 
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = LEVEC
c          No longer referenced. See remark 2 below.
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used. 
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsaupd  for the 
c          five modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), dsaupd  returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          6 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.        
c
c  IPNTR   Integer array of length 11.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
c          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
c          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZ in WORKL.
c          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
c
c          Note: IPNTR(8:10) is only referenced by dseupd . See Remark 2.
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd  if RVEC = .TRUE. See Remarks.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision  work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
c          subroutine dseupd  uses this output.
c          See Data Distribution Note below.  
c
c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV .
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iterations allowed
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array WORKL is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informatinal error from LAPACK routine dsteqr .
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -13: NEV and WHICH = 'BE' are incompatable.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization. The user is advised to check that
c                   enough workspace and array storage has been allocated.
c
c
c\Remarks
c  1. The converged Ritz values are always returned in ascending 
c     algebraic order.  The computed Ritz values are approximate
c     eigenvalues of OP.  The selection of WHICH should be made
c     with this in mind when Mode = 3,4,5.  After convergence, 
c     approximate eigenvalues of the original problem may be obtained 
c     with the ARPACK subroutine dseupd . 
c
c  2. If the Ritz vectors corresponding to the converged Ritz values
c     are needed, the user must call dseupd  immediately following completion
c     of dsaupd . This is new starting with version 2.1 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular 
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requrement is that NCV > NEV.
c     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.   The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically.
c
c  5. If IPARAM(7) = 2 then in the Reverse commuication interface the user
c     must do the following. When IDO = 1, Y = OP * X is to be computed.
c     When IPARAM(7) = 2 OP = inv(B)*A. After computing A*X the user
c     must overwrite X with A*X. Y is then the solution to the linear set
c     of equations B*Y = A*X.
c
c  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
c     NP = IPARAM(8) shifts in locations: 
c     1   WORKL(IPNTR(11))           
c     2   WORKL(IPNTR(11)+1)         
c                        .           
c                        .           
c                        .      
c     NP  WORKL(IPNTR(11)+NP-1). 
c
c     The eigenvalues of the current tridiagonal matrix are located in 
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
c     order defined by WHICH. The associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note:
c
c  Fortran-D syntax:
c  ================
c  REAL       RESID(N), V(LDV,NCV), WORKD(3*N), WORKL(LWORKL)
c  DECOMPOSE  D1(N), D2(N,NCV)
c  ALIGN      RESID(I) with D1(I)
c  ALIGN      V(I,J)   with D2(I,J)
c  ALIGN      WORKD(I) with D1(I)     range (1:N)
c  ALIGN      WORKD(I) with D1(I-N)   range (N+1:2*N)
c  ALIGN      WORKD(I) with D1(I-2*N) range (2*N+1:3*N)
c  DISTRIBUTE D1(BLOCK), D2(BLOCK,:)
c  REPLICATED WORKL(LWORKL)
c
c  Cray MPP syntax:
c  ===============
c  REAL       RESID(N), V(LDV,NCV), WORKD(N,3), WORKL(LWORKL)
c  SHARED     RESID(BLOCK), V(BLOCK,:), WORKD(BLOCK,:)
c  REPLICATED WORKL(LWORKL)
c  
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c  8. R.B. Lehoucq, D.C. Sorensen, "Implementation of Some Spectral
c     Transformations in a k-Step Arnoldi Method". In Preparation.
c
c\Routines called:
c     dsaup2   ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     dlamch   LAPACK routine that determines machine constants.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/15/93: Version ' 2.4' 
c
c\SCCS Information: @(#) 
c FILE: saupd.F   SID: 2.8   DATE OF SID: 04/10/01   RELEASE: 2 
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsaupd 
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &     ipntr, workd, workl, lworkl, info )
c
c     %----------------------------------------------------%
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision 
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(11)
      Double precision 
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision 
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw, 
     &           ldh, ldq, mxiter, mode, nb,
     &           nev0, next, np, ritz, j
      save       bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, mxiter, mode, nb,
     &           nev0, next, np, ritz
!$omp threadprivate(bounds, ierr, ih, iq, ishift, iupd, iw)
!$omp threadprivate(ldh, ldq, mxiter, mode, nb)
!$omp threadprivate(nev0, next, np, ritz)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dsaup2
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dlamch 
      external   dlamch 
 
      NOPX=0; NBX=0; NRORTH=0
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
      if (ido .eq. 0) then
c
         ierr   = 0
         ishift = iparam(1)
         mxiter = iparam(3)
c         nb     = iparam(4)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev .or.  ncv .gt. n) then
            ierr = -3
         end if
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        %----------------------------------------------%
c
         np     = ncv - nev
c 
         if (mxiter .le. 0)                     ierr = -4
         if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LA' .and.
     &       which .ne. 'SA' .and.
     &       which .ne. 'BE')                   ierr = -5
         if (bmat .ne. 'I' .and. bmat .ne. 'G') ierr = -6
c
         if (lworkl .lt. ncv**2 + 8*ncv)        ierr = -7
         if (mode .lt. 1 .or. mode .gt. 5) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         else if (nev .eq. 1 .and. which .eq. 'BE') then
                                                ierr = -13
         end if
c 
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c 
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)                         nb = 1
         if (tol .le. zero)                     tol = dlamch ('EpsMach')
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev 
c 
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, ncv**2 + 8*ncv
            workl(j) = zero
 10      continue
c 
c        %-------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
c        | etc... and the remaining workspace.                   |
c        | Also update pointer to be used on output.             |
c        | Memory is laid out as follows:                        |
c        | workl(1:2*ncv) := generated tridiagonal matrix        |
c        | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
c        | workl(3*ncv+1:3*ncv+ncv) := computed error bounds     |
c        | workl(4*ncv+1:4*ncv+ncv*ncv) := rotation matrix Q     |
c        | workl(4*ncv+ncv*ncv+1:7*ncv+ncv*ncv) := workspace     |
c        %-------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritz   = ih     + 2*ldh
         bounds = ritz   + ncv
         iq     = bounds + ncv
         iw     = iq     + ncv**2
         next   = iw     + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritz
         ipntr(7) = bounds
         ipntr(11) = iw
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Lanczos Iteration. |
c     %-------------------------------------------------------%
c
      call dsaup2  
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritz),
     &     workl(bounds), workl(iq), ldq, workl(iw), ipntr, workd,
     &     info )
c
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c 
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within dsaup2 .               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
 9000 continue
c 
      return
c
c     %---------------%
c     | End of dsaupd  |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsconv
c
c\Description: 
c  Convergence testing for the symmetric Arnoldi eigenvalue routine.
c
c\Usage:
c  call dsconv
c     ( N, RITZ, BOUNDS, TOL, NCONV )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.
c
c  RITZ    Double precision array of length N.  (INPUT)
c          The Ritz values to be checked for convergence.
c
c  BOUNDS  Double precision array of length N.  (INPUT)
c          Ritz estimates associated with the Ritz values in RITZ.
c
c  TOL     Double precision scalar.  (INPUT)
c          Desired relative accuracy for a Ritz value to be considered
c          "converged".
c
c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines called:
c     dlamch  LAPACK routine that determines machine constants. 
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\SCCS Information: @(#) 
c FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
c
c\Remarks
c     1. Starting with version 2.4, this routine no longer uses the
c        Parlett strategy using the gap conditions. 
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsconv (n, ritz, bounds, tol, nconv)
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    n, nconv
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           ritz(n), bounds(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i
      Double precision
     &           temp, eps23
c
c     %-------------------%
c     | External routines |
c     %-------------------%
c
      Double precision
     &           dlamch
      external   dlamch
!
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      eps23 = dlamch('Epsilon-Machine') 
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
      nconv  = 0
      do 10 i = 1, n
c
c        %-----------------------------------------------------%
c        | The i-th Ritz value is considered "converged"       |
c        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
c        %-----------------------------------------------------%
c
         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) .le. tol*temp ) then
            nconv = nconv + 1
         end if
c
   10 continue
c 
      return
c
c     %---------------%
c     | End of dsconv |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dseigt
c
c\Description: 
c  Compute the eigenvalues of the current symmetric tridiagonal matrix
c  and the corresponding error bounds given the current residual norm.
c
c\Usage:
c  call dseigt
c     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          RNORM contains the residual norm corresponding to the current
c          symmetric tridiagonal matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the symmetric tridiagonal matrix H.
c
c  H       Double precision N by 2 array.  (INPUT)
c          H contains the symmetric tridiagonal matrix with the 
c          subdiagonal in the first column starting at H(2,1) and the 
c          main diagonal in 2nd column.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  EIG     Double precision array of length N.  (OUTPUT)
c          On output, EIG contains the N eigenvalues of H possibly 
c          unsorted.  The BOUNDS arrays are returned in the
c          same sorted order as EIG.
c
c  BOUNDS  Double precision array of length N.  (OUTPUT)
c          On output, BOUNDS contains the error estimates corresponding
c          to the eigenvalues EIG.  This is equal to RNORM times the
c          last components of the eigenvectors corresponding to the
c          eigenvalues in EIG.
c
c  WORKL   Double precision work array of length 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from dstqrb.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dstqrb  ARPACK routine that computes the eigenvalues and the
c             last components of the eigenvectors of a symmetric
c             and tridiagonal matrix.
c     dcopy   Level 1 BLAS that copies one vector to another.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dseigt 
     &   ( rnorm, n, h, ldh, eig, bounds, workl, ierr )
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    ierr, ldh, n
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           eig(n), bounds(n), h(ldh,2), workl(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           zero
      parameter (zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, k
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dstqrb
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------% 
c
      call dcopy  (n, h(1,2), 1, eig, 1)
      call dcopy  (n-1, h(2,1), 1, workl, 1)
      call dstqrb (n, eig, workl, bounds, workl(n+1), ierr)
      if (ierr .ne. 0) go to 9000
c     %-----------------------------------------------%
c     | Finally determine the error bounds associated |
c     | with the n Ritz values of H.                  |
c     %-----------------------------------------------%
c
      do 30 k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
   30 continue
c 
 9000 continue
      return
c
c     %---------------%
c     | End of dseigt |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsesrt
c
c\Description:
c  Sort the array X in the order specified by WHICH and optionally 
c  apply the permutation to the columns of the matrix A.
c
c\Usage:
c  call dsesrt
c     ( WHICH, APPLY, N, X, NA, A, LDA)
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X is sorted into increasing order of magnitude.
c          'SM' -> X is sorted into decreasing order of magnitude.
c          'LA' -> X is sorted into increasing order of algebraic.
c          'SA' -> X is sorted into decreasing order of algebraic.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to A.
c          APPLY = .FALSE. -> do not apply the sorted order to A.
c
c  N       Integer.  (INPUT)
c          Dimension of the array X.
c
c  X      Double precision array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.
c
c  NA      Integer.  (INPUT)
c          Number of rows of the matrix A.
c
c  A      Double precision array of length NA by N.  (INPUT/OUTPUT)
c         
c  LDA     Integer.  (INPUT)
c          Leading dimension of A.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines
c     dswap  Level 1 BLAS that swaps the contents of two vectors.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     12/15/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO and 
c               the ARPACK code dsortr
c
c\SCCS Information: @(#) 
c FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsesrt (which, apply, n, x, na, a, lda)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    lda, n, na
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           x(0:n-1), a(lda, 0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Double precision
     &           temp
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dswap
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c 
      if (which .eq. 'SA') then
c
c        X is sorted into decreasing order of algebraic.
c
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            if (x(j).lt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        X is sorted into decreasing order of magnitude.
c
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j.lt.0) go to 60
c
            if (abs(x(j)).lt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c
      else if (which .eq. 'LA') then
c
c        X is sorted into increasing order of algebraic.
c
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c           
            if (x(j).gt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c 
      else if (which .eq. 'LM') then
c
c        X is sorted into increasing order of magnitude.
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (abs(x(j)).gt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
c
 9000 continue
      return
c
c     %---------------%
c     | End of dsesrt |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: dseupd 
c
c\Description: 
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) the corresponding approximate eigenvectors,
c
c      (2) an orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace,
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost 
c  of n*nev if both are requested (in this case a separate array Z must be 
c  supplied).
c
c  These quantities are obtained from the Lanczos factorization computed
c  by DSAUPD  for the linear operator OP prescribed by the MODE selection
c  (see IPARAM(7) in DSAUPD  documentation.)  DSAUPD  must be called before
c  this routine is called. These approximate eigenvalues and vectors are 
c  commonly called Ritz values and Ritz vectors respectively.  They are 
c  referred to as such in the comments that follow.   The computed orthonormal 
c  basis for the invariant subspace corresponding to these Ritz values is 
c  referred to as a Lanczos basis.
c
c  See documentation in the header of the subroutine DSAUPD  for a definition 
c  of OP as well as other terms and the relation of computed Ritz values 
c  and vectors of OP with respect to the given problem  A*z = lambda*B*z.  
c
c  The approximate eigenvalues of the original problem are returned in
c  ascending algebraic order.  The user may elect to call this routine
c  once for each desired Ritz vector and store it peripherally if desired.
c  There is also the option of computing a selected set of these vectors
c  with a single call.
c
c\Usage:
c  call dseupd  
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
c       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether Ritz vectors corresponding to the Ritz value 
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors.
c
c  HOWMNY  Character*1  (INPUT) 
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute NEV Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' , SELECT is used as a workspace for
c          reordering the Ritz values.
c
c  D       Double precision  array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd  transformed to
c          those of the original eigensystem A*z = lambda*B*z. If 
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
c          as the those of A*z = lambda*B*z.
c
c  Z       Double precision  N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Arnoldi/Lanczos basis array V computed by DSAUPD .
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   Double precision   (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DSAUPD  that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DSEUPD  following the last call
c         to DSAUPD .  These arguments MUST NOT BE MODIFIED between
c         the the last call to DSAUPD  and the call to DSEUPD .
c
c  Two of these parameters (WORKL, INFO) are also output parameters:
c
c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:4*ncv) contains information obtained in
c          dsaupd .  They are not changed by dseupd .
c          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
c          untransformed Ritz values, the computed error estimates,
c          and the associated eigenvector matrix of H.
c
c          Note: IPNTR(8:10) contains the pointer into WORKL for addresses
c          of the above information computed by dseupd .
c          -------------------------------------------------------------
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd  if RVEC = .TRUE. See Remarks.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Information error from LAPACK routine dsteqr .
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -14: DSAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
c          = -16: HOWMNY = 'S' not yet implemented
c          = -17: DSEUPD  got a different count of the number of converged
c                 Ritz values than DSAUPD  got.  This indicates the user
c                 probably made an error in passing data from DSAUPD  to
c                 DSEUPD  or that the data was modified before entering 
c                 DSEUPD .
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c
c\Remarks
c  1. The converged Ritz values are always returned in increasing 
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it. 
c
c\Routines called:
c     dsesrt   ARPACK routine that sorts an array X, and applies the
c             corresponding permutation to a matrix A.
c     dsortr   dsortr   ARPACK sorting routine.
c     dgeqr2   LAPACK routine that computes the QR factorization of
c             a matrix.
c     dlacpy   LAPACK matrix copy routine.
c     dlamch   LAPACK routine that determines machine constants.
c     dorm2r   LAPACK routine that applies an orthogonal matrix in
c             factored form.
c     dsteqr   LAPACK routine that computes eigenvalues and eigenvectors
c             of a tridiagonal matrix.
c     dger     Level 2 BLAS rank one update to a matrix.
c     dcopy    Level 1 BLAS that copies one vector to another .
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dscal    Level 1 BLAS that scales a vector.
c     dswap    Level 1 BLAS that swaps the contents of two vectors.
!
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas
c     Dept. of Computational & 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/15/93: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: seupd.F   SID: 2.11   DATE OF SID: 04/10/01   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine dseupd (rvec  , howmny, select, d    ,
     &                   z     , ldz   , sigma , bmat ,
     &                   n     , which , nev   , tol  ,
     &                   resid , ncv   , v     , ldv  ,
     &                   iparam, ipntr , workd , workl,
     &                   lworkl, info )
c
c     %----------------------------------------------------%
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision      
     &           sigma, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Double precision 
     &           d(nev)     , resid(n)  , v(ldv,ncv),
     &           z(ldz, nev), workd(2*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision 
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds , ierr   , ih    , ihb   , ihd   ,
     &           iq     , iw     , j     , k     , ldh   ,
     &           ldq    , mode   , nconv , next  ,
     &           ritz   , irz    , ibd   , np    , ishift,
     &           leftptr, rghtptr, numcnv, jj
      Double precision 
     &           bnorm2 , rnorm, temp, temp1, eps23
      logical    reord
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dger   , dgeqr2 , dlacpy , dorm2r , dscal , 
     &           dsesrt , dsteqr , dswap  , dsortr 
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dnrm2 , dlamch 
      external   dnrm2 , dlamch 
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c     %--------------%
c     | Quick return |
c     %--------------%
c
      if (nconv .eq. 0) go to 9000
      ierr = 0
c
      if (nconv .le. 0)                        ierr = -14 
      if (n .le. 0)                            ierr = -1
      if (nev .le. 0)                          ierr = -2
      if (ncv .le. nev .or.  ncv .gt. n)       ierr = -3
      if (which .ne. 'LM' .and.
     &    which .ne. 'SM' .and.
     &    which .ne. 'LA' .and.
     &    which .ne. 'SA' .and.
     &    which .ne. 'BE')                     ierr = -5
      if (bmat .ne. 'I' .and. bmat .ne. 'G')   ierr = -6
      if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) 
     &                                         ierr = -15
      if (rvec .and. howmny .eq. 'S')           ierr = -16
c
      if (rvec .and. lworkl .lt. ncv**2+8*ncv) ierr = -7
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else if (mode .eq. 4 ) then
         type = 'BUCKLE'
      else if (mode .eq. 5 ) then
         type = 'CAYLEY'
      else 
                                               ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')     ierr = -11
      if (nev .eq. 1 .and. which .eq. 'BE')    ierr = -12
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c     
c     %-------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
c     | etc... and the remaining workspace.                   |
c     | Also update pointer to be used on output.             |
c     | Memory is laid out as follows:                        |
c     | workl(1:2*ncv) := generated tridiagonal matrix H      |
c     |       The subdiagonal is stored in workl(2:ncv).      |
c     |       The dead spot is workl(1) but upon exiting      |
c     |       dsaupd  stores the B-norm of the last residual   |
c     |       vector in workl(1). We use this !!!             |
c     | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
c     |       The wanted values are in the first NCONV spots. |
c     | workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates   |
c     |       The wanted values are in the first NCONV spots. |
c     | NOTE: workl(1:4*ncv) is set by dsaupd  and is not      |
c     |       modified by dseupd .                             |
c     %-------------------------------------------------------%
c
c     %-------------------------------------------------------%
c     | The following is used and set by dseupd .              |
c     | workl(4*ncv+1:4*ncv+ncv) := used as workspace during  |
c     |       computation of the eigenvectors of H. Stores    |
c     |       the diagonal of H. Upon EXIT contains the NCV   |
c     |       Ritz values of the original system. The first   |
c     |       NCONV spots have the wanted values. If MODE =   |
c     |       1 or 2 then will equal workl(2*ncv+1:3*ncv).    |
c     | workl(5*ncv+1:5*ncv+ncv) := used as workspace during  |
c     |       computation of the eigenvectors of H. Stores    |
c     |       the subdiagonal of H. Upon EXIT contains the    |
c     |       NCV corresponding Ritz estimates of the         |
c     |       original system. The first NCONV spots have the |
c     |       wanted values. If MODE = 1,2 then will equal    |
c     |       workl(3*ncv+1:4*ncv).                           |
c     | workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is  |
c     |       the eigenvector matrix for H as returned by     |
c     |       dsteqr . Not referenced if RVEC = .False.        |
c     |       Ordering follows that of workl(4*ncv+1:5*ncv)   |
c     | workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=         |
c     |       Workspace. Needed by dsteqr  and by dseupd .      |
c     | GRAND total of NCV*(NCV+8) locations.                 |
c     %-------------------------------------------------------%
c
c
      ih     = ipntr(5)
      ritz   = ipntr(6)
      bounds = ipntr(7)
      ldh    = ncv
      ldq    = ncv
      ihd    = bounds + ldh
      ihb    = ihd    + ldh
      iq     = ihb    + ldh
      iw     = iq     + ldh*ncv
      next   = iw     + 2*ncv
      ipntr(4)  = next
      ipntr(8)  = ihd
      ipntr(9)  = ihb
      ipntr(10) = iq
c
c     %----------------------------------------%
c     | irz points to the Ritz values computed |
c     |     by _seigt before exiting _saup2.   |
c     | ibd points to the Ritz estimates       |
c     |     computed by _seigt before exiting  |
c     |     _saup2.                            |
c     %----------------------------------------%
c
      irz = ipntr(11)+ncv
      ibd = irz+ncv
c
c
c     %---------------------------------%
c     | Set machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch ('Epsilon-Machine') 
      eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c     %---------------------------------------%
c     | RNORM is B-norm of the RESID(1:N).    |
c     | BNORM2 is the 2 norm of B*RESID(1:N). |
c     | Upon exit of dsaupd  WORKD(1:N) has    |
c     | B*RESID(1:N).                         |
c     %---------------------------------------%
c
      rnorm = workl(ih)
      if (bmat .eq. 'I') then
         bnorm2 = rnorm
      else if (bmat .eq. 'G') then
         bnorm2 = dnrm2 (n, workd, 1)
      end if
c
      if (rvec) then
c
         reord = .false.
c
c        %---------------------------------------------------%
c        | Use the temporary bounds array to store indices   |
c        | These will be used to mark the select array later |
c        %---------------------------------------------------%
c
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
c
c        %-------------------------------------%
c        | Select the wanted Ritz values.      |
c        | Sort the Ritz values so that the    |
c        | wanted ones appear at the tailing   |
c        | NEV positions of workl(irr) and     |
c        | workl(iri).  Move the corresponding |
c        | error estimates in workl(bound)     |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call dsgets (ishift, which       , nev          ,
     &                np    , workl(irz)  , workl(bounds),
     &                workl)
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23, abs(workl(irz+ncv-j)) )
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          workl(ibd+jj-1) .le. tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by _saupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the _saupd data.           |
c        %-----------------------------------------------------------%
c
         if (numcnv .ne. nconv) then
            info = -17
            go to 9000
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine _steqr to compute the eigenvalues and |
c        | eigenvectors of the final symmetric tridiagonal matrix H. |
c        | Initialize the eigenvector matrix Q to the identity.      |
c        %-----------------------------------------------------------%
c
         call dcopy (ncv-1, workl(ih+1), 1, workl(ihb), 1)
         call dcopy (ncv, workl(ih+ldh), 1, workl(ihd), 1)
c
         call dsteqr ('Identity', ncv, workl(ihd), workl(ihb),
     &                workl(iq) , ldq, workl(iw), ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c
         if (reord) then
c
c           %---------------------------------------------%
c           | Reordered the eigenvalues and eigenvectors  |
c           | computed by _steqr so that the "converged"  |
c           | eigenvalues appear in the first NCONV       |
c           | positions of workl(ihd), and the associated |
c           | eigenvectors appear in the first NCONV      |
c           | columns.                                    |
c           %---------------------------------------------%
c
            leftptr = 1
            rghtptr = ncv
c
            if (ncv .eq. 1) go to 30
c
 20         if (select(leftptr)) then
c
c              %-------------------------------------------%
c              | Search, from the left, for the first Ritz |
c              | value that has not converged.             |
c              %-------------------------------------------%
c
               leftptr = leftptr + 1
c
            else if ( .not. select(rghtptr)) then
c
c              %----------------------------------------------%
c              | Search, from the right, the first Ritz value |
c              | that has converged.                          |
c              %----------------------------------------------%
c
               rghtptr = rghtptr - 1
c
            else
c
c              %----------------------------------------------%
c              | Swap the Ritz value on the left that has not |
c              | converged with the Ritz value on the right   |
c              | that has converged.  Swap the associated     |
c              | eigenvector of the tridiagonal matrix H as   |
c              | well.                                        |
c              %----------------------------------------------%
c
               temp = workl(ihd+leftptr-1)
               workl(ihd+leftptr-1) = workl(ihd+rghtptr-1)
               workl(ihd+rghtptr-1) = temp
               call dcopy (ncv, workl(iq+ncv*(leftptr-1)), 1,
     &                    workl(iw), 1)
               call dcopy (ncv, workl(iq+ncv*(rghtptr-1)), 1,
     &                    workl(iq+ncv*(leftptr-1)), 1)
               call dcopy (ncv, workl(iw), 1,
     &                    workl(iq+ncv*(rghtptr-1)), 1)
               leftptr = leftptr + 1
               rghtptr = rghtptr - 1
c
            end if
c
            if (leftptr .lt. rghtptr) go to 20
c
 30      end if
c
c        %----------------------------------------%
c        | Load the converged Ritz values into D. |
c        %----------------------------------------%
c
         call dcopy (nconv, workl(ihd), 1, d, 1)
c
      else
c
c        %-----------------------------------------------------%
c        | Ritz vectors not required. Load Ritz values into D. |
c        %-----------------------------------------------------%
c
         call dcopy (nconv, workl(ritz), 1, d, 1)
         call dcopy (ncv, workl(ritz), 1, workl(ihd), 1)
c
      end if
c
c     %------------------------------------------------------------------%
c     | Transform the Ritz values and possibly vectors and corresponding |
c     | Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values |
c     | (and corresponding data) are returned in ascending order.        |
c     %------------------------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
c        %---------------------------------------------------------%
c        | Ascending sort of wanted Ritz values, vectors and error |
c        | bounds. Not necessary if only Ritz values are desired.  |
c        %---------------------------------------------------------%
c
         if (rvec) then
            call dsesrt ('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call dcopy (ncv, workl(bounds), 1, workl(ihb), 1)
         end if
c
      else 
c 
c        %-------------------------------------------------------------%
c        | *  Make a copy of all the Ritz values.                      |
c        | *  Transform the Ritz values back to the original system.   |
c        |    For TYPE = 'SHIFTI' the transformation is                |
c        |             lambda = 1/theta + sigma                        |
c        |    For TYPE = 'BUCKLE' the transformation is                |
c        |             lambda = sigma * theta / ( theta - 1 )          |
c        |    For TYPE = 'CAYLEY' the transformation is                |
c        |             lambda = sigma * (theta + 1) / (theta - 1 )     |
c        |    where the theta are the Ritz values returned by dsaupd .  |
c        | NOTES:                                                      |
c        | *The Ritz vectors are not affected by the transformation.   |
c        |  They are only reordered.                                   |
c        %-------------------------------------------------------------%
c
         call dcopy  (ncv, workl(ihd), 1, workl(iw), 1)
         if (type .eq. 'SHIFTI') then 
            do 40 k=1, ncv
               workl(ihd+k-1) = one / workl(ihd+k-1) + sigma
  40        continue
         else if (type .eq. 'BUCKLE') then
            do 50 k=1, ncv
               workl(ihd+k-1) = sigma * workl(ihd+k-1) / 
     &                          (workl(ihd+k-1) - one)
  50        continue
         else if (type .eq. 'CAYLEY') then
            do 60 k=1, ncv
               workl(ihd+k-1) = sigma * (workl(ihd+k-1) + one) /
     &                          (workl(ihd+k-1) - one)
  60        continue
         end if
c 
c        %-------------------------------------------------------------%
c        | *  Store the wanted NCONV lambda values into D.             |
c        | *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)   |
c        |    into ascending order and apply sort to the NCONV theta   |
c        |    values in the transformed system. We will need this to   |
c        |    compute Ritz estimates in the original system.           |
c        | *  Finally sort the lambda`s into ascending order and apply |
c        |    to Ritz vectors if wanted. Else just sort lambda`s into  |
c        |    ascending order.                                         |
c        | NOTES:                                                      |
c        | *workl(iw:iw+ncv-1) contain the theta ordered so that they  |
c        |  match the ordering of the lambda. We`ll use them again for |
c        |  Ritz vector purification.                                  |
c        %-------------------------------------------------------------%
c
         call dcopy (nconv, workl(ihd), 1, d, 1)
         call dsortr ('LA', .true., nconv, workl(ihd), workl(iw))
         if (rvec) then
            call dsesrt ('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call dcopy (ncv, workl(bounds), 1, workl(ihb), 1)
            call dscal (ncv, bnorm2/rnorm, workl(ihb), 1)
            call dsortr ('LA', .true., nconv, d, workl(ihb))
         end if
c
      end if 
c 
c     %------------------------------------------------%
c     | Compute the Ritz vectors. Transform the wanted |
c     | eigenvectors of the symmetric tridiagonal H by |
c     | the Lanczos basis matrix V.                    |
c     %------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A') then
c    
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(iq,ldq).                                |
c        %----------------------------------------------------------%
c     
         call dgeqr2 (ncv, nconv        , workl(iq) ,
     &                ldq, workl(iw+ncv), workl(ihb),
     &                ierr)
c
c        %--------------------------------------------------------%
c        | * Postmultiply V by Q.                                 |   
c        | * Copy the first NCONV columns of VQ into Z.           |
c        | The N by NCONV matrix Z is now a matrix representation |
c        | of the approximate invariant subspace associated with  |
c        | the Ritz values in workl(ihd).                         |
c        %--------------------------------------------------------%
c     
         call dorm2r ('Right', 'Notranspose', n        ,
     &                ncv    , nconv        , workl(iq),
     &                ldq    , workl(iw+ncv), v        ,
     &                ldv    , workd(n+1)   , ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
c
c        %-----------------------------------------------------%
c        | In order to compute the Ritz estimates for the Ritz |
c        | values in both systems, need the last row of the    |
c        | eigenvector matrix. Remember, it`s in factored form |
c        %-----------------------------------------------------%
c
         do 65 j = 1, ncv-1
            workl(ihb+j-1) = zero 
  65     continue
         workl(ihb+ncv-1) = one
         call dorm2r ('Left', 'Transpose'  , ncv       ,
     &                1     , nconv        , workl(iq) ,
     &                ldq   , workl(iw+ncv), workl(ihb),
     &                ncv   , temp         , ierr)
c
      else if (rvec .and. howmny .eq. 'S') then
c
c     Not yet implemented. See remark 2 above.
c
      end if
c
      if (type .eq. 'REGULR' .and. rvec) then
c
            do 70 j=1, ncv
               workl(ihb+j-1) = rnorm * abs( workl(ihb+j-1) )
 70         continue
c
      else if (type .ne. 'REGULR' .and. rvec) then
c
c        %-------------------------------------------------%
c        | *  Determine Ritz estimates of the theta.       |
c        |    If RVEC = .true. then compute Ritz estimates |
c        |               of the theta.                     |
c        |    If RVEC = .false. then copy Ritz estimates   |
c        |              as computed by dsaupd .             |
c        | *  Determine Ritz estimates of the lambda.      |
c        %-------------------------------------------------%
c
         call dscal  (ncv, bnorm2, workl(ihb), 1)
         if (type .eq. 'SHIFTI') then 
c
            do 80 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) ) 
     &                        / workl(iw+k-1)**2
 80         continue
c
         else if (type .eq. 'BUCKLE') then
c
            do 90 k=1, ncv
               workl(ihb+k-1) = sigma * abs( workl(ihb+k-1) )
     &                        / (workl(iw+k-1)-one )**2
 90         continue
c
         else if (type .eq. 'CAYLEY') then
c
            do 100 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1)
     &                        / workl(iw+k-1)*(workl(iw+k-1)-one) )
 100        continue
c
         end if
c
      end if
c
c     %-------------------------------------------------%
c     | Ritz vector purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 3,4,5. See reference 7               |
c     %-------------------------------------------------%
c
      if (rvec .and. (type .eq. 'SHIFTI' .or. type .eq. 'CAYLEY')) then
c
         do 110 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)
     &                  / workl(iw+k)
 110     continue
c
      else if (rvec .and. type .eq. 'BUCKLE') then
c
         do 120 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)
     &                  / (workl(iw+k)-one)
 120     continue
c
      end if 
c
      if (type .ne. 'REGULR')
     &   call dger  (n, nconv, one, resid, 1, workl(iw), 1, z, ldz)
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of dseupd |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsgets
c
c\Description: 
c  Given the eigenvalues of the symmetric tridiagonal matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of 
c  degree NP which filters out components of the unwanted eigenvectors 
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: This is called even in the case of user specified shifts in 
c  order to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call dsgets
c     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> KEV eigenvalues of largest magnitude are retained.
c          'SM' -> KEV eigenvalues of smallest magnitude are retained.
c          'LA' -> KEV eigenvalues of largest value are retained.
c          'SA' -> KEV eigenvalues of smallest value are retained.
c          'BE' -> KEV eigenvalues, half from each end of the spectrum.
c                  If KEV is odd, compute one more from the high end.
c
c  KEV      Integer.  (INPUT)
c          KEV+NP is the size of the matrix H.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be computed.
c
c  RITZ    Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          On INPUT, RITZ contains the eigenvalues of H.
c          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues 
c          are in the first NP locations and the wanted part is in 
c          the last KEV locations.  When exact shifts are selected, the
c          unwanted part corresponds to the shifts to be applied.
c
c  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTS  Double precision array of length NP.  (INPUT/OUTPUT)
c          On INPUT:  contains the user specified shifts if ISHIFT = 0.
c          On OUTPUT: contains the shifts sorted into decreasing order 
c          of magnitude with respect to the Ritz estimates contained in
c          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dsortr  ARPACK utility sorting routine.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     xx/xx/93: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsgets ( ishift, which, kev, np, ritz, bounds, shifts )
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      integer    ishift, kev, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           bounds(kev+np), ritz(kev+np), shifts(np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    kevd2
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dswap, dcopy, dsortr
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    max, min
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------%
c
      if (which .eq. 'BE') then
c
c        %-----------------------------------------------------%
c        | Both ends of the spectrum are requested.            |
c        | Sort the eigenvalues into algebraically increasing  |
c        | order first then swap high end of the spectrum next |
c        | to low end in appropriate locations.                |
c        | NOTE: when np < floor(kev/2) be careful not to swap |
c        | overlapping locations.                              |
c        %-----------------------------------------------------%
c
         call dsortr ('LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2 
         if ( kev .gt. 1 ) then
            call dswap ( min(kevd2,np), ritz, 1, 
     &                   ritz( max(kevd2,np)+1 ), 1)
            call dswap ( min(kevd2,np), bounds, 1, 
     &                   bounds( max(kevd2,np)+1 ), 1)
         end if
c
      else
c
c        %----------------------------------------------------%
c        | LM, SM, LA, SA case.                               |
c        | Sort the eigenvalues of H into the desired order   |
c        | and apply the resulting order to BOUNDS.           |
c        | The eigenvalues are sorted so that the wanted part |
c        | are always in the last KEV locations.               |
c        %----------------------------------------------------%
c
         call dsortr (which, .true., kev+np, ritz, bounds)
      end if
c
      if (ishift .eq. 1 .and. np .gt. 0) then
c     
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first.       |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when the shifts  |
c        | are applied in subroutine dsapps.                     |
c        %-------------------------------------------------------%
c     
         call dsortr ('SM', .true., np, bounds, ritz)
         call dcopy (np, ritz, 1, shifts, 1)
      end if
c 
      return
c
c     %---------------%
c     | End of dsgets |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsortc
c
c\Description:
c  Sorts the complex array in XREAL and XIMAG into the order 
c  specified by WHICH and optionally applies the permutation to the
c  real array Y. It is assumed that if an element of XIMAG is
c  nonzero, then its negative is also an element. In other words,
c  both members of a complex conjugate pair are to be sorted and the
c  pairs are kept adjacent to each other.
c
c\Usage:
c  call dsortc
c     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
c          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
c          'LR' -> sort XREAL into increasing order of algebraic.
c          'SR' -> sort XREAL into decreasing order of algebraic.
c          'LI' -> sort XIMAG into increasing order of magnitude.
c          'SI' -> sort XIMAG into decreasing order of magnitude.
c          NOTE: If an element of XIMAG is non-zero, then its negative
c                is also an element.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to array Y.
c          APPLY = .FALSE. -> do not apply the sorted order to array Y.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  XREAL,  Double precision array of length N.  (INPUT/OUTPUT)
c  XIMAG   Real and imaginary part of the array to be sorted.
c
c  Y       Double precision array of length N.  (INPUT/OUTPUT)
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#) 
c FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsortc (which, apply, n, xreal, ximag, y)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision     
     &           xreal(0:n-1), ximag(0:n-1), y(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Double precision     
     &           temp, temp1, temp2
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision     
     &           dlapy2
      external   dlapy2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c 
      if (which .eq. 'LM') then
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into increasing order of magnitude. |
c        %------------------------------------------------------%
c
   10    continue
         if (igap .eq. 0) go to 9000
c
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            temp1 = dlapy2(xreal(j),ximag(j))
            temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
c
            if (temp1.gt.temp2) then
                temp = xreal(j)
                xreal(j) = xreal(j+igap)
                xreal(j+igap) = temp
c
                temp = ximag(j)
                ximag(j) = ximag(j+igap)
                ximag(j+igap) = temp
c
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                go to 30
            end if
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------------%
c
   40    continue
         if (igap .eq. 0) go to 9000
c
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j .lt. 0) go to 60
c
            temp1 = dlapy2(xreal(j),ximag(j))
            temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
c
            if (temp1.lt.temp2) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c 
      else if (which .eq. 'LR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into increasing order of algebraic. |
c        %------------------------------------------------%
c
   70    continue
         if (igap .eq. 0) go to 9000
c
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (xreal(j).gt.xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c 
      else if (which .eq. 'SR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into decreasing order of algebraic. |
c        %------------------------------------------------%
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (xreal(j).lt.xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
c 
      else if (which .eq. 'LI') then
c
c        %------------------------------------------------%
c        | Sort XIMAG into increasing order of magnitude. |
c        %------------------------------------------------%
c
  130    continue
         if (igap .eq. 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
c
            if (j.lt.0) go to 150
c
            if (abs(ximag(j)).gt.abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 150
            endif
            j = j-igap
            go to 140
  150    continue
         igap = igap / 2
         go to 130
c 
      else if (which .eq. 'SI') then
c
c        %------------------------------------------------%
c        | Sort XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------%
c
  160    continue
         if (igap .eq. 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
c
            if (j.lt.0) go to 180
c
            if (abs(ximag(j)).lt.abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 180
            endif
            j = j-igap
            go to 170
  180    continue
         igap = igap / 2
         go to 160
      end if
c 
 9000 continue
      return
c
c     %---------------%
c     | End of dsortc |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsortr
c
c\Description:
c  Sort the array X1 in the order specified by WHICH and optionally 
c  applies the permutation to the array X2.
c
c\Usage:
c  call dsortr
c     ( WHICH, APPLY, N, X1, X2 )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X1 is sorted into increasing order of magnitude.
c          'SM' -> X1 is sorted into decreasing order of magnitude.
c          'LA' -> X1 is sorted into increasing order of algebraic.
c          'SA' -> X1 is sorted into decreasing order of algebraic.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to X2.
c          APPLY = .FALSE. -> do not apply the sorted order to X2.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  X1      Double precision array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.
c
c  X2      Double precision array of length N.  (INPUT/OUTPUT)
c          Only referenced if APPLY = .TRUE.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     12/16/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#) 
c FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsortr (which, apply, n, x1, x2)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           x1(0:n-1), x2(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Double precision
     &           temp
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c 
      if (which .eq. 'SA') then
c
c        X1 is sorted into decreasing order of algebraic.
c
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            if (x1(j).lt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        X1 is sorted into decreasing order of magnitude.
c
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j.lt.0) go to 60
c
            if (abs(x1(j)).lt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c
      else if (which .eq. 'LA') then
c
c        X1 is sorted into increasing order of algebraic.
c
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c           
            if (x1(j).gt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c 
      else if (which .eq. 'LM') then
c
c        X1 is sorted into increasing order of magnitude.
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (abs(x1(j)).gt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
c
 9000 continue
      return
c
c     %---------------%
c     | End of dsortr |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dstqrb
c
c\Description:
c  Computes all eigenvalues and the last component of the eigenvectors
c  of a symmetric tridiagonal matrix using the implicit QL or QR method.
c
c  This is mostly a modification of the LAPACK routine dsteqr.
c  See Remarks.
c
c\Usage:
c  call dstqrb
c     ( N, D, E, Z, WORK, INFO )
c
c\Arguments
c  N       Integer.  (INPUT)
c          The number of rows and columns in the matrix.  N >= 0.
c
c  D       Double precision array, dimension (N).  (INPUT/OUTPUT)
c          On entry, D contains the diagonal elements of the
c          tridiagonal matrix.
c          On exit, D contains the eigenvalues, in ascending order.
c          If an error exit is made, the eigenvalues are correct
c          for indices 1,2,...,INFO-1, but they are unordered and
c          may not be the smallest eigenvalues of the matrix.
c
c  E       Double precision array, dimension (N-1).  (INPUT/OUTPUT)
c          On entry, E contains the subdiagonal elements of the
c          tridiagonal matrix in positions 1 through N-1.
c          On exit, E has been destroyed.
c
c  Z       Double precision array, dimension (N).  (OUTPUT)
c          On exit, Z contains the last row of the orthonormal 
c          eigenvector matrix of the symmetric tridiagonal matrix.  
c          If an error exit is made, Z contains the last row of the
c          eigenvector matrix associated with the stored eigenvalues.
c
c  WORK    Double precision array, dimension (max(1,2*N-2)).  (WORKSPACE)
c          Workspace used in accumulating the transformation for 
c          computing the last components of the eigenvectors.
c
c  INFO    Integer.  (OUTPUT)
c          = 0:  normal return.
c          < 0:  if INFO = -i, the i-th argument had an illegal value.
c          > 0:  if INFO = +i, the i-th eigenvalue has not converged
c                              after a total of  30*N  iterations.
c
c\Remarks
c  1. None.
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.
c     lsame   LAPACK character comparison routine.
c     dlae2   LAPACK routine that computes the eigenvalues of a 2-by-2 
c             symmetric matrix.
c     dlaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric 
c             matrix.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanst  LAPACK routine that computes the norm of a matrix.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dlaset  LAPACK matrix initialization routine.
c     dlasr   LAPACK routine that applies an orthogonal transformation to 
c             a matrix.
c     dlasrt  LAPACK sorting routine.
c     dsteqr  LAPACK routine that computes eigenvalues and eigenvectors
c             of a symmetric tridiagonal matrix.
c     xerbla  LAPACK error handler routine.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\SCCS Information: @(#) 
c FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     1. Starting with version 2.5, this routine is a modified version
c        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
c        only commeted out and new lines inserted.
c        All lines commented out have "c$$$" at the beginning.
c        Note that the LAPACK version 1.0 subroutine SSTEQR contained
c        bugs. 
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dstqrb ( n, d, e, z, work, info )
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    info, n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           d( n ), e( n-1 ), z( n ), work( 2*n-2 )
c
c     .. parameters ..
      Double precision               
     &                   zero, one, two, three
      parameter          ( zero = 0.0D+0, one = 1.0D+0, 
     &                     two = 2.0D+0, three = 3.0D+0 )
      integer            maxit
      parameter          ( maxit = 30 )
c     ..
c     .. local scalars ..
      integer            i, icompz, ii, iscale, j, jtot, k, l, l1, lend,
     &                   lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1,
     &                   nm1, nmaxit
      Double precision               
     &                   anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2,
     &                   s, safmax, safmin, ssfmax, ssfmin, tst
c     ..
c     .. external functions ..
      logical            lsame
      Double precision
     &                   dlamch, dlanst, dlapy2
      external           lsame, dlamch, dlanst, dlapy2
c     ..
c     .. external subroutines ..
      external           dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr,
     &                   dlasrt, dswap, xerbla
c     ..
c     .. intrinsic functions ..
      intrinsic          abs, max, sign, sqrt
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
c
c$$$      IF( LSAME( COMPZ, 'N' ) ) THEN
c$$$         ICOMPZ = 0
c$$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
c$$$         ICOMPZ = 1
c$$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
c$$$         ICOMPZ = 2
c$$$      ELSE
c$$$         ICOMPZ = -1
c$$$      END IF
c$$$      IF( ICOMPZ.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
c$$$     $         N ) ) ) THEN
c$$$         INFO = -6
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'SSTEQR', -INFO )
c$$$         RETURN
c$$$      END IF
c
c    *** New starting with version 2.5 ***
c
      icompz = 2
c    *************************************
c
c     quick return if possible
c
      if( n.eq.0 )
     $   return
c
      if( n.eq.1 ) then
         if( icompz.eq.2 )  z( 1 ) = one
         return
      end if
c
c     determine the unit roundoff and over/underflow thresholds.
c
      eps = dlamch( 'e' )
      eps2 = eps**2
      safmin = dlamch( 's' )
      safmax = one / safmin
      ssfmax = sqrt( safmax ) / three
      ssfmin = sqrt( safmin ) / eps2
c
c     compute the eigenvalues and eigenvectors of the tridiagonal
c     matrix.
c
c$$      if( icompz.eq.2 )
c$$$     $   call dlaset( 'full', n, n, zero, one, z, ldz )
c
c     *** New starting with version 2.5 ***
c
      if ( icompz .eq. 2 ) then
         do 5 j = 1, n-1
            z(j) = zero
  5      continue
         z( n ) = one
      end if
c     *************************************
c
      nmaxit = n*maxit
      jtot = 0
c
c     determine where the matrix splits and choose ql or qr iteration
c     for each block, according to whether top or bottom diagonal
c     element is smaller.
c
      l1 = 1
      nm1 = n - 1
c
   10 continue
      if( l1.gt.n )
     $   go to 160
      if( l1.gt.1 )
     $   e( l1-1 ) = zero
      if( l1.le.nm1 ) then
         do 20 m = l1, nm1
            tst = abs( e( m ) )
            if( tst.eq.zero )
     $         go to 30
            if( tst.le.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+
     $          1 ) ) ) )*eps ) then
               e( m ) = zero
               go to 30
            end if
   20    continue
      end if
      m = n
c
   30 continue
      l = l1
      lsv = l
      lend = m
      lendsv = lend
      l1 = m + 1
      if( lend.eq.l )
     $   go to 10
c
c     scale submatrix in rows and columns l to lend
c
      anorm = dlanst( 'i', lend-l+1, d( l ), e( l ) )
      iscale = 0
      if( anorm.eq.zero )
     $   go to 10
      if( anorm.gt.ssfmax ) then
         iscale = 1
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n,
     $                info )
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n,
     $                info )
      else if( anorm.lt.ssfmin ) then
         iscale = 2
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n,
     $                info )
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n,
     $                info )
      end if
c
c     choose between ql and qr iteration
c
      if( abs( d( lend ) ).lt.abs( d( l ) ) ) then
         lend = lsv
         l = lendsv
      end if
c
      if( lend.gt.l ) then
c
c        ql iteration
c
c        look for small subdiagonal element.
c
   40    continue
         if( l.ne.lend ) then
            lendm1 = lend - 1
            do 50 m = l, lendm1
               tst = abs( e( m ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m+1 ) )+
     $             safmin )go to 60
   50       continue
         end if
c
         m = lend
c
   60    continue
         if( m.lt.lend )
     $      e( m ) = zero
         p = d( l )
         if( m.eq.l )
     $      go to 80
c
c        if remaining matrix is 2-by-2, use dlae2 or dlaev2
c        to compute its eigensystem.
c
         if( m.eq.l+1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
               work( l ) = c
               work( n-1+l ) = s
c$$$               call dlasr( 'r', 'v', 'b', n, 2, work( l ),
c$$$     $                     work( n-1+l ), z( 1, l ), ldz )
c
c              *** New starting with version 2.5 ***
c
               tst      = z(l+1)
               z(l+1) = c*tst - s*z(l)
               z(l)   = s*tst + c*z(l)
c              *************************************
            else
               call dlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
            end if
            d( l ) = rt1
            d( l+1 ) = rt2
            e( l ) = zero
            l = l + 2
            if( l.le.lend )
     $         go to 40
            go to 140
         end if
c
         if( jtot.eq.nmaxit )
     $      go to 140
         jtot = jtot + 1
c
c        form shift.
c
         g = ( d( l+1 )-p ) / ( two*e( l ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )
c
         s = one
         c = one
         p = zero
c
c        inner loop
c
         mm1 = m - 1
         do 70 i = mm1, l, -1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m-1 )
     $         e( i+1 ) = r
            g = d( i+1 ) - p
            r = ( d( i )-g )*s + two*c*b
            p = s*r
            d( i+1 ) = g + p
            g = c*r - b
c
c           if eigenvectors are desired, then save rotations.
c
            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = -s
            end if
c
   70    continue
c
c        if eigenvectors are desired, then apply saved rotations.
c
         if( icompz.gt.0 ) then
            mm = m - l + 1
c$$$            call dlasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ),
c$$$     $                  z( 1, l ), ldz )
c
c             *** New starting with version 2.5 ***
c
              call dlasr( 'r', 'v', 'b', 1, mm, work( l ), 
     &                    work( n-1+l ), z( l ), 1 )
c             *************************************                             
         end if
c
         d( l ) = d( l ) - p
         e( l ) = g
         go to 40
c
c        eigenvalue found.
c
   80    continue
         d( l ) = p
c
         l = l + 1
         if( l.le.lend )
     $      go to 40
         go to 140
c
      else
c
c        qr iteration
c
c        look for small superdiagonal element.
c
   90    continue
         if( l.ne.lend ) then
            lendp1 = lend + 1
            do 100 m = l, lendp1, -1
               tst = abs( e( m-1 ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m-1 ) )+
     $             safmin )go to 110
  100       continue
         end if
c
         m = lend
c
  110    continue
         if( m.gt.lend )
     $      e( m-1 ) = zero
         p = d( l )
         if( m.eq.l )
     $      go to 130
c
c        if remaining matrix is 2-by-2, use dlae2 or dlaev2
c        to compute its eigensystem.
c
         if( m.eq.l-1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
c$$$               work( m ) = c
c$$$               work( n-1+m ) = s
c$$$               call dlasr( 'r', 'v', 'f', n, 2, work( m ),
c$$$     $                     work( n-1+m ), z( 1, l-1 ), ldz )
c
c               *** New starting with version 2.5 ***
c
                tst      = z(l)
                z(l)   = c*tst - s*z(l-1)
                z(l-1) = s*tst + c*z(l-1)
c               ************************************* 
            else
               call dlae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
            end if
            d( l-1 ) = rt1
            d( l ) = rt2
            e( l-1 ) = zero
            l = l - 2
            if( l.ge.lend )
     $         go to 90
            go to 140
         end if
c
         if( jtot.eq.nmaxit )
     $      go to 140
         jtot = jtot + 1
c
c        form shift.
c
         g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l-1 ) / ( g+sign( r, g ) ) )
c
         s = one
         c = one
         p = zero
c
c        inner loop
c
         lm1 = l - 1
         do 120 i = m, lm1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m )
     $         e( i-1 ) = r
            g = d( i ) - p
            r = ( d( i+1 )-g )*s + two*c*b
            p = s*r
            d( i ) = g + p
            g = c*r - b
c
c           if eigenvectors are desired, then save rotations.
c
            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = s
            end if
c
  120    continue
c
c        if eigenvectors are desired, then apply saved rotations.
c
         if( icompz.gt.0 ) then
            mm = l - m + 1
c$$$            call dlasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ),
c$$$     $                  z( 1, m ), ldz )
c
c           *** New starting with version 2.5 ***
c
            call dlasr( 'r', 'v', 'f', 1, mm, work( m ), work( n-1+m ),
     &                  z( m ), 1 )
c           *************************************                             
         end if
c
         d( l ) = d( l ) - p
         e( lm1 ) = g
         go to 90
c
c        eigenvalue found.
c
  130    continue
         d( l ) = p
c
         l = l - 1
         if( l.ge.lend )
     $      go to 90
         go to 140
c
      end if
c
c     undo scaling if necessary
c
  140 continue
      if( iscale.eq.1 ) then
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1,
     $                d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv ),
     $                n, info )
      else if( iscale.eq.2 ) then
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1,
     $                d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv ),
     $                n, info )
      end if
c
c     check for no convergence to an eigenvalue after a total
c     of n*maxit iterations.
c
      if( jtot.lt.nmaxit )
     $   go to 10
      do 150 i = 1, n - 1
         if( e( i ).ne.zero )
     $      info = info + 1
  150 continue
      go to 190
c
c     order eigenvalues and eigenvectors.
c
  160 continue
      if( icompz.eq.0 ) then
c
c        use quick sort
c
         call dlasrt( 'i', n, d, info )
c
      else
c
c        use selection sort to minimize swaps of eigenvectors
c
         do 180 ii = 2, n
            i = ii - 1
            k = i
            p = d( i )
            do 170 j = ii, n
               if( d( j ).lt.p ) then
                  k = j
                  p = d( j )
               end if
  170       continue
            if( k.ne.i ) then
               d( k ) = d( i )
               d( i ) = p
c$$$               call dswap( n, z( 1, i ), 1, z( 1, k ), 1 )
c           *** New starting with version 2.5 ***
c
               p    = z(k)
               z(k) = z(i)
               z(i) = p
c           *************************************
            end if
  180    continue
      end if
c
  190 continue
      return
c
c     %---------------%
c     | End of dstqrb |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: zgetv0
c
c\Description: 
c  Generate a random initial residual vector for the Arnoldi process.
c  Force the residual vector to be in the range of the operator OP.  
c
c\Usage:
c  call zgetv0
c     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, 
c       IPNTR, WORKD, IERR )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to zgetv0.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B in the (generalized)
c          eigenvalue problem A*x = lambda*B*x.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  ITRY    Integer.  (INPUT)
c          ITRY counts the number of times that zgetv0 is called.  
c          It should be set to 1 on the initial call to zgetv0.
c
c  INITV   Logical variable.  (INPUT)
c          .TRUE.  => the initial residual vector is given in RESID.
c          .FALSE. => generate a random initial residual vector.
c
c  N       Integer.  (INPUT)
c          Dimension of the problem.
c
c  J       Integer.  (INPUT)
c          Index of the residual vector to be generated, with respect to
c          the Arnoldi process.  J > 1 in case of a "restart".
c
c  V       Complex*16 N by J array.  (INPUT)
c          The first J-1 columns of V contain the current Arnoldi basis
c          if this is a "restart".
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          Initial residual vector to be generated.  If RESID is 
c          provided, force RESID into the range of the operator OP.
c
c  RNORM   Double precision scalar.  (OUTPUT)
c          B-norm of the generated residual.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c
c  WORKD   Complex*16 work array of length 2*N.  (REVERSE COMMUNICATION).
c          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
c
c  IERR    Integer.  (OUTPUT)
c          =  0: Normal exit.
c          = -1: Cannot generate a nontrivial restarted residual vector
c                in the range of the operator OP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     zlarnv  LAPACK routine for generating a random vector. 
c     zgemv   Level 2 BLAS routine for matrix vector multiplication.
c     zcopy   Level 1 BLAS that copies one vector to another.
c     zdotc   Level 1 BLAS that computes the scalar product of two vectors.
c     dznrm2  Level 1 BLAS that computes the norm of a vector. 
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas            
c
c\SCCS Information: @(#)
c FILE: getv0.F   SID: 2.3   DATE OF SID: 08/27/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zgetv0 
     &   ( ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, 
     &     ipntr, workd, ierr )
c 
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      logical    initv
      integer    ido, ierr, itry, j, ldv, n
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Complex*16
     &           resid(n), v(ldv,j), workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      Double precision
     &           rzero
      parameter  (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0),
     &            rzero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    first, inits, orth
      integer    idist, iseed(4), iter, jj
      Double precision
     &           rnorm0
      Complex*16
     &           cnorm
      save       first, iseed, inits, iter, orth, rnorm0
!$omp threadprivate(first, iseed, inits, iter, orth, rnorm0)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zcopy, zgemv, zlarnv
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dznrm2, dlapy2
      Complex*16
     &           zdotc
      external   zdotc, dznrm2, dlapy2
c
c     %-----------------%
c     | Data Statements |
c     %-----------------%
c
      data       inits /.true./
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
      NOPX=0; NBX=0
c     %-----------------------------------%
c     | Initialize the seed of the LAPACK |
c     | random number generator           |
c     %-----------------------------------%
c
      if (inits) then
          iseed(1) = 1
          iseed(2) = 3
          iseed(3) = 5
          iseed(4) = 7
          inits = .false.
      end if
c
      if (ido .eq.  0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        %-------------------------------%
c
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
c
c        %-----------------------------------------------------%
c        | Possibly generate a random starting vector in RESID |
c        | Use a LAPACK random number generator used by the    |
c        | matrix generation routines.                         |
c        |    idist = 1: uniform (0,1)  distribution;          |
c        |    idist = 2: uniform (-1,1) distribution;          |
c        |    idist = 3: normal  (0,1)  distribution;          |
c        %-----------------------------------------------------%
c
         if (.not.initv) then
            idist = 2
            call zlarnv (idist, iseed, n, resid)
         end if
c 
c        %----------------------------------------------------------%
c        | Force the starting vector into the range of OP to handle |
c        | the generalized problem when B is possibly (singular).   |
c        %----------------------------------------------------------%
c
         if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call zcopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         end if
      end if
c 
c     %----------------------------------------%
c     | Back from computing B*(initial-vector) |
c     %----------------------------------------%
c
      if (first) go to 20
c
c     %-----------------------------------------------%
c     | Back from computing B*(orthogonalized-vector) |
c     %-----------------------------------------------%
c
      if (orth)  go to 40
c 
c     %------------------------------------------------------%
c     | Starting vector is now in the range of OP; r = OP*r; |
c     | Compute B-norm of starting vector.                   |
c     %------------------------------------------------------%
c
      first = .TRUE.
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call zcopy (n, workd(n+1), 1, resid, 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call zcopy (n, resid, 1, workd, 1)
      end if
c 
   20 continue
c
      first = .FALSE.
      if (bmat .eq. 'G') then
          cnorm  = zdotc (n, resid, 1, workd, 1)
          rnorm0 = sqrt(dlapy2(dble(cnorm),dimag(cnorm)))
      else if (bmat .eq. 'I') then
           rnorm0 = dznrm2(n, resid, 1)
      end if
      rnorm  = rnorm0
c
c     %---------------------------------------------%
c     | Exit if this is the very first Arnoldi step |
c     %---------------------------------------------%
c
      if (j .eq. 1) go to 50
c 
c     %----------------------------------------------------------------
c     | Otherwise need to B-orthogonalize the starting vector against |
c     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
c     | This is the case where an invariant subspace is encountered   |
c     | in the middle of the Arnoldi factorization.                   |
c     |                                                               |
c     |       s = V^{T}*B*r;   r = r - V*s;                           |
c     |                                                               |
c     | Stopping criteria used for iter. ref. is discussed in         |
c     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
c     %---------------------------------------------------------------%
c
      orth = .TRUE.
   30 continue
c
      call zgemv ('C', n, j-1, one, v, ldv, workd, 1, 
     &            zero, workd(n+1), 1)
      call zgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1, 
     &            one, resid, 1)
c 
c     %----------------------------------------------------------%
c     | Compute the B-norm of the orthogonalized starting vector |
c     %----------------------------------------------------------%
c
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call zcopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call zcopy (n, resid, 1, workd, 1)
      end if
c 
   40 continue
c
      if (bmat .eq. 'G') then
         cnorm = zdotc (n, resid, 1, workd, 1)
         rnorm = sqrt(dlapy2(dble(cnorm),dimag(cnorm)))
      else if (bmat .eq. 'I') then
         rnorm = dznrm2(n, resid, 1)
      end if
c
c     %--------------------------------------%
c     | Check for further orthogonalization. |
c     %--------------------------------------%
c
      if (rnorm .gt. 0.717*rnorm0) go to 50
c 
      iter = iter + 1
      if (iter .le. 1) then
c
c        %-----------------------------------%
c        | Perform iterative refinement step |
c        %-----------------------------------%
c
         rnorm0 = rnorm
         go to 30
      else
c
c        %------------------------------------%
c        | Iterative refinement step "failed" |
c        %------------------------------------%
c
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = rzero
         ierr = -1
      end if
c 
   50 continue
c
      ido = 99
c 
 9000 continue
      return
c
c     %---------------%
c     | End of zgetv0 |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: znaitr
c
c\Description: 
c  Reverse communication interface for applying NP additional steps to 
c  a K step nonsymmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in znaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call znaitr
c     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and do not need to be
c          recomputed in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.  See znaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current size of V and H.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  NB      Integer.  (INPUT)
c          Blocksize to be used in the recurrence.          
c          Only work for NB = 1 right now.  The goal is to have a 
c          program that implement both the block and non-block method.
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Double precision scalar.  (INPUT/OUTPUT)
c          B-norm of the starting residual on input.
c          B-norm of the updated residual r_{k+p} on output.
c
c  V       Complex*16 N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K 
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Complex*16 (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
c          H is used to store the generated upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not 
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On input, WORKD(1:N) = B*RESID and is used to save some 
c          computation at the first step.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of the spanning invariant subspace of OP found.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     zgetv0  ARPACK routine to generate the initial vector.
c     zlanhs  LAPACK routine that computes various norms of a matrix.
c     zlascl  LAPACK routine for careful scaling of a matrix.
c     dlabad  LAPACK routine for defining the underflow and overflow
c             limits.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     zgemv   Level 2 BLAS routine for matrix vector multiplication.
c     zaxpy   Level 1 BLAS that computes a vector triad.
c     zcopy   Level 1 BLAS that copies one vector to another .
c     zdotc   Level 1 BLAS that computes the scalar product of two vectors. 
c     zscal   Level 1 BLAS that scales a vector.
c     zdscal  Level 1 BLAS that scales a complex vector by a real number. 
c     dznrm2  Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c 
c\SCCS Information: @(#)
c FILE: naitr.F   SID: 2.3   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c  
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        ( At present tol is zero )
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in znaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        H(:,j) = w_{j};
c        H(j,j-1) = rnorm
c        rnorm = || r_(j) ||
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};   
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf 
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znaitr
     &   (ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh, 
     &    ipntr, workd, info)
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, nb, np
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Complex*16
     &           h(ldh,k+np), resid(n), v(ldv,k+np), workd(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      Double precision
     &           rone, rzero
      parameter (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0), 
     &           rone = 1.0D+0, rzero = 0.0D+0)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      Double precision
     &           rtemp(2)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    first, orth1, orth2, rstart, step3, step4
      integer    ierr, i, infol, ipj, irj, ivj, iter, itry, j, 
     &           jj
      Double precision            
     &           ovfl, smlnum, tst1, ulp, unfl, betaj,
     &           temp1, rnorm1, wnorm
      Complex*16
     &           cnorm
c
      save       first, orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, ovfl,
     &           betaj, rnorm1, smlnum, ulp, unfl, wnorm
!$omp threadprivate(first, orth1, orth2, rstart, step3, step4)
!$omp threadprivate(ierr, ipj, irj, ivj, iter, itry, j, ovfl)
!$omp threadprivate(betaj, rnorm1, smlnum, ulp, unfl, wnorm)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zaxpy, zcopy, zscal, zdscal, zgemv, zgetv0, 
     &           dlabad
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Complex*16
     &           zdotc 
      Double precision            
     &           dlamch,  dznrm2, zlanhs, dlapy2
      external   zdotc, dznrm2, zlanhs, dlamch, dlapy2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic  dimag, dble, max, sqrt 
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
      data       first / .true. /
c
      NOPX=0; NBX=0; NRORTH=0
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------%
c        | Set machine-dependent constants for the |
c        | the splitting and deflation criterion.  |
c        | If norm(H) <= sqrt(OVFL),               |
c        | overflow should not occur.              |
c        | REFERENCE: LAPACK subroutine zlahqr     |
c        %-----------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = dble(one / unfl)
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        %-------------------------------%
c
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
         j      = k + 1
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
c 
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true. when ....                        |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         zgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %-----------------------------%
c     | Else this is the first step |
c     %-----------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
!
 1000 continue
c
c        %---------------------------------------------------%
c        | STEP 1: Check if the B norm of j-th residual      |
c        | vector is zero. Equivalent to determine whether   |
c        | an exact j-step Arnoldi factorization is present. |
c        %---------------------------------------------------%
c
         betaj = rnorm
         if (rnorm .gt. rzero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted.         |
c           %---------------------------------------------%
c 
            betaj  = rzero
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
            call zgetv0 (ido, bmat, itry, .false., n, j, v, ldv, 
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               ido = 99
               go to 9000
            end if
c 
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call zcopy (n, resid, 1, v(1,j), 1)
         if ( rnorm .ge. unfl) then
             temp1 = rone / rnorm
             call zdscal (n, temp1, v(1,j), 1)
             call zdscal (n, temp1, workd(ipj), 1)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine zlascl               |
c            %-----------------------------------------%
c
             call zlascl ('General', i, i, rnorm, rone,
     &                    n, 1, v(1,j), n, infol)
             call zlascl ('General', i, i, rnorm, rone,  
     &                    n, 1, workd(ipj), n, infol)
         end if
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         nopx  = nopx + 1
         call zcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c 
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c 
         go to 9000 
   50    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
c        | if step3 = .true.                |
c        %----------------------------------%
c
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call zcopy (n, workd(irj), 1, resid, 1)
c 
c        %---------------------------------------%
c        | STEP 4:  Finish extending the Arnoldi |
c        |          factorization to length j.   |
c        %---------------------------------------%
c
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call zcopy (n, resid, 1, workd(ipj), 1)
         end if
   60    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
c        | if step4 = .true.                |
c        %----------------------------------%
c
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
         if (bmat .eq. 'G') then  
             cnorm = zdotc (n, resid, 1, workd(ipj), 1)
             wnorm = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
             wnorm = dznrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c 
         call zgemv ('C', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, h(1,j), 1)
c
c        %--------------------------------------%
c        | Orthogonalize r_{j} against V_{j}.   |
c        | RESID contains OP*v_{j}. See STEP 3. | 
c        %--------------------------------------%
c
         call zgemv ('N', n, j, -one, v, ldv, h(1,j), 1,
     &               one, resid, 1)
c
         if (j .gt. 1) h(j,j-1) = dcmplx(betaj, rzero)
c
         orth1 = .true.
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call zcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call zcopy (n, resid, 1, workd(ipj), 1)
         end if 
   70    continue
c 
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then         
            cnorm = zdotc (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
            rnorm = dznrm2(n, resid, 1)
         end if
c 
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        | The following test determines whether the sine of the     |
c        | angle between  OP*x and the computed residual is less     |
c        | than or equal to 0.717.                                   |
c        %-----------------------------------------------------------%
c
         if ( rnorm .gt. 0.717*wnorm ) go to 100
c
         iter  = 0
         nrorth = nrorth + 1
c 
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c 
   80    continue
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call zgemv ('C', n, j, one, v, ldv, workd(ipj), 1, 
     &               zero, workd(irj), 1)
c
c        %---------------------------------------------%
c        | Compute the correction to the residual:     |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
c        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
c        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
c        %---------------------------------------------%
c
         call zgemv ('N', n, j, -one, v, ldv, workd(irj), 1, 
     &               one, resid, 1)
         call zaxpy (j, one, workd(irj), 1, h(1,j), 1)
c 
         orth2 = .true.
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call zcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call zcopy (n, resid, 1, workd(ipj), 1)
         end if 
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c 
         if (bmat .eq. 'G') then         
             cnorm  = zdotc (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
             rnorm1 = dznrm2(n, resid, 1)
         end if
c 
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if ( rnorm1 .gt. 0.717*rnorm ) then
c
c           %---------------------------------------%
c           | No need for further refinement.       |
c           | The cosine of the angle between the   |
c           | corrected residual vector and the old |
c           | residual vector is greater than 0.717 |
c           | In other words the corrected residual |
c           | and the old residual vector share an  |
c           | angle of less than arcCOS(0.717)      |
c           %---------------------------------------%
c
            rnorm = rnorm1
c 
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required.     |
c           %-------------------------------------------%
c
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue 
            rnorm = rzero
         end if
c 
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c 
  100    continue
c 
         rstart = .false.
         orth2  = .false.
c 
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            ido = 99
            do 110 i = max(1,k), k+np-1
c     
c              %--------------------------------------------%
c              | Check for splitting and deflation.         |
c              | Use a standard test as in the QR algorithm |
c              | REFERENCE: LAPACK subroutine zlahqr        |
c              %--------------------------------------------%
c     
               tst1 = dlapy2(dble(h(i,i)),dimag(h(i,i)))
     &              + dlapy2(dble(h(i+1,i+1)), dimag(h(i+1,i+1)))
               if( tst1.eq.dble(zero) )
     &              tst1 = zlanhs( '1', k+np, h, ldh, workd(n+1) )
               if( dlapy2(dble(h(i+1,i)),dimag(h(i+1,i))) .le. 
     &                    max( ulp*tst1, smlnum ) ) 
     &             h(i+1,i) = zero
 110        continue
c     
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c 
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of znaitr |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: znapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP implicit shifts resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix which is the product of rotations
c  and reflections resulting from the NP bulge change sweeps.
c  The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call znapps
c     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, 
c       WORKL, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. size of matrix A.
c
c  KEV     Integer.  (INPUT/OUTPUT)
c          KEV+NP is the size of the input matrix H.
c          KEV is the size of the updated matrix HNEW. 
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFT   Complex*16 array of length NP.  (INPUT)
c          The shifts to be applied.
c
c  V       Complex*16 N by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, V contains the current KEV+NP Arnoldi vectors.
c          On OUTPUT, V contains the updated KEV Arnoldi vectors
c          in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Complex*16 (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, H contains the current KEV+NP by KEV+NP upper 
c          Hessenberg matrix of the Arnoldi factorization.
c          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
c          matrix in the KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          On INPUT, RESID contains the the residual vector r_{k+p}.
c          On OUTPUT, RESID is the update residual vector rnew_{k} 
c          in the first KEV locations.
c
c  Q       Complex*16 KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations and reflections
c          during the bulge chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex*16 work array of length (KEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  WORKD   Complex*16 work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     zlacpy  LAPACK matrix copy routine.
c     zlanhs  LAPACK routine that computes various norms of a matrix.
c     zlartg  LAPACK Givens rotation construction routine.
c     zlaset  LAPACK matrix initialization routine.
c     dlabad  LAPACK routine for defining the underflow and overflow
c             limits.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     zgemv   Level 2 BLAS routine for matrix vector multiplication.
c     zaxpy   Level 1 BLAS that computes a vector triad.
c     zcopy   Level 1 BLAS that copies one vector to another.
c     zscal   Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c\SCCS Information: @(#)
c FILE: napps.F   SID: 2.3   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the sublocks of
c     the Hessenberg matrix H and not just to the submatrix that it
c     comes from. Deflation as in LAPACK routine zlahqr (QR algorithm
c     for upper Hessenberg matrices ) is used.
c     Upon output, the subdiagonals of H are enforced to be non-negative
c     real numbers.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znapps
     &   ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, 
     &     workl, workd )
c
c     %----------------------------------------------------%
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16
     &           h(ldh,kev+np), resid(n), shift(np), 
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      Double precision
     &           rzero
      parameter (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0),
     &           rzero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, iend, istart, j, jj, kplusp
      logical    first
      Complex*16
     &           cdum, f, g, h11, h21, r, s, sigma, t
      Double precision             
     &           c,  ovfl, smlnum, ulp, unfl, tst1
      save       first, ovfl, smlnum, ulp, unfl 
!$omp threadprivate(first, ovfl, smlnum, ulp, unfl)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zaxpy, zcopy, zgemv, zscal, zlacpy, zlartg, 
     &           zlaset, dlabad
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision                 
     &           zlanhs, dlamch, dlapy2
      external   zlanhs, dlamch, dlapy2
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs, dimag, conjg, dcmplx, max, min, dble
c
c     %---------------------%
c     | Statement Functions |
c     %---------------------%
c
      Double precision     
     &           zabs1
      zabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine zlahqr           |
c        %-----------------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = dble(one / unfl)
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------%
c
      kplusp = kev + np 
c 
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
      call zlaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
      do 110 jj = 1, np
         sigma = shift(jj)
         istart = 1
   20    continue
c
         do 30 i = istart, kplusp-1
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine zlahqr    |
c           %----------------------------------------%
c
            tst1 = zabs1( h( i, i ) ) + zabs1( h( i+1, i+1 ) )
            if( tst1.eq.rzero )
     &         tst1 = zlanhs( '1', kplusp-jj+1, h, ldh, workl )
            if ( abs(dble(h(i+1,i))) 
     &           .le. max(ulp*tst1, smlnum) )  then
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        | or if the current block starts after the point |
c        | of compression since we'll discard this stuff  |
c        %------------------------------------------------%
c
         if ( istart .eq. iend .or. istart .gt. kev) go to 100
c
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         f = h11 - sigma
         g = h21
c 
         do 80 i = istart, iend-1
c
c           %------------------------------------------------------%
c           | Construct the plane rotation G to zero out the bulge |
c           %------------------------------------------------------%
c
            call zlartg (f, g, c, s, r)
            if (i .gt. istart) then
               h(i,i-1) = r
               h(i+1,i-1) = zero
            end if
c
c           %---------------------------------------------%
c           | Apply rotation to the left of H;  H <- G'*H |
c           %---------------------------------------------%
c
            do 50 j = i, kplusp
               t        =  c*h(i,j) + s*h(i+1,j)
               h(i+1,j) = -conjg(s)*h(i,j) + c*h(i+1,j)
               h(i,j)   = t   
   50       continue
c
c           %---------------------------------------------%
c           | Apply rotation to the right of H;  H <- H*G |
c           %---------------------------------------------%
c
            do 60 j = 1, min(i+2,iend)
               t        =  c*h(j,i) + conjg(s)*h(j,i+1)
               h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
               h(j,i)   = t   
   60       continue
c
c           %-----------------------------------------------------%
c           | Accumulate the rotation in the matrix Q;  Q <- Q*G' |
c           %-----------------------------------------------------%
c
            do 70 j = 1, min(i+jj, kplusp)
               t        =   c*q(j,i) + conjg(s)*q(j,i+1)
               q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
               q(j,i)   = t   
   70       continue
c
c           %---------------------------%
c           | Prepare for next rotation |
c           %---------------------------%
c
            if (i .lt. iend-1) then
               f = h(i+1,i)
               g = h(i+2,i)
            end if
   80    continue
c
c        %-------------------------------%
c        | Finished applying the shift.  |
c        %-------------------------------%
c 
  100    continue
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
         istart = iend + 1
         if (iend .lt. kplusp) go to 20
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
  110 continue
c
c     %---------------------------------------------------%
c     | Perform a similarity transformation that makes    |
c     | sure that the compressed H will have non-negative |
c     | real subdiagonal elements.                        |
c     %---------------------------------------------------%
c
      do 120 j=1,kev
         if ( dble( h(j+1,j) ) .lt. rzero .or.
     &        dimag( h(j+1,j) ) .ne. rzero ) then
            t = h(j+1,j) / dlapy2(dble(h(j+1,j)),dimag(h(j+1,j)))
            call zscal( kplusp-j+1, conjg(t), h(j+1,j), ldh )
            call zscal( min(j+2, kplusp), t, h(1,j+1), 1 )
            call zscal( min(j+np+1,kplusp), t, q(1,j+1), 1 )
            h(j+1,j) = dcmplx( dble( h(j+1,j) ), rzero )
         end if
  120 continue
c
      do 130 i = 1, kev
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine zlahqr.       |
c        | Note: Since the subdiagonals of the        |
c        | compressed H are nonnegative real numbers, |
c        | we take advantage of this.                 |
c        %--------------------------------------------%
c
         tst1 = zabs1( h( i, i ) ) + zabs1( h( i+1, i+1 ) )
         if( tst1 .eq. rzero )
     &       tst1 = zlanhs( '1', kev, h, ldh, workl )
         if( dble( h( i+1,i ) ) .le. max( ulp*tst1, smlnum ) ) 
     &       h(i+1,i) = zero
 130  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero, 
     &                workd(n+1), 1)
c 
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
      do 140 i = 1, kev
         call zgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call zcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call zlacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
c 
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zcopy (n, workd(n+1), 1, v(1,kev+1), 1)
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call zscal (n, q(kplusp,kev), resid, 1)
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zaxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
c
 9000 continue
      return
c
c     %---------------%
c     | End of znapps |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: znaup2 
c
c\Description: 
c  Intermediate level interface called by znaupd .
c
c\Usage:
c  call znaup2 
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, 
c       Q, LDQ, WORKL, IPNTR, WORKD, RWORK, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in znaupd .
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in znaupd .
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during
c          each Arnoldi iteration.
c          If ISHIFT=1, NP is adjusted dynamically at each iteration
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV since a leading block of the current
c          upper Hessenberg matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Complex*16  N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV 
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Complex*16  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZ    Complex*16  array of length NEV+NP.  (OUTPUT)
c          RITZ(1:NEV)  contains the computed Ritz values of OP.
c
c  BOUNDS  Complex*16  array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to 
c          the computed Ritz values.
c          
c  Q       Complex*16  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex*16  work array of length at least 
c          (NEV+NP)**2 + 3*(NEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Complex*16  work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in ZNAUPD .
c
c  RWORK   Double precision    work array of length  NEV+NP ( WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.  
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16 
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     zgetv0   ARPACK initial vector generation routine. 
c     znaitr   ARPACK Arnoldi factorization routine.
c     znapps   ARPACK application of implicit shifts routine.
c     zneigh   ARPACK compute Ritz values and error bounds routine. 
c     zngets   ARPACK reorder Ritz values and error bounds routine.
c     zsortc   ARPACK sorting routine.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     zcopy    Level 1 BLAS that copies one vector to another .
c     zdotc    Level 1 BLAS that computes the scalar product of two vectors. 
c     zswap    Level 1 BLAS that swaps two vectors.
c     dznrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice Universitya
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c 
c\SCCS Information: @(#)
c FILE: naup2.F   SID: 2.6   DATE OF SID: 06/01/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znaup2 
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, 
     &     ishift, mxiter, v, ldv, h, ldh, ritz, bounds, 
     &     q, ldq, workl, ipntr, workd, rwork, info )
c
c     %----------------------------------------------------%
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision   
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Complex*16 
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), 
     &           resid(n), ritz(nev+np),  v(ldv,nev+np), 
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
       Double precision   
     &           rwork(nev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16 
     &           one, zero
      Double precision 
     &           rzero
      parameter (one = (1.0D+0, 0.0D+0) , zero = (0.0D+0, 0.0D+0) ,
     &           rzero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    cnorm , getv0, initv , update, ushift
      integer    ierr  , iter , kplusp, nconv, 
     &           nevbef, nev0 , np0   , nptemp, i    ,
     &           j    
      Complex*16 
     &           cmpnorm
      Double precision 
     &           rnorm , eps23, rtemp
      character  wprime*2
c
      save       cnorm,  getv0, initv , update, ushift, 
     &           rnorm,  iter , kplusp, nconv ,
     &           nevbef, nev0 , np0   , eps23
!$omp threadprivate(cnorm,  getv0, initv , update, ushift)
!$omp threadprivate(rnorm,  iter , kplusp, nconv)
!$omp threadprivate(nevbef, nev0 , np0   , eps23)
c
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(3)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zcopy , zgetv0 , znaitr , zneigh , zngets , znapps ,
     &           zsortc , zswap 
c
c     %--------------------%
c     | External functions |
c     %--------------------%
c
      Complex*16 
     &           zdotc 
      Double precision   
     &           dznrm2 , dlamch , dlapy2 
      external   zdotc , dznrm2 , dlamch , dlapy2 
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic  dimag , dble , min, max
c
      NBX=0
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c 
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvalues.     |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 0
c 
c        %---------------------------------%
c        | Get machine dependent constant. |
c        %---------------------------------%
c
         eps23 = dlamch ('Epsilon-Machine')
         eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c 
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call zgetv0  (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. rzero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c 
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call znaitr  (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, 
     &             h, ldh, ipntr, workd, info)
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c 
 1000 continue
c
         iter = iter + 1
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine znapps .                  |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call znaitr (ido, bmat, n, nev, np,    mode,  resid, rnorm,
     &               v  , ldv , h, ldh, ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         call zneigh  (rnorm, kplusp, h, ldh, ritz, bounds,
     &                q, ldq, workl, rwork,  ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZ,    |
c        | and BOUNDS respectively.                          | 
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
c
c        %--------------------------------------------------%
c        | Make a copy of Ritz values and the corresponding |
c        | Ritz estimates obtained from zneigh .             |
c        %--------------------------------------------------%
c
         call zcopy (kplusp,ritz,1,workl(kplusp**2+1),1)
         call zcopy (kplusp,bounds,1,workl(kplusp**2+kplusp+1),1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | bounds are in the last NEV loc. of RITZ           |
c        | BOUNDS respectively.                              |
c        %---------------------------------------------------%
c
         call zngets  (ishift, which, nev, np, ritz, bounds)
c 
c        %------------------------------------------------------------%
c        | Convergence test: currently we use the following criteria. |
c        | The relative accuracy of a Ritz value is considered        |
c        | acceptable if:                                             |
c        |                                                            |
c        | error_bounds(i) .le. tol*max(eps23, magnitude_of_ritz(i)). |
c        |                                                            |
c        %------------------------------------------------------------%
c
         nconv  = 0
c
         do 25 i = 1, nev
            rtemp = max( eps23, dlapy2 ( dble (ritz(np+i)),
     &                                  dimag (ritz(np+i)) ) ) 
            if ( dlapy2 (dble (bounds(np+i)),dimag (bounds(np+i))) 
     &                 .le. tol*rtemp ) then
               nconv = nconv + 1
            end if
   25    continue
c 
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c     
         if ( (nconv .ge. NEV0_FIX) .or.
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to zneupd  if needed               |
c           %------------------------------------------%
!
            h(3,1) = dcmplx (rnorm,rzero)
c
c           %----------------------------------------------%
c           | Sort Ritz values so that converged Ritz      |
c           | values appear within the first NEV locations |
c           | of ritz and bounds, and the most desired one |
c           | appears at the front.                        |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call zsortc (wprime, .true., kplusp, ritz, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23, magnitude of the Ritz value).  |
c           %--------------------------------------------------%
c
            do 35 j = 1, nev0 
                rtemp = max( eps23, dlapy2 ( dble (ritz(j)),
     &                                       dimag (ritz(j)) ) )
                bounds(j) = bounds(j)/rtemp
 35         continue
c
c           %---------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz |
c           | estimates.  This will push all the converged ones |
c           | towards the front of ritz, bounds (in the case    |
c           | when NCONV < NEV.)                                |
c           %---------------------------------------------------%
c
            wprime = 'LM'
            call zsortc (wprime, .true., nev0, bounds, ritz)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, nev0
                rtemp = max( eps23, dlapy2 ( dble (ritz(j)),
     &                                       dimag (ritz(j)) ) )
                bounds(j) = bounds(j)*rtemp
 40         continue
c
c           %-----------------------------------------------%
c           | Sort the converged Ritz values again so that  |
c           | the "threshold" value appears at the front of |
c           | ritz and bound.                               |
c           %-----------------------------------------------%
c
            call zsortc (which, .true., nconv, ritz, bounds)
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. nev0) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. nev0)  info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. nev0) .and. (ishift .eq. 1) ) then
c     
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call zngets  (ishift, which, nev, np, ritz, bounds)
c
         end if              
         if (ishift .eq. 0) then
c
c           %-------------------------------------------------------%
c           | User specified shifts: pop back out to get the shifts |
c           | and return them in the first 2*NP locations of WORKL. |
c           %-------------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
   50    continue
         ushift = .false.
c
         if ( ishift .ne. 1 ) then
c 
c            %----------------------------------%
c            | Move the NP shifts from WORKL to |
c            | RITZ, to free up WORKL           |
c            | for non-exact shift case.        |
c            %----------------------------------%
c
             call zcopy  (np, workl, 1, ritz, 1)
         end if
c
c        %---------------------------------------------------------%
c        | Apply the NP implicit shifts by QR bulge chasing.       |
c        | Each shift is applied to the whole upper Hessenberg     |
c        | matrix H.                                               |
c        | The first 2*N locations of WORKD are used as workspace. |
c        %---------------------------------------------------------%
c
         call znapps  (n, nev, np, ritz, v, ldv, 
     &                h, ldh, resid, q, ldq, workl, workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to znaitr .  |
c        %---------------------------------------------%
c
         cnorm = .true.
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call zcopy  (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call zcopy  (n, resid, 1, workd, 1)
         end if
c 
  100    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then         
            cmpnorm = zdotc  (n, resid, 1, workd, 1)
            rnorm = sqrt(dlapy2 (dble (cmpnorm),dimag (cmpnorm)))
         else if (bmat .eq. 'I') then
            rnorm = dznrm2 (n, resid, 1)
         end if
         cnorm = .false.
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = nconv
c     
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
 9000 continue
c
c     %---------------%
c     | End of znaup2  |
c     %---------------%
c
      return
      end
c\BeginDoc
c
c\Name: znaupd
c
c\Description: 
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This is intended to be used to find a few eigenpairs of a 
c  complex linear operator OP with respect to a semi-inner product defined 
c  by a hermitian positive semi-definite real matrix B. B may be the identity 
c  matrix.  NOTE: if both OP and B are real, then dsaupd or dnaupd should
c  be used.
c
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  znaupd is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M hermitian positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M hermitian semi-definite
c           ===> OP =  inv[A - sigma*M]*M   and  B = M. 
c           ===> shift-and-invert mode 
c           If OP*x = amu*x, then lambda = sigma + 1/amu.
c  
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call znaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to znaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          znaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = M * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute and return the shifts in the first 
c                    NP locations of WORKL.
c          IDO = 99: done
c          -------------------------------------------------------------
c          After the initialization phase, when the routine is used in 
c          the "shift-and-invert" mode, the vector M * X is already 
c          available and does not need to be recomputed in forming OP*X.
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Double precision  scalar.  (INPUT)
c          Stopping criteria: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = dlamch('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine dlamch).
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 1 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Arnoldi vectors are generated, the algorithm generates 
c          approximately NCV-NEV Arnoldi vectors at each subsequent update 
c          iteration. Most of the cost in generating each Arnoldi vector is 
c          in the matrix-vector operation OP*x. (See remark 4 below.)
c
c  V       Complex*16 array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to filter out
c          the components of the unwanted eigenvector.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are to be provided by the user via
c                      reverse communication.  The NCV eigenvalues of 
c                      the Hessenberg matrix H are returned in the part
c                      of WORKL array corresponding to RITZ.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to 
c                      restarting the iteration from the beginning 
c                      after updating the starting vector with a linear
c                      combination of Ritz vectors associated with the 
c                      "wanted" eigenvalues.
c          ISHIFT = 2: other choice of internal shift to be defined.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced 
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.  
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3; See under \Description of znaupd for the 
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), _naupd returns NP, the number
c          of shifts the user is to provide. 0 < NP < NCV-NEV.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg
c                    matrix H in WORKL.
c          IPNTR(6): pointer to the  ritz value array  RITZ
c          IPNTR(7): pointer to the (projected) ritz vector array Q
c          IPNTR(8): pointer to the error BOUNDS array in WORKL.
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by zneupd. See Remark 2 below.
c
c          IPNTR(9): pointer to the NCV RITZ values of the 
c                    original system.
c          IPNTR(10): Not Used
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     zneupd if RVEC = .TRUE. See Remark 2 below.
c
c          -------------------------------------------------------------
c          
c  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note below.  
c
c  WORKL   Complex*16 work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 5*NCV.
c
c  RWORK   Double precision  work array of length NCV (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 1 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   User input error highly likely.  Please
c                   check actual array dimensions and layout.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when using
c     Mode = 3.  When operating in Mode = 3 setting WHICH = 'LM' will
c     compute the NEV eigenvalues of the original problem that are
c     closest to the shift SIGMA . After convergence, approximate eigenvalues 
c     of the original problem may be obtained with the ARPACK subroutine zneupd.
c
c  2. If a basis for the invariant subspace corresponding to the converged Ritz 
c     values is needed, the user must call zneupd immediately following 
c     completion of znaupd. This is new starting with release 2 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular 
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requirement is that NCV > NEV + 1.
c     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.  The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically. 
c     See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
c     NP = IPARAM(8) complex shifts in locations
c     WORKL(IPNTR(14)), WORKL(IPNTR(14)+1), ... , WORKL(IPNTR(14)+NP).
c     Eigenvalues of the current upper Hessenberg matrix are located in
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are ordered
c     according to the order defined by WHICH.  The associated Ritz estimates
c     are located in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... ,
c     WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note: 
c
c  Fortran-D syntax:
c  ================
c  Complex*16 resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Complex*16 resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c  
c  CM2/CM5 syntax:
c  ==============
c  
c-----------------------------------------------------------------------
c
c     include   'ex-nonsym.doc'
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "_Complex_ Shift and Invert Strategies for
c     Double precision Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     znaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c 
c\SCCS Information: @(#)
c FILE: naupd.F   SID: 2.9   DATE OF SID: 07/21/02   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &     ipntr, workd, workl, lworkl, rwork, info )
c
c     %----------------------------------------------------%
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision 
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      Complex*16
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      Double precision  
     &           rwork(ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      parameter (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0))
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw, 
     &           ldh, ldq, levec, mode, mxiter, nb,
     &           nev0, next, np, ritz, j
      save       bounds, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, levec, mode, mxiter, nb,
     &           nev0, next, np, ritz
!$omp threadprivate(bounds, ih, iq, ishift, iupd, iw)
!$omp threadprivate(ldh, ldq, levec, mode, mxiter, nb)
!$omp threadprivate(nev0, next, np, ritz)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   znaup2
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dlamch
      external   dlamch
c
      NOPX=0; NBX=0; NRORTH=0
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        %-------------------------------%
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         ierr   = 0
         ishift = iparam(1)
c         levec  = iparam(2)
         mxiter = iparam(3)
c         nb     = iparam(4)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
         if (n .le. 0) then
             ierr = -1
         else if (nev .le. 0) then
             ierr = -2
         else if (ncv .le. nev .or.  ncv .gt. n) then
             ierr = -3
         else if (mxiter .le. 0) then
             ierr = -4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 5*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 3) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         end if
c 
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c 
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)				nb = 1
         if (tol .le. 0.0D+0 )			tol = dlamch('EpsMach')
         if (ishift .ne. 0  .and.  
     &       ishift .ne. 1  .and.
     &       ishift .ne. 2) 			ishift = 1
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev 
c 
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, 3*ncv**2 + 5*ncv
            workl(j) = zero
  10     continue
c 
c        %-------------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
c        | etc... and the remaining workspace.                         |
c        | Also update pointer to be used on output.                   |
c        | Memory is laid out as follows:                              |
c        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
c        | workl(ncv*ncv+1:ncv*ncv+ncv) := the ritz values             |
c        | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv)   := error bounds        |
c        | workl(ncv*ncv+2*ncv+1:2*ncv*ncv+2*ncv) := rotation matrix Q |
c        | workl(2*ncv*ncv+2*ncv+1:3*ncv*ncv+5*ncv) := workspace       |
c        | The final workspace is needed by subroutine zneigh called   |
c        | by znaup2. Subroutine zneigh calls LAPACK routines for      |
c        | calculating eigenvalues and the last row of the eigenvector |
c        | matrix.                                                     |
c        %-------------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritz   = ih     + ldh*ncv
         bounds = ritz   + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritz
         ipntr(7) = iq
         ipntr(8) = bounds
         ipntr(14) = iw
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
      call znaup2 
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritz), 
     &     workl(bounds), workl(iq), ldq, workl(iw), 
     &     ipntr, workd, rwork, info )
c 
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP.              |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c 
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within znaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of znaupd |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: zneigh
c
c\Description:
c  Compute the eigenvalues of the current upper Hessenberg matrix
c  and the corresponding Ritz estimates given the current residual norm.
c
c\Usage:
c  call zneigh
c     ( RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          Residual norm corresponding to the current upper Hessenberg 
c          matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the matrix H.
c
c  H       Complex*16 N by N array.  (INPUT)
c          H contains the current upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZ    Complex*16 array of length N.  (OUTPUT)
c          On output, RITZ(1:N) contains the eigenvalues of H.
c
c  BOUNDS  Complex*16 array of length N.  (OUTPUT)
c          On output, BOUNDS contains the Ritz estimates associated with
c          the eigenvalues held in RITZ.  This is equal to RNORM 
c          times the last components of the eigenvectors corresponding 
c          to the eigenvalues in RITZ.
c
c  Q       Complex*16 N by N array.  (WORKSPACE)
c          Workspace needed to store the eigenvectors of H.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex*16 work array of length N**2 + 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  This is needed to keep the full Schur form
c          of H and also in the calculation of the eigenvectors of H.
c
c  RWORK   Double precision  work array of length N (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end. 
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from zlahqr or ztrevc.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\Routines called:
c     zlacpy  LAPACK matrix copy routine.
c     zlahqr  LAPACK routine to compute the Schur form of an
c             upper Hessenberg matrix.
c     zlaset  LAPACK matrix initialization routine.
c     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper triangular form
c     zcopy   Level 1 BLAS that copies one vector to another. 
c     zdscal  Level 1 BLAS that scales a complex vector by a real number.
c     dznrm2  Level 1 BLAS that computes the norm of a vector.
c     
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c\SCCS Information: @(#)
c FILE: neigh.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zneigh (rnorm, n, h, ldh, ritz, bounds, 
     &                   q, ldq, workl, rwork, ierr)
c
c     %----------------------------------------------------%
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    ierr, n, ldh, ldq
      Double precision     
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16     
     &           bounds(n), h(ldh,n), q(ldq,n), ritz(n),
     &           workl(n*(n+3)) 
      Double precision 
     &           rwork(n)
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16     
     &           one, zero
      Double precision
     &           rone
      parameter  (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0),
     &           rone = 1.0D+0)
c 
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    select(1)
      integer    j
      Complex*16     
     &           vl(1)
      Double precision
     &           temp
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zlacpy, zlahqr, ztrevc, zcopy, 
     &           zdscal
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dznrm2
      external   dznrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------%
c
c     %----------------------------------------------------------%
c     | 1. Compute the eigenvalues, the last components of the   |
c     |    corresponding Schur vectors and the full Schur form T |
c     |    of the current upper Hessenberg matrix H.             |
c     |    zlahqr returns the full Schur form of H               | 
c     |    in WORKL(1:N**2), and the Schur vectors in q.         |
c     %----------------------------------------------------------%
c
      call zlacpy ('All', n, n, h, ldh, workl, n)
      call zlaset ('All', n, n, zero, one, q, ldq)
      call zlahqr (.true., .true., n, 1, n, workl, ldh, ritz,
     &             1, n, q, ldq, ierr)
      if (ierr .ne. 0) go to 9000
c
      call zcopy (n, q(n-1,1), ldq, bounds, 1)
c     %----------------------------------------------------------%
c     | 2. Compute the eigenvectors of the full Schur form T and |
c     |    apply the Schur vectors to get the corresponding      |
c     |    eigenvectors.                                         |
c     %----------------------------------------------------------%
c
      call ztrevc ('Right', 'Back', select, n, workl, n, vl, n, q, 
     &             ldq, n, n, workl(n*n+1), rwork, ierr)
c
      if (ierr .ne. 0) go to 9000
c
c     %------------------------------------------------%
c     | Scale the returning eigenvectors so that their |
c     | Euclidean norms are all one. LAPACK subroutine |
c     | ztrevc returns each eigenvector normalized so  |
c     | that the element of largest magnitude has      |
c     | magnitude 1; here the magnitude of a complex   |
c     | number (x,y) is taken to be |x| + |y|.         |
c     %------------------------------------------------%
c
      do 10 j=1, n
            temp = dznrm2( n, q(1,j), 1 )
            call zdscal ( n, rone / temp, q(1,j), 1 )
   10 continue
c
c     %----------------------------%
c     | Compute the Ritz estimates |
c     %----------------------------%
c
      call zcopy(n, q(n,1), n, bounds, 1)
      call zdscal(n, rnorm, bounds, 1)
c
 9000 continue
      return
c
c     %---------------%
c     | End of zneigh |
c     %---------------%
c
      end
c\BeginDoc
c 
c\Name: zneupd 
c 
c\Description: 
c  This subroutine returns the converged approximations to eigenvalues 
c  of A*z = lambda*B*z and (optionally): 
c 
c      (1) The corresponding approximate eigenvectors; 
c 
c      (2) An orthonormal basis for the associated approximate 
c          invariant subspace; 
c 
c      (3) Both.  
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal 
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied). 
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to ZNAUPD.  ZNAUPD must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such 
c  in the comments that follow.   The computed orthonormal basis for the 
c  invariant subspace corresponding to these Ritz values is referred to as a 
c  Schur basis. 
c 
c  The definition of OP as well as other terms and the relation of computed
c  Ritz values and vectors of OP with respect to the given problem
c  A*z = lambda*B*z may be found in the header of ZNAUPD.  For a brief 
c  description, see definitions of IPARAM(7), MODE and WHICH in the
c  documentation of ZNAUPD.
c
c\Usage:
c  call zneupd 
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT, 
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, 
c       WORKL, LWORKL, RWORK, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT)
c          Specifies whether a basis for the invariant subspace corresponding
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
c                                See Remarks below.
c
c  HOWMNY  Character*1  (INPUT)
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors;
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the  Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT need not be initialized 
c          but it is used as internal workspace.
c
c  D       Complex*16 array of dimension NEV+1.  (OUTPUT)
c          On exit, D contains the  Ritz  approximations 
c          to the eigenvalues lambda for A*z = lambda*B*z.
c
c  Z       Complex*16 N by NEV array    (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
c          Z represents approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z.
c
c          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required, 
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi 
c          basis array V computed by ZNAUPD.  In this case the Arnoldi basis 
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ) is required.  
c          In any case,  LDZ .ge. 1 is required.
c
c  SIGMA   Complex*16  (INPUT)
c          If IPARAM(7) = 3 then SIGMA represents the shift. 
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  WORKEV  Complex*16 work array of dimension 2*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to ZNAUPD that was just completed.               ****
c
c  NOTE: The remaining arguments 
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, 
c           WORKD, WORKL, LWORKL, RWORK, INFO 
c
c         must be passed directly to ZNEUPD following the last call 
c         to ZNAUPD.  These arguments MUST NOT BE MODIFIED between
c         the the last call to ZNAUPD and the call to ZNEUPD.
c
c  Three of these parameters (V, WORKL and INFO) are also output parameters:
c
c  V       Complex*16 N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by ZNAUPD .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+2*ncv) contains information obtained in
c          znaupd.  They are not changed by zneupd.
c          WORKL(ncv*ncv+2*ncv+1:3*ncv*ncv+4*ncv) holds the
c          untransformed Ritz values, the untransformed error estimates of 
c          the Ritz values, the upper triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by zneupd.
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the NCV RITZ values of the
c                     original system.
c          IPNTR(10): Not used
c          IPNTR(11): pointer to the NCV corresponding error estimates.
c          IPNTR(12): pointer to the NCV by NCV upper triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     zneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine csheqr
c                could not be reordered by LAPACK routine ztrsen.
c                Re-enter subroutine zneupd with IPARAM(5)=NCV and
c                increase the size of the array D to have
c                dimension at least dimension NCV and allocate at least NCV
c                columns for Z. NOTE: Not necessary if Z and V share
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 1 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation.
c                This should never happened.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine ztrevc.
c          = -10: IPARAM(7) must be 1,2,3
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: ZNAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: ZNEUPD got a different count of the number of converged
c                 Ritz values than ZNAUPD got.  This indicates the user
c                 probably made an error in passing data from ZNAUPD to
c                 ZNEUPD or that the data was modified before entering
c                 ZNEUPD
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
c     "How to Implement the Spectral Transformation", Math Comp.,
c     Vol. 48, No. 178, April, 1987 pp. 664-673. 
c
c\Routines called:
c     zgeqr2  LAPACK routine that computes the QR factorization of 
c             a matrix.
c     zlacpy  LAPACK matrix copy routine.
c     zlahqr  LAPACK routine that computes the Schur form of a
c             upper Hessenberg matrix.
c     zlaset  LAPACK matrix initialization routine.
c     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper triangular form.
c     ztrsen  LAPACK routine that re-orders the Schur form.
c     zunm2r  LAPACK routine that applies an orthogonal matrix in 
c             factored form.
c     dlamch  LAPACK routine that determines machine constants.
c     ztrmm   Level 3 BLAS matrix times an upper triangular matrix.
c     zgeru   Level 2 BLAS rank one update to a matrix.
c     zcopy   Level 1 BLAS that copies one vector to another .
c     zscal   Level 1 BLAS that scales a vector.
c     zdscal  Level 1 BLAS that scales a complex vector by a real number.
c     dznrm2  Level 1 BLAS that computes the norm of a complex vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented. 
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .true. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c       transpose( V(:,1:IPARAM(5)) ) * V(:,1:IPARAM(5)) = I
c     are approximately satisfied.
c     Here T is the leading submatrix of order IPARAM(5) of the 
c     upper triangular matrix stored workl(ipntr(12)). 
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas 
c     Dept. of Computational & 
c     Applied Mathematics 
c     Rice University 
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: neupd.F   SID: 2.8   DATE OF SID: 07/21/02   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine zneupd(rvec , howmny, select, d     ,
     &                   z    , ldz   , sigma , workev,
     &                   bmat , n     , which , nev   ,
     &                   tol  , resid , ncv   , v     ,
     &                   ldv  , iparam, ipntr , workd ,
     &                   workl, lworkl, rwork , info  )
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Complex*16     
     &           sigma
      Double precision 
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           rwork(ncv)
      Complex*16
     &           d(nev)     , resid(n)     , v(ldv,ncv),
     &           z(ldz, nev), 
     &           workd(3*n) , workl(lworkl), workev(2*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      parameter  (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0))
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds, ierr  , ih    , ihbds, iheig , nconv ,
     &           invsub, iuptri, iwev  , j    , ldh   , ldq   ,
     &           mode  , ritz  , wr   , k     , irz   ,
     &           ibd   , outncv, iq    , np   , numcnv, jj    ,
     &           ishift
      Complex*16
     &           rnorm, temp, vl(1)
      Double precision
     &           conds, sep, rtemp, eps23
      logical    reord
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zcopy , zgeru, zgeqr2, zlacpy, zmout,
     &           zunm2r, ztrmm, 
     &           zlahqr
c  
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dznrm2, dlamch, dlapy2
      external   dznrm2, dlamch, dlapy2
c
      Complex*16
     &           zdotc
      external   zdotc
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
c     %-------------------------------%
c     | Quick return                  |
c     | Check for incompatible input  |
c     %-------------------------------%
c
      ierr = 0
c
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 4*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else 
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')    ierr = -11
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c 
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, WORKEV, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+ncv) := ritz values            |
c     | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv) := error bounds     |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by ZNEUPD.                 |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := The untransformed |
c     |                                      Ritz values.         |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                                      error bounds of      |
c     |                                      the Ritz values      |
c     | workl(ncv*ncv+4*ncv+1:2*ncv*ncv+4*ncv) := Holds the upper |
c     |                                      triangular matrix    |
c     |                                      for H.               |
c     | workl(2*ncv*ncv+4*ncv+1: 3*ncv*ncv+4*ncv) := Holds the    |
c     |                                      associated matrix    |
c     |                                      representation of    |
c     |                                      the invariant        |
c     |                                      subspace for H.      |
c     | GRAND total of NCV * ( 3 * NCV + 4 ) locations.           |
c     %-----------------------------------------------------------%
c     
      ih     = ipntr(5)
      ritz   = ipntr(6)
      iq     = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheig  = bounds + ldh
      ihbds  = iheig  + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheig
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wr = 1
      iwev = wr + ncv
c
c     %-----------------------------------------%
c     | irz points to the Ritz values computed  |
c     |     by _neigh before exiting _naup2.    |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irz = ipntr(14) + ncv*ncv
      ibd = irz + ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c
      if (rvec) then
c
         reord = .false.
c
c        %---------------------------------------------------%
c        | Use the temporary bounds array to store indices   |
c        | These will be used to mark the select array later |
c        %---------------------------------------------------%
c
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
c
c        %-------------------------------------%
c        | Select the wanted Ritz values.      |
c        | Sort the Ritz values so that the    |
c        | wanted ones appear at the tailing   |
c        | NEV positions of workl(irr) and     |
c        | workl(iri).  Move the corresponding |
c        | error estimates in workl(ibd)       |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call zngets(ishift, which     , nev          ,
     &                np    , workl(irz), workl(bounds))
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            rtemp = max(eps23,
     &                 dlapy2 ( dble(workl(irz+ncv-j)),
     &                          dimag(workl(irz+ncv-j)) ))
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          dlapy2( dble(workl(ibd+jj-1)),
     &          dimag(workl(ibd+jj-1)) )
     &          .le. tol*rtemp) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by dnaupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the dnaupd data.           |
c        %-----------------------------------------------------------%
c
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
c
c        %-------------------------------------------------------%
c        | Call LAPACK routine zlahqr to compute the Schur form |
c        | of the upper Hessenberg matrix returned by ZNAUPD.   |
c        | Make a copy of the upper Hessenberg matrix.           |
c        | Initialize the Schur vector matrix Q to the identity. |
c        %-------------------------------------------------------%
c
         call zcopy(ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call zlaset('All', ncv, ncv          , 
     &                zero , one, workl(invsub),
     &                ldq)
         call zlahqr(.true., .true.       , ncv          , 
     &                1     , ncv          , workl(iuptri),
     &                ldh   , workl(iheig) , 1            ,
     &                ncv   , workl(invsub), ldq          ,
     &                ierr)
         call zcopy(ncv         , workl(invsub+ncv-1), ldq,
     &               workl(ihbds), 1)
c
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c
         if (reord) then
c
c           %-----------------------------------------------%
c           | Reorder the computed upper triangular matrix. |
c           %-----------------------------------------------%
c
            call ztrsen('None'       , 'V'          , select      ,
     &                   ncv          , workl(iuptri), ldh         ,
     &                   workl(invsub), ldq          , workl(iheig),
     &                   nconv        , conds        , sep         , 
     &                   workev       , ncv          , ierr)
c
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
         end if
c
c        %---------------------------------------------%
c        | Copy the last row of the Schur basis matrix |
c        | to workl(ihbds).  This vector will be used  |
c        | to compute the Ritz estimates of converged  |
c        | Ritz values.                                |
c        %---------------------------------------------%
c
         call zcopy(ncv         , workl(invsub+ncv-1), ldq,
     &               workl(ihbds), 1)
c 
c        %--------------------------------------------%
c        | Place the computed eigenvalues of H into D |
c        | if a spectral transformation was not used. |
c        %--------------------------------------------%
c
         if (type .eq. 'REGULR') then
            call zcopy(nconv, workl(iheig), 1, d, 1)
         end if
c
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(invsub,ldq).                            |
c        %----------------------------------------------------------%
c
         call zgeqr2(ncv , nconv , workl(invsub),
     &                ldq , workev, workev(ncv+1),
     &                ierr)
c
c        %--------------------------------------------------------%
c        | * Postmultiply V by Q using zunm2r.                    |
c        | * Copy the first NCONV columns of VQ into Z.           |
c        | * Postmultiply Z by R.                                 |
c        | The N by NCONV matrix Z is now a matrix representation |
c        | of the approximate invariant subspace associated with  |
c        | the Ritz values in workl(iheig). The first NCONV       | 
c        | columns of V are now approximate Schur vectors         |
c        | associated with the upper triangular matrix of order   |
c        | NCONV in workl(iuptri).                                |
c        %--------------------------------------------------------%
c
         call zunm2r('Right', 'Notranspose', n            ,
     &                ncv    , nconv        , workl(invsub),
     &                ldq    , workev       , v            ,
     &                ldv    , workd(n+1)   , ierr)
         call zlacpy('All', n, nconv, v, ldv, z, ldz)
c
         do 20 j=1, nconv
c
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | triangular form of workl(iuptri,ldq).             |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones.          |
c           %---------------------------------------------------%
c
            if ( dble( workl(invsub+(j-1)*ldq+j-1) ) .lt. 
     &                  dble(zero) ) then
               call zscal(nconv, -one, workl(iuptri+j-1), ldq)
               call zscal(nconv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
c
 20      continue
c
         if (howmny .eq. 'A') then
c
c           %--------------------------------------------%
c           | Compute the NCONV wanted eigenvectors of T |
c           | located in workl(iuptri,ldq).              |
c           %--------------------------------------------%
c
            do 30 j=1, ncv
               if (j .le. nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
c
            call ztrevc('Right', 'Select'     , select       ,
     &                   ncv    , workl(iuptri), ldq          ,
     &                   vl     , 1            , workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev , rwork        , ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | ztrevc returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1.                                   |
c           %------------------------------------------------%
c
            do 40 j=1, nconv
                  rtemp = dznrm2(ncv, workl(invsub+(j-1)*ldq), 1)
                  rtemp = dble(one) / rtemp
                  call zdscal ( ncv, rtemp,
     &                 workl(invsub+(j-1)*ldq), 1 )
c
c                 %------------------------------------------%
c                 | Ritz estimates can be obtained by taking |
c                 | the inner product of the last row of the |
c                 | Schur basis of H with eigenvectors of T. |
c                 | Note that the eigenvector matrix of T is |
c                 | upper triangular, thus the length of the |
c                 | inner product can be set to j.           |
c                 %------------------------------------------%
c 
                  workev(j) = zdotc(j, workl(ihbds), 1,
     &                        workl(invsub+(j-1)*ldq), 1)
 40         continue
c
c           %---------------------------------------%
c           | Copy Ritz estimates into workl(ihbds) |
c           %---------------------------------------%
c 
            call zcopy(nconv, workev, 1, workl(ihbds), 1)
c
c           %----------------------------------------------%
c           | The eigenvector matrix Q of T is triangular. |
c           | Form Z*Q.                                    |
c           %----------------------------------------------%
c
            call ztrmm('Right'   , 'Upper'      , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
         end if 
c
      else
c
c        %--------------------------------------------------%
c        | An approximate invariant subspace is not needed. |
c        | Place the Ritz values computed ZNAUPD into D.    |
c        %--------------------------------------------------%
c
         call zcopy(nconv, workl(ritz), 1, d, 1)
         call zcopy(nconv, workl(ritz), 1, workl(iheig), 1)
         call zcopy(nconv, workl(bounds), 1, workl(ihbds), 1)
c
      end if
c
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
         if (rvec) 
     &      call zscal(ncv, rnorm, workl(ihbds), 1)
c      
      else
c     
c        %---------------------------------------%
c        |   A spectral transformation was used. |
c        | * Determine the Ritz estimates of the |
c        |   Ritz values in the original system. |
c        %---------------------------------------%
c
         if (rvec) 
     &      call zscal(ncv, rnorm, workl(ihbds), 1)
c    
         do 50 k=1, ncv
            temp = workl(iheig+k-1)
            workl(ihbds+k-1) = workl(ihbds+k-1) / temp / temp
  50     continue
c  
      end if
c
c     %-----------------------------------------------------------%
c     | *  Transform the Ritz values back to the original system. |
c     |    For TYPE = 'SHIFTI' the transformation is              |
c     |             lambda = 1/theta + sigma                      |
c     | NOTES:                                                    |
c     | *The Ritz vectors are not affected by the transformation. |
c     %-----------------------------------------------------------%
c    
      if (type .eq. 'SHIFTI') then
         do 60 k=1, nconv
            d(k) = one / workl(iheig+k-1) + sigma
  60     continue
      end if
c
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 3. See reference 3.                  |
c     %-------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta.                           |
c        %------------------------------------------------%
c
         do 100 j=1, nconv
            if (workl(iheig+j-1) .ne. zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) /
     &                      workl(iheig+j-1)
            endif
 100     continue
!
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call zgeru (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
      return
c     
c     %---------------%
c     | End of zneupd|
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: zngets
c
c\Description: 
c  Given the eigenvalues of the upper Hessenberg matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of 
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: call this even in the case of user specified shifts in order
c  to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call zngets
c      ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> want the KEV eigenvalues of largest magnitude.
c          'SM' -> want the KEV eigenvalues of smallest magnitude.
c          'LR' -> want the KEV eigenvalues of largest REAL part.
c          'SR' -> want the KEV eigenvalues of smallest REAL part.
c          'LI' -> want the KEV eigenvalues of largest imaginary part.
c          'SI' -> want the KEV eigenvalues of smallest imaginary part.
c
c  KEV     Integer.  (INPUT)
c          The number of desired eigenvalues.
c
c  NP      Integer.  (INPUT)
c          The number of shifts to compute.
c
c  RITZ    Complex*16 array of length KEV+NP.  (INPUT/OUTPUT)
c          On INPUT, RITZ contains the the eigenvalues of H.
c          On OUTPUT, RITZ are sorted so that the unwanted
c          eigenvalues are in the first NP locations and the wanted
c          portion is in the last KEV locations.  When exact shifts are 
c          selected, the unwanted part corresponds to the shifts to 
c          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
c          are further sorted so that the ones with largest Ritz values
c          are first.
c
c  BOUNDS  Complex*16 array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\Routines called:
c     zsortc  ARPACK sorting routine.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c\SCCS Information: @(#)
c FILE: ngets.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. This routine does not keep complex conjugate pairs of
c        eigenvalues together.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zngets ( ishift, which, kev, np, ritz, bounds)
c
c     %----------------------------------------------------%
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      integer    ishift, kev, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16
     &           bounds(kev+np), ritz(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      parameter (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0))
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     %-------------------------------%
c 
      call zsortc (which, .true., kev+np, ritz, bounds)
c     
      if ( ishift .eq. 1 ) then
c     
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when the shifts  |
c        | are applied in subroutine znapps.                     |
c        | Be careful and use 'SM' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c     
         call zsortc ( 'SM', .true., np, bounds, ritz )
c
      end if
c     
      return
c     
c     %---------------%
c     | End of zngets |
c     %---------------%
c     
      end
c\BeginDoc
c
c\Name: zsortc
c
c\Description:
c  Sorts the Complex*16 array in X into the order 
c  specified by WHICH and optionally applies the permutation to the
c  Double precision  array Y. 
c
c\Usage:
c  call zsortc
c     ( WHICH, APPLY, N, X, Y )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> sort X into increasing order of magnitude.
c          'SM' -> sort X into decreasing order of magnitude.
c          'LR' -> sort X with real(X) in increasing algebraic order 
c          'SR' -> sort X with real(X) in decreasing algebraic order
c          'LI' -> sort X with imag(X) in increasing algebraic order
c          'SI' -> sort X with imag(X) in decreasing algebraic order
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to array Y.
c          APPLY = .FALSE. -> do not apply the sorted order to array Y.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  X       Complex*16 array of length N.  (INPUT/OUTPUT)
c          This is the array to be sorted.
c
c  Y       Complex*16 array of length N.  (INPUT/OUTPUT)
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines called:
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c     Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#)
c FILE: sortc.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zsortc (which, apply, n, x, y)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16     
     &           x(0:n-1), y(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Complex*16     
     &           temp
      Double precision 
     &           temp1, temp2
c
c     %--------------------%
c     | External functions |
c     %--------------------%
c
      Double precision
     &           dlapy2
c
c     %--------------------%
c     | Intrinsic Functions |
c     %--------------------%
       Intrinsic
     &           dble, dimag
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c 
      if (which .eq. 'LM') then
c
c        %--------------------------------------------%
c        | Sort X into increasing order of magnitude. |
c        %--------------------------------------------%
c
   10    continue
         if (igap .eq. 0) go to 9000
c
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            temp1 = dlapy2(dble(x(j)),dimag(x(j)))
            temp2 = dlapy2(dble(x(j+igap)),dimag(x(j+igap)))
c
            if (temp1.gt.temp2) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
c
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                go to 30
            end if
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        %--------------------------------------------%
c        | Sort X into decreasing order of magnitude. |
c        %--------------------------------------------%
c
   40    continue
         if (igap .eq. 0) go to 9000
c
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j .lt. 0) go to 60
c
            temp1 = dlapy2(dble(x(j)),dimag(x(j)))
            temp2 = dlapy2(dble(x(j+igap)),dimag(x(j+igap)))
c
            if (temp1.lt.temp2) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c 
      else if (which .eq. 'LR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into increasing order of algebraic. |
c        %------------------------------------------------%
c
   70    continue
         if (igap .eq. 0) go to 9000
c
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (dble(x(j)).gt.dble(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c 
      else if (which .eq. 'SR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into decreasing order of algebraic. |
c        %------------------------------------------------%
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (dble(x(j)).lt.dble(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
c 
      else if (which .eq. 'LI') then
c
c        %--------------------------------------------%
c        | Sort XIMAG into increasing algebraic order |
c        %--------------------------------------------%
c
  130    continue
         if (igap .eq. 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
c
            if (j.lt.0) go to 150
c
            if (dimag(x(j)).gt.dimag(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 150
            endif
            j = j-igap
            go to 140
  150    continue
         igap = igap / 2
         go to 130
c 
      else if (which .eq. 'SI') then
c
c        %---------------------------------------------%
c        | Sort XIMAG into decreasing algebraic order  |
c        %---------------------------------------------%
c
  160    continue
         if (igap .eq. 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
c
            if (j.lt.0) go to 180
c
            if (dimag(x(j)).lt.dimag(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 180
            endif
            j = j-igap
            go to 170
  180    continue
         igap = igap / 2
         go to 160
      end if
c 
 9000 continue
      return
c
c     %---------------%
c     | End of zsortc |
c     %---------------%
c
      end
