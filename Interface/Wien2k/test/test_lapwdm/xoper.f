      subroutine Xoper(L,S2,X)
!     Xoper calculates matrices
!     <ML' MS'|L*S|ML  MS> and also others spin and orbital operators
!     including spin dipolar and orbital dipolar parts of the hyperfine field
!     this version is dimensioned for S=1/2 and L smaller than 4. It may
!     be used for any S, L if dimensions are increased
!     List of operators calculated
!     kls=1   Unit operator
!     kls=2   Sdzeta (dzeta||M)
!     kls=3   Ldzeta (dzeta||M)
!     kls=4   spin-orbit coupling (L.S)
!     kls=5   dipolar contribution to hyperfine field
!             Bdip=coef*[L(L+1)Sdzeta-(3/2){(L.S)Ldzeta+Ldzeta(L.S)}]
!     kls=6   Ldzeta**2
!     kls=7   WIEN2k_05 approx. to Bdip:
!             Bdip=coef*[L(L+1)-3Ldzeta**2]Sdzeta
!             P.Novak 2006 
!     kls=8   T=(L.S)Ldzeta+Ldzeta(L.S)    P.Novak 2007
      USE param
      IMPLICIT REAL*8(A-H,O-Z)
      integer S2
      complex*16 sksi(2,2),seta(2,2),sdzeta(2,2)   &
                ,sx(2,2),sy(2,2),sz(2,2)           &
                ,dsx(14,14),dsy(14,14),dsz(14,14)  &
                ,dsksi(14,14),dseta(14,14),dsdzeta(14,14)  &
                ,dLx(14,14),dLy(14,14),dLz(14,14)  &
                ,dLksi(14,14),dLeta(14,14),dLdzeta(14,14)  &
                ,X(0:LMAX2,-LMAX2:LMAX2,-LMAX2:LMAX2,2,2)  &
                ,Xop(14,14),soc(14,14),Xhelp(14,14)  
      common /aver/ krad,kls,cx(-20:20,20),iprx
      common /spinorb/ dsx,dsy,dsz,dsksi,dseta,dsdzeta  &
                ,dLx,dLy,dLz,dLksi,dLeta,dLdzeta  
!      abksi is parameter ksi for single electron equivalent hamiltonian
!      see Abragam, Bleaney, Electron Paramagnetic Resonance .. eq. (17.44)
!      note that hyperfine field of Abragam, Bleaney is negative of WIEN Bhf
      abksi=2.d0/dfloat((2*L-1)*(2*L+3))
      n=(2*L+1)*(S2+1) 
  ! spin-orbit coupling. X that needs s-o: s-o, Bdip, T
       do i=1,n
        do j=1,n
          soc(i,j)=(0.d0,0.d0)
          do k=1,n
           soc(i,j)=soc(i,j)+   &
!          dlx(i,k)*dsx(k,j)+dly(i,k)*dsy(k,j)+dlz(i,k)*dsz(k,j)
           dlksi(i,k)*dsksi(k,j)+dleta(i,k)*dseta(k,j)+dldzeta(i,k)*dsdzeta(k,j)
          enddo
        enddo
       enddo
! First calculate Xop(n,n)
write(76,*)kls,' kls'
      do i=1,n
       do j=1,n
        Xop(i,j)=(0.,0.)
       enddo
      enddo
      if(kls.eq.1)then   ! Unit operator
       do i=1,n
        Xop(i,i)=(1.,0.)
       enddo
      else if(kls.eq.2)then  !Sdzeta
       do i=1,n
        Xop(i,i)=dsdzeta(i,i)
       enddo
      else if(kls.eq.3)then  !Ldzeta
       do i=1,n
        do j=1,n
         Xop(i,j)=dLdzeta(i,j)
        enddo
       enddo 
      else if(kls.eq.4)then
        do i=1,n
         do j=1,n
          Xop(i,j)=soc(i,j)
         enddo
        enddo 
      else if(kls.eq.5)then    ! Bdip
        do i=1,n
         Xop(i,i)=dfloat(L*(L+1))*dsdzeta(i,i)
         do j=1,n
          do k=1,n
           Xop(i,j)=Xop(i,j)-  &
                    (3./2.)*(dldzeta(i,k)*soc(k,j)+soc(i,k)*dldzeta(k,j))
          enddo
         enddo
        enddo
       do i=1,n
        do j=1,n
          Xop(i,j)=abksi*Xop(i,j)
        enddo
       enddo
      else if(kls.eq.6)then    !Ldzeta**2
       do i=1,n
        do j=1,n
         do k=1,n
          Xop(i,j)=Xop(i,j)+dldzeta(i,k)*dldzeta(k,j)
         enddo
        enddo
       enddo
      else if(kls.eq.7)then    !WIEN2k_05 approx. of Bdip
       do i=1,n
        Xhelp(i,i)=dfloat(L*(L+1))    ! Xhelp=-3*Ldzeta**2+L*(L+1)
        do j=1,n
         do k=1,n
          Xhelp(i,j)=Xhelp(i,j)-3.*dldzeta(i,k)*dldzeta(k,j)
         enddo
        enddo
       enddo
       do i=1,n                ! Xop=abksi*[3*Ldzeta**2-L*(L+1)]*Sdzeta
        do j=1,n
         do k=1,n
          Xop(i,j)=Xop(i,j)+Xhelp(i,k)*dsdzeta(k,j)
         enddo
         Xop(i,j)=abksi*Xop(i,j)
        enddo
       enddo
      else if(kls.eq.8)then    ! (L.S)Ldzeta+Ldzeta(L.S)
       do i=1,n
        do j=1,n
         do k=1,n
          Xop(i,j)=Xop(i,j)+  &
                    (dldzeta(i,k)*soc(k,j)+soc(i,k)*dldzeta(k,j))
         enddo
        enddo
       enddo
      endif
      nL=2*L+1
   !rewrite Xop(i,j) to X(L,ML1,ML2,MS1,MS2)
   ! invert the spin indices of X, MS=1=up, MS=-1=dn
   ! because of opposite dn, up sequence of Xop and Xqtl
       M1=0
       DO 16 ML1=-L,L
        M1=M1+1
        M2=0          
        DO 16 ML2=-L,L
         M2=M2+1
         X(L,ML1,ML2,1,1)=Xop(M1+NL,M2+NL)   !Xupup
         X(L,ML1,ML2,2,1)=Xop(M1,M2+NL)      !Xdnup
         X(L,ML1,ML2,1,2)=Xop(M1+NL,M2)      !Xupdn
         X(L,ML1,ML2,2,2)=Xop(M1,M2)         !Xdndn
16    continue
      if(iprx.gt.0)then
       write(6,112)L
112    format(' X operator, S=1/2, L=',i2,' sequence')
       write(6,*)'             -1/2                     1/2'
       write(6,*)'     -L, -L+1 ......L           -L, -L+1 ......L'
       call printx(n,Xop)
      endif
      return
      end
