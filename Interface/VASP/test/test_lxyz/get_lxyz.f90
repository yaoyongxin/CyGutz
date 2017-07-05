      SUBROUTINE SETUP_L(L)
      IMPLICIT NONE
      
      INTEGER L,M,M_,I,J,K
      integer,parameter::q=8, iu=9
      REAL(q) C_UP,C_DW
      
      COMPLEX(q) U_C2R(2*L+1,2*L+1),U_R2C(2*L+1,2*L+1),TMP(2*L+1,2*L+1)
      COMPLEX(q) L_OP_C(2*L+1,2*L+1,3),L_OP_R(2*L+1,2*L+1,3)

      open(iu, file='v_lxyz.txt', status='replace')

! set up L operator (in units of h_bar) for complex spherical harmonics y_lm     
!
!   |y_lm1> L_k <y_lm2| = |y_lm1> L_OP_C(m1,m2,k) <y_lm2| , where k=x,y,z
!
      L_OP_C=(0._q,0._q)
      
      DO M=1,2*L+1
         M_=M-L-1   
         C_UP=SQRT(REAL((L-M_)*(L+M_+1),q))/2.0_q
         C_DW=SQRT(REAL((L+M_)*(L-M_+1),q))/2.0_q
         ! fill x-component
         IF ((M_+1)<= L) L_OP_C(M+1,M,1)=C_UP
         IF ((M_-1)>=-L) L_OP_C(M-1,M,1)=C_DW
         ! fill y-component
         IF ((M_+1)<= L) L_OP_C(M+1,M,2)=-CMPLX(0,C_UP,q)
         IF ((M_-1)>=-L) L_OP_C(M-1,M,2)= CMPLX(0,C_DW,q)
         ! fill z-component
         L_OP_C(M,M,3)=M_
      ENDDO

      DO I=1,3
        WRITE(iu,*) 'component in comp_sph_harm',I 
        call write_matrix(L_OP_C(:,:,i),2*L+1,iu)
      ENDDO
      
! set up transformation matrix real->complex spherical harmonics
!
!  <y_lm1|Y_lm2> = U_R2C(m1,m2) 
! 
! where y_lm and Y_lm are, respectively, the complex and real 
! spherical harmonics
!
      U_R2C=(0._q,0._q)
          
      DO M=1,2*L+1
         M_=M-L-1
         IF (M_>0) THEN
            U_R2C( M_+L+1,M)=(-1)**M_/SQRT(2._q)
            U_R2C(-M_+L+1,M)=1/SQRT(2._q)
         ENDIF
         IF (M_==0) THEN
            U_R2C(L+1,L+1)=1
         ENDIF
         IF (M_<0) THEN
            U_R2C( M_+L+1,M)= CMPLX(0,1/SQRT(2._q),q)
            U_R2C(-M_+L+1,M)=-CMPLX(0,(-1)**M_/SQRT(2._q),q)
         ENDIF
      ENDDO

! set up transformation matrix complex->real spherical harmonics
!
!  <Y_lm1|y_lm2> = U_C2R(m1,m2)
! 
! where y_lm and Y_lm are, respectively, the complex and real 
! spherical harmonics
!
      U_C2R=(0._q,0._q)
      
      DO M=1,2*L+1
         M_=M-L-1
         IF (M_>0) THEN
            U_C2R( M_+L+1,M)=(-1)**M_/SQRT(2._q)
            U_C2R(-M_+L+1,M)=CMPLX(0,(-1)**M_/SQRT(2._q),q)
         ENDIF
         IF (M_==0) THEN
            U_C2R(L+1,L+1)=1
         ENDIF
         IF (M_<0) THEN
            U_C2R( M_+L+1,M)=-CMPLX(0,1/SQRT(2._q),q)
            U_C2R(-M_+L+1,M)=1/SQRT(2._q)
         ENDIF
      ENDDO


      WRITE(iu,*) '<comp_sph_harm | real_sph_harm>'
      call write_matrix(U_R2C,2*L+1,iu)


      TMP=(0._q,0._q)
      DO M=1,2*L+1
      DO M_=1,2*L+1
         DO I=1,2*L+1
            TMP(M,M_)=TMP(M,M_)+U_C2R(M,I)*U_R2C(I,M_)
         ENDDO
      ENDDO
      ENDDO

      WRITE(iu,*) ' identity check'
      call write_matrix(tmp,2*L+1,iu)

      write(*,*)maxval(abs(U_C2R-transpose(conjg(U_R2C)))),'unitary check!'

! Calculate L operator (in units of h_bar) with respect to 
! the real spherical harmonics Y_lm
!
!    |Y_lm1> L_k <Y_lm2| = |Y_lm1> L_OP_R(m1,m2,k) <Y_lm2| , where k=x,y,z
!
! n.b. L_OP_R(m1,m2,k)= \sum_ij U_C2R(m1,i) L_OP_C(i,j,k) U_R2C(j,m2)
!
      L_OP_R=(0._q,0._q)

      DO M=1,2*L+1
      DO M_=1,2*L+1
         DO I=1,2*L+1
         DO J=1,2*L+1
            L_OP_R(M,M_,:)=L_OP_R(M,M_,:)+U_C2R(M,I)*L_OP_C(I,J,:)*U_R2C(J,M_)
         ENDDO
         ENDDO      
      ENDDO
      ENDDO
      
      DO I=1,3
        WRITE(iu,*) 'component in real_sph_harm',I 
        call write_matrix(L_OP_R(:,:,i),2*L+1,iu)
      ENDDO

      close(iu)

      end SUBROUTINE SETUP_L



    subroutine write_matrix(a,n,iu)
    integer n,iu
    complex(8) a(n,n)

    integer m1,m2

    do m1=1,n
        do m2=1,n
            write(iu,'(f12.6)',advance='no')real(a(m1,m2))
        enddo
        write(iu,*)
    enddo
    write(iu,*)
    do m1=1,n
        do m2=1,n
            write(iu,'(f12.6)',advance='no')aimag(a(m1,m2))
        enddo
        write(iu,*)
    enddo
    write(iu,*)
    return

    end subroutine write_matrix      
