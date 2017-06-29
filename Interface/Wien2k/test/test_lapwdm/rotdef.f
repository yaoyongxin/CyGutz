SUBROUTINE ROTDEF
  !                                                                       
  !     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN      
  !     ATOM TO ITS CORRESPONDING POSITION OF AN EQUIVALENT ATOM.
  !               
        USE param
        USE struct
        USE rotat
      IMPLICIT REAL*8 (A-H,O-Z)


  CHARACTER*67    ERRMSG
  REAL(8)   x(3),x1(3),toler,toler2,one

  DATA TOLER/1.D-7/,ONE/1.D0/
!
! Additions for averaging
        integer jtimes,nthis
        real*8 X00, Y00, Z00, X0, Y0, Z0,XA,YA,ZA
        real*8 XB, YB, ZB, TT
!
  toler2=1.5d0*toler            
!
!---------------------------------------------------------------------  
!  Patch up symmetry issues with positions in limited precision
      INDEX=0
      DO JATOM=1,NAT
        INDEX=INDEX+MULT(JATOM)
      ENDDO
      DO JTIMES=1,2
!      write(88,*)'Pass ',Jtimes
      DO JATOM=1,index                                                             
         NTHIS=0
         X00=POS(1,JATOM)
         Y00=POS(2,JATOM)
         Z00=POS(3,JATOM)
         X0=0
         Y0=0
         Z0=0                                           
            DO I=1,IORD                                              
            XA=TAU(1,I)+1.D0
            YA=TAU(2,I)+1.D0
            ZA=TAU(3,I)+1.D0                                                      
            DO J=1,3 
                XA=XA+dble(IZ(J,1,I))*POS(J,JATOM)                               
                YA=YA+dble(IZ(J,2,I))*POS(J,JATOM) 
                ZA=ZA+dble(IZ(J,3,I))*POS(J,JATOM)
            ENDDO
            XA = mod(XA,1.D0)
            YA = mod(YA,1.D0)
            ZA = mod(ZA,1.D0)
            if(xa.gt.0.999999)XA=XA-1.D0
            if(ya.gt.0.999999)YA=YA-1.D0
            if(Za.gt.0.999999)ZA=ZA-1.D0
            XB=ABS(XA-X00)*aa                                      
            YB=ABS(YA-Y00)*bb                 
            ZB=ABS(ZA-Z00)*cc                                      
                                      
            IF(XB.LT.1d-2.AND.YB.LT.1d-2.AND.ZB.LT.1d-2) THEN        
                NTHIS=NTHIS+1                                             
                X0=X0+XA
                Y0=Y0+YA
                Z0=Z0+ZA
            ENDIF
            ENDDO
       if(nthis.gt.1)then
       TT=1.D0/dble(nthis)                         
       POS(1,JATOM)=X0*TT
       POS(2,JATOM)=Y0*TT
       POS(3,JATOM)=Z0*TT
!      write(6,88)'OLD     ',X00,Y00,Z00
!      write(6,88)'Average ',X0*TT,Y0*TT,Z0*TT,NTHIS,JATOM
       endif
 88    format(a,3f16.13,2i3)
      ENDDO
      ENDDO
!---------------------------------------------------------------------  
!
  INDEX=0                                                           
  DO JATOM=1,NAT                                                 
     INDEX1=INDEX+1
     NCOUNT=0                                                 
     DO M=1,MULT(JATOM)                                          
        INDEX=INDEX+1
        DO I=1,IORD                                              
           x(1:3)=0.d0
           x=MATMUL(TRANSPOSE(iz(1:3,1:3,i)),pos(1:3,index1))
           x(1:3)=x(1:3)+tau(1:3,i)
           x1(1:3)=MOD(ABS(X(1:3)-POS(1:3,INDEX))+toler,one)-toler
!          WRITE(6,*) 'JATOM,INDEX,I',JATOM,INDEX,I                   
!          WRITE(6,*) ABS(X1(1:3)),toler
           IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
              NCOUNT=NCOUNT+1                                             
              itr(i,jatom)=M
           END IF
           !....check positions for centered lattices
           if(lattic(1:1).eq.'B') then
              x1(1:3)=mod(x1(1:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 itr(i,jatom)=M
              END IF
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1(1:2)=mod(x1(1:2)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 itr(i,jatom)=M
              END IF
              x1(1:2)=mod(x1(1:2)+0.5d0,one)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1(1)=mod(x1(1)+0.5d0+toler,one)-toler
              x1(3)=mod(x1(3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 itr(i,jatom)=M
              END IF
              x1(1)=mod(x1(1)+0.5d0,one)
              x1(3)=mod(x1(3)+0.5d0,one)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              x1(2:3)=mod(x1(2:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 itr(i,jatom)=M
              END IF
           end if
        ENDDO
        !
        !           Error: no symmetry operation found
        !
     ENDDO
!    Debug
!    write(6,*)'ITR array ',jatom,itr(1:iord,jatom)
     if(NCOUNT .eq. 0)GOTO 900

  ENDDO
  RETURN     

  !        Error messages
  !
900 CALL OUTERR('ROTDEF','no symmetry operation found.')
  WRITE(ERRMSG,9000) jatom, index          
  CALL OUTERR('ROTDEF',ERRMSG)
  WRITE(ERRMSG,9010) (POS(I1,index1),I1=1,3) 
  CALL OUTERR('ROTDEF',ERRMSG)
  WRITE(ERRMSG,9020) (POS(I1,INDEX),I1=1,3) 
  CALL OUTERR('ROTDEF',ERRMSG)
  STOP 'ROTDEF - Error'
  !
9000 FORMAT ('for jatom, index',I2,i2)
9010 FORMAT ('atomposition of jatom',3F12.7)
9020 FORMAT ('atomposition of index',3F12.7)
9030 FORMAT ('NCOUNT=',I2)
END SUBROUTINE ROTDEF
