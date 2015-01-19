SUBROUTINE ROTDEF (iz,tau,iord,nat,pos,ndif,rotij,tauij,mult,lattic)
  !                                                                       
  !     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN      
  !     ATOM TO ITS CORRESPONDING POSITION OF AN EQUIVALENT ATOM.
  !               
  IMPLICIT NONE
  INTEGER,INTENT(IN)  :: iord,nat,ndif
  INTEGER,INTENT(IN)  :: iz(3,3,iord),mult(nat)
  CHARACTER*4,INTENT(IN) :: lattic
  REAL(8),INTENT(IN)  :: pos(3,ndif),tau(3,iord)
  REAL(8),INTENT(OUT) :: rotij(3,3,ndif),tauij(3,ndif)
  CHARACTER*67    ERRMSG
  REAL(8)   x(3),x1(3),toler,toler2,one
  INTEGER   i,i1,j,k,l,m,index,index1,jatom,ncount
  DATA TOLER/1.D-7/,ONE/1.D0/
  toler2=1.5d0*toler            
  INDEX=0                                                           
  NCOUNT=0                                                          
  DO JATOM=1,NAT                                                 
     INDEX1=INDEX+1                                                 
     DO M=1,MULT(JATOM)                                          
        INDEX=INDEX+1                                               
        DO I=1,IORD                                              
           x(1:3)=0.d0
           x=MATMUL(TRANSPOSE(iz(1:3,1:3,i)),pos(1:3,index1))
           x(1:3)=x(1:3)+tau(1:3,i)
           x1(1:3)=MOD(ABS(X(1:3)-POS(1:3,INDEX))+toler,one)-toler
!           WRITE(*,*) 'JATOM,INDEX,I',JATOM,INDEX,I                   
!           WRITE(*,*) ABS(X1(1:3)),toler
           IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
              NCOUNT=NCOUNT+1                                             
              TAUIJ(1:3,INDEX)=TAU(1:3,I)
              ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
              GOTO 30                                                     
           END IF
           !....check positions for centered lattices
           if(lattic(1:1).eq.'B') then
              x1(1:3)=mod(x1(1:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,INDEX)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                 GOTO 30                                                     
              END IF
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1(1:2)=mod(x1(1:2)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,INDEX)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                 GOTO 30                                                     
              END IF
              x1(1:2)=mod(x1(1:2)+0.5d0,one)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1(1)=mod(x1(1)+0.5d0+toler,one)-toler
              x1(3)=mod(x1(3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,INDEX)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                 GOTO 30                                                     
              END IF
              x1(1)=mod(x1(1)+0.5d0,one)
              x1(3)=mod(x1(3)+0.5d0,one)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              x1(2:3)=mod(x1(2:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,INDEX)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                 GOTO 30                                                     
              END IF
           end if
        ENDDO
        !
        !           Error: no symmetry operation found
        !
        GOTO 900
30      CONTINUE                                                       
     ENDDO
  ENDDO
  IF(NCOUNT.NE.INDEX) GOTO 910
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
910 CALL OUTERR('ROTDEF','ROTIJ not defined for all atoms of basis.')
  WRITE (ERRMSG,9030) NCOUNT
  CALL OUTERR('ROTDEF',ERRMSG)
  STOP 'ROTDEF - Error'
  !
9000 FORMAT ('for jatom, index',I2,i2)
9010 FORMAT ('atomposition of jatom',3F12.7)
9020 FORMAT ('atomposition of index',3F12.7)
9030 FORMAT ('NCOUNT=',I2)
END SUBROUTINE ROTDEF
