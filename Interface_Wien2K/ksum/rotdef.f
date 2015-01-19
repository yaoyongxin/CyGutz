SUBROUTINE ROTDEF
  !                                                                       
  !     ROTDEF GENERATES THE ROTATION-MATRICES ROTIJ(3,3,ntot) TAUIJ(3,ntot) FOR      
  !     NONSYMMORPHIC STRUCTURES. THESE MATRICES PERFORM ROTATIONS        
  !     FROM THE GENERAL COORDINATION-SYSTEM INTO LOCAL SYSTEMS WITH      
  !     SPECIFIED POINTSYMMETRY.                                          
  !     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN      
  !     ATOM TO IT'S CORRESPONDING POSITION OF AN EQUIVALENT              
  !     ATOM.                                                             
  !                                                                       
  USE param
  USE defs
  USE struct
  USE sym2
  IMPLICIT NONE !IMPLICIT REAL*8 (A-H,O-Z)
  !
  CHARACTER*67  :: ERRMSG
  REAL*8        :: TOLER, ONE, X, Y, Z, X1, Y1, Z1
  INTEGER       :: ntot, INDEX, NCOUNT, JATOM, M, INDEX1, I, J, L, I1
  LOGICAL       :: FOUND
  DATA TOLER /1.D-4/,ONE /1.D0/                                      
  !                                                                       
  ntot = sum(mult) ! all atoms
  
  ALLOCATE (rotij(3,3,ntot),tauij(3,ntot))
  
  WRITE(6,*) 'Inside ROTDEF'
  INDEX=0                                                           
  NCOUNT=0                                                          
  DO JATOM=1,NAT          ! 20
     INDEX1=INDEX+1                                                 
     DO M=1,MULT(JATOM)   ! 30
        INDEX=INDEX+1                                               
        FOUND = .FALSE.
        DO I=1,IORD       ! 25
           X=0.D0                                                      
           DO J=1,3       ! 26
              X=X+iz(J,1,I)*POS(J,INDEX1)                               
           ENDDO          ! 26
           X=X+TAU(1,I) + 1.D0                                           
           X= MOD (X,ONE)                                              
           Y=0.D0                                                      
           DO  J=1,3      ! 27
              Y=Y+IZ(J,2,I)*POS(J,INDEX1)                               
           ENDDO          ! 27
           Y=Y+TAU(2,I) + 1.D0                                           
           Y= MOD (Y,ONE)                                              
           Z=0.D0                                                      
           DO  J=1,3      ! 28
              Z=Z+IZ(J,3,I)*POS(J,INDEX1)                               
           ENDDO          ! 28
           Z=Z+TAU(3,I) + 1.D0                                           
           Z= MOD (Z,ONE)                                              
           X1=ABS(X-POS(1,INDEX))                                      
           Y1=ABS(Y-POS(2,INDEX))                                      
           Z1=ABS(Z-POS(3,INDEX))                                      
           !WRITE(*,'(I3,1x,I3,1x,I3)') JATOM, M, I
           !WRITE(*, '(f8.4,2x,f8.4,2x,f8.4,10x,f8.4,2x,f8.4,2x,f8.4,2x)') X, Y, Z, POS(1,INDEX), POS(2,INDEX), POS(3,INDEX)
           !WRITE(*, '(f8.4,2x,f8.4,2x,f8.4)') X1, Y1, Z1
           !            WRITE(*,*) 'JATOM,INDEX,I',JATOM,INDEX,I                   
           !            WRITE(*,*) X1,Y1,Z1                                        
           IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN
              
              NCOUNT=NCOUNT+1                                             
              DO J=1,3      ! 29
                 TAUIJ(J,INDEX)=TAU(J,I)                                    
                 DO  L=1,3  ! 29
                    ROTIJ(J,L,INDEX)=IZ(J,L,I)                                
                 ENDDO      ! 29
              ENDDO         ! 29
              
              FOUND = .TRUE.
              EXIT !GOTO 30                                                     
           END IF
           !....check positions for centered lattices
           if(lattic(1:1).eq.'B') then
              x1=mod(x1+0.500000001d0,one)
              y1=mod(y1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
                 NCOUNT=NCOUNT+1                                             
                 DO J=1,3       ! 129
                    TAUIJ(J,INDEX)=TAU(J,I)                                    
                    DO  L=1,3   ! 129
                       ROTIJ(J,L,INDEX)=IZ(J,L,I)                                
                    ENDDO       ! 129
                 ENDDO          ! 129
                 FOUND = .TRUE.
                 EXIT !GOTO 30                                                     
              END IF
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1=mod(x1+0.500000001d0,one)
              y1=mod(y1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
                 NCOUNT=NCOUNT+1                                             
                 DO J=1,3     ! 128
                    TAUIJ(J,INDEX)=TAU(J,I)                                    
                    DO L=1,3  ! 128                                                 
                       ROTIJ(J,L,INDEX)=IZ(J,L,I)                                
                    ENDDO     ! 128
                 ENDDO        ! 128
                 FOUND = .TRUE.
                 EXIT !GOTO 30                                                     
              END IF
              x1=mod(x1+0.5d0,one)
              y1=mod(y1+0.5d0,one)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1=mod(x1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
                 NCOUNT=NCOUNT+1                                             
                 DO  J=1,3        ! 127                                               
                    TAUIJ(J,INDEX)=TAU(J,I)                                    
                    DO  L=1,3     ! 127
                       ROTIJ(J,L,INDEX)=IZ(J,L,I)                                
                    ENDDO         ! 127
                 ENDDO            ! 127
                 FOUND = .TRUE.
                 EXIT  !GOTO 30                                                     
              END IF
              x1=mod(x1+0.5d0,one)
              z1=mod(z1+0.5d0,one)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              y1=mod(y1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
                 NCOUNT=NCOUNT+1                                             
                 DO J=1,3      ! 126
                    TAUIJ(J,INDEX)=TAU(J,I)                                    
                    DO L=1,3   ! 126
                       ROTIJ(J,L,INDEX)=IZ(J,L,I)                                
                    ENDDO      ! 126
                 ENDDO         ! 126
                 FOUND = .TRUE.
                 EXIT  !GOTO 30                                                     
              END IF
           end if

        ENDDO !25         CONTINUE

        IF (.NOT.FOUND) THEN
           !
           ! Error: no symmetry operation found
           !
           CALL OUTERR('ROTDEF','no symmetry operation found.')
           WRITE(ERRMSG,9000) jatom, index          
           CALL OUTERR('ROTDEF',ERRMSG)
           WRITE(ERRMSG,9010) (POS(I1,JATOM),I1=1,3) 
           CALL OUTERR('ROTDEF',ERRMSG)
           WRITE(ERRMSG,9020) (POS(I1,INDEX),I1=1,3) 
           CALL OUTERR('ROTDEF',ERRMSG)
           STOP 'ROTDEF - Error'
        ELSE
           WRITE(6,'(A,I3,1x,A,I3,1x,A,I3,1x)') 'Operation', I, 'transforms atom ', index, 'from first atom ', index1
           !WRITE(6,*) 'rotij='
           !do j=1,3
           !   do l=1,3
           !      WRITE(6,'(f14.8,1x)',advance='no') ROTIJ(j,l,INDEX)
           !   enddo
           !   WRITE(6,*)
           !enddo
           !WRITE(6,'(A,1x,f14.8,1x,f14.8,1x,f14.8)') 'tauij=', TAUIJ(:,INDEX)
           !WRITE(6,*)
        ENDIF
     ENDDO       !  30  
  ENDDO          !  20
  IF(NCOUNT.NE.INDEX) THEN

     CALL OUTERR('ROTDEF','ROTIJ not defined for all atoms of basis.')
     WRITE (ERRMSG,9030) NCOUNT
     CALL OUTERR('ROTDEF',ERRMSG)
     STOP 'ROTDEF - Error'
  ENDIF

  RETURN                                                            

  !        Error messages
  !
  !
9000 FORMAT ('for jatom, index',I2,i2)
9010 FORMAT ('atomposition of jatom',3F12.7)
9020 FORMAT ('atomposition of index',3F12.7)
9030 FORMAT ('NCOUNT=',I2)
END SUBROUTINE ROTDEF
