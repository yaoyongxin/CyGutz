      SUBROUTINE ROTDEF
!                                                                       
!     ROTDEF GENERATES THE ROTATION-MATRICES ROTLOC(3,3,JATOM) FOR      
!     NONSYMMORPHIC STRUCTURES. THESE MATRICES PERFORM ROTATIONS        
!     FROM THE GENERAL COORDINATION-SYSTEM INTO LOCAL SYSTEMS WITH      
!     SPECIFIED POINTSYMMETRY.                                          
!     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN      
!     ATOM TO IT'S CORRESPONDING POSITION OF AN EQUIVALENT              
!     ATOM.                                                             
        USE param
        USE struct
        USE rotat
      IMPLICIT REAL*8 (A-H,O-Z)
!                                                                       
!
      CHARACTER*67       ERRMSG
!
      DATA TOLER/1.D-4/,ONE/1.D0/                                      
      toler2=toler/2.d0
!-----------------------------------------------------------------------
!                                                                       

        IF (LATTIC(1:1).EQ.'B') THEN
        GOTO 1000
	ELSE IF (LATTIC(1:1).EQ.'F') THEN
        GOTO 2000
        ELSE IF (LATTIC(1:3).EQ.'CXY') THEN
        GOTO 3000
        ELSE IF (LATTIC(1:3).EQ.'CYZ') THEN
        GOTO 4000
        ELSE IF (LATTIC(1:3).EQ.'CXZ') THEN
        GOTO 5000
        ELSE
        GOTO 6000
        END IF 

	RETURN

 1000  CONTINUE
	LFIRST=1
        DO 1200 JATOM=1,NAT
	if(jatom.gt.1)LFIRST=LFIRST+MULT(JATOM-1)
            DO 1200 I=1,IORD
            DO 1300 MU=1,MULT(JATOM)
            INDEX=LFIRST+MU-1
            X=pos(1,index)
            Y=pos(2,index)
            Z=pos(3,index)
            X1=0.D0
            DO J=1,3
            X1=X1+IZ(J,1,I)*POS(J,lfirst)
            END DO
            X1=X1+TAU(1,I) + 5.D0
            X1= MOD (X1+toler2,ONE)-toler2
            Y1=0.D0
            DO J=1,3
            Y1=Y1+IZ(J,2,I)*POS(J,lfirst)
            END DO
            Y1=Y1+TAU(2,I) + 5.D0
            Y1= MOD (Y1+toler2,ONE)-toler2
            Z1=0.D0
            DO J=1,3
            Z1=Z1+IZ(J,3,I)*POS(J,lfirst)
            END DO
            Z1=Z1+TAU(3,I) + 5.D0
            Z1= MOD (Z1+toler2,ONE)-toler2
            X2= MOD (X1+0.5,ONE)
            Y2= MOD (Y1+0.5,ONE)
            Z2= MOD (Z1+0.5,ONE)

            tt=ABS(X-X1) &
             +ABS(Y-Y1) &
             +ABS(Z-Z1)
	    if (tt.lt.toler) then
	    itr(i,jatom)=MU
            goto 1200
	    end if

            tt=ABS(X-X2) &
             +ABS(Y-Y2) &
             +ABS(Z-Z2)
            if (tt.lt.toler) then
            itr(i,jatom)=MU
            goto 1200
            end if

 1300       CONTINUE
 1200       CONTINUE
 	    RETURN

 2000  CONTINUE
        LFIRST=1
        DO 2200 JATOM=1,NAT
        if(jatom.gt.1)LFIRST=LFIRST+MULT(JATOM-1)
            DO 2200 I=1,IORD
            ITR(I,JATOM)=0
            DO 2300 MU=1,MULT(JATOM)
            INDEX=LFIRST+MU-1
            X=pos(1,index)
            Y=pos(2,index)
            Z=pos(3,index)
            X1=0.D0
            DO J=1,3
            X1=X1+IZ(J,1,I)*POS(J,lfirst)
            END DO
            X1=X1+TAU(1,I) + 5.D0
            X1= MOD (X1+toler2,ONE)-toler2
            Y1=0.D0
            DO J=1,3
            Y1=Y1+IZ(J,2,I)*POS(J,lfirst)
            END DO
            Y1=Y1+TAU(2,I) + 5.D0
            Y1= MOD (Y1+toler2,ONE)-toler2
            Z1=0.D0
            DO J=1,3
            Z1=Z1+IZ(J,3,I)*POS(J,lfirst)
            END DO
            Z1=Z1+TAU(3,I) + 5.D0
            Z1= MOD (Z1+toler2,ONE)-toler2

            X2= MOD (X1+0.5,ONE)
            Y2= MOD (Y1+0.5,ONE)
            Z2= MOD (Z1+0.0,ONE)

            X3= MOD (X1+0.0,ONE)
            Y3= MOD (Y1+0.5,ONE)
            Z3= MOD (Z1+0.5,ONE)

            X4= MOD (X1+0.5,ONE)
            Y4= MOD (Y1+0.0,ONE)
            Z4= MOD (Z1+0.5,ONE)

            tt=ABS(X-X1) &
             +ABS(Y-Y1) &
             +ABS(Z-Z1)
            if (tt.lt.toler) then
            itr(i,jatom)=MU
            goto 2200
            end if

            tt=ABS(X-X2) &
             +ABS(Y-Y2) &
             +ABS(Z-Z2)
            if (tt.lt.toler) then
            itr(i,jatom)=MU
            goto 2200
            end if

            tt=ABS(X-X3) &
             +ABS(Y-Y3) &
             +ABS(Z-Z3)
            if (tt.lt.toler) then
            itr(i,jatom)=MU
            goto 2200
            end if

            tt=ABS(X-X4) &
             +ABS(Y-Y4) &
             +ABS(Z-Z4)
            if (tt.lt.toler) then
            itr(i,jatom)=MU
            goto 2200
            end if

 2300       CONTINUE
 2200       CONTINUE
	    RETURN

 3000  CONTINUE
        LFIRST=1
        DO 3200 JATOM=1,NAT
        if (jatom.gt.1)LFIRST=LFIRST+MULT(JATOM-1)
            II=0
            DO 3200 I=1,IORD
            DO 3300 MU=1,MULT(jatom)
            INDEX=LFIRST+MU-1
            X=pos(1,index)
            Y=pos(2,index)
            Z=pos(3,index)
            X1=0.D0
            DO J=1,3
            X1=X1+IZ(J,1,I)*POS(J,lfirst)
            END DO
            X1=X1+TAU(1,I) + 5.D0
            X1= MOD (X1+toler2,ONE)-toler2
            Y1=0.D0
            DO J=1,3
            Y1=Y1+IZ(J,2,I)*POS(J,lfirst)
            END DO
            Y1=Y1+TAU(2,I) + 5.D0
            Y1= MOD (Y1+toler2,ONE)-toler2
            Z1=0.D0
            DO J=1,3
            Z1=Z1+IZ(J,3,I)*POS(J,lfirst)
            END DO
            Z1=Z1+TAU(3,I) + 5.D0
            Z1= MOD (Z1+toler2,ONE)-toler2

            X2= MOD (X1+0.5,ONE)
            Y2= MOD (Y1+0.5,ONE)
            Z2= MOD (Z1+0.0,ONE)

            tt=ABS(X-X1) &
             +ABS(Y-Y1) &
             +ABS(Z-Z1)
            if (tt.lt.toler) then
            itr(i,jatom)=mu
            goto 3200
            end if

            tt=ABS(X-X2) &
             +ABS(Y-Y2) &
             +ABS(Z-Z2)
            if (tt.lt.toler) then
            itr(i,jatom)=mu
            goto 3200
            end if

 3300       CONTINUE
 3200       CONTINUE
	    RETURN

 4000  CONTINUE
        LFIRST=1
        DO 4200 JATOM=1,NAT
        if (jatom.gt.1) LFIRST=LFIRST+MULT(JATOM-1)
            DO 4200 I=1,IORD
            DO 4300 MU=1,mult(jatom)
            INDEX=LFIRST+MU-1
            X=pos(1,index)
            Y=pos(2,index)
            Z=pos(3,index)
            X1=0.D0
            DO J=1,3
            X1=X1+IZ(J,1,I)*POS(J,lfirst)
            END DO
            X1=X1+TAU(1,I) + 5.D0
            X1= MOD (X1+toler2,ONE)-toler2
            Y1=0.D0
            DO J=1,3
            Y1=Y1+IZ(J,2,I)*POS(J,lfirst)
            END DO
            Y1=Y1+TAU(2,I) + 5.D0
            Y1= MOD (Y1+toler2,ONE)-toler2
            Z1=0.D0
            DO J=1,3
            Z1=Z1+IZ(J,3,I)*POS(J,lfirst)
            END DO
            Z1=Z1+TAU(3,I) + 5.D0
            Z1= MOD (Z1+toler2,ONE)-toler2

            X2= MOD (X1+0.0,ONE)
            Y2= MOD (Y1+0.5,ONE)
            Z2= MOD (Z1+0.5,ONE)

            tt=ABS(X-X1) &
             +ABS(Y-Y1) &
             +ABS(Z-Z1)
            if (tt.lt.toler) then
            itr(i,jatom)=mu
            goto 4200
            end if

            tt=ABS(X-X2) &
             +ABS(Y-Y2) &
             +ABS(Z-Z2)
            if (tt.lt.toler) then
            itr(i,jatom)=mu
            goto 4200
            end if
 4300       CONTINUE
 4200       CONTINUE
	    RETURN

 5000  CONTINUE
        LFIRST=1
        DO 5200 JATOM=1,NAT
        if(jatom.gt.1)LFIRST=LFIRST+MULT(JATOM-1)
            II=0
            DO 5200 I=1,IORD
            DO 5300 MU=1,mult(jatom)
            index=lfirst+mu-1
            X=pos(1,index)
            Y=pos(2,index)
            Z=pos(3,index)
            X1=0.D0
            DO J=1,3
            X1=X1+IZ(J,1,I)*POS(J,lfirst)
            END DO
            X1=X1+TAU(1,I) + 5.D0
            X1= MOD (X1+toler2,ONE)-toler2
            Y1=0.D0
            DO J=1,3
            Y1=Y1+IZ(J,2,I)*POS(J,lfirst)
            END DO
            Y1=Y1+TAU(2,I) + 5.D0
            Y1= MOD (Y1+toler2,ONE)-toler2
            Z1=0.D0
            DO J=1,3
            Z1=Z1+IZ(J,3,I)*POS(J,lfirst)
            END DO
            Z1=Z1+TAU(3,I) + 5.D0
            Z1= MOD (Z1+toler2,ONE)-toler2

            X2= MOD (X1+0.5,ONE)
            Y2= MOD (Y1+0.0,ONE)
            Z2= MOD (Z1+0.5,ONE)

            tt=ABS(X-X1) &
             +ABS(Y-Y1) &
             +ABS(Z-Z1)
            if (tt.lt.toler) then
            itr(i,jatom)=mu
            goto 5200
            end if

            tt=ABS(X-X2) &
             +ABS(Y-Y2) &
             +ABS(Z-Z2)
            if (tt.lt.toler) then
            itr(i,jatom)=mu
            goto 5200
            end if

 5300       CONTINUE
 5200       CONTINUE
	    RETURN

 6000  CONTINUE
        LFIRST=1
        DO 6200 JATOM=1,NAT
        if(jatom.gt.1)LFIRST=LFIRST+MULT(JATOM-1)
            DO 6200 I=1,IORD
            DO 6300 MU=1,MULT(jatom)
            index=lfirst+mu-1
            X=pos(1,index)
            Y=pos(2,index)
            Z=pos(3,index)
            X1=0.D0
            DO J=1,3
            X1=X1+IZ(J,1,I)*POS(J,lfirst)
            END DO
            X1=X1+TAU(1,I) + 5.D0
            X1= MOD (X1+toler2,ONE)-toler2
            Y1=0.D0
            DO J=1,3
            Y1=Y1+IZ(J,2,I)*POS(J,lfirst)
            END DO
            Y1=Y1+TAU(2,I) + 5.D0
            Y1= MOD (Y1+toler2,ONE)-toler2
            Z1=0.D0
            DO J=1,3
            Z1=Z1+IZ(J,3,I)*POS(J,lfirst)
            END DO
            Z1=Z1+TAU(3,I) + 5.D0
            Z1= MOD (Z1+toler2,ONE)-toler2

            tt=ABS(X-X1) &
             +ABS(Y-Y1) &
             +ABS(Z-Z1)
            if (tt.lt.toler) then
            itr(i,jatom)=mu
            end if

 6300       CONTINUE
 6200       CONTINUE
            RETURN                                                            

 9000 FORMAT ('for jatom, index',I2,i2)
 9010 FORMAT ('atomposition of jatom',3F12.7)
 9020 FORMAT ('atomposition of index',3F12.7)
 9030 FORMAT ('NCOUNT=',I2)
      END                                                               
