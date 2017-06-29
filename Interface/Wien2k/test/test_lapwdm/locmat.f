SUBROUTINE LOCMAT(ICASE,LL,XGL,XLO)
  USE param
  USE struct
  USE case
	IMPLICIT REAL*8 (A-H,O-Z)

	COMPLEX*16 XGL,XLO,TLOC,SUM,CZERO
	DIMENSION XGL(2*LXDOS+1,2*LXDOS+1), &
                  XLO(2*LXDOS+1,2*LXDOS+1), &
                  TLOC(2*LXDOS+1,2*LXDOS+1)

	DATA DD /1D0/, CZERO /(0D0,0D0)/
	
	pi=acos(-1.d0)
	CALL EULER(ROTLOC(1,1,IATOM(ICASE)),a,b,c)
        dd=1.d0
	CALL DMAT(LL,a,b,c,DD,TLOC)
         write(6,3)a*180/pi,b*180/pi,c*180/pi
 3   format('euler-locrot:',3f7.1)
	N=2*LL+1
	DO I=1,N
	DO J=1,N
	SUM=CZERO
	DO K=1,N
	DO L=1,N
        SUM=SUM+TLOC(k,i)*XGL(K,L)*CONJG(TLOC(l,j))
	END DO
	END DO
	XLO(I,J)=SUM
	END DO
	END DO

	END
