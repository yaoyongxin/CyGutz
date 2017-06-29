SUBROUTINE ORB(LL,XGL,tx,ty,tz)
USE param
        IMPLICIT REAL*8 (A-H,O-Z)

        COMPLEX*16 XGL,TLOC
        DIMENSION XGL(2*LXDOS+1,2*LXDOS+1), &
                  TLOC(2*LXDOS+1,2*LXDOS+1)


        DATA DD /1D0/,ZERO /0D0/
	PI=ACOS(-1.D0)

        N=2*LL+1
	a=0.0
	b=pi/2
	c=pi
        CALL DMAT(LL,a,b,c,DD,TLOC)
        SUM=ZERO
	DO I=1,N
        DO K=1,N
        DO L=1,N
!        SUM=SUM+REAL(TLOC(l,i)*XGL(K,L)*CONJG(TLOC(k,i)))*(I-LL-1)
           SUM=SUM+DBLE(TLOC(K,i)*XGL(K,L)*CONJG(TLOC(L,i)))*(I-LL-1)
        END DO
        END DO
        END DO
        TX=SUM

        a=0.0
        b=pi/2
        c=pi/2
        CALL DMAT(LL,a,b,c,DD,TLOC)
        SUM=ZERO
        DO I=1,N
        DO K=1,N
        DO L=1,N
!        SUM=SUM+REAL(TLOC(l,i)*XGL(K,L)*CONJG(TLOC(k,i)))*(I-LL-1)
           SUM=SUM+DBLE(TLOC(K,i)*XGL(K,L)*CONJG(TLOC(L,i)))*(I-LL-1)
        END DO
        END DO
        END DO
        TY=SUM

        SUM=ZERO
        DO I=1,N
        SUM=SUM+DBLE(XGL(I,I))*(I-LL-1)
        END DO
        TZ=SUM

        write(6,55)TX,TY,TZ
!        write(21,55)TX,TY,TZ
55      format(3f11.6,' Lx,Ly,Lz in global orthog. system')
	END




