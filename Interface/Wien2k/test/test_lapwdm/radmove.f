USE param
SUBROUTINE RADMOVE(is,L)
	IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /ATSPDT/ EL(0:LMAX2),P(0:LMAX2),DP(0:LMAX2),PE(0:LMAX2), &
                      DPE(0:LMAX2),PEI(0:LMAX2) 
      COMMON /TMP/    PX(2),DPX(2),PEX(2),DPEX(2),ALOX(2),BLOX(2) &
                      ,CLOX(2)
      COMMON /RADFU/   RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2), &
                       RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)
      common /loabc/   alo(0:lomax),blo(0:lomax), &
                       clo(0:lomax), &
                       elo(0:lomax),plo(0:lomax), &
                       dplo(0:lomax),pelo(0:lomax), &
                       dpelo(0:lomax),peilo(0:lomax), &
                       pi12lo(0:lomax),pe12lo(0:lomax), &
                       a1lo(nrad,0:lomax),b1lo(nrad,0:lomax)
      COMMON /RFIS/    A(NRAD,2),B(NRAD,2), &
                       C(NRAD,2), &
                       AA(NRAD,2),BB(NRAD,2), &
                       CC(NRAD,2)

      DO IR=1,NRAD
	  A(IR,IS)=RRAD1(IR,L)
          AA(IR,IS)=RRAD2(IR,L)
          B(IR,IS)=RADE1(IR,L)
          BB(IR,IS)=RADE2(IR,L)
          C(IR,IS)=a1lo(IR,L)
          CC(IR,IS)=b1lo(IR,L)
      END DO
          PX(IS)=P(L)
	  DPX(IS)=DP(L)
 	  PEX(IS)=PE(L)
	  DPEX(IS)=DPE(L)
          IF (L.LE.LOMAX) THEN
	  ALOX(IS)=ALO(L)
	  BLOX(IS)=BLO(L)
	  CLOX(IS)=CLO(L)
          END IF
      END

