SUBROUTINE ATPAR (JATOM,itape,jtape,is,ISPIN,QoffDiag)
  USE param
  USE com_mpi
  USE struct
  IMPLICIT REAL*8 (A-H,O-Z)
  logical  QoffDiag  ! off-diagonal l needed
  LOGICAL    loor(0:lomax),lapw(0:lmax2),rlo(1:nloat,0:lomax)                                       
  DIMENSION  emist(0:lomax,nloat),E(0:LMAX2),elo(0:LOMAX,nloat),pei(0:lmax2),ilo(0:lomax),vr(nrad) 
  ! Output are the following quatities:
  !   /lolog/ , used later in LoMain
  !   /loabc/ , used later in LoMain
  !   ri_mat from /RINTEG/
  ! 
  !   P, DP from    /ATSPDT/
  !   RF1,RF2 from  /RADFU/
  !   
  COMMON /UHELP/   A(NRAD),B(NRAD),AE(NRAD),BE(NRAD)
  COMMON /ATSPDT/  P(0:LMAX2,2,nrf),DP(0:LMAX2,2,nrf)
  COMMON /RADFU/   RF1(NRAD,0:LMAX2,2,nrf),RF2(NRAD,0:LMAX2,2,nrf)
  COMMON /RINTEG/  RI_MAT(0:lmax2,0:lmax2,nrf,nrf,2,2)
  common /loabc/   alo(0:lomax,2,nloat,nrf)
  common /lolog/   nlo,nlov,nlon,loor,ilo,lapw  
!---------------------------------------------------------------------  
!.....READ TOTAL SPHERICAL POTENTIAL V(0,0) OF TAPEjtape=VSP               
!     NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)     
! itape points to the vector file, jtape to the potential file                        
  READ(jtape,1980)
  READ(jtape,2000)IDUMMY
  READ(jtape,2031)
  READ(jtape,2022)(VR(J),J=1,JRJ(JATOM))
  READ(jtape,2031)
  READ(jtape,2030)
  DO J=1,JRJ(JATOM)
     VR(J)=VR(J)/2.0D0  
  ENDDO
  if (ispin.eq.2 .and. (myrank.EQ.master .OR. fastFilesystem)) then
     if(is.eq.1)then
	write(6,*) 'ATPAR for spin up'
     else
        write(6,*) 'ATPAR for spin down'  
     endif
  end if
  nlo=0
  nlov=0
  nlon=0
  do i=0,lomax
     ilo(i)=0
  end do
  DO I=1,JATOM      ! 15
     READ(itape)E   ! linearization energies
     READ(itape)elo ! local
     IF(i.EQ.jatom) THEN
        DO l=0,lmax2
           lapw(l)=.TRUE. ! local
           IF(e(l).GT.150.) THEN
              e(l)=e(l)-200.d+0
              lapw(l)=.FALSE.
           ENDIF
        ENDDO
     ENDIF
     rlo(:,:)=.FALSE. ! local
     do l = 0,lomax ! 15
        loor(l)=.FALSE. ! local
        DO k=1,nloat
           IF (i.EQ.jatom) THEN
              IF (elo(l,k).LT.(995.d+0)) THEN
                 ilo(l)=ilo(l)+1
                 nlo=nlo+((2*l+1))*mult(i)
                 IF(.NOT.lapw(l).AND.k.EQ.1) CYCLE !GOTO 666
                 IF(k.EQ.nloat) THEN
                    rlo(ilo(l),l)=.TRUE.
                    CYCLE
                    !GOTO 666
                 ENDIF
                 loor(l)=.TRUE.
                 !666                 CONTINUE
              ENDIF
           ELSE
              IF (elo(l,k).LT.(995.d+0)) nlov=nlov+((2*l+1))*mult(i)
           ENDIF
        ENDDO
     enddo!15    continue
  ENDDO
  IF(JATOM.NE.NAT) THEN !GOTO 30
     DO I=JATOM+1,NAT   ! 20
        READ(itape) EMIST                                                    
        READ(itape) EMIST                                                    
        do l=0,lomax    ! 20
           do k=1,nloat ! 20
              if (emist(l,k).lt.(995.0D+0)) nlon=nlon+((2*l+1))*mult(i)
           enddo !20  continue
        enddo
     ENDDO
  ENDIF !30  CONTINUE     

  if (myrank.EQ.master .OR. fastFilesystem) then
     WRITE(6,7) ANAME(JATOM)
     WRITE(6,5) E
     WRITE(6,14)
  endif
!                                                                       
  DO l=0,LMAX2  ! 70
     DELE=2.0D-3
     DELEI=0.25D0/DELE
     FL=L
     EI=E(l)/2.0d0
     !     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE                  
     !     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES          
     !                                                                       
     E1=EI-DELE
     ! OUTPUT IS: UVB,DUVB,NODEL
     CALL OUTWIN(REL,VR,r0(JATOM),DX(JATOM),jrj(JATOM),E1,FL,UVB,DUVB,NODEL,zz(jatom))
     ! OUTPUT IS OVLP
     CALL rint13(A,B,A,B,OVLP,JATOM)
     TRX=1.0D0/SQRT(OVLP)
     IMAX=jrj(JATOM)
     DO M=1,IMAX
        AE(M)=TRX*A(M)
        BE(M)=TRX*B(M)
     ENDDO !45 CONTINUE                                                          
     UVB=TRX*UVB                                                       
     DUVB=TRX*DUVB                                                     
     E1=EI+DELE                                                        
     ! OUTPUT IS: UVE,DUVE,NODE
     CALL OUTWIN(REL,VR,r0(JATOM),DX(JATOM),jrj(JATOM),E1,FL,UVE,DUVE,NODE,zz(jatom))
     CALL rint13(A,B,A,B,OVLP,JATOM)
     TRX=1.0d0/SQRT(OVLP)
     UVE=DELEI*(TRX*UVE-UVB)
     DUVE=DELEI*(TRX*DUVE-DUVB)
     IMAX=jrj(JATOM)
     DO M=1,IMAX ! 50
        AE(M)=DELEI*(TRX*A(M)-AE(M))
        BE(M)=DELEI*(TRX*B(M)-BE(M))
     ENDDO !50 CONTINUE
     ! AE and BE set! 
!
!     CALCULATE FUNCTION AT EI
!
     ! OUTPUT IS: UV,DUV,NODES
     CALL OUTWIN(REL,VR,r0(JATOM),DX(JATOM),jrj(JATOM),EI,FL,UV,DUV,NODES,zz(jatom))
     CALL rint13(A,B,A,B,OVLP,JATOM)
     TRX=1.0d0/SQRT(OVLP)
     P(l,is,1)=TRX*UV
     DP(l,is,1)=TRX*DUV
     IMAX=jrj(JATOM)
     DO M=1,IMAX ! 60
        A(M)=TRX*A(M)
        B(M)=TRX*B(M)
     ENDDO ! 60
     ! A and B set
!                                                                       
!     INSURE ORTHOGONALIZATION                                          
!                                                                       
     CALL rint13(A,B,AE,BE,CROSS,JATOM)
     TRY=-CROSS
     IMAX=jrj(JATOM)
     DO M=1,IMAX ! 55
        AE(M)=(AE(M)+TRY*A(M))
        BE(M)=(BE(M)+TRY*B(M))
     ENDDO ! 55
     IMAX=jrj(JATOM)
     DO I=1,IMAX  ! 80
        RF1(I,l,is,1)=A(I)                                                   
        RF2(I,l,is,1)=B(I)                                                   
        RF1(I,l,is,2)=AE(I)                                                  
        RF2(I,l,is,2)=BE(I)                                                  
     ENDDO   !80    CONTINUE                                                          
     P(l,is,2)=UVE+TRY*P(l,is,1)                                                
     DP(l,is,2)=DUVE+TRY*DP(l,is,1) 
     CALL rint13(AE,BE,AE,BE,PEI(l),JATOM) ! PEI: What is it good for
     if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,8) L,P(l,is,1),DP(l,is,1),P(l,is,2),DP(l,is,2)
  ENDDO ! 70
! Results saved for P(l,is,irf),DP(l,is,irf)
!                         
! nun fur lo
!
  DO l=0,lomax ! 170
     irf=2
     do jlo=1,ilo(l) ! 170
        if (.not.((.not.lapw(l)).and.(jlo.eq.1))) then ! goto 180
           irf=irf+1
           DELE=2.0D-3                                                       
           DELEI=0.25D0/DELE                                                 
           FL=L                                                              
           EI=elo(l,jlo)/2.d0                                                         
!                                                                       
!     CALCULATE FUNCTION AT EI 
           IF(rlo(jlo,l)) THEN
              ei=elo(l,nloat)/2.d0
              kappa=l
              ! output: uv,duv,nodes
              CALL diracout(rel,vr,r0(jatom),dx(jatom),jrj(jatom),ei,fl,kappa,uv,duv,nodes,zz(jatom))
              ! output: b
              CALL dergl(a,b,r0(jatom),dx(jatom),jrj(jatom))
              DO m = 1, jrj(jatom)
                 r_m = r0(jatom)*exp(dx(jatom)*(m-1))
                 b(m) = b(m)*r_m/(2.d0*clight+(elo(l,jlo)-2.d0*vr(m)/r_m)/(2.d0*clight))
                 b(m)=b(m)*clight
              ENDDO
           ELSE
              CALL outwin(rel,vr,r0(jatom),dx(jatom),jrj(jatom),ei,fl,uv,duv,nodes,zz(jatom))
           ENDIF
!                                                                       
           CALL rint13(A,B,A,B,OVLP,JATOM)
           TRX=1.0d0/SQRT(OVLP)
           P(l,is,irf)=TRX*UV
           DP(l,is,irf)=TRX*DUV
           IMAX=jrj(JATOM)
           DO M=1,IMAX  ! 160
              rf1(M,l,is,irf)=TRX*A(M)
              rf2(M,l,is,irf)=TRX*B(M)
           ENDDO !   160
                                                    
           CALL rint13(rf1(1,l,is,1),rf2(1,l,is,1),rf1(1,l,is,irf),rf2(1,l,is,irf),pi12lo,JATOM)
           CALL rint13(rf1(1,l,is,2),rf2(1,l,is,2),rf1(1,l,is,irf),rf2(1,l,is,irf),pe12lo,JATOM)
           !180 continue         
        endif
        call abcd(l,jatom,pei(l),pi12lo,pe12lo,is,jlo,lapw(l))
        if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,*)
     enddo
  ENDDO ! 170  continue
!... CALCULATION OF RADIAL INTEGRALS
  IF (.NOT.(ISPIN.EQ.2.AND.IS.EQ.1)) THEN !  GOTO 469
     CALL RADINT(JATOM,ISPIN,QoffDiag) ! computes ri_mat from rf1 and rf2
  ENDIF
  !469 
  RETURN                                                            
                                                                       
2022 FORMAT(3X,4E19.12)
2 FORMAT(5E14.7)                                                    
3 FORMAT(16I5)                                                      
4 FORMAT(I4,E16.7)                                                  
5 FORMAT(10X,' ENERGY PARAMETERS ARE',7F7.2)                        
6 FORMAT(10X,'E(',I2,2H)=,F10.4)                                    
7 FORMAT(/10X,'ATOMIC PARAMETERS FOR ',A10/)                        
14 FORMAT(/11X,1HL,5X,4HU(R),10X,5HU'(R),9X,5HDU/DE,8X,6HDU'/DE,6X,7HNORM-U')
8 FORMAT(10X,I2,5E14.6,5X,3I2)
9 FORMAT(10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)
11 FORMAT(7F10.5)
13 FORMAT(////,':POS',i2.2,':',1x,'AT.NR.',I3,2X,'POSITION =',3F8.5,2X,'MULTIPLICITY =',I3)
1040 FORMAT(8I2)
1050 FORMAT(A10,I5,5X,2F10.9,I5,F5.2)
1060 FORMAT(//,3X,'NOT EQUIV ATOM ',A10,'  LOCAL ROTATION MATRIX')
1070 FORMAT(30X,3F10.5)
1080 FORMAT(13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)
1090 FORMAT(30X,3F10.5,5X,F10.5)
1980 FORMAT(3X)
2000 FORMAT(16X,I2//)
2030 FORMAT(///)
2031 FORMAT(/)
END SUBROUTINE ATPAR
