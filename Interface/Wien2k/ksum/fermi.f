subroutine eord(e,eb,nemax,k,isospin) 
  USE param
  IMPLICIT REAL*8 (A-H,O-Z)                                        
  dimension eb(nume,nkpt,isospin),e(nemax*isospin,k)
  do k1=1,k
     ne2=0
     do ispin=1,isospin
        do ne1=1,nemax
           ne2=ne2+1
           e(ne2,k1)=eb(ne1,k1,ispin)
        enddo
     enddo
  enddo
  return

end subroutine eord


subroutine eweigh(ef,weight,nemax,k,jspin,ne,eb,nbmax)       
  USE param
  USE com_mpi
  USE bandm
  IMPLICIT REAL*8 (A-H,O-Z)                                        
  dimension eb(nume,nkpt,2)
  dimension weight(nemax*jspin,k),ne(k)
  common /correct/ cordeg,icor
  TEST1=2.D-8                                                      
  nbmax=0  
  emax=-10.d0                                                      
  DO 8 KK=1,K                                                      
     NNN=NE(KK)             
     NNN=NEmax    
     !
     ! ### 2007-08-31, Clas Persson, Juergen Spitaler, and Claudia Ambrosch-Draxl
     ! We suggest to change on page 90 in the usersguide.pdf [note, now abs(eval)]:
     !
     !                                     ...    abs(eval).gt.100 specifies
     ! the use of the standard tetrahedron method instead of the modified one
     ! (see above). Using eval.lt.0 in combination with TETRA forces the
     ! occupancy of degenerate states to be equal in order to avoid incorrect
     ! split of these states, which otherwise may occur for partially filled
     ! states. Degeneracy is here determined by a tolerance of abs(eval) for
     ! -1.gt.eval.lt.0, and of 1e-6 Ry for eval.le.-1.
     !
     if(cordeg.lt.0.d0) then
        if(kk.eq.1) ifcpt=1
        do jspin1=1,jspin
           nn=1
           do while(nn.le.nnn)
              ifcp=0
              ndeg=1
              wecp=weight(nn+(jspin1-1)*nnn,kk)
              !
              ! check degeneracy
              !        do while((nn+ndeg).le.nnn.and. &
              !                 abs(eb(nn+ndeg,kk,jspin1)-eb(nn,kk,jspin1)).lt.abs(cordeg))
              do while((nn+ndeg).le.nnn)
                 if(abs(eb(nn+ndeg,kk,jspin1)-eb(nn,kk,jspin1)).gt.abs(cordeg)) exit
                 if(abs(weight(nn+ndeg+(jspin1-1)*nnn,kk)-weight(nn+(jspin1-1)*nnn,kk)).ge.1e-6) ifcp=1
                 wecp=wecp+weight((nn+ndeg)+(jspin1-1)*nnn,kk)
                 ndeg=ndeg+1
              enddo
              !
              ! equalizes occupancy and outputs to case.output2
              if(ifcp.eq.1) then
                 if(ifcpt.eq.1) then
                    if (myrank.EQ.master .OR. fastFilesystem) write(6,'("k-pnt spin band  energy          old/new occ.    ")')
                    ifcpt=0
                 endif
                 do icp=0,ndeg-1
                    if (myrank.EQ.master .OR. fastFilesystem) write(6,50)kk,jspin1,nn,eb(nn+icp,kk,jspin1),weight((nn+icp)+(jspin1-1)*nnn,kk)
                    weight((nn+icp)+(jspin1-1)*nnn,kk) = wecp/dble(ndeg)
                    if (myrank.EQ.master .OR. fastFilesystem) write(6,'(f9.6)') weight((nn+icp)+(jspin1-1)*nnn,kk)
                 enddo
              endif
              nn=nn+ndeg
           enddo
        enddo
     endif
50   format(3i5,f10.6,2x,f9.6,"/",$)
! ### END equalizes occupancy of degenerate states

!
     DO 8 NN=1,NNN   
        do jspin1=1,jspin
           if(abs(weight(nn+(jspin1-1)*nnn,kk)).gt.test1) then
              nbmax=max(nbmax,nn)
              emax=max(emax,eb(nn,kk,jspin1))                      
           end if
        enddo
   8 CONTINUE
     if (myrank.EQ.master .OR. fastFilesystem) write(6,*) '  number of occupied bands:',nbmax
     if (myrank.EQ.master .OR. fastFilesystem) write(6,*) '  highest energy:',emax
     if(ef.gt.emax.and.abs(emax-ebmin(min(nbmax+1,nume))).gt.1.d-3) then
        if (myrank.EQ.master .OR. fastFilesystem) write(6,*) 'insulator !'
        if((ef-emax).lt.1.d-4.or.(ef-emax).ge.1.d-4) then
           if (myrank.EQ.master .OR. fastFilesystem) write(6,*) 'EF-inconsistency corrected'
           if (myrank.eq.master) write(21,888) 
888        format('       Insulator, EF-inconsistency corrected')
           ef=emax
        endif
     endif
     return

end subroutine eweigh


  !     -----------------------------------------------------------------
  !     -----------------------------------------------------------------
  !     ----                                                         ----
  !     ----  BLOCK KINTR                                            ----
  !     ----  CALCULATION OF SAMPLING WEIGHTS                        ----
  !     ----                                                         ----
  !     -----------------------------------------------------------------
  !     -----------------------------------------------------------------
  !     .....................................................DOS.........
SUBROUTINE DOS(NB,NKP,EB,WGHT,RNTOT,D,EF,W,NWX)                  
  USE com_mpi
  USE param
  !USE parallel
  !     **                                                              *
  !     **  CALCULATES THE SAMPLING WEIGHTS FROM TETRAHEDRON INTEGRATION*
  !     **                                                              *
  !     **  INPUT :                                                     *
  !     **    NB          NUMBER OF BANDS                               *
  !     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                *
  !     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )          *
  !     **    RNTOT       NUMBER OF OCCUPIED STATES                     *
  !     **    W           (INTEGER) WORK ARRAY                          *
  !     **    NWX         LENGTH OF WROK ARRAY W                        *
  !     **  OUTPUT :                                                    *
  !     **    WGHT        SAMPLING WEIGHTS                              *
  !     **    EF          FERMI LEVEL                                   *
  !     **                                                              *
  !     **  AUTHOR : PETER E. BLOECHL                                   *
  !     **                                                              *
  !     **  SUBROUTINES USED:                                           *
  !     **  EFERMI,SAMFAC,DEF0,TETR0,INITDR,TETR1,TOTNOS,EFI,WEIGHT
  !     **  DRVAL
  !     **                                                           **KI
  IMPLICIT DOUBLE PRECISION (A-H,P-Z)                              
  IMPLICIT INTEGER (O)                                             
  DIMENSION EB(NB*NKP),WGHT(NB*NKP)                                
  CHARACTER *67    ERRMSG
  INTEGER W(NWX)                                                   
  DATA MWRIT/100/ICHECK/1/TOLMAX/1.D-5/                            
  ! a larger TOLMAX may not lead to the correct number of electrons!
  ! eventually one could issue just a warning, but continue if normalization
  ! is only slightly violated.
  !     -----------------------------------------------------------------
  !     --  CALCULATE FERMI LEVEL                                       -
  !     -----------------------------------------------------------------
  !      write(*,*) NB,NKP,RNTOT,NWX
  !      write(*,*) eb(1)
  !           write(6,*) 'before efermi in dos'
  CALL EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,W,NWX)                     
                                                                       
  if (myrank.EQ.master .OR. fastFilesystem) write(6,*) '  FERMI ENERGY AT ',EF                                   
  !     -----------------------------------------------------------------
  !     --  CALCULATE WEIGHTS                                           -
  !     -----------------------------------------------------------------
  CALL SAMFAC(NB,NKP,EB,EF,WGHT,W,NWX)                             
                                                                       
  IF(ICHECK.EQ.0) RETURN                                           
  !     -----------------------------------------------------------------
  !     --  CHECK WHETHER SUMRULE IS FULLFILLED                         -
  !     -----------------------------------------------------------------
  SUM=0.D0                                                         
  DO I=1,NB*NKP                                                
     SUM=SUM+WGHT(I)                                                  
  ENDDO
! YYX begin 
  if (myrank.EQ.master .OR. fastFilesystem) WRITE (6,9000) SUM,RNTOT
  IF(DABS(SUM-RNTOT).GT.1.E-9) THEN
!  IF(DABS(SUM-RNTOT).GT.TOLMAX) THEN                               
!     if (myrank.EQ.master .OR. fastFilesystem) WRITE (6,9000) SUM,RNTOT
     if (myrank.eq.master) write(21,9001) SUM,RNTOT
  ELSE                                                             
     if (myrank.EQ.master .OR. fastFilesystem) write(6,*) '  SUM RULE OF TETRAHEDRON INTEGRATION CHECKED '
  END IF
! YYX end
  RETURN                                                           
900 CALL OUTERR('FERMI',' INTEGRATION FAILED.....STOP IN DOS')
  WRITE (ERRMSG,9000) SUM,RNTOT
  CALL OUTERR('FERMI',ERRMSG)
  STOP  'FERMI - Error'
9000 FORMAT(' RESULT OF INTEGRATION: ',f20.12, '; SHOULD BE: ',f20.12)
9001 FORMAT(':WARN : RESULT OF INTEGRATION: ',f20.12, '; SHOULD BE: ',f20.12)
END SUBROUTINE DOS

!     .....................................................EFERMI......
SUBROUTINE EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,W,NWX)               
  !     **                                                              *
  !     **  CALCUALTES THE FERMILEVEL BY INTEGRATION OF THE TOTAL       *
  !     **  DENSITY OF STATES                                           *
  !     **                                                              *
  !     **  INPUT :                                                     *
  !     **    RNTOT       NUMBER OF OCCUPIED STATES                     *
  !     **    NB          NUMBER OF BANDS                               *
  !     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                *
  !     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )          *
  !     **    W           (INTEGER) WORK ARRAY                          *
  !     **    NWX         LENGTH OF WROK ARRAY W                        *
  !     **    TOLMAX      TOLERANCE IN THE NUMBER OF STATES AT EF       *
  !     **  OUTPUT :                                                    *
  !     **    EF          FERMI LEVEL                                   *
  !     **                                                              *
  USE com_mpi
  USE param
  IMPLICIT DOUBLE PRECISION (A-H,P-Z)                              
  IMPLICIT INTEGER (O)                                             
  DIMENSION EB(NB,NKP),E(4),IKP(4)                                 
  INTEGER W(NWX)                                                   
  CHARACTER*67 ERRMSG
  DATA NP/1000/                                                    
  CALL DEF0(NWX)                                                   
  !     -----------------------------------------------------------------
  !     --  FIND EMIN EMAX     (ENERGYBANDS ARE ASSUMED                 -
  !     --                      TO BE ORDERED WITH RESPECT TO SYMMETRY  -
  !     -----------------------------------------------------------------
      EMIN=EB(1,1)                                                     
      EMAX=EB(NB,1)                                                    
      DO 100 IK=1,NKP                                                  
      DO 100 IB=1,NB                                                   
      EMIN=DMIN1(EMIN,EB(IB,IK))
      EMAX=DMAX1(EMAX,EB(IB,IK))
100   CONTINUE      
      DE=(EMAX-EMIN)/DBLE(NP-1)
      EMAX=EMAX+DE  
      EMIN=EMIN-DE 
      CALL DEFDR(ONOS,NP)
      CALL DEFDR(OSOS,NP)
      CALL TETR0(VOL1,NKP,NTET,MWRIT) 
      CALL DEFI(OWORK,5*MWRIT)                                         
      ILOOP=0                                                          
1000  CONTINUE                                                         
      ILOOP=ILOOP+1                                                    
!      write(*,*) 'in efermi',iloop
      CALL INITDR(W(ONOS),NP)                                          
      CALL INITDR(W(OSOS),NP)                                          
      INIT=1                                                           
      SUM=0.D0                                                         
      DO 200 ITET=1,NTET                                               
      CALL TETR1(INIT,ITET,IWGHT,IKP,MWRIT,W(OWORK))                   
      VOL=VOL1*DBLE(IWGHT)                                             
      SUM=SUM+VOL                                                      
      DO 220 IB=1,NB                                                   
      DO 210 I=1,4                                                     
      E(I)=EB(IB,IKP(I))                                               
210   CONTINUE                                                         
!c      if(iloop.eq.4)write(*,*)'totnos',itet,ntet,ib,nb,(e(i),i=1,4)
      CALL TOTNOS(VOL,E,EMIN,EMAX,NP,W(ONOS),W(OSOS))                  
220   CONTINUE                                                         
200   CONTINUE                                                         
      IF(DABS(SUM-1.D0).GT.1.D-5) GOTO 900
!     -----------------------------------------------------------------
!     --  GET FERMI LEVEL                                             -
!     -----------------------------------------------------------------
      tol=tolmax
!         write(*,*) 'in efi'
      CALL EFI(RNTOT,TOL,EF,EMIN,EMAX,NP,W(ONOS),W(OSOS))              
!         write(*,*) 'after efi',tol,tolmax,emax,emin,np,rntot
                                                                       
!     -----------------------------------------------------------------
!     --  CHECK ACCURACY AND RESTART IF NECCESARY                     -
!     -----------------------------------------------------------------
      IF(TOL.GT.TOLMAX) THEN                                           
        ESTEP=(EMAX-EMIN)/DBLE(NP-1)     
        IP=1+(EF-EMIN)/ESTEP                                           
!          write(*,*) ip,DRVAL(W(ONOS),IP),DRVAL(W(ONOS),IP+1)
        EMIN=EMIN+ESTEP*DBLE(IP-1)                                     
        EMAX=EMIN+ESTEP                                                

!         if(estep.lt.1.d-8) then
         if(estep.lt.1.d-10) then

         if (myrank.eq.master) write(*,*) 'WARNING: EF not accurate, new emin,emax,NE-min,', &
         'NE-max',emin,emax,DRVAL(W(ONOS),IP),DRVAL(W(ONOS),IP+1)
          ef=(emin+emax)/2.d0
         goto 2000
         endif
        IF(RNTOT-DRVAL(W(ONOS),IP).LE.TOLMAX) THEN                     
          EF=EMIN                                                      
          GOTO 2000                                                    
        ELSE IF(DRVAL(W(ONOS),IP+1)-RNTOT.LE.TOLMAX) THEN              
          EF=EMAX                                                      
          GOTO 2000                                                    
        END IF                                                         
        IF(ILOOP.GT.5) GOTO 910
        if (myrank.EQ.master .OR. fastFilesystem) write(6,*) 'ILOOP ',ILOOP                                          
        GOTO 1000                                                      
      END IF                                                           
2000  CONTINUE                                                         
      CALL RLSE(ONOS)                                                  
      RETURN                                                           
  900 CALL OUTERR('FERMI',' TETRAHEDRA DO NOT FILL VOLUME')
      CALL OUTERR('FERMI',' STOP IN EFERMI')
      WRITE (ERRMSG,9000) SUM
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
  910 CALL OUTERR('FERMI',' CANNOT FIND FERMI LEVEL')
      CALL OUTERR('FERMI',' STOP IN EFERMI')
      WRITE (ERRMSG,9010) TOL
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9020) EMIN,EMAX
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9030) DRVAL(W(ONOS),IP),DRVAL(W(ONOS),IP+1)
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 9000 FORMAT(' SUM ',f10.5, ' SHOULD BE 1.')
 9010 FORMAT(' TOL ',f10.5)
 9020 FORMAT(' EMIN ',f10.5, ' EMAX ',f10.5)
 9030 FORMAT(' NOS(EMIN) ',f10.5, ' NOS(EMAX) ',f10.5)
 END SUBROUTINE EFERMI


                                                             
!     .....................................................EFI.........
 SUBROUTINE EFI(RNTOT,TOL,EFERMI,EMIN,EMAX,NP,NOS,SOS)            
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
   DOUBLE PRECISION NOS(NP),NOSUP,NOSLOW,NOSIP                      
   CHARACTER*67 ERRMSG
   DIMENSION SOS(NP)                                                

      ADD=0.D0                                                         
      DO 100 I=1,NP                                                    
      ADD=ADD+SOS(I)                                                   
      NOS(I)=NOS(I)+ADD                                                
100   CONTINUE                                                         
      IF(NOS(1).GT.RNTOT+.5D0*TOL.OR.NOS(NP).LT.RNTOT-.5D0*TOL) GOTO 900
      IPUP=NP                                                          
      IPLOW=1                                                          
      nosup=nos(ipup)
      noslow=nos(iplow)
      DO 200 IFIND=1,NP                                                
      IP=IPLOW+0.5*(IPUP-IPLOW)                                        
      NOSIP=NOS(IP)                                                    
      IF(RNTOT-NOSIP.GT.0.d0) THEN                                        
        IPLOW=IP                                                       
        NOSLOW=NOSIP                                                   
      ELSE                                                             
        IPUP=IP                                                        
        NOSUP=NOSIP                                                    
      END IF                                                           
      IF(IPUP.EQ.IPLOW+1) GOTO 300                                     
200   CONTINUE                                                         
      GOTO 910
300   CONTINUE                                                         
      TOL=NOSUP-NOSLOW                                                 
      ESTEP=(EMAX-EMIN)/DBLE(NP-1)                                     
      ELOW=EMIN+DBLE(IPLOW-1)*ESTEP                                    
      DNOS=NOSUP-NOSLOW                                                
      IF(DNOS.NE.0.D0) THEN                                            
        EFERMI=ELOW+(RNTOT-NOSLOW)/(NOSUP-NOSLOW)*ESTEP                
      ELSE                                                             
        EFERMI=ELOW                                                    
      END IF                                                           
      IF(EFERMI-ELOW.LT.0) PRINT*,'ERROR IN EFI '                      
      RETURN                                                           
  900 CALL OUTERR('FERMI','EFERMI OUT OF ENERGY RANGE')
      CALL OUTERR('FERMI','STOP IN EFI')
      WRITE (ERRMSG,9000) EMIN
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9010) NOS(1)
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9020) EMAX
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9030) NOS(NP)
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9040) ADD
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9050) (SOS(I),I=100,1000,100)
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9060) (NOS(I),I=100,1000,100)
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
  910 CALL OUTERR('FERMI','EFERMI NOT FOUND')
      CALL OUTERR('FERMI','STOP IN EFI')
      STOP  'FERMI - Error'
 9000 FORMAT('ENERGY OF LOWER BOUND                 :',f10.5)
 9010 FORMAT('NUMBER OF STATES AT THE LOWER BOUND   :',f10.5)
 9020 FORMAT('ENERGY OF UPPER BOUND                 :',f10.5)
 9030 FORMAT('NUMBER OF STATES AT THE UPPER BOUND   :',f10.5)
 9040 FORMAT('ADD ',f10.5)
 9050 FORMAT('SOS ',10f5.3)
 9060 FORMAT('NOS ',10f5.5)
END SUBROUTINE EFI
                                                             
!     .....................................................TOTNOS......
SUBROUTINE TOTNOS(VOL,E,EMIN,EMAX,NP,NOS,SOS)                    
  !     **                                                              *
  !     **  CALCULATES THE INTEGRATED DOS                               *
  !     **  FOR ONE TETRAHEDRON ON THE ENERGYMESH                       *
  !     **  INPUT :                                                     *
  !     **    VOL         WEIGHT OF THIS TETRAHEDRON                    *
  !     **    E           ENERGIES AT THE EDGEPOINTS                    *
  !     **    EMIN        MINIMUM VALUE OF ENERGY MESH FOR NOS AND SOS  *
  !     **    EMAX        MAXIMUM VALUE OF ENERGY MESH FOR NOS AND SOS  *
  !     **    NP          NUMBER OF POINTS ON THE ENERGY MESH           *
  !     **  OUTPUT:                                                     *
  !     **    SOS         > NOS(E)+ SUM OVER E: SOS(E) =                *
  !     **    NOS         > NUMBER OF STATES BELOW E                    *
  !     **                                                              *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DOUBLE PRECISION NOS                                             
      DIMENSION E(4)                                                   
      DIMENSION NOS(NP),SOS(NP)                                        
!     -----------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            -
!     -----------------------------------------------------------------
      X=DMIN1(E(1),E(2),E(3),E(4))                                     
      IF(X.GE.EMAX) THEN                                               
        RETURN                                                         
      END IF                                                           
      X=DMAX1(E(1),E(2),E(3),E(4))                                     
      IF(X.LE.EMIN) THEN                                               
        SOS(1)=SOS(1)+VOL                                              
        RETURN                                                         
      END IF                                                           
!     -----------------------------------------------------------------
!     --  ORDER ENERGIES                                              -
!     -----------------------------------------------------------------
      DO 100 I=1,3                                                     
      DO 100 J=I+1,4                                                   
      SVAR1=DMIN1(E(I),E(J))                                           
      SVAR2=DMAX1(E(I),E(J))                                           
      E(I)=SVAR1                                                       
      E(J)=SVAR2                                                       
100   CONTINUE                                                         
!     -----------------------------------------------------------------
!     --  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 -
!     -----------------------------------------------------------------
      E21=E(2)-E(1)                                                    
      if(e21.lt.1.d-10) e21=1.d-10
      E31=E(3)-E(1)                                                    
      E41=E(4)-E(1)                                                    
      E32=E(3)-E(2)                                                    
      E42=E(4)-E(2)                                                    
      E43=E(4)-E(3)                                                    
      ESTEP=(EMAX-EMIN)/DBLE(NP-1)                                     
      IMIN=IDINT(2.D0+(E(1)-EMIN)/ESTEP)                               
      IMIN=MAX0(1,IMIN)                                                
      IMAX=IDINT(1.D0+(E(2)-EMIN)/ESTEP)                               
      IMAX=MIN0(NP,IMAX)                                               
      EN=EMIN+ESTEP*(IMIN-1)                                           
      IF(IMAX.GE.IMIN) THEN                                            
        A=VOL/(E21*E31*E41)                                            
        DO 200 I=IMIN,IMAX                                             
        NOS(I)=NOS(I)+A*(EN-E(1))**3                                   
        EN=EN+ESTEP                                                    
200     CONTINUE                                                       
      END IF                                                           
      IMIN=MAX0(1,IMAX+1)                                              
      IMAX=INT(1.D0+(E(3)-EMIN)/ESTEP)                                 
      IMAX=MIN0(NP,IMAX)                                               
      IF(IMAX.GE.IMIN) THEN                                            
        A=VOL*E21**2/(E31*E41)                                         
        B=3.D0*VOL*E21/(E31*E41)                                       
        C=3.D0*VOL/(E31*E41)                                           
        D=-VOL/(E32*E41*E31*E42)*(E31+E42)                             
        DO 300 I=IMIN,IMAX                                             
        DE=EN-E(2)                                                     
!       NOS(I)=NOS(I)+A+B*DE+C*DE**2+D*DE**3                           
        NOS(I)=NOS(I)+A+DE*(B+DE*(C+D*DE))                             
        EN=EN+ESTEP                                                    
300     CONTINUE                                                       
      END IF                                                           
      IMIN=MAX0(1,IMAX+1)                                              
      IMAX=INT(1.D0+(E(4)-EMIN)/ESTEP)                                 
      IMAX=MIN0(NP,IMAX)                                               
      IF(E43.GT.0.D0) THEN                                             
        A=VOL                                                          
        D=VOL/(E41*E42*E43)                                            
        DO 400 I=IMIN,IMAX                                             
        NOS(I)=NOS(I)+A+D*(EN-E(4))**3                                 
        EN=EN+ESTEP                                                    
400     CONTINUE                                                       
      END IF                                                           
      IMIN=MAX0(1,IMAX+1)                                              
      IF(IMIN.GT.NP) RETURN                                            
      SOS(IMIN)=SOS(IMIN)+VOL                                          
      RETURN                                                           
      END                                                              
!     .....................................................SAMFAC......
      SUBROUTINE SAMFAC(NB,NKP,EB,EF,WGHT,W,NWX)                       
!     **                                                              *
!     **  CALCULATES SAMPLING WEIGHTS                                 *
!     **  INPUT :                                                     *
!     **    NB          NUMBER OF BANDS                               *
!     **    NKP         NUMBER OF K-POINTS                            *
!     **    EF          FERMI LEVEL                                   *
!     **    W           INTEGER WORK ARRAY                            *
!     **    NWX         LENTH OF WORK ARRAY                           *
!     **  OUTPUT :                                                    *
!     **    WGHT        SAMPLING WEIGHTS                              *
!     **                                                              *
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                              
      IMPLICIT INTEGER (O)                                             
      DIMENSION EB(NB,NKP),WGHT(NB,NKP)                                
      DIMENSION E(4),WGHT0(4),IKP(4)                                   
      INTEGER W(NWX)                                                   
      CALL DEF0(NWX)                                                   
      CALL INITDR(WGHT,NB*NKP)                                         
      CALL TETR0(VOL0,NKP,NTET,MWRIT) 
      CALL DEFI(OWORK,5*MWRIT)                                         
      INIT=1                                                           
      DO 100 ITET=1,NTET                                               
      CALL TETR1(INIT,ITET,IWGHT,IKP,MWRIT,W(OWORK))                   
      VOL=VOL0*DBLE(IWGHT)                                             
      DO 100 IB=1,NB                                                   
      DO 110 I=1,4                                                     
      E(I)=EB(IB,IKP(I))                                               
110   CONTINUE                                                         
      CALL INITDR(WGHT0,4)                                             
      CALL WEIGHT(VOL,E,EF,WGHT0)                                      
      DO 120 I=1,4                                                     
      WGHT(IB,IKP(I))=WGHT(IB,IKP(I))+WGHT0(I)                         
120   CONTINUE                                                         
100   CONTINUE                                                         
      CALL RLSE(OWORK)                                                 
      RETURN                                                           
      END                                                              
!     .....................................................WHEIGT......
      SUBROUTINE WEIGHT(VOL,E,EF,WGHT)                                 
!     **                                                              *
!     **  CALCULATES THE WEIGHTS FOR TETRAHEDRON-SAMPLING             *
!     **  CORRESPONDING TO INTEGRATION OVER ONE TETRAHEDRON           *
!     **                                                              *
!     **  CORRECTION FOR THE NONLINEAR SHAPE INCLUDED IF ICOR=1       *
!     **                                                              *
!     **  AUTHOR : P.BLOECHL                                          *
!     **                                                              *
!     **    VOL.........VOLUME OF THIS TETRAHEDRON                    *
!     **    EF..........FERMI ENERGY                                  *
!     **    D...........KT (NOT USED)                                 *
!     **    E...........ENERGIES AT THE EDGEPOINTS                    *
!     **                                                              *
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                              
      DIMENSION E(4),WGHT(4)                                           
      DIMENSION FA(4),FB(4),INDEX(4)                                   
      common /correct/ cordeg,icor
!     -----------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            -
!     -----------------------------------------------------------------
      X=DMIN1(E(1),E(2),E(3),E(4))                                     
      IF(X.GE.EF) THEN                                                 
        CALL INITDR(WGHT,4)                                            
        RETURN                                                         
      END IF                                                           
      X=DMAX1(E(1),E(2),E(3),E(4))                                     
      IF(X.LE.EF) THEN                                                 
        VPRIME=.25D0*VOL                                               
        DO 10 I=1,4                                                    
        WGHT(I)=VPRIME                                                 
10      CONTINUE                                                       
        RETURN                                                         
      END IF                                                           
!     -----------------------------------------------------------------
!     --  ORDER ENERGIES                                              -
!     -----------------------------------------------------------------
!     -- INDEX HOLDS THE ORIGINAL POSITION OF THE ENERGIES AND WEIGHTS 
      DO 120 I=1,4                                                     
      INDEX(I)=I                                                       
120   CONTINUE                                                         
      DO 100 I=1,3                                                     
      IP=I                                                             
      DO 110 J=I+1,4                                                   
      IF(E(IP).GT.E(J)) IP=J                                           
110   CONTINUE                                                         
      IF(IP.GT.I) THEN                                                 
        X=E(IP)                                                        
        E(IP)=E(I)                                                     
        E(I)=X                                                         
        K=INDEX(IP)                                                    
        INDEX(IP)=INDEX(I)                                             
        INDEX(I)=K                                                     
      END IF                                                           
100   CONTINUE                                                         
!     -----------------------------------------------------------------
!     --  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 -
!     -----------------------------------------------------------------
      E21=E(2)-E(1)                                                    
      E31=E(3)-E(1)                                                    
      E41=E(4)-E(1)                                                    
      E32=E(3)-E(2)                                                    
      E42=E(4)-E(2)                                                    
      E43=E(4)-E(3)                                                    
      CALL INITDR(WGHT,4)                                              
      IF(EF.GT.E(1).AND.EF.LE.E(2)) THEN                               
        DE=EF-E(1)                                                     
        VPRIME=.25D0*VOL*DE**3/(E21*E31*E41)                           
        WGHT(1)=VPRIME*(4.D0-DE/E21-DE/E31-DE/E41)                     
        WGHT(2)=VPRIME*DE/E21                                          
        WGHT(3)=VPRIME*DE/E31                                          
        WGHT(4)=VPRIME*DE/E41                                          
!       ------  PARAMETERS FOR CORRECION                               
        DOS=3.D0*VPRIME*4.D0/(EF-E(1))                                 
      ELSE IF(EF.GT.E(2).AND.EF.LT.E(3)) THEN                          
        DE1=EF-E(1)                                                    
        DE2=EF-E(2)                                                    
        DE3=E(3)-EF                                                    
        DE4=E(4)-EF                                                    
!       ------  TETRAHEDRON X1,X2,X13',X14'                            
        VPRIME=VOL*DE1**2/(E41*E31)*.25D0                              
        WGHT(2)=VPRIME                                                 
        WGHT(3)=VPRIME*(DE1/E31)                                       
        WGHT(4)=VPRIME*(DE1/E41)                                       
        WGHT(1)=VPRIME*(3.D0-DE1/E41-DE1/E31)                          
!       ------  TETRAHEDRON X2,X13',X23',X14'                          
        VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                     
        WGHT(1)=WGHT(1)+VPRIME*(2.D0-DE1/E31-DE1/E41)                  
        WGHT(2)=WGHT(2)+VPRIME*(2.D0-DE2/E32)                          
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32+DE1/E31)                       
        WGHT(4)=WGHT(4)+VPRIME*(DE1/E41)                               
!       ------  TETRAHEDRON X2,X23',X24',X14'                          
        VPRIME=.25D0*VOL*DE2**2*DE4/(E42*E32*E41)                      
        WGHT(1)=WGHT(1)+VPRIME*(1.D0-DE1/E41)                          
        WGHT(2)=WGHT(2)+VPRIME*(3.D0-DE2/E32-DE2/E42)                  
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32)                               
        WGHT(4)=WGHT(4)+VPRIME*(DE2/E42+DE1/E41)                       
!       ------  DOS=A+B*(EF-E2)+C*(EF-E2)**2                           
        DA=3.D0*VOL*E21/(E31*E41)                                      
        DB=6.D0*VOL/(E31*E41)                                          
        DC=-3.D0*VOL/(E32*E41*E31*E42)*(E31+E42)                       
        DOS=DA+DB*DE2+DC*DE2**2                                        
      ELSE IF(EF.GE.E(3).AND.EF.LT.E(4)) THEN                          
        DE=E(4)-EF                                                     
        VPRIME=.25D0*VOL*DE**3/(E41*E42*E43)                           
        VOL14=.25D0*VOL                                                
        WGHT(1)=VOL14-VPRIME*DE/E41                                    
        WGHT(2)=VOL14-VPRIME*DE/E42                                    
        WGHT(3)=VOL14-VPRIME*DE/E43                                    
        WGHT(4)=VOL14-VPRIME*(4.D0-DE/E41-DE/E42-DE/E43)               
!       ------  PARAMETERS FOR CORRECION                               
        DOS=3.D0*VPRIME*4.D0/(E(4)-EF)                                 
      ELSE                                                             
        GOTO 900
      END IF                                                           
!     -----------------------------------------------------------------
!     --  ADD CORRECTION FOR QUADRATIC DEVIATION                      -
!     -----------------------------------------------------------------
      IF(ICOR.EQ.1) THEN                                               
        DO 500 M=1,4                                                   
        DO 510 N=1,4                                                   
        WGHT(M)=WGHT(M)+.25D0*(E(N)-E(M))*DOS*.1D0                     
510     CONTINUE                                                       
500     CONTINUE                                                       
      END IF                                                           
!     -----------------------------------------------------------------
!     --  REORDER WEIGHTS                                             -
!     -----------------------------------------------------------------
      DO 600 I=1,4                                                     
      FA(INDEX(I))=WGHT(I)                                             
      FB(INDEX(I))=E(I)                                                
600   CONTINUE                                                         
      DO 610 I=1,4                                                     
      WGHT(I)=FA(I)                                                    
      E(I)=FB(I)                                                       
610   CONTINUE                                                         
      RETURN                                                           
  900 CALL OUTERR('FERMI','ERROR IN TETINT')
      STOP  'FERMI - Error'
      END                                                              
!     .....................................................TETR0.......
      SUBROUTINE TETR0(VOL,NKP,NTET,MWRIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      REWIND 14                                                        
      READ(14,1234)NKP1,NTET,VOL,MWRIT,NREC
      if(nkp1.ne.nkp) goto 900
 1234 format(2i10,e20.12,2i10) 
      RETURN                                                           
 900  CALL OUTERR('FERMI', &
                  'number of k-points inconsistent when reading kgen')
      CALL OUTERR('FERMI','check IN1 and KGEN files!')
      STOP  'FERMI - Error'
      END                                                              
!     .....................................................TETR1.......
      SUBROUTINE TETR1(INIT,ITET,IWGHT,IKP,MWRIT,IWORK)                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION IKP(4),IWORK(5*MWRIT)                                  
      SAVE IPOS,ntet                                                   
      IF(INIT.EQ.1) THEN                                               
!       READ(15,REC=1)NKP,NTET,V,MWRIT,NREC                            
        REWIND 14                                                      
        READ(14,1234)NKP,NTET,V,MWRIT,NREC                             
 1234 format(2i10,e20.12,2i10) 
        INIT=0                                                         
        IPOS=0                                                         
      END IF                                                           
      IREC=(ITET-1)/MWRIT+1                                            
      IF(ITET.GT.NTET) GOTO 900
      IF(IREC.NE.IPOS) THEN                                            
!       READ(15,REC=IREC+1)IWORK                                       
        READ(14,1235)IWORK                                             
 1235 format(6i10) 
!       ---------------------------------------------------------------
        NLEFT=NTET-(IREC-1)*MWRIT                                      
        NLEFT=5*MIN0(MWRIT,NLEFT)                                      
!       ---------------------------------------------------------------
        IPOS=IREC                                                      
      END IF                                                           
      IP=5*(ITET-1-(IPOS-1)*MWRIT)                                     
      IWGHT=IWORK(IP+1)                                                
      DO 100 I=1,4                                                     
      IKP(I)=IWORK(IP+1+I)                                             
100   CONTINUE                                                         
      RETURN                                                           
  900 CALL OUTERR('FERMI','ASK FOR NONEXISTING TETRAHEDRON')
      CALL OUTERR('FERMI','STOP IN TETR1')
      STOP  'FERMI - Error'
      END                                                              
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     ----  BLOCK MIXPIC                                            ---
!     ----  MIXED PICKLES                                           ---
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     .................................................................
      SUBROUTINE DEF0(NMAX)                                            
!     **                                                              *
!     **  ALLOCATES SPACE ON A (INTEGER WORK ARRAY)                   *
!     **  USE DEF0 BEFORE USE TO GIVE AVAILABLE SPACE ON WORK ARRAY   *
!     **  USE DEFI TO ALLOCATE SPACE FOR INTEGER ARRAYS               *
!     **  USE DEFDR TO ALLOCATE SPACE FOR INTEGER ARRAYS              *
!     **  USE RLSE TO RELEASE SPACE                                   *
!     **                                                              *
!     **  INPUT:                                                      *
!     **    NMAX        LENGTH OF INTEGER WORK ARRAY                  *
!     **    LENG        NUMBER OF ELEMENTS IN THE ARRAY TO BE         *
!     **                MAPPED ONTO THE WORK ARRAY (ENTRY: DEFI,DEFDR)*
!     **    ONAME       ALL ARRAYS FROM POINTER ONAME ARE DROPPED     *
!     **                (ENTRY: RLSE)                                 *
!     **  OUTPUT :                                                    *
!     **    ONAME       POINTER OF ARRAY TO BE ALLOCATED              *
!     **                                          (ENTRY: DEFI,DEFDR) *
!     **  REMARKS :                                                   *
!     **    AN INTEGER NUMBER IS ASSUMED TO HAVE 4 BYTES              *
!     **    A DOUBLE PRECISION NUMBER IS ASSUMED TO HAVE 8 BYTES      *
!     **                                                              *
      IMPLICIT INTEGER (O)                                             
      SAVE OMAX,OMAXX                                                  
      OMAX=1                                                           
      OMAXX=NMAX                                                       
      RETURN                                                           
!     =================================================================
      ENTRY DEFI(ONAME,LENG)                                           
      ONAME=OMAX                                                       
      OMAX=OMAX+LENG                                                   
      IF(OMAX.GT.OMAXX) GOTO 9999                                      
      RETURN                                                           
!     =================================================================
      ENTRY DEFDR(ONAME,LENG)                                          
      ONAME=OMAX                                                       
      OMAX=OMAX+LENG*2                                                 
      IF(OMAX.GT.OMAXX) GOTO 9999                                      
      RETURN                                                           
!     =================================================================
      ENTRY RLSE(ONAME)                                                
      IF(ONAME.LE.0.OR.ONAME.GT.OMAX) GOTO 9999                        
      OMAX=ONAME                                                       
      RETURN                                                           
9999  CONTINUE                                                         
      CALL OUTERR('FERMI','ERROR IN DEF0')
      STOP  'FERMI - Error'
      END                                                              
!     .................................................................
      FUNCTION DRVAL(NAME,INDEX)                                       
      DOUBLE PRECISION NAME(INDEX),drval                               
      DRVAL=NAME(INDEX)                                                
      RETURN                                                           
      END                                                              
!     .................................................................
      SUBROUTINE INITI(NAME,LENG)                                      
      DIMENSION NAME(LENG)                                             
      DO 100 I=1,LENG                                                  
      NAME(I)=0                                                        
100   CONTINUE                                                         
      RETURN                                                           
      END                                                              
!     .................................................................
      SUBROUTINE INITDR(NAME,LENG)                                     
      DOUBLE PRECISION NAME(LENG)                                      
      DO 100 I=1,LENG                                                  
      NAME(I)=0.D0                                                     
100   CONTINUE                                                         
      RETURN                                                           
      END                                                              
!

double precision function erfc(xx)
  !
  !
  !     sandia mathematical program library
  !     applied mathematics division 2613
  !     sandia laboratories
  !     albuquerque, new mexico  87185
  !     control data 6600/7600  version 7.2  may 1978
  !
  !
  !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
  !  * the primary document for the library of which this routine is     
  !  * part is sand77-1441.                                              
  !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
  !
  !     written by j.e. vogel from approximations derived by w.j. cody .
  !
  !     abstract
  !
  !         erfc(x) computes 2.0/sqrt(pi) times the integral from x to
  !         infinity of exp(-x**2). this is done using rational approx-
  !         imations.  eleven correct significant figures are provided.
  !
  !     description of parameters
  !
  !         x may be any real value
  !
  !     erfc is documented completely in sc-m-70-275.
  !
  !     .. Scalar Arguments ..
  double precision xx
  !     ..
  !     .. Local Scalars ..
  double precision a, r, sqpi, x, x2, xi2
  integer i
  !     ..
  !     .. Local Arrays ..
  double precision p1(4), p2(6), p3(4), q1(4), q2(6), q3(4)
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic abs, exp
  !     ..
  !     .. Data statements ..
  data (p1(i),i=1,4)/242.6679552305318D0, 21.97926161829415D0, &
       6.996383488619136D0, -3.560984370181539D-02/, &
       (q1(i),i=1,4)/215.0588758698612D0, 91.16490540451490D0, &
       15.08279763040779D0, 1.0D0/, (p2(i),i=1,6)/22.898992851659D0, &
       26.094746956075D0, 14.571898596926D0, 4.2677201070898D0, &
       .56437160686381D0, -6.0858151959688D-06/, &
       (q2(i),i=1,6)/22.898985749891D0, 51.933570687552D0, &
       50.273202863803D0, 26.288795758761D0, 7.5688482293618D0, &
       1.0D0/, (p3(i),i=1,4)/-1.21308276389978D-2, &
       -.1199039552681460D0, -.243911029488626D0, &
       -3.24319519277746D-2/, (q3(i),i=1,4)/4.30026643452770D-02, &
       .489552441961437D0, 1.43771227937118D0, 1.0D0/
  data sqpi/.564189583547756D0/
  !
  x = abs(xx)
  x2 = x*x
  if (xx .lt. -6.0D0) then
     erfc = 2.0D0
  else if (x .gt. 25.8D0) then
     erfc = 0.0D0
  else if (x .gt. 4.0D0) then
     xi2 = 1.D0/x2
     r = xi2*(p3(1)+xi2*(p3(2)+xi2*(p3(3)+xi2*p3(4))))/(q3(1)+xi2*(q3(2)+xi2*(q3(3)+xi2*q3(4))))
     a = exp(-x2)*(sqpi+r)/x
     if (xx .lt. 0.0D0) then
        erfc = 2.0D0 - a
     else
        erfc = a
     end if
  else if (x .gt. .46875D0) then
     a = exp(-x2)*(p2(1)+x*(p2(2)+x*(p2(3)+x*(p2(4)+x*(p2(5)+x*p2(6))))))
     a = a/(q2(1)+x*(q2(2)+x*(q2(3)+x*(q2(4)+x*(q2(5)+x*q2(6))))))
     if (xx .le. 0.0D0) then
        erfc = 2.0D0 - a
     else
        erfc = a
     end if
  else
     a = x*(p1(1)+x2*(p1(2)+x2*(p1(3)+x2*p1(4))))
     a = a/(q1(1)+x2*(q1(2)+x2*(q1(3)+x2*q1(4))))
     if (xx .lt. 0.D0) a = -a
     erfc = 1.0D0 - a
  end if
  return
end function erfc
