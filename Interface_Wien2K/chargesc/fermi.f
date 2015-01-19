SUBROUTINE FERMI(nbmax,sumw0)
  !                                                                      
  !     NELEC IS THE NUMBER OF ELECTRONS IN THIS COMPOUND                
  !     EF IS THE FERMI-ENERGY TO BE CALCULATED                          
  !                                                                      
  USE param
  USE com_mpi
  USE bandm
  USE char
  USE kpp1
  USE xa2
  USE com,only : weigh,ef,elecn,xwt,nat,nband,rel,nk,nb,minwav,maxwav,init_com,GDELTA
  use struk
  USE GMPI
  IMPLICIT REAL*8 (A-H,O-Z)
  !
  REAL*8, intent(out) :: sumw0
  INTEGER, intent(out):: nbmax
  !
  CHARACTER*10    KNAME
  CHARACTER*67    ERRMSG                                          
  CHARACTER*180   FNAME
  REAL*8,ALLOCATABLE  :: e_local(:,:)
  common /mkef/delef,ts2
  !--------------------------------------------------------------------- 
  !find nkpt, nmat and nume in energy file
  k=0  
  sumw0=0

  DO ivector=1,nvector
     itape=30
     if (vector_para) then
        open(itape,FILE=fvectors(ivector, 3),STATUS='old',FORM='formatted')
     endif
     DO I=1,NAT
        READ(itape,'(f9.5)') EMIST
        READ(itape,'(f9.5)') EMIST
     ENDDO
     ios=0
     DO WHILE (ios == 0)
        READ(itape,'(3e19.12,a10,2i6,F5.1)',IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
        IF (ios /= 0) CYCLE
        k=k+1
        nmat=MAX(n,nmat)
        nume=MAX(nen,nume)
        DO ii=1,NEn
           READ(itape,*) NUM,E1
        ENDDO
        sumw0=sumw0+wgh
     ENDDO
     if (vector_para) then
        close(itape)
     else
        REWIND(itape)
     endif
  END DO
  nkpt=k 
  IF(vector_para)THEN
    CALL IMAX1_ALL_MPI(nmat)
    CALL IMAX1_ALL_MPI(nume)
    CALL ISUM1_ALL_MPI(nkpt)
    CALL DSUM1_ALL_MPI(sumw0)
  ENDIF
  nkpt=nkpt+1
  if (myrank.EQ.master .OR. FastFilesystem) WRITE(6,*) 'nume=', nume, 'nkpt=', nkpt 

  allocate(e_local(2*nkpt,nume))
  CALL init_bandm(nume)
  CALL init_com(nkpt,nume)
  CALL init_kpp(nkpt)
  CALL init_xa2(nume,nkpt)
  If(efmod.eq.'TETRA') then	
     !     tetraeder method of bloechl (for ef.lt.100  corr is used!!)
     call fermi_tetra(nbmax)
     ts2=0.d0
     return	
  ELSE ! Not used actually, just for estimation
     call fermi_GAUSS(nbmax)
     return
  end if
  STOP ' ERROR: efmod not implemented.'
  !
  !        Error messages
  !
900 WRITE (ERRMSG,9000) NKPT
  CALL OUTERR('FERMI',ERRMSG)
  STOP  'FERMI - Error'
920 WRITE (ERRMSG,9020) nume,ne(k) 
  CALL OUTERR('FERMI',ERRMSG)
  STOP  'FERMI - Error'
930 CALL OUTERR('FERMI',' SUMUP AND SUMDN IN FERMI NOT EQUAL')
  CALL OUTERR('FERMI',' CHECK UP AND DOWN INPUTS OF LAPW1 ')
  WRITE (ERRMSG,9030) SUMUP,SUMW,DIF             
  CALL OUTERR('FERMI',ERRMSG)
  STOP  'FERMI - Error'
940 WRITE (ERRMSG,9041) ELECN
  CALL OUTERR('FERMI',ERRMSG)
  WRITE (ERRMSG,9040) ELECN,ELN,EMIN,INDEX
  CALL OUTERR('FERMI',ERRMSG)
  CALL OUTERR('FERMI','    INCREASE ENERGY WINDOW IN CASE.IN1')
  CALL OUTERR('FERMI',' OR INCREASE PARAMETER NUME ')
  CALL OUTERR('FERMI',' OR DECREASE NE IN CASE.IN2 ')
  STOP  'FERMI - Error'
950 WRITE (ERRMSG,9060) IUNIT
  CALL OUTERR('LAPW2',ERRMSG)
  WRITE (ERRMSG,9070) FNAME
  CALL OUTERR('LAPW2',ERRMSG)
  STOP 'FERMI - Error'
  !
9000 FORMAT ('NUMBER OF K-POINTS .GT. NKPT =',i6)
9020 FORMAT ('nume,ne(k)',2i6)
9030 FORMAT ('SUMUP,SUMW,DIF', 3f10.7)
9040 FORMAT ('ELECN,ELN,EMIN,INDEX', 3f9.5,i5)
9041 FORMAT (' NOT ENOUGH EIGENVALUES FOR',f4.0,' ELECTRONS')
9050 FORMAT (8F10.7) 
9060 FORMAT('can''t open unit: ',I2)
9070 FORMAT('       filename: ',A50)
9080 FORMAT('         status: ',A,'  form: ',A)
5001 FORMAT(3e19.12,a10,2i6,F5.1)
199 FORMAT(3F10.5,3X,A10,2I5,F5.2)                                   
170 FORMAT(F5.2)                                                     
201 FORMAT(8(3I3))                                                   
202 FORMAT(I5,F10.5)                                                 
203 FORMAT(8F9.6)                                                    
END SUBROUTINE FERMI

SUBROUTINE fermi_GAUSS(nbmax)
  !                                                                      
  !     NELEC IS THE NUMBER OF ELECTRONS IN THIS COMPOUND                
  !     EF IS THE FERMI-ENERGY TO BE CALCULATED                          
  !                                                                      
  USE param
  USE com_mpi
  USE bandm
  USE kpp1
  USE xa2,only : e,ne
  USE com,only : weigh,ef,elecn,xwt,nat,nband,rel,nk,nb,minwav,maxwav,JSPIN,GISO,GDELTA
  use struk
  USE GMPI
  IMPLICIT REAL*8 (A-H,O-Z)                                        
  PARAMETER (nw=250000)                                            
  CHARACTER *10    KNAME                                
  CHARACTER *67    ERRMSG                                          
  CHARACTER *180   FNAME
  LOGICAL          EFERM
  INTEGER          :: iw(nw),RSPO
  REAL*8           :: eb(nume,nkpt,2),WE(nume,nkpt),WTK(NKPT)
  INTEGER          NEHELP(NKPT,2)
  
  MINWAV=18000                                                     
  MAXWAV=0                                                         
  nemax=0
  TEST=8.9D0
  NEHELP=0                                      
  !     for babi test1 changed to lower value, otherwise no degeneracy fo
  TEST1=2.D-8                                                      
  !                                                                      
  !.....MORE DIMENS.REPRESENTATION,IF DELTA E LESS THEN TEST1            
  ELECN=ELECN-1.D-10                                               
  !para begin
  ! ensure to mimick the proper vector file
  ! if running on multiple processors
  k1=0
  ispin=1
  WTK=0
  !.....READ FROM TAPE 10,WRITTEN BY LAPW1                               
  DO ivector=1,nvector
     itape=30
     if (vector_para) then
        open(itape,FILE=fvectors(ivector, 3),STATUS='old',FORM='formatted')
     endif
     DO I=1,NAT
        READ(itape,'(f9.5)') EMIST
        READ(itape,'(f9.5)') EMIST
     ENDDO
     ios=0; K=0
     DO WHILE (ios == 0)
        READ(itape,'(3e19.12,a10,2i6,F5.1)',IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
        IF (ios /= 0) CYCLE
        k=k+1
        if (vector_para) then
           ikp = vectors(ivector,3)+K  ! successive index in k-point table from case.klist
        else
           ikp = K                     ! successive index in k-point table from case.klist
        endif
        NEHELP(IKP,ispin)=NEn
        WTK(IKP)=wgh
        if(nehelp(ikp,ispin).gt.nume) then
          WRITE (ERRMSG,9020) nume,nehelp(k,ispin)
          CALL OUTERR('FERMI',ERRMSG)
          STOP  'FERMI - Error'
        endif
        IF(N.GT.MAXWAV) MAXWAV=N
        IF(N.LT.MINWAV) MINWAV=N
        DO ii=1,NEn
           READ(itape,*) NUM,E1
           Eb(num,IKP,ispin)=E1
        ENDDO
     ENDDO
     if (vector_para) then
        close(itape)
     else
        REWIND(itape)
     endif
  ENDDO

  IF(vector_para)THEN
    CALL ISUM_ALL_MPI(NEHELP,NKPT*2)
    CALL IMAX1_ALL_MPI(MAXWAV); CALL IMIN1_ALL_MPI(MINWAV)
    CALL DSUM_ALL_MPI(Eb,nume*nkpt*2)
    CALL DSUM_ALL_MPI(WTK,nkpt)
  ENDIF
  SUMW=SUM(WTK)
  WTK=WTK/SUMW

  NE=RESHAPE(NEHELP,(/NKPT*2/))
  NEMAX=MAXVAL(nehelp)
  DO IKP=1,NKPT; DO II=1,NE(IKP)
    E1=Eb(II,IKP,1)
    IF(E1.GT.ebmax(II))ebmax(II)=E1
    IF(E1.LT.ebmin(II))ebmin(II)=E1
  ENDDO; ENDDO

  do i=1,nkpt
     do j=NE(I)+1,nume
        eb(j,i,1)=4.d0
     enddo
  enddo
  eb(:,:,2)=4.d0

  RSPO=3-MAX(GISO,JSPIN); WE=0
  CALL FERMI_GS(ELECN,GDELTA,RSPO,JSPIN,Eb(:,:,1),WE,NUME,NKPT,WTK,NKPT,3.9D0,1.D-4,EF,TS2,NBMAX)
  DO i=1,nkpt; DO J=1,NE(I)
    WEIGH(I,J)=WE(J,I)
  ENDDO; ENDDO 
  RETURN                                             
9020 FORMAT ('nume,ne',2i6)


END SUBROUTINE fermi_GAUSS

SUBROUTINE fermi_tetra(nbmax)
  !                                                                      
  !     NELEC IS THE NUMBER OF ELECTRONS IN THIS COMPOUND                
  !     EF IS THE FERMI-ENERGY TO BE CALCULATED                          
  !                                                                      
  USE param
  USE com_mpi
  USE bandm
  USE kpp1
  USE xa2,only : weight,e,ne
  USE com,only : weigh,ef,elecn,xwt,jspin=>nspin,GISO,nat,nband,rel,nk,nb,minwav,maxwav
  use struk
  USE GMPI
  IMPLICIT REAL*8 (A-H,O-Z)                                        
  PARAMETER (nw=250000)                                            
  CHARACTER *10    KNAME                                
  CHARACTER *67    ERRMSG                                          
  CHARACTER *180   FNAME
  LOGICAL          EFERM
  INTEGER          :: iw(nw)
  REAL*8           :: eb(nume,nkpt,2)
  INTEGER          NEHELP(NKPT,2)
  common /correct/ cordeg,icor
  !--------------------------------------------------------------------- 
  !  
  !Clas0
  !     icor switches non-linear correction  on/off
  icor=1
  if(abs(ef).ge.100.d0) icor=0
  !     cordeg is a switch to correct occupancy of degenerate states
  cordeg=-1.d-6
  if(ef.gt.0.d0)  cordeg=-ef   
  if(ef-100.d0.gt.0.d0) cordeg=-ef+100.d0
  if(ef.lt.0.d0)        cordeg=ef
  if(ef.lt.-100.d0)     cordeg=ef+100.d0
  if(abs(cordeg).gt.0.01d0)  cordeg=-1.d-6
  !
  if (myrank.EQ.master .OR. fastFilesystem) write(6,'(" BZ-integration with TETRA-program.   icor=:",I2  )') icor
  if(cordeg.lt.0.d0 .and. (myrank.EQ.master .OR. fastFilesystem) ) write(6,'(" Equal occupancy of degenerate states, tol=:",E8.1)') -cordeg
  
  do j=1,nume*2*nkpt
     e(j)=3.0d0
  enddo
  eb=0

  MINWAV=18000                                                     
  MAXWAV=0                                                         
  nemax=0
  JSPIN=GISO
  TEST=8.9D0
  NEHELP=0                                      
  !     for babi test1 changed to lower value, otherwise no degeneracy fo
  TEST1=2.D-8                                                      
  !                                                                      
  !.....MORE DIMENS.REPRESENTATION,IF DELTA E LESS THEN TEST1            
  ELECN=ELECN-1.D-10                                               
  !para begin
  ! ensure to mimick the proper vector file
  ! if running on multiple processors
  k1=0
  ispin=1
  !.....READ FROM TAPE 10,WRITTEN BY LAPW1                               
  DO ivector=1,nvector
     itape=30
     if (vector_para) then
        open(itape,FILE=fvectors(ivector, 3),STATUS='old',FORM='formatted')
     endif
     DO I=1,NAT
        READ(itape,'(f9.5)') EMIST
        READ(itape,'(f9.5)') EMIST
     ENDDO
     ios=0; K=0
     DO WHILE (ios == 0)
        READ(itape,'(3e19.12,a10,2i6,F5.1)',IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
        IF (ios /= 0) CYCLE
        k=k+1
        if (vector_para) then
           ikp = vectors(ivector,3)+K  ! successive index in k-point table from case.klist
        else
           ikp = K                     ! successive index in k-point table from case.klist
        endif
        NEHELP(IKP,ispin)=NEn
        if(nehelp(ikp,ispin).gt.nume) then
          WRITE (ERRMSG,9020) nume,nehelp(k,ispin)
          CALL OUTERR('FERMI',ERRMSG)
          STOP  'FERMI - Error'
        endif
        IF(N.GT.MAXWAV) MAXWAV=N
        IF(N.LT.MINWAV) MINWAV=N
        DO ii=1,NEn
           READ(itape,*) NUM,E1
           Eb(num,IKP,ispin)=E1
        ENDDO
     ENDDO
     if (vector_para) then
        close(itape)
     else
        REWIND(itape)
     endif
  ENDDO

  IF(vector_para)THEN
    CALL ISUM_ALL_MPI(NEHELP,NKPT*2)
    CALL IMAX1_ALL_MPI(MAXWAV); CALL IMIN1_ALL_MPI(MINWAV)
    CALL DSUM_ALL_MPI(Eb,nume*nkpt*2)
  ENDIF

  NE=RESHAPE(NEHELP,(/NKPT*2/))
  NEMAX=MAXVAL(nehelp)
  DO IKP=1,NKPT; DO II=1,NE(IKP)
    E1=Eb(II,IKP,1)
    IF(E1.GT.ebmax(II))ebmax(II)=E1
    IF(E1.LT.ebmin(II))ebmin(II)=E1
  ENDDO; ENDDO

  do i=1,nkpt
     do j=NE(I)+1,nume
        eb(j,i,1)=3.0d0
     enddo
  enddo
  eb(:,:,2)=3.d0

  K=NKPT-1
1234 format(20(10f8.4,/))
  if (myrank.EQ.master .OR. fastFilesystem) write(6,*)'call eord...' 
  call eord(e,eb,nemax,k,jspin) 
  if (myrank.EQ.master .OR. fastFilesystem) write(6,*)'call dos...' 
  call dos(nemax*jspin,k,e,weight,elecn/2.d0*jspin,d,ef,iw,nw) 
  if (myrank.EQ.master .OR. fastFilesystem) write(6,*)'call eweigh...'
  call eweigh(ef,weigh,weight,nemax,k,jspin,nehelp,eb,nbmax)
  !     REWIND TAPES AND RETURN                             
  RETURN                                             
  !
  !        Error messages
  !
900 WRITE (ERRMSG,9000) NKPT
  CALL OUTERR('FERMI',ERRMSG)
  STOP  'FERMI - Error'
920 WRITE (ERRMSG,9020) nume,nehelp(k,ispin) 
  CALL OUTERR('FERMI',ERRMSG)
  STOP  'FERMI - Error'
930 CALL OUTERR('FERMI',' # of k-points in up and down not equal:')
  WRITE (ERRMSG,9030) k1,k
  CALL OUTERR('FERMI',ERRMSG)
  STOP  'FERMI - Error'
950 WRITE (ERRMSG,9060) 
  CALL OUTERR('LAPW2',ERRMSG)
  WRITE (ERRMSG,9070) FNAME
  CALL OUTERR('LAPW2',ERRMSG)
  STOP 'FERMI - Error'
!
!#ifdef Extended
!5001 format(3D27.20,a10,2i6,F5.1)
!#else
5001 FORMAT(3e19.12,a10,2i6,F5.1)
!#endif
  !
9000 FORMAT ('NUMBER OF K-POINTS .GT. NKPT =',i6)
9020 FORMAT ('nume,ne',2i6)
9030 FORMAT ('k1, k', 2i6,' check INPUTS OF LAPW1 ')
9060 FORMAT('can''t open unit: ')
9070 FORMAT('       filename: ',A50)
!                                                                      
END SUBROUTINE fermi_tetra

subroutine eord(e,eb,nemax,k,jspin) 
  USE param
  IMPLICIT REAL*8 (A-H,O-Z)                                        
  dimension eb(nume,nkpt,2),e(nemax*jspin,k)
  do k1=1,k
     ne2=0
     do ispin=1,jspin
        do ne1=1,nemax
           ne2=ne2+1
           e(ne2,k1)=eb(ne1,k1,ispin)
        enddo
     enddo
  enddo
  return
end subroutine eord

subroutine eweigh(ef,weigh,weight,nemax,k,jspin,ne,eb,nbmax)       
  USE param
  !USE parallel
  USE com_mpi
  USE bandm
  IMPLICIT REAL*8 (A-H,O-Z)                                        
  dimension eb(nume,nkpt,2)
  dimension weigh(2*nkpt,nume),weight(nemax*jspin,k),ne(k)
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
   8  WEIGH(KK,NN)=WEIGHT(nn,KK) *2.0d0 / jspin                        
        
     if (myrank.EQ.master .OR. fastFilesystem) write(6,*) '  number of occupied bands:',nbmax
     if (myrank.EQ.master .OR. fastFilesystem) write(6,*) '  highest energy:',emax
     if(ef.gt.emax.and.abs(emax-ebmin(nbmax+1)).gt.1.d-3) then
        if (myrank.EQ.master .OR. fastFilesystem) write(6,*) 'insulator !'
        if((ef-emax).lt.1.d-4.or.(ef-emax).ge.1.d-4) then
           if (myrank.EQ.master .OR. fastFilesystem) write(6,*) 'EF-inconsistency corrected'
           if (myrank.eq.master) write(21,888) 
888        format('       Insulator, EF-inconsistency corrected')
           ef=emax
        endif
     endif
     ef=ef+0.5d0            
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
 
  IF(DABS(SUM-RNTOT).GT.TOLMAX) THEN                               
     if (myrank.EQ.master .OR. fastFilesystem) WRITE (6,9000) SUM,RNTOT
     if (myrank.eq.master) write(21,9001) SUM,RNTOT
  ELSE                                                             
     if (myrank.EQ.master .OR. fastFilesystem) write(6,*) '  SUM RULE OF TETRAHEDRON INTEGRATION CHECKED '
  END IF
  RETURN                                                           
900 CALL OUTERR('FERMI',' INTEGRATION FAILED.....STOP IN DOS')
  WRITE (ERRMSG,9000) SUM,RNTOT
  CALL OUTERR('FERMI',ERRMSG)
  STOP  'FERMI - Error'
9000 FORMAT(' RESULT OF INTEGRATION: ',f10.5, '; SHOULD BE: ',f10.5)
9001 FORMAT(':WARN : RESULT OF INTEGRATION: ',f10.5, '; SHOULD BE: ',f10.5)
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
         if(estep.lt.1.d-8) then
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
!      DATA ICOR/1/                                                    
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
      !COMMON /IPROC/ IPROC                       
!     READ(15,REC=1)NKP,NTET,VOL,MWRIT,NREC                            
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

!***********************************************************************
      SUBROUTINE FERMI_GS(ZTOT,DELTA,RSPO,JSPIN,E,WE,NBMAX,NKST,WTK,NKPT,&
                         &ECUT,ESTEP,EF,TS2,NBTOP)
      IMPLICIT NONE
      INTEGER,PARAMETER :: gq=8 ! Precision
      REAL(gq)   ,PARAMETER :: PI =3.141592653589793238_gq
      INTEGER,INTENT(IN)::NBMAX,NKST,NKPT,RSPO,JSPIN
      REAL(gq),INTENT(IN)::ZTOT,DELTA,ECUT,ESTEP,E(NBMAX,NKST),WTK(NKPT)
      REAL(gq),INTENT(OUT)::WE(NBMAX,NKST),EF,TS2
      INTEGER,INTENT(OUT)::NBTOP
! LOCAL
      REAL(gq) EFERMI,EMIN,EMAX,NEL,NEL0,DNNT,EF0,EF1,DE,WT,FAC,E1,ETA
      INTEGER IKS,IB,IKP,NNT,NC,NB_TOP
      REAL(gq),PARAMETER::TOL=1.E-6_gq
      INTEGER,PARAMETER::NCMAX=10000

      NB_TOP=0
      NNT=ZTOT/RSPO; DNNT=ZTOT/RSPO-NNT
      IF(NNT==0)THEN
        EF0=E(1,1)
      ELSE
        EF0=E(NNT,1)+DNNT*(E(NNT+1,1)-E(NNT,1))
      endif
      EFERMI=EF0

      DO NC=1,NCMAX
      NEL=0; WE=0
      EMAX=EF0-1000._gq; EMIN=EF0+1000_gq
      DO IKS=1,NKST
      IF(IKS<=NKPT)THEN; IKP=IKS; ELSE; IKP=IKS-NKPT; ENDIF
      DO IB=1,NBMAX
        IF(E(IB,IKS)>ECUT)CYCLE
        E1=E(IB,IKS)
        EMIN=MIN(EMIN,E1)
        DE=(E1-EFERMI)/DELTA
        IF(DE<-3._gq)THEN
          WT=2._gq
        ELSEIF(DE<0._gq)THEN
          WT=2-ERFC(-DE)
        ELSEIF(DE<3._gq)THEN
          WT=ERFC(DE)
        ELSE
          WT=0._gq
        ENDIF
        WT=WT/2*RSPO
        NEL=NEL+WT*WTK(IKP)
        WE(IB,IKS)=WT*WTK(IKP)
        IF(WT>1.E-5_gq)THEN
          EMAX=MAX(E1,EMAX); NB_TOP=MAX(NB_TOP,IB)
        ENDIF
      ENDDO; ENDDO

! write(*,*)ABS(ZTOT-NEL),ZTOT,NEL
! special treatment of the first iteration step
      IF(NC==1)THEN
        EFERMI=EFERMI+ESTEP
        NEL0=NEL
        CYCLE
      ENDIF
! all other steps
      IF(ABS(ZTOT-NEL)<TOL)GOTO 100 ! LOCATED EF
      EF1=EFERMI
      FAC=(NEL0-NEL)/(ZTOT-NEL)
      IF(ABS(FAC)>.1_gq)THEN
        EFERMI=EFERMI+(EF0-EFERMI)/FAC
        EF0=EF1; NEL0=NEL
      ELSE
        EF0=EF1; NEL0=NEL
        IF((ZTOT-NEL)<0)THEN
          EFERMI=EFERMI-ESTEP
        ELSE
          EFERMI=EFERMI+ESTEP
        ENDIF
      ENDIF
      ENDDO ! NC

      WRITE(0,*) 'ERROR IN FERMI_GS: FERMILEVEL NOT CONVERGED!'; STOP
100   CONTINUE
      EF=EFERMI
      EMAX=EMAX+0.0001_gq
      EMIN=EMIN-0.0001_gq
      NBTOP=NB_TOP
! Caluculate energycorrection caused by gaussian smearing
! Start loop over k-points (for spinpol.systems divide by 2
      ETA=0
      DO IKS=1,NKST; 
      IF(IKS<=NKPT)THEN; IKP=IKS; ELSE; IKP=IKS-NKPT; ENDIF
      DO IB=1,NBMAX
        IF(E(IB,IKS)>ECUT)CYCLE
        DE=(E(IB,IKS)-EFERMI)/DELTA
        DE=DE*DE
        IF(DE<15._gq)ETA=ETA+0.5*DELTA*EXP(-DE)*WTK(IKP)
      ENDDO; ENDDO
      ETA=-ETA*2._gq/SQRT(PI)/JSPIN
      TS2=ETA/2._gq
      RETURN

      END SUBROUTINE FERMI_GS
