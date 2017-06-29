      PROGRAM LAPWDM                                                      
! LAPW2DM calculates the density matrix
! for new calculation of <X> and for rotationally invariant LDA+U method
! LAPW2DM is modified LAPW2 packages
! last change P.Novak 08.010.2001 novak fzu.cz
        USE param
        USE struct
        USE rotat
        USE sporb
        USE case
        USE matpdf
        USE symop
        USE xxa
        USE xdos
        USE com
        USE xa
        USE xa3
        USE ams
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*5      COORD 
      CHARACTER*4     MAG                                      
      CHARACTER*10    KNAME
      CHARACTER*11     STATUS,FORM                                      
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN
      CHARACTER*180     FNAME,VECFN  

      LOGICAL         SO,MAGN
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
      COMMON/ANGL/     THETA,PHI
      common /aver/ krad,kls,cx(-20:20,20),iprx
      character*60 xrad(-20:20),comment(-20:20,20)
      character*70 xls(12)
!                                                                       
      DIMENSION  S(3),l(1:lmax2+1)
      DATA             SO/.false./,MAGN/.false./                           
      
!-----------------------------------------------------------------------  
!                                                                       
      CALL init_ams
      CALL GTFNAM(DEFFN,ERRFN,IPROC)
      CALL ERRFLG(ERRFN,'Error in LAPW2DM')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
   10 CONTINUE
         READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         if(iunit.eq.9) then
         VECFN=FNAME
           do i=180,5,-1
          if((fname(i-4:i).eq.'rsoup').or.(fname(i-4:i).eq.'rsodn'))then
            SO=.true.
            MAG='MAGN'
            goto 10
          else if(fname(i-4:i).eq.'ector')then
            SO=.false.
            MAG='NMAG'
            goto 10
          else if(fname(i-4:i).eq.'torso')then
            SO=.true. 
            MAG='NMAG'
            goto 10
          else if((fname(i-4:i).eq.'torup').or.(fname(i-4:i).eq.'tordn'))then
            SO=.false.
            MAG='MAGN'
            goto 10
          endif
           enddo
         endif
      GOTO 10
   20 CONTINUE
      CLOSE (1)

      iprx=0
      write(6,*)
      if (mag.eq.'MAGN') then 
       if(so)then
       magn=.true.
        write(6,*)'MAGNETIC SYSTEM WITH SPIN-ORBIT COUPLING'
       else
        write(6,*)'MAGNETIC SYSTEM WITHOUT SPIN-ORBIT COUPLING'
       endif
      else
       if(so)then
        write(6,*)'NONMAGNETIC SYSTEM WITH SPIN-ORBIT COUPLING'
       else
        write(6,*)'NONMAGNETIC SYSTEM WITHOUT SPIN-ORBIT COUPLING'
       endif
      end if
      CALL CPUTIM(Tstart)

!.....READ STRUCT                                                       
      CALL init_struct

      CALL init_rotat(nsym,nato)
      CALL init_case(nato,lmax2)
      CALL init_sporb(nsym)
      CALL init_matpdf(nsym)
      CALL init_symop(nsym)
      CALL init_xdos(lxdos,ndif,lmax2)
!....Find nume, nmat and nkpt
     k=0  
      DO I=1,NAT 
         READ(50,'(f9.5)') EMIST
         READ(50,'(f9.5)') EMIST
      ENDDO
      DO I=1,NAT                                                  
         READ(51,'(f9.5)',iostat=ios) EMIST
         READ(51,'(f9.5)',iostat=ios) EMIST
      ENDDO
      IF(ios.eq.0) JSPIN=2
      DO
         READ(50,'(3e19.12,a10,2i6)',IOSTAT=ios) Sxx,Txx,Zxx,KNAME,N,NEn
         IF (ios /= 0) EXIT
         k=k+1
         nmat=MAX(n,nmat)
         nume=MAX(nen,nume)
         DO ii=1,nen
            READ(50,*) NUM,E1
         ENDDO
         IF(jspin.EQ.2) THEN
            READ(51,'(3e19.12,a10,2i6)',IOSTAT=ios) Sxx,Txx,Zxx,KNAME,N,NEn
            nmat=MAX(n,nmat)
            nume=MAX(nen,nume)
            DO ii=1,nen
               READ(51,*) NUM,E1
            ENDDO
         ENDIF
      ENDDO
      nkpt=k+1
      REWIND(50)
      REWIND(51)

      CALL init_com(nume,nkpt)
      CALL init_xa(LMAX2,NMAT,NRAD,NUME)
      CALL init_xa3(nume,nmat)
      CALL init_xxa(lmax2,nume,ndif,nrf)
!.....READ INPUT AND POTE  
      do i=1,nato
	do j=0,lmax2
	  lcase(i,j)=.false.
	end do
      end do                                             
      write(6,*)
      write(6,*)'ATOMS AND ORB. NUMBERS TO BE CONSIDERED:'
! if more information needed, put IPRINT>0
      IPRINT=0
      READ(5,*)EMIN
      READ(5,*)NATOM
	DO I=1,NATOM
	READ(5,*)IATOM(I),n,(l(j),j=1,n)
        do j=1,n
	lcase(I,l(j))=.true.
	end do
	write(6,*)'ATOM TYPE',IATOM(I),'  L:',(l(j),j=1,n)
	END DO
 1234 FORMAT(//,1A)
      WRITE(6,800)                                                      
      WRITE(6,805)  TITLE                                               
      WRITE(6,810)  LATTIC                                              
      WRITE(6,820)  AA,BB,CC                                            
      WRITE(6,840)  NAT                                                 
      WRITE(6,850)  IREL                                                
!      WRITE(6,870)  COORD  
      CALL LATGEN
      write(6,*)' alpha test',(alpha(i),i=1,3)
      CALL SYMOPER
      write(6,*)'SO=',SO
      IF (SO) THEN
      do i=1,3
      read(4,*)
      end do
      read(4,*)s(1),s(2),s(3)
      write(6,121)(S(i),i=1,3)
      write(21,121)(S(i),i=1,3)
121   format(' Spin-polarized + s-o calculation, M||',3f7.3)
      write(6,*)' alpha test',(alpha(i),i=1,3)
!      CALL ANGLE(S,THETA,PHI)   
      CALL ANGLE(S)   
      write(6,*)' alpha test',(alpha(i),i=1,3)
      pi=3.141592653
      th=theta*180/pi
      ph=phi*180/pi
111   format(2f6.1,' angle (M,z), angle (M,x) deg')
      write(6,111)th,ph 
      CALL SYM(MAGN)
      ELSE
      DO I=1,IORD
      det(i)=1.
      END DO  
      END IF

!....calculation of sym. op. in p,d,f basis
      CALL LMTMAT    

!....reading weights calculated by lapw2
      CALL READW
      WRITE(6,*)'WEIGHTS READ'
! read input for X operator
      krad=0
      kls=0
      read(5,*,end=11)krad,kls
      if(krad.ne.0)then
      do i=-20,20
       do j=1,20
        cx(i,j)=1.
       enddo
      enddo
      xrad(0) = ' Density matrix calculation '
      xrad(1) = ' I '
      xrad(2) = ' (1/r**3)  (large component only) '
      xrad(3) = ' (1/r**3)S (large component only) '
      xrad(4) = ' (1/r**3)S (large component only without rel. mass enh.) '
       xls(1)=' I '
       xls(2)=' S(dzeta) '
       xls(3)=' L(dzeta) '
       xls(4)=' s-o coupling (L.S) '
       xls(5)=' <L||alpha||L>*[-L*(L+1)*Sdzeta+(3/2){(L.S)Ldzeta+Ldzeta(L.S)}] '
       xls(6)=' L(dzeta)**2'
       xls(7)=' <L||alpha||L>*[-3*L(dzeta)**2+L*(L+1)]*S(dzeta) '
       comment(3,3)=' Bhf(orb) in T' 
       comment(4,3)=' Bhf(orb) in T, no rel. mass enhancement' 
       comment(3,5)=' Bhf(dip) in T' 
       comment(4,5)=' Bhf(dip) in T, no rel. mass enhancement' 
       comment(3,7)=' Bhf(dip) in T, diag. terms only' 
       comment(4,7)=' Bhf(dip) in T, diag. terms only, no rel. mass enhancement' 
      if(krad.gt.10)then
       write(6,108)krad-10
       write(21,108)krad-10
108    format(' Calculation of <X>, X=r**',i1)
      else if(krad.lt.-10)then
       write(6,109)abs(krad)-10
       write(21,109)abs(krad)-10
109    format(' Calculation of <X>, X=1/r**',i1)
      else
       write(6,*)' Calculation of <X>, X=c*Xr(r)*Xls(l,s)'
       write(21,*)' Calculation of <X>, X=c*Xr(r)*Xls(l,s)'
       write(6,102)xrad(krad)
       write(21,102)xrad(krad)
       write(6,101)xls(kls)
       write(21,101)xls(kls)
       if(((krad.eq.3).or.(krad.eq.4)).and. &
           ((kls.eq.3).or.(kls.eq.5).or.(kls.eq.7)))then
        cx(krad,kls)=12.5169
        write(21,103)cx(krad,kls),comment(krad,kls)
        write(6,103)cx(krad,kls),comment(krad,kls)
103     format('  c=',f9.5,a60)
       else
        write(21,105)cx(krad,kls)
        write(6,105)cx(krad,kls)
105     format('  c=',f9.5)
       endif
      endif
      write(6,*) &
      ' atom   L        up          dn         total'             
      write(21,*) &
      ' atom   L        up          dn         total'             
 101  format('  Xls(l,s) =',a70)
102   format('  Xr(r)    =',a70)
      endif
!      
11    continue
                                                                       
!.....CALCULATE CHARGE DENSITY CLM(R) IN SPHERES,  PARTIAL CHARGES      

      itape=10
      jtape=18
      call l2main(SO,magn)
      CALL CPUTIM(TTIME)
!      TCLM=TTIME                                          
!                                                                       
! 1213 CALL CPUTIM(TTIME)
      TFOUR=TTIME                                          
!                                                                       
!.....CALCULATE CPUTIME REQUIRED                                        
      TTOTAL=TFOUR-TSTART                                               
!      TFOUR=TFOUR-TCLM                                                  
!      PFOUR=TFOUR/TTOTAL*100.                                           
!      TCLM=TCLM-TFERMI                                                  
!      PCLM=TCLM/TTOTAL*100.                                             
!      TFERMI=TFERMI-TSTART                                              
!      PFERMI=TFERMI/TTOTAL*100.                                         
      WRITE(6,2000)                                                     
      WRITE(6,2010) TTOTAL,100.0                                        
!      WRITE(6,2020) TFERMI,PFERMI                                       
!      WRITE(6,2030) TCLM,PCLM                                           
!      WRITE(6,2040) TFOUR,PFOUR                                         
      CALL ERRCLR(ERRFN)
      STOP 'LAPWDM END'                                                 
!
!        error handling
!
  910 INFO = 1
!
!        'lapw2.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('lapwdm',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('lapwdm',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('lapwdm',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('lapwdm',ERRMSG)
      GOTO 999
  930 INFO = 3
!
!        illegal number of equivalent atoms
!
      CALL OUTERR('lapwdm','MULT .EQ. 0')
      GOTO 999
  950 INFO = 5
      CALL OUTERR('lapwdm','Too many atoms (NATO too small).')
      GOTO 999
  955 INFO = 6
      CALL OUTERR('lapwdm','Too many atoms (NDIF too small).')
      GOTO 999
  956 INFO = 56
      CALL OUTERR('lapwdm','LXDOS must be 3 for ISPLIT=999.')
      GOTO 999
  960 INFO = 7
!
!        Error reading file 'lapw2.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('lapwdm',ERRMSG)
      GOTO 999
  999 STOP 'lapwdm - Error'
!                                                                       
!                                                                       
  43  FORMAT(3X,A77)                                                    
  44  FORMAT(I3,A77)                                                    
 700  FORMAT(I3,A77)                                                    
 800  FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ',            &
             'I N F O R M A T I O N',/,30X,50(1H-),//)                  
 805  FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)                             
 810  FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
 820  FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
 830  FORMAT(3X,'SYMMETRY ATTITUDE IS',9X,'= ',A4)                      
 840  FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
 850  FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)                    
 860  FORMAT(3X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)                 
 870  FORMAT(3X,'TYPE OF COORDINATES IN DSPLIT= ',A5)                   
 1000 FORMAT(A80)                                                       
 1010 FORMAT(A4,24X,I2,1x,a4,/,13X,A4,18X,A4)                                 
 1020 FORMAT(6F10.7,10X,F10.7)                                          
 1030 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1031 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1040 FORMAT(///,3X,'ERROR IN lapwdm : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1002 FORMAT(3F10.5,I5)                                                 
 1003 FORMAT(A5)                                                        
 1004 FORMAT(A5,f10.5)                                                        
 1060 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
             //   ,':FER  :',1X,'F E R M I - ENERGY',11X,'= ',F9.5)            
 1061 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(TETRAH.M.)','= ',F9.5)           
 1062 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(GAUSS-.M.)','= ',F9.5)           
 1063 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(FERMI-SM.)','= ',F9.5)           
 2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)                       
 2010 FORMAT(12X,'TOTAL       : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2020 FORMAT(12X,'PART FERMI  : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2030 FORMAT(12X,'PART CLM    : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2040 FORMAT(12X,'PART FOURIR : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 6000 FORMAT(///,3X,'ERROR IN lapwdm : MULT(JATOM)=0 ...', &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
      END                                                               
