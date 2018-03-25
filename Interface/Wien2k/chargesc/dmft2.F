PROGRAM DMFT2
  USE param
  USE defs
  USE reallocate
  USE ams
  USE normal
  USE kpp1
  USE struk
  USE xa2
  USE com
  USE char
  USE dmf
  USE com_mpi
  USE GMPI
  use ghdf5_base
  IMPLICIT REAL*8 (A-H,O-Z)

#ifdef MPI
  include 'mpif.h'
#endif

  REAL*8, allocatable       ::  ZZ(:)
  REAL*8,PARAMETER       :: Ry2eV= 13.60569193
  REAL*8           GMAX
  CHARACTER*4      IREL,RELA,cform,RCFILE
  CHARACTER*5      COORD, CHARS
  CHARACTER*10    KNAME
  CHARACTER*11     STATUS,FORM                                      
  CHARACTER*67       ERRMSG
  CHARACTER*80       DEFFN, ERRFN
  CHARACTER*77     RMARGN                                           
  CHARACTER*180     FNAME,FNAME1
  CHARACTER*200 :: vec_up, vec_dn, ene_up, ene_dn
  INTEGER :: latom, i, icix, wicix, wndim, size, minsigind
  REAL*8  :: xx(3), xz(3), thetaz, phiz, ddd, check, sw(3)
  INTEGER :: INDEX, INDEX1, MI, JR, JC
  REAL*8 :: POST(3)
  DATA             RELA/'RELA'/                           
  save TTOTAL,TCLM,TFOUR,tstart,errfn                          

  CALL start_MPI()

  Qcomplex=.FALSE.
  nlo=0
  RCFILE='NOFI'
!-----------------------------------------------------------------------  
!                                                                       
  CALL init_ams

  CALL GTFNAM(DEFFN,ERRFN)

  if (myrank.EQ.master) CALL ERRFLG(ERRFN,'Error in DMFT2')

  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)

  DO
     READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL

     lngth = Len_Trim(FNAME)
     if ( FNAME(lngth-1:lngth).EQ.'_x' ) then
        FNAME = FNAME(:lngth-2)
        vector_para = .True.
     else if (IUNIT.EQ.6) then
        if (myrank.EQ.master .OR. fastFilesystem/=0) then
           CALL FilenameMPI(FNAME)    ! Each processor will have different information file
           OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        endif
     else if (iunit.eq.8.or.iunit.eq.21) then
        if (myrank.eq.master) OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
     elseif (iunit.eq.13)then
        fname1=fname
     else
        OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
     endif

     if(iunit.eq.18) then
        do i=180,5,-1
           if(fname(i-4:i).eq.'vspdn')then
               sspin="2"
               exit
           endif
        enddo
     endif
     if(iunit.eq.3) then
        read(iunit,'(A5)',end=12) CHARS
        Qcomplex=.true. ! switches Qcomplex  on if case.in1c exists and is nonzero
12      continue
     end if
  ENDDO
20 CONTINUE
  CLOSE (1)

  if (nprocs.EQ.1) then
     write(6,*)'Running DMFT2 in single processor mode'
     write(6,*)
  else if (myrank.EQ.master .OR. fastFilesystem/=0) then
     write(6,*)'Running DMFT2 in mpi mode'
     write(6,*)
  endif
  
  !                                                                       
  !.....READ STRUCT                                                       
  READ(20,1000) TITLE                                               
  READ(20,1010) LATTIC,NAT,cform,IREL    
  allocate (rmt(nat),v(nat),iatnr(nat),mult(nat),isplit(nat))
  allocate (rotloc(3,3,nat))
  allocate (r0(nat),dx(nat),jri(nat))
  allocate (aname(nat))
  allocate (ZZ(nat))
  allocate ( pos(3,48*nat) )
  !.....READ IN LATTICE CONSTANTS                                         
  READ(20,1020) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
  IF(ABS(ALPHA(1)).LT.TEST) ALPHA(1)=NINETY
  IF(ABS(ALPHA(2)).LT.TEST) ALPHA(2)=NINETY
  IF(ABS(ALPHA(3)).LT.TEST) ALPHA(3)=NINETY
  IF(IREL.EQ.'RELA') REL=.TRUE.                                     
  IF(IREL.EQ.'NREL') REL=.FALSE.                                    
  !
  !  read crystal-structure (atompositions, symmetry-parameters,
  !                                muffin-tin radius, ...)
  !  'INDEX' counts all atoms in the unit cell,
  !  'JATOM' counts only the non-equivalent atoms
  INDEX=0                                                           
  DO JATOM = 1,NAT                                               
     INDEX=INDEX+1
     READ(20,1030) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM),ISPLIT(JATOM)
     if(isplit(jatom).eq.99)    lxdos=   3           
     if(isplit(jatom).eq.88)    lxdos=   3          
     IF (MULT(JATOM) .EQ. 0) THEN
        !...illegal number of equivalent atoms
        if (myrank.EQ.master .OR. fastFilesystem/=0) WRITE (6,6000) JATOM, INDEX, MULT(JATOM)
        GOTO 930
     ENDIF
     DO M=1,MULT(JATOM)-1                                     
        INDEX=INDEX+1                                            
        READ(20,1031) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
     ENDDO
     READ(20,1050) ANAME(JATOM),JRI(JATOM),R0(JATOM),RMT(JATOM),ZZ(JATOM)
     DX(JATOM)=LOG(RMT(JATOM)/R0(JATOM)) / (JRI(JATOM)-1)           
     RMT(JATOM)=R0(JATOM)*EXP( DX(JATOM)*(JRI(JATOM)-1) )           
     READ(20,1051) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)                
  ENDDO

  allocate (pu1u1(0:lxdos,0:lxdos),pu1ue(0:lxdos,0:lxdos),pu1u2(0:lxdos,0:lxdos,nloat),pueue(0:lxdos,0:lxdos),pueu2(0:lxdos,0:lxdos,nloat),pu2u2(0:lxdos,nloat,0:lxdos,nloat))
  ALLOCATE (rotij(3,3,index),tauij(3,index))
  ndif=index
  CALL doreallocate(pos, 3, ndif)
  !.....READ INPUT AND POTE                                               
  READ(5,1003)  MODUS,MODUS1,coord                                        
  if (myrank.EQ.master .OR. fastFilesystem/=0) write(6,*)  ' Modus: ', MODUS
  esepermin=0.5d0
  eseper0=0.05d0
  READ(5,*,err=4716)  EMIN,ELECN,esepermin,eseper0
  if (myrank.EQ.master .OR. fastFilesystem/=0) write(6,*) 'no read error'
4716 continue

  READ(5,*) efmod,GDELTA
  DO I=1,NAT
     READ(5,*) IDUMMY
  ENDDO
  READ(5,*,END=4715) GMAX

  READ(5,'(A4)',END=4715) RCFILE
4715 CONTINUE
  
  REWIND 5

  READ(5,1234) DUMMY
  if (myrank.EQ.master) write(21,720)gmax
  
  if (myrank.EQ.master .OR. fastFilesystem/=0) WRITE(6,*) 'RECPR LIST: ',RCFILE
     
  if (myrank.EQ.master .OR. fastFilesystem/=0) then
     WRITE(6,800)                                                      
     WRITE(6,805)  TITLE                                               
     WRITE(6,810)  LATTIC                                              
     WRITE(6,820)  AA,BB,CC                                            
     WRITE(6,840)  NAT                                                 
     WRITE(6,850)  IREL                                                
     WRITE(6,870)  COORD
  endif
  REL=.FALSE.                                                       
  IF(IREL.EQ.RELA) REL=.TRUE.                                       

  CALL LATGEN(NAT)
  CALL CPUTIM(TTIME)
  TSTART=TTIME
  call walltim(TTIME)
  TSTART_w=TTIME
  
  CALL RECFIL(FNAME1,GMAX,RCFILE,KXMAX,KYMAX,KZMAX,NWAVE)

  if (vector_para) then
     FNAME = '_processes_'
     CALL FilenameMPI2(FNAME)
     nvector=0
     vectors=0
     open(999,FILE=FNAME,STATUS='OLD',ERR=88,IOSTAT=ios) ! opens _processes_x
     nkp=0
     DO
        READ (999,*,END=88,ERR=970) ivector,inkp,ik_start,vec_up,vec_dn,ene_dn,ene_up  ! can jump to 20
        nvector = nvector+1
        nkp = nkp+inkp
        vectors(nvector,1) = ivector
        vectors(nvector,2) = inkp
        vectors(nvector,3) = ik_start
        fvectors(nvector,1) = vec_up
        fvectors(nvector,2) = vec_dn
        fvectors(nvector,3) = ene_dn
        fvectors(nvector,4) = ene_up
        if (myrank.EQ.master .OR. FastFilesystem/=0) WRITE(6,*) 'FNAMES=', fvectors(nvector, :)
     ENDDO
88   CONTINUE
     close(999)
  else
     nvector=1
  endif

  call gh5_init()
#ifdef MPI
  CALL INI_GMPI(MYRANK,NPROCS,MASTER,VECTOR_PARA,NVECTOR,VECTORS)
#endif


  if(myrank.EQ.master)then
      call gh5_open_r('GPARAMBANDS.h5', f_id)
      call gh5_read(nkpt, "/kptdim", f_id)
      call gh5_read(nmat, "/nkmax", f_id)
      call gh5_read(nume, "/nbmax", f_id)
      call gh5_read(iso, "/iso", f_id)
      call gh5_read(ispin_pol, "/ispin", f_id)
      call gh5_close(f_id)
#ifdef MPI
      call mpi_bcast(nkpt,1,mpi_integer,master,mpi_comm_world,ierr)
      call mpi_bcast(nmat,1,mpi_integer,master,mpi_comm_world,ierr)
      call mpi_bcast(nume,1,mpi_integer,master,mpi_comm_world,ierr)
      call mpi_bcast(iso,1,mpi_integer,master,mpi_comm_world,ierr)
      call mpi_bcast(ispin_pol,1,mpi_integer,master,mpi_comm_world,ierr)
#endif
  endif

  iso_orig = iso
  iso = max(1,iso/ispin_pol)

  if (.not.vector_para) then
     vectors(1,2)=nkpt
  endif

  if(iso_orig==2)then
      Qcomplex=.true.
  endif

  !.....CALCULATE CHARGE DENSITY CLM(R) IN SPHERES,  PARTIAL CHARGES      
  call l2main(cform,nlo,coord,zz,nwave,TCLM,TCLM_w,TFOUR,TFOUR_w)
  
  !.....CALCULATE CPUTIME REQUIRED                                        
  TTOTAL=TFOUR-TSTART                                               
  TFOUR=TFOUR-TCLM                                                  
  PFOUR=TFOUR/TTOTAL*100.                                           
  TCLM=TCLM-TSTART                                                  
  PCLM=TCLM/TTOTAL*100.                                             
  if (myrank.EQ.master .OR. fastFilesystem/=0) then
     WRITE(6,2000) 
     WRITE(6,2010) TTOTAL,100.0
     WRITE(6,2030) TCLM,PCLM
     WRITE(6,2040) TFOUR,PFOUR
  endif

  TTOTAL_w=TFOUR_w-TSTART_w                                               
  TFOUR_w=TFOUR_w-TCLM_w                                                  
  PFOUR_w=TFOUR_w/TTOTAL_w*100.                                           
  TCLM_w=TCLM_w-TSTART_w                                                  
  PCLM_w=TCLM_w/TTOTAL_w*100.                                             
  if (myrank.EQ.master .OR. fastFilesystem/=0) then
     WRITE(6,2001)
     WRITE(6,2010) TTOTAL_w,100.0
     WRITE(6,2030) TCLM_w,PCLM_w
     WRITE(6,2040) TFOUR_w,PFOUR_w
  endif
  if (myrank.EQ.master) CALL ERRCLR(ERRFN)

  call gh5_end()

  CALL stop_MPI()
  
  STOP ! ' DMFT2 END'                                                 
!
!        error handling
!
910 INFO = 1
!
!        'lapw2.def' couldn't be opened
!
  WRITE (ERRMSG,9000) FNAME
  CALL OUTERR('DMFT2',ERRMSG)
  GOTO 999
970 INFO = 8
!
! Error in parallel: could not find .processes.x
!
  WRITE (ERRMSG,'(A,A)')  'file open error:', trim(ADJUSTL(FNAME))
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999
920 INFO = 2
!
!        file FNAME couldn't be opened
!
  WRITE (ERRMSG,9010) IUNIT
  CALL OUTERR('DMFT2',ERRMSG)
  WRITE (ERRMSG,9020) FNAME
  CALL OUTERR('DMFT2',ERRMSG)
  WRITE (ERRMSG,9030) STATUS, FORM
  CALL OUTERR('DMFT2',ERRMSG)
  GOTO 999
930 INFO = 3
!
!        illegal number of equivalent atoms
!
  CALL OUTERR('DMFT2','MULT .EQ. 0')
  GOTO 999
956 INFO = 56
  CALL OUTERR('DMFT2','LXDOS must be 3 for ISPLIT=99 or 88')
  CALL OUTERR('DMFT2','LXDOS must be 3 for ISPLIT=99 or 88')
  GOTO 999
960 INFO = 7
!
!        Error reading file 'lapw2.def'
!
  WRITE (ERRMSG,9040) FNAME
  CALL OUTERR('DMFT2',ERRMSG)
  GOTO 999
999 STOP 'DMFT2 - Error. Check file lapw2.error'
!                                                                       
!                                                                       
121 format('        Energy to separate low and high energy','states: ',f10.5)
720 FORMAT(':GMA  :',' POTENTIAL AND CHARGE CUT-OFF',f7.2,' Ry**.5')
800 FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ','I N F O R M A T I O N',/,30X,50(1H-),//)                  
805 FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)                             
810 FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
820 FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
830 FORMAT(3X,'SYMMETRY ATTITUDE IS',9X,'= ',A4)                      
840 FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
850 FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)                    
860 FORMAT(3X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)                 
870 FORMAT(3X,'TYPE OF COORDINATES IN DSPLIT= ',A5)                   
901 format(a,i8,2F16.10,D18.6)
1000 FORMAT(A80)                                                       
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
1040 FORMAT(///,3X,'ERROR IN DMFT2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F10.8)                                             
1002 FORMAT(3F10.5,I5)                                                 
1003 FORMAT(A5,a5,a5)                                                        
1234 FORMAT(//,1A)
2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)                       
2001 FORMAT(//,3X,'=====>>> WALL TIME SUMMARY',/)                       
2010 FORMAT(12X,'TOTAL       : ',F8.1,5X,'... ',F4.0,' PERCENT')       
2030 FORMAT(12X,'PART CLM    : ',F8.1,5X,'... ',F4.0,' PERCENT')       
2040 FORMAT(12X,'PART FOURIR : ',F8.1,5X,'... ',F4.0,' PERCENT')       
2061 FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')    !EECE
6000 FORMAT(///,3X,'ERROR IN DMFT2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
9000 FORMAT('can''t open definition file ',A40)
9010 FORMAT('can''t open unit: ',I2)
9020 FORMAT('       filename: ',A50)
9030 FORMAT('         status: ',A,'  form: ',A)
9040 FORMAT('Error reading file: ',A47)
877 FORMAT(' Energy window:',f10.4,' < E < ',f10.3)
887 FORMAT(' Greens function calculation for',i3,' atoms')
120 FORMAT(' New z axis || ',3f9.4)
1013 FORMAT('LOCAL ROT MATRIX:   ',3f10.7,/,20x,3f10.7,/,20x,3f10.7)
126 FORMAT(' Magnetic system with s-o coupling; M theta, phi:',2f8.4) 
END PROGRAM DMFT2
