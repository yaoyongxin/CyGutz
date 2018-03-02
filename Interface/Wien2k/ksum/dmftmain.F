PROGRAM DMFTMAIN
  USE param
  USE struct
  USE case
  USE sym2
  USE com
  USE abc
  USE kpts, only : numkpt, MKNAME, MIPGR, mweight, tweight, &
        &MSX, MSY, MSZ, allocate_kpts, deallocate_kpts
  USE com_mpi
  use ghdf5_base

  IMPLICIT REAL*8 (A-H,O-Z)
  integer,external::get_sym_ie

  complex*16, allocatable :: trmat(:,:)
  real*8, allocatable     :: rcf(:,:)
  complex*16      qtl1(24), IMAG
  dimension       xx(3),xz(3),acf(14),bcf(14),qtle(34)
  CHARACTER*4     adum
  CHARACTER*5     MODSYM,CHAR
  CHARACTER*10    KNAME
  CHARACTER*11    STATUS,FORM                                      
  CHARACTER*67    ERRMSG
  CHARACTER*80    DEFFN, ERRFN
  CHARACTER*200   FNAME, vec_up, vec_dn, ene_up, ene_dn
  CHARACTER*161   TEXT
  INTEGER,ALLOCATABLE :: ieband(:)        
  !INTEGER         IPROC
  LOGICAL         Qcomplex, newform
  INTEGER         kindex, nsymop
  INTEGER         wndim, icix, wicix, size, imatsubara
  REAL*8       :: EF_LDA, Ry2eV
  REAL*8       :: POST(3), BKRLOC(3,3)
  REAL*8,ALLOCATABLE:: wtmp(:)
  INTEGER      :: projector
  INTEGER      :: strlen, locrot, latom, natm, wfirst, iat, im, wat
  COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
  COMMON /RADFU/  RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2), RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)
  DIMENSION  S(3)
  DATA Qcomplex /.false./, IMAG/(0.0D0,1.0D0)/, Ry2eV/13.60569193D0/
  DATA strlen/200/


  CALL start_MPI()
  CALL GTFNAM(DEFFN,ERRFN)
  call set_gfile(deffn)

  if (myrank.EQ.master) CALL ERRFLG(ERRFN,'Error in DMFT')
  
  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910,IOSTAT=ios) ! opens gwien1.def file

  vector_para = .False.  
  ! The following few lines read gwien1.def file 
  ! containing all filename definitions
  DO 
     READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL  ! can jump to 20
     if(iunit.eq.6 .AND. (myrank.EQ.master .OR. fastFilesystem)) &
            &CALL FilenameMPI(FNAME)    
     ! Each processor will have different information file 

     lngth = Len_Trim(FNAME)
     if ( FNAME(lngth-1:lngth).EQ.'_x' ) then
        FNAME = FNAME(:lngth-2)
        vector_para = .True.
     else
        if ( iunit.lt.180  .or. myrank.eq.master .or. &
            &(iunit.eq.6.and.fastFilesystem)) then            
        ! output files (iunit>=100) should be open on master node only.
           OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        endif
     endif

     if(iunit.eq.9)  VECFN(1)=FNAME
     if(iunit.eq.10) VECFN(2)=FNAME
     if(iunit.eq.59) VECFN(3)=FNAME
     if(iunit.eq.60) VECFN(4)=FNAME
     if(iunit.eq.18) then
        do i=strlen,5,-1
           if(fname(i:i).ne.' ') then
              if(fname(i-2:i).eq.'pup') ispin_pol=2
              if(fname(i-2:i).eq.'pdn') ispin_pol=2 
           endif
        enddo
     endif
     if(iunit.eq.7) then
        read(iunit,123,end=12)CHAR
        Qcomplex=.true. ! switches Qcomplex  on if case.in1c exists and is nonzero
12      continue
     end if
     if(iunit.eq.9) then
        iloop1: do i=strlen,8,-1
            if(fname(i-7:i).eq.'vectorso') then
               iso=2
               Qcomplex=.true.
               exit iloop1
            endif
        enddo iloop1
     endif
  ENDDO ! 10
20 CONTINUE
  CLOSE (1)
  ! finished gwien1.def

  isospin=max(1,ispin_pol/iso)

  !.....READING STRUCT 
  ! It reads case.struct file and saves data into module named "param". These variables are available in dmftmain.f.
  CALL init_struct                                                      
  if (myrank.EQ.master .OR. fastFilesystem) then
     write(6,'(A)')'******Program DMFT: Computes Gloc and Delta for Dynamical Mean Field Theory calculation ******'
     if (Qcomplex) write(6,*)'CALCULATION WITH COMPLEX VECTORS'
     
     if (iso==2)then
        if(ispin_pol.eq.2)then
            write(6,*)'CALCULATION WITH SOC-UP AND DN VECTORS REQUIRED'
        else
            write(6,*)'CALCULATION WITH SOC VECTORS REQUIRED'
        endif
     endif
  endif

  if (nprocs.EQ.1) then
     write(6,*)'Running DMFT in single processor mode'
     write(6,*)' '
  else
     write(6,*)'Running DMFT in mpi mode'
     write(6,*)' '
  endif

  ! read klist and associated weights from file 66==case.klist
  ! first read how many k-points exists and the old/new file type
  ITAPE = 13
  numkpt = 0
  newform = .TRUE.
  DO
     READ (ITAPE, '(A20)', IOSTAT=ios) TEXT
     IF (numkpt==0 .AND. TEXT(15:16) .NE. ' ') THEN
        newform=.FALSE.
     ENDIF
     KNAME = TEXT(1:10)
     IF (KNAME .EQ. 'END       ' .OR. ios.ne.0) EXIT
     numkpt = numkpt + 1
  ENDDO
  NKPT=numkpt
  REWIND(ITAPE)

  if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,*) 'numkpt=', numkpt
 
  CALL allocate_kpts(numkpt)
  MWEIGHT=0
  DO kindex=1,numkpt
     IF(newform) THEN
        READ (ITAPE,5101) MKNAME(KINDEX), ISX, ISY, ISZ, IDV, MWEIGHT(KINDEX), E1, E2, MIPGR(KINDEX)
     ELSE
        READ (ITAPE,5100) MKNAME(KINDEX), ISX, ISY, ISZ, IDV, MWEIGHT(KINDEX), E1, E2, MIPGR(KINDEX)
     ENDIF
     MSX(KINDEX) = DBLE(ISX)/DBLE(IDV)
     MSY(KINDEX) = DBLE(ISY)/DBLE(IDV)
     MSZ(KINDEX) = DBLE(ISZ)/DBLE(IDV)
     if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,'(4f10.5,1x)') MSX(kindex), MSY(kindex), MSZ(kindex), MWEIGHT(kindex)
  END DO
  TWEIGHT = sum(MWEIGHT)
  MWEIGHT=MWEIGHT/TWEIGHT

  ! iord : number of symmetry operations from case.struct
  CALL init_sym2(iord) ! allocates memory for: idet(iord),opimat(3,3,iord),phase(iord),tmat(3,3,iord)
  
  if (myrank.EQ.master .OR. fastFilesystem) then
     WRITE(6,800) 
     WRITE(6,805)   TITLE                                               
     WRITE(6,810)   LATTIC                                              
     WRITE(6,820)   AA,BB,CC                                            
     WRITE(6,840)   NAT                                                 
     WRITE(6,850)   IREL                             
     write(6,*)'IORD=',IORD
  endif
 
  if (vector_para) then
     FNAME = '_processes_'
     CALL FilenameMPI2(FNAME)
     !WRITE(6,*) 'Trying to open', FNAME
     nvector=0
     vectors=0
     open(999,FILE=FNAME,STATUS='OLD',ERR=88,IOSTAT=ios) ! opens _processes_x                                                                                 
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
        if (myrank.EQ.master .OR. FastFilesystem) WRITE(6,*) 'FNAMES=', fvectors(nvector, :)
     ENDDO
88   CONTINUE
     close(999)
  else
     nvector=1
  endif
!     LATGEN GENERATES TWO BRAVAIS MATRICES, DEFINES THE VOLUME OF      
!     THE UNIT CELL AND CALLS ROTDEF                                    
!     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS      
!                 GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN      
!                 COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN      
!                 WFTAPE) INTO CARTESIAN SYSTEM                         
!     BR2(3,3) :  TRANSFORMS A RECIPROCAL LATTICE VECTOR OF A SPE-      
!                 CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )       
!                 TO CARTESIAN SYSTEM                                   
!      
  CALL LATGEN ! for results look at case.outputq  


! Reads current input file case.indmfl
! 5 is connected to current input file case.indmfl
  !print *, 'Start reading 5'
  READ(5,*) EMIN,EMAX,projector
  EMIN = EMIN/Ry2eV; EMAX = EMAX/Ry2eV
  read(5,*) natom

  if (myrank.EQ.master .OR. fastFilesystem) then
     write(6, 877) Emin,Emax  ! Energy window
     write(6, 887) natom      ! Projected density of states calculation for xxx atoms
  endif
  
  ALLOCATE(nl(natom))
  allocate(ll(natom,4), qsplit(natom,4), cix(natom,4), iatom(natom))
  ALLOCATE(crotloc(3,3,natom))
  crotloc=0

  natm = sum(mult)
  ALLOCATE( isort(natm), ifirst(natm) )
  wfirst  = 1          ! first atom of this particular sort
  do iat=1,nat         ! over all sorts
     do im=1,MULT(iat) ! over all atoms of this sort
        wat = wfirst + im-1  ! atom number
        isort(wat)=iat       ! sort of each atom
        ifirst(wat)=wfirst
     enddo
     wfirst = wfirst + MULT(iat)
  enddo

  DO i=1,natom
     READ(5,*) latom,nL(i) ! read from case.inq
     iatom(i) = latom       ! The succesive number of atom (all atoms counted)
     jatom = isort(latom)   ! The sort of the atom
     if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,'(A,I2,A)'), '--------- atom ', i, ' -------------'
     do j=1,nL(i)
        read(5,*) LL(i,j), qsplit(i,j), cix(i,j) ! LL contains all L-quantum numbers which need to be computed.
     enddo

     if (myrank.EQ.master .OR. fastFilesystem) then
        do j=1,nL(i)
           write(6, '(A,I3,2x,A,I3)') 'l=', LL(i,j), 'qsplit=', qsplit(i,j)!, 'nrfmax=', nrfmax(i,j)
        enddo
        write(6,*)'Symmetrization over eq. k-points is performed'
     endif

     crotloc(:,:,i) = rotloc(:,:,jatom)  ! Rotation for this type from the struct file

  enddo

  READ(5,*) ! comment: Next few lines contain instructions (transformation,index) for all correlated orbitals
  READ(5,*) ncix, maxdim ! number of independent cix blocks

  if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,*) '********** Start Reading Cix file *************'
  ALLOCATE(CF(maxdim,maxdim,ncix))      ! transformation matrix
  ALLOCATE(wtmp(2*maxdim))              ! temp
  
  CF=0
  do icix=1,ncix
     READ(5,*) wicix, wndim  ! icix, size-of-matrix, L
     if (wicix.ne.icix) then
        print *, 'Something wrong reading case.indmfl file. Boilig out...'
        goto 999
     endif
     READ(5,*) ! Comment: Transformation matrix follows
     do i=1,wndim
        READ(5,*) (wtmp(j),j=1,2*wndim)
        do j=1,wndim
           CF(i,j,icix) = dcmplx(wtmp(2*j-1),wtmp(2*j))
        enddo
     enddo

     if (myrank.EQ.master .OR. fastFilesystem) then
        write(6,*)' Correlated block number', icix
        write(6,*)' Real part of unitary matrix:'
        do m1=1,wndim
           do m2=1,wndim
              write(6,'(F8.4)',advance='no') dble(cf(m1,m2,icix))
           enddo
           write(6,*)
        enddo
        write(6,*)' Imaginary part of unitary matrix:'
        do m1=1,wndim
           do m2=1,wndim
              write(6,'(F8.4)',advance='no') aimag(cf(m1,m2,icix))
           enddo
           write(6,*)
        enddo
     endif
  enddo
  DEALLOCATE(wtmp)

  ! For spin-orbit coupling crates transformation "rotloc"      
  !....reading *.inso
  if (iso==2) then
     do i=1,3
        read(4,*)
     end do
     read(4,*)s(1),s(2),s(3)
     call angle(s,theta,fi)
     if (myrank.EQ.master .OR. fastFilesystem) write(6,126)theta,fi
     !write(21,126)theta,fi
     !.... test of collinearity of the z-axis with spin quantization axis
     cost=cos(theta)
     cosf=cos(fi)
     sint=sin(theta)
     sinf=sin(fi)
     xz(1)=sint*cosf
     xz(2)=sint*sinf
     xz(3)=cost
     do i=1,natom
        dsum=crotloc(3,1,i)*xz(1)+crotloc(3,2,i)*xz(2)+crotloc(3,3,i)*xz(3)
        if(abs(dsum-1.).gt.0.00001)then  ! old z not along M, new rotloc
           crotloc(1,1,i)=cosf*cost
           crotloc(2,1,i)=-sinf
           crotloc(3,1,i)=cosf*sint
           crotloc(1,2,i)=sinf*cost
           crotloc(2,2,i)=cosf
           crotloc(3,2,i)=sinf*sint
           crotloc(1,3,i)=-sint
           crotloc(2,3,i)=0
           crotloc(3,3,i)=cost
           if (myrank.EQ.master .OR. fastFilesystem) then
              write(6,*)' New local rotation matrix with z || M'
              write(6,1013)((crotloc(j,j1,i),j1=1,3),j=1,3)  ! written as in .struct
           endif
        endif
     end do
  end if  ! spin-orbit
  
  !---------------------------------------------------------------------------------!
  !--- WRITE Output for the BR1 matrix and the ROTIJ Matrix in file case.rotlm:  ---!
  !---  they are read by LectureROTIJ in lecture.f in SRC_elnes / or SRC_mdff/.  ---!
  !---  These lines are in SRC_lapw2lm/l2main.frc: surch *PH*                    ---!
  !---  The number of the output file is 22. his extention must be .rotlm        ---!
  !---------------------------------------------------------------------------------!
  CALL ROTDEF 
  if (myrank.EQ.master) then

     WRITE(22,*) 'BR1'
     DO JR=1, 3
        WRITE(22, 9101) BR1(1,JR), BR1(2,JR), BR1(3,JR) !--- Writting BR1 ----!
     ENDDO
     WRITE(22,*) 'BR2'
     DO JR=1, 3
        WRITE(22, 9101) BR2(1,JR), BR2(2,JR), BR2(3,JR) !--- Writting BR1 ----!
     ENDDO
 
     INDEX = 0
     DO LATOM=1, NAT
        INDEX1 = INDEX+1
        DO LATEQ=1, MULT(LATOM)
           INDEX = INDEX + 1
           WRITE(22, 9100) LATOM, LATEQ, INDEX !---- inequivalent atomnumber  ---!
           DO JR=1, 3
              WRITE(22, 9101) (ROTIJ(JC,JR,INDEX),JC=1,3) !--- Writting ROTIJ ---!
           ENDDO
           ! Rotates first atom by rotij, the transformation which makes the
           ! first and current atom equivalent.
           POST = matmul(POS(:,index1), ROTIJ(:,:,index))
           POST = POST + TAUIJ(:,index)
           write(6,'(A,I2,A)',advance='no') 'Actual position of atom ', INDEX, ' is:'
           write(6,'(f10.6,2x,f10.6,2x,f10.6)') POST
        ENDDO
     ENDDO
     write(6,*)
     write(6,'(A)') 'Combined transformation (acting on k-point) for correlated atoms'
     write(6,'(A)') 'including users rotation and internal local rotation'
     DO icase=1,natom
        latom = iatom(icase)
        write(6, '(A,I3,1x,A,1x,I3)') 'catom', icase, 'atom', latom
        BKRLOC = matmul(crotloc(:,:,icase),matmul(BR1,rotij(:,:,latom)))
        DO JR=1,3
           WRITE(6, '(3F10.5)') (BKRLOC(JR,JC),JC=1,3)
        ENDDO
        WRITE(6,*)
     ENDDO
     write(6,*)
  endif

  ! If SO-coupling and symmetrization over all k-points
  ! not trivial
  if (iso==2) then
     CALL SYM(THETA,FI)
  endif

  !*****  finds nmat and nume in energy file *****
  ! 59 -> case.energydn
  ! 60 -> case.energy
  ! using case.energy, we will find out what is:
  !  nmat -- maximum number of reciprocal vectors
  !  nume -- maximum number of bands kept
  k=0
  nmat=0
  nume=0
  DO ivector=1,nvector

     if (vector_para) then
        DO is=1,iso    !------ over up/dn ---------------------!
           FNAME = fvectors(ivector, is+2)
           !FNAME = TRIM(VECFN(is+2))//'_'
           !CALL FilenameMPI3(FNAME, vectors(ivector,1))
           itape=61-is
           WRITE(6,*) 'FNAME=', FNAME, itape
           open(itape,FILE=FNAME,STATUS='old',FORM='formatted')
        ENDDO
     endif

     do is=1,isospin
        itape = 61-is
        DO I=1,NAT
           READ(itape,'(f9.5)',END=1005) EMIST
           READ(itape,'(f9.5)',END=1005) EMIST
        ENDDO
        ios=0
        DO WHILE (ios == 0)
           READ(itape,'(3e19.12,a10,2i6)',IOSTAT=ios) SS,T,Z,KNAME,N,NEn
           IF (ios /= 0) CYCLE
           k=k+1
           nmat=MAX(n,nmat)
           nume=MAX(nen,nume)
           DO ii=1,nen
              READ(itape,*) NUM,E1
           ENDDO
        ENDDO
     end do
1005 CONTINUE

     DO is=1,isospin    !------ over up/dn ---------------------!
        itape=61-is
        if (vector_para) then
           close(itape)
        else
           REWIND(itape)
        endif
     END DO
  END DO

#ifdef MPI
  CALL IMAX1_ALL_MPI(NUME)
#endif

  CALL init_abc(nume,nmat,lmax2,ndim2,nrf) ! Allocates variables:

  READ(3,*); READ(3,*)TMP,ELECN; READ(3,*)efmod,delta

  call gh5_init()
  sym_ie=get_sym_ie(iord,iz,tau)
  IF(MYRANK.eq.MASTER)THEN
      if(gfile==2)then
          call gh5_open_rw('GPARAMBANDS.h5', log_file_id)
      else
          call gh5_open_w('GPARAMBANDS.h5', log_file_id)
      endif
      if(gfile==1)then
          call gh5_write(iso,'/iso/',log_file_id)
          call gh5_write(ispin_pol,'/ispin',log_file_id)
          call gh5_write(nume,'/nbmax',log_file_id)
          call gh5_write(iord,'/symnop',log_file_id)
          call gh5_write(sym_ie,'/symie',log_file_id)
          call gh5_write(nmat,'/nkmax',log_file_id)
      endif
  ENDIF

#ifdef MPI
  call gmpi_wrt()
#endif

  CALL GUTZ_SET_EF_LDA(NUME,NUMKPT,projector)
  if (myrank.EQ.master .OR. fastFilesystem) then
     WRITE(6 ,1060)  ELECN,EF
  endif
  nsymop = iord

  call l2main(Qcomplex,nsymop,projector)  !!!! This calls the main subroutine

  IF(MYRANK.eq.MASTER)THEN
      call gh5_close(log_file_id)
  ENDIF
  call gh5_end()

  DEALLOCATE(CF)
  DEALLOCATE(isort, ifirst)
  DEALLOCATE(crotloc)
  CALL deallocate_kpts


  CALL CPUTIM(TTIME)
                                                                       
!.....CALCULATE CPUTIME REQUIRED                                        
  if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,2000)                                                     
  if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,2010) TTIME                                        
  if (myrank.EQ.master) CALL ERRCLR(ERRFN)

  CALL stop_MPI()

  STOP !' DMFT END'

  
!   error handling
!
910 INFO = 1
!
!   'DMFTx.def' couldn't be opened
!910
  print *, 'Definition file error', TRIM(DEFFN), 'could not be opened!'
  print *, 'ios=', ios
  WRITE (ERRMSG,9000) DEFFN
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999

920 INFO = 2
!
!  file FNAME couldn't be opened
!
  WRITE (ERRMSG,9010) IUNIT
  CALL OUTERR('DMFT',ERRMSG)
  WRITE (ERRMSG,9020) FNAME
  CALL OUTERR('DMFT',ERRMSG)
  WRITE (ERRMSG,9030) STATUS, FORM
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999
951 INFO = 5
  write(6,*) 'error:',FNAME
  CALL OUTERR('DMFT','file open error')
  GOTO 999
960 INFO = 7
!
!        Error reading file 'dmft1.def'
!
970 INFO = 8
!
! Error in parallel: could not find .processes.x
!
  WRITE (ERRMSG,'(A,A)')  'file open error:', trim(ADJUSTL(FNAME))
  WRITE (ERRMSG,9040) FNAME
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999
999 STOP 'DMFT - Error'
!                                                                       
!                                                                       

9100 FORMAT('inequivalent atomnumber ',I3,' number ',I2,' total ',I4)
9101 FORMAT(3F10.5)

123 FORMAT(A5)
815 FORMAT(' Compound:',a60)
666 FORMAT(a4)
667 FORMAT(38x,f10.5)
887 FORMAT(' Projected density of states calculation for',i3,' atoms')
877 FORMAT(' Energy window:',f10.4,' < E < ',f10.3)
4030 FORMAT(' Ordering of DOS in DMFT file for: ',a60,/)
120 FORMAT(' New z axis || ',3f9.4)
121 FORMAT(' New x axis || ',3f9.4)
1013 FORMAT('users (crotloc) matrix: ',3f10.7,/,24x,3f10.7,/,24x,3f10.7)
101 FORMAT(' L=',i2,'. Unitary transFORMATion read from input')
555 FORMAT(7f9.4)
111 FORMAT(' Population matrix diagonal in L for L=',4i3)
110 FORMAT(' Complete population matrix with <L1|L2> terms for L=',4i3)
112 FORMAT(' Number of cases for TETRA=',i4)
776 FORMAT(4f12.6,' Emin,dE,Emax,dgauss')
777 FORMAT(i5,'                         number of cases for tetra')
126 FORMAT(' Magnetic system with s-o coupling; M theta, phi:',2f8.4) 
4031 FORMAT(' atom',i4,' ordering of projected DOS ')
987 FORMAT(a5)              
1052 FORMAT(' Data for interstital DOS correspond to atom index ',i4)
1050 FORMAT(F10.5,I3,I5,F8.5,3X,39F8.5)
1051 FORMAT(F10.5,I3,F8.5,3X,39F8.5)
3000 FORMAT(A80,/)                                                       
3010 FORMAT(1X,'LATTICE CONST.=',3F8.4,3X,'FERMI ENERGY=',F10.5)       
3020 FORMAT(I5,' < NMAT <',I5,3X,'SPIN=',I1,3X,'NAT=',I3, 6x,'SO',i2,' KLmax',i3)
3031 FORMAT(1X,'JATOM',I3,2X,'MULT=',I2,60x)                                 
3030 FORMAT(1X,'JATOM',I3,2X,'MULT=',I2,2x,'ISPLIT=',i2,2x,a)              
3040 FORMAT(1X,'BAND:',I4)                                          
800 FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ', 'I N F O R M A T I O N',/,30X,50(1H-),//)
805 FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)
810 FORMAT(3X,'LATTICE',22X,'= ',A4)
820 FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)
840 FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)
850 FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)
1003 FORMAT(A4)
1004 FORMAT(2A5)
2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)
2010 FORMAT(12X,'TOTAL       : ',F8.1)
9000 FORMAT('can''t open definition file ',A40)
9010 FORMAT('can''t open unit: ',I2)
9020 FORMAT('       filename: ',A50)
9030 FORMAT('         status: ',A,'  form: ',A)
9040 FORMAT('Error reading file: ',A47)
1060 FORMAT(//,':NLDA  :',1X,'NUMBER OF ELECTRONS',10X,'=',F8.3,//   ,':FLDA  :',1X,'F E R M I - ENERGY',11X,'= ',F9.5)
5100 FORMAT(A10,4I5,3F5.2,A3)
5101 FORMAT(A10,4I10,3F5.2,A3)

END PROGRAM DMFTMAIN
