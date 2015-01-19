      MODULE TB_MODULE
      USE prec; USE gconstant
      IMPLICIT NONE

! TYPES FOR TB
      TYPE latt
         REAL(gq) :: SCALE
         REAL(gq) :: A(3,3),B(3,3)
         REAL(gq) :: ANORM(3),BNORM(3)
         REAL(gq) :: OMEGA
      END TYPE latt

      TYPE type_info
!only T_INFO
        CHARACTER*40 SZNAM2           ! name of poscar file
        INTEGER NTYP                  ! number of types
        INTEGER NIONS,NIONSL          ! actual number of ions
        LOGICAL LDIRCO                ! positions in direct/recproc. lattice
        REAL(gq), POINTER :: POSION(:,:)  ! positions usually same as DYN%POSION
        INTEGER, POINTER :: ITYP(:)   ! type for each ion
        INTEGER, POINTER :: NITYP(:)  ! number of ions for each type
      END TYPE type_info

      TYPE kpoints_struct
!only  KPOINTS,KPT_BAND,KPT_FS
        INTEGER NKPTS,NKPT(3)        ! actual number of k-points
        REAL(gq) UNIT(3,3)            ! with 2pi/scale
        REAL(gq),POINTER:: VKPT(:,:)  ! coordinate of k-point
        REAL(gq),POINTER:: WTKPT(:)   ! symmetry weight-factor for each k-point
        REAL(gq)        :: VOLWGT     ! volume weight for each tetrahedron
        INTEGER,POINTER:: IDTET(:,:) ! link list for tetrahedron
        INTEGER     :: NTET          ! number of tetrahedrons
        INTEGER,POINTER:: IKPT(:,:,:) ! temporary work array for tetrahedron method
        LOGICAL     :: LTET          ! use tetrahedron method ?
        LOGICAL     :: LUNIT         ! use additional units?
        INTEGER ISMEAR               ! type of smearing
        REAL(gq) SIGMA                ! type of smearing
        INTEGER NEDOS
        REAL(gq) EMIN                 ! minimal E for DOS
        REAL(gq) EMAX                 ! maximal E for DOS
        REAL(gq),POINTER :: DOS(:,:),DOSI(:,:)
        CHARACTER*40  SZNAMK         ! name of k-points file
      END TYPE kpoints_struct

      TYPE wavespin
!only W
        INTEGER NBANDS,ISPIN
        REAL(gq) EFERMI,NET,RSPIN
        COMPLEX(gq),POINTER:: CPTWFP(:,:,:,:) ! wavefunction
        REAL(gq),   POINTER:: FERWE(:,:,:)    ! fermi-weight for each band
        REAL(gq),   POINTER:: CELEN(:,:,:)    ! eigenvalues
      END TYPE wavespin

      TYPE QO_TYPE
        LOGICAL LORTHIN,LORTHOUT
        INTEGER NSP  ! ISPIN FOR TBH/S
        INTEGER DIM,CRNG(2),CDIM ! CORRELATED ORBITAL RANGE
        INTEGER RDIM ! REDUCED DIMENSION DUE TO FROZEN ORBITALS
        INTEGER,POINTER :: IDX(:) ! RE-INDEX OF ORBITALS DUE TO FROZEN ORBITALS
        COMPLEX(gq),POINTER :: HAO(:,:,:),PREROT(:,:)
        COMPLEX(gq),POINTER :: SR(:,:,:,:)
        COMPLEX(gq),POINTER :: HR(:,:,:,:)
        COMPLEX(gq),POINTER :: SK(:,:,:,:)
        COMPLEX(gq),POINTER :: HK(:,:,:,:)
      END TYPE QO_TYPE

      TYPE CELL_TYPE
        INTEGER DIM,ORIG
        INTEGER,ALLOCATABLE :: R(:,:),NDEG(:)
      END TYPE CELL_TYPE

      TYPE TB_YZ
        TYPE(QO_TYPE) QO
      END TYPE TB_YZ

! VARIABLES
      TYPE (latt)          ,SAVE :: LATT_CUR
      TYPE (type_info)     ,SAVE :: T_INFO
      TYPE (kpoints_struct),SAVE :: KPOINTS,KPT_BAND,KPT_FS
      TYPE (wavespin)      ,SAVE :: W
      TYPE (CELL_TYPE)     ,SAVE :: CELL
      TYPE (QO_TYPE)       ,SAVE :: QB
      TYPE (TB_YZ),ALLOCATABLE,SAVE :: TYZ(:)

      CONTAINS

! PUBLIC SUBROUTINES
!***********************************************************************
      SUBROUTINE INI_TB(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO

      CALL RD_POSCAR(IU,IO)
      CALL RD_KPOINTS(KPOINTS,'KPOINTS',IU,IO)
      ALLOCATE(TYZ(T_INFO%NIONS))
      CALL READ_TBHS(IU,IO)
      CALL ALLOC_gqB()
      CALL READ_HAO(IU)
      CALL TB_PAR_ANA_TRANS(IO)
      CALL READ_PREROT(IU,IO)
      CALL ROT_TBHS()
      CALL SET_TBHSK()
!      CALL OUT_TBHSK(IU)
      CALL SET_BAND()
      RETURN

      END SUBROUTINE INI_TB

! PRIVATE SUBROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_POSCAR(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      CHARACTER*255  INPLIN,INPWRK
      INTEGER, EXTERNAL :: NITEMS
      CHARACTER*1  CSEL
      INTEGER I,NT,NI,NSCALE
      REAL(gq) SCALEX,SCALEY,SCALEZ

      WRITE(IO,'(" READING POSCAR")')
      OPEN(IU,FILE='POSCAR',STATUS='OLD')
      READ(IU,'(A40)') T_INFO%SZNAM2
      WRITE(IO,*)T_INFO%SZNAM2
! one scaling parameter or three
      READ(IU,'(A)') INPLIN
      NSCALE=NITEMS(INPLIN,INPWRK,.TRUE.,'F')
      IF (NSCALE==1) THEN
        READ(INPLIN,*) LATT_CUR%SCALE
        SCALEX=1; SCALEY=1; SCALEZ=1
      ELSE IF (NSCALE==3) THEN
        LATT_CUR%SCALE=1
        READ(INPLIN,*) SCALEX,SCALEY,SCALEZ
      ELSE
        WRITE(0,*)'ERROR: there must be 1 or 3 items on line 2 of POSCAR'
        STOP
      ENDIF
      DO I=1,3
        READ(IU,*) LATT_CUR%A(1,I),LATT_CUR%A(2,I),LATT_CUR%A(3,I)
      ENDDO
      IF (LATT_CUR%SCALE<0._gq) THEN
!----alternatively give a volume (=abs(scale)) and adjust the lengths of
!----the three lattice vectors to get the correct desired volume ... :
         CALL LATTIC(LATT_CUR)
         LATT_CUR%SCALE=(ABS(LATT_CUR%SCALE)  &
     &                 / ABS(LATT_CUR%OMEGA))**(1._gq/3._gq)
      ENDIF

      LATT_CUR%A(1,:) =LATT_CUR%A(1,:)*SCALEX*LATT_CUR%SCALE
      LATT_CUR%A(2,:) =LATT_CUR%A(2,:)*SCALEY*LATT_CUR%SCALE
      LATT_CUR%A(3,:) =LATT_CUR%A(3,:)*SCALEZ*LATT_CUR%SCALE

      CALL LATTIC(LATT_CUR)
      WRITE(IO,'(" LATT_CUR%SCALE=",F12.4)')LATT_CUR%SCALE
      WRITE(IO,'(" LATT_CUR%A:")')
      WRITE(IO,'(3F12.4)')LATT_CUR%A
      WRITE(IO,'(" LATT_CUR%B:")')
      WRITE(IO,'(3F12.4)')LATT_CUR%B

      IF (LATT_CUR%OMEGA<0) THEN
        WRITE(0,*)'ERROR: the triple product of the basis vectors ', &
     &     'is negative exchange two basis vectors'
        STOP
      ENDIF

! we are mainly interested in this (6th) line ...
      READ(IU,'(A)') INPLIN
! how many words/data items? --> number of ion types on file POSCAR!
      T_INFO%NTYP=NITEMS(INPLIN,INPWRK,.TRUE.,'I')
      ALLOCATE(T_INFO%NITYP(T_INFO%NTYP))
!-----number of atoms per type
      READ(INPLIN,*) (T_INFO%NITYP(NT),NT=1,T_INFO%NTYP)
! how many ions do we have on file POSCAR ... ?
      T_INFO%NIONS=0
      DO NI=1,T_INFO%NTYP
         T_INFO%NIONS=T_INFO%NIONS+T_INFO%NITYP(NI)
      END DO
      ALLOCATE(T_INFO%ITYP(T_INFO%NIONS))
!---- Set up the table from which we get type of each ion
      NI=1
      DO NT=1,T_INFO%NTYP
      DO NI=NI,T_INFO%NITYP(NT)+NI-1
        T_INFO%ITYP(NI)=NT
      ENDDO
      ENDDO
! posion
      READ(IU,'(A1)') CSEL
      IF (CSEL=='K'.OR.CSEL=='k'.OR. &
     &    CSEL=='C'.OR.CSEL=='c') THEN
        CSEL='K'
        WRITE(IO,*)'Positions in cartesian coordinates'
        T_INFO%LDIRCO=.FALSE.
      ELSE
        WRITE(IO,*)'Positions in direct lattice'
        T_INFO%LDIRCO=.TRUE.
      ENDIF
      ALLOCATE(T_INFO%POSION(3,T_INFO%NIONS)); T_INFO%POSION=0
      DO NI=1,T_INFO%NIONS
        READ(IU,*,ERR=400,END=400) T_INFO%POSION(:,NI)
      ENDDO

      IF (CSEL=='K') THEN
        T_INFO%POSION(1,:)=LATT_CUR%SCALE*T_INFO%POSION(1,:)*SCALEX
        T_INFO%POSION(2,:)=LATT_CUR%SCALE*T_INFO%POSION(2,:)*SCALEY
        T_INFO%POSION(3,:)=LATT_CUR%SCALE*T_INFO%POSION(3,:)*SCALEZ
        CALL KARDIR(T_INFO%NIONS,T_INFO%POSION,LATT_CUR%B)
      ENDIF
      CALL TOPRIM(T_INFO%NIONS,T_INFO%POSION)

      CLOSE(IU)
      RETURN

 400 CONTINUE
      WRITE(0,*)' No initial positions read in'
      STOP

      END SUBROUTINE RD_POSCAR

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_KPOINTS(KPOINTS,FNAME,IU,IO)
      USE prec
      IMPLICIT NONE
      TYPE (kpoints_struct) KPOINTS
      INTEGER IO,IU
      CHARACTER(*) FNAME
      INTEGER NKP,INDX,NINTER,NKPX,NKPY,NKPZ,I,IX,IY,IZ,IRED
      INTEGER,PARAMETER :: NKDIM=100
      CHARACTER*1   CSEL,CLINE
      REAL(gq) SHIFT(3),NKPIF(2,3),X(3)
! required for reallocation
      INTEGER IERR,N
      LOGICAL LTRS
      REAL(gq),ALLOCATABLE   :: VKPT2(:,:)

      LTRS=.FALSE.
      KPOINTS%NKPT=0
      WRITE(IO,'(" READING KPOINTS")')
      OPEN(IU,FILE=FNAME,STATUS='OLD')
      READ(IU,'(A40)') KPOINTS%SZNAMK
      WRITE(IO,*) KPOINTS%SZNAMK
      READ(IU,*) NINTER
      READ(IU,'(A1)') CSEL
      IF(CSEL=='U'.OR.CSEL=='u')THEN
        WRITE(IO,'(" ADDITIONAL UNIT USED FOR GENERATING KPOINTS:")')
        KPOINTS%LUNIT=.TRUE.
        READ(IU,*)KPOINTS%UNIT
        WRITE(IO,'(5X,3F10.4)')KPOINTS%UNIT
        READ(IU,'(A1)') CSEL
      ELSE
        KPOINTS%LUNIT=.FALSE.
        KPOINTS%UNIT=LATT_CUR%B
      ENDIF
      IF (CSEL=='L'.OR.CSEL=='l') THEN
         CLINE='L'
         READ(IU,'(A1)') CSEL
         WRITE(IO,*)'Interpolating k-points between supplied coordinates'
      ELSE
         CLINE=" "
      ENDIF
      IF (CSEL=='K'.OR.CSEL=='k'.OR. &
     &    CSEL=='C'.OR.CSEL=='c') THEN
        CSEL='K'
        WRITE(IO,*)'k-points in cartesian coordinates'
      ELSE
        WRITE(IO,*)'k-points in reciprocal lattice'
      ENDIF
!=======================================================================
! read in a set of k-points and interpolate NKPTS between each
!=======================================================================
      IF (NINTER>0) THEN
      kr: IF (CLINE=='L') THEN
        ALLOCATE(VKPT2(3,NKDIM)); VKPT2=0
        NKP=0  ! counter for the number of k-points already read in
        DO
          NKP=NKP+1
          IF (NKP>NKDIM) THEN
            WRITE(0,*)'ERROR: MAIN: increase NKDIM'
            STOP
          ENDIF
          READ(IU,*,IOSTAT=IERR) VKPT2(:,NKP)
          IF (IERR/=0) EXIT
        ENDDO
        NKP=NKP-1
        KPOINTS%NKPTS=(NKP/2)*NINTER
        ALLOCATE(KPOINTS%VKPT(3,KPOINTS%NKPTS),KPOINTS%WTKPT(KPOINTS%NKPTS))
        IF (CSEL=='K') THEN
          VKPT2(:,1:NKP)=VKPT2(:,1:NKP)/LATT_CUR%SCALE
          CALL KARDIR(NKP,VKPT2,LATT_CUR%A)
        ELSEIF(KPOINTS%LUNIT)THEN
          DO I=1,NKP; VKPT2(:,I)=MATMUL(KPOINTS%UNIT,VKPT2(:,I)); ENDDO
          VKPT2(:,1:NKP)=VKPT2(:,1:NKP)/LATT_CUR%SCALE
          CALL KARDIR(NKP,VKPT2,LATT_CUR%A)
        ENDIF

        INDX=0
        ! make NKPTS even
        NKP=(NKP/2)*2
        DO NKP=1,NKP-1,2
          SHIFT=(VKPT2(:,NKP+1)-VKPT2(:,NKP))/(NINTER-1)
          DO N=0,NINTER-1
          INDX=INDX+1
          KPOINTS%VKPT(:,INDX)=VKPT2(:,NKP)+SHIFT*N
          KPOINTS%WTKPT(INDX)=1._gq/KPOINTS%NKPTS
          ENDDO
        ENDDO
      ELSE kr
        WRITE(0,'(" ERROR: INITIAL NKPTS>0 while not line mode! Not supported!")')
        STOP
      ENDIF kr
      ELSE
!=======================================================================
! Automatic generation of a mesh if NINTER<=0:
!=======================================================================
        WRITE(IO,*) 'Automatic generation of k-mesh.'
! k-lattice basis vectors in cartesian or reciprocal coordinates?
! Always use gamma center k-mesh.
        IF ((CSEL=='G').OR.(CSEL=='g')) THEN
! G: time reversal symmetry; g: no time reversal symmetry
! Here we give the Monkhorst-Pack conventions ... :
          READ(IU,*) NKPX,NKPY,NKPZ
          IF(FNAME.NE.'KPOINTS_FS')THEN
! For simplicity, adjust NKPX/Y/Z to odd
          NKPIF(:,1)=NKPX/2; NKPIF(:,2)=NKPY/2; NKPIF(:,3)=NKPZ/2
          NKPX=NKPX/2*2+1  ; NKPY=NKPY/2*2+1  ; NKPZ=NKPZ/2*2+1
          KPOINTS%NKPT(1)=NKPX; KPOINTS%NKPT(2)=NKPY; KPOINTS%NKPT(3)=NKPZ
          WRITE(IO,'(" Adjusted k-mesh:",3I3)')NKPX,NKPY,NKPZ
! Set up k-points
          LTRS=CSEL=='G'
          IRED=1
          IF(CSEL=='G')THEN
          DO I=3,1,-1
            IF(NKPIF(1,I).EQ.0)CYCLE
            NKPIF(1,I)=0; IRED=I; EXIT
          ENDDO
          ENDIF
          ELSE ! KPOINTS_FS
          KPOINTS%NKPT(1)=NKPX; KPOINTS%NKPT(2)=NKPY; KPOINTS%NKPT(3)=NKPZ
          NKPIF(1,:)=0; NKPIF(2,:)=KPOINTS%NKPT(:)
          ENDIF

          KPOINTS%NKPTS=1
          DO I=1,3
            KPOINTS%NKPTS=KPOINTS%NKPTS*(NKPIF(2,I)+NKPIF(1,I)+1)
          ENDDO
          WRITE(IO,'(" Total number of k-points:",I4)')KPOINTS%NKPTS
          ALLOCATE(KPOINTS%VKPT(3,KPOINTS%NKPTS),KPOINTS%WTKPT(KPOINTS%NKPTS))
          NKP=0
! FOR XCRYSDEN
          DO IX=-NKPIF(1,1),NKPIF(2,1); X(1)=REAL(IX,KIND=gq)
          DO IY=-NKPIF(1,2),NKPIF(2,2); X(2)=REAL(IY,KIND=gq)
          DO IZ=-NKPIF(1,3),NKPIF(2,3); X(3)=REAL(IZ,KIND=gq)
            NKP=NKP+1
            KPOINTS%VKPT(1,NKP)=X(1)/NKPX
            KPOINTS%VKPT(2,NKP)=X(2)/NKPY
            KPOINTS%VKPT(3,NKP)=X(3)/NKPZ
            IF(X(IRED).EQ.0._gq.OR.CSEL=='g')THEN
              KPOINTS%WTKPT(NKP)=1._gq/(NKPX*NKPY*NKPZ)
            ELSE
              KPOINTS%WTKPT(NKP)=2._gq/(NKPX*NKPY*NKPZ)
            ENDIF
          ENDDO; ENDDO; ENDDO
          IF(KPOINTS%LUNIT)THEN
            DO I=1,KPOINTS%NKPTS; KPOINTS%VKPT(:,I)=MATMUL(KPOINTS%UNIT,KPOINTS%VKPT(:,I)); ENDDO
            KPOINTS%VKPT(:,1:NKP)=KPOINTS%VKPT(:,1:NKP)/LATT_CUR%SCALE
            CALL KARDIR(NKP,KPOINTS%VKPT,LATT_CUR%A)
          ENDIF
        ELSE
          WRITE(IO,'(" CSEL \= G,g: Not supported!")')
          STOP
        ENDIF
      ENDIF
! Determine LTET
      IF(FNAME.EQ.'KPOINTS'.AND.((KPOINTS%NKPT(1).GT.1.AND.KPOINTS%NKPT(2).GT.1).OR. &
                                &(KPOINTS%NKPT(1).GT.1.AND.KPOINTS%NKPT(3).GT.1).OR. &
                                &(KPOINTS%NKPT(2).GT.1.AND.KPOINTS%NKPT(3).GT.1)))THEN
!        KPOINTS%LTET=.TRUE.
!        KPOINTS%ISMEAR=-5
        KPOINTS%LTET=.FALSE.
        KPOINTS%ISMEAR=0
        KPOINTS%SIGMA =0.01_gq
      ELSE
        KPOINTS%LTET=.FALSE.
        KPOINTS%ISMEAR=0
        KPOINTS%SIGMA =0.01_gq
      ENDIF
      CALL SET_UP_TET(KPOINTS,LTRS,IO)
      CLOSE(IU)
      WRITE(IO,'(" KPOINTS (50,Direct):")')
      WRITE(IO,'(" NK   K    WT")')
      DO NKP=1,MIN(KPOINTS%NKPTS,5)
        WRITE(IO,'(I4,3F8.3,F10.5)')NKP,KPOINTS%VKPT(:,NKP),KPOINTS%WTKPT(NKP)
      ENDDO
      WRITE(IO,'(" KPOINTS (50,Cartesian):")')
      WRITE(IO,'(" NK   K    WT")')
      DO NKP=1,MIN(KPOINTS%NKPTS,5)
        WRITE(IO,'(I4,3F8.3,F10.5)')NKP,MATMUL(LATT_CUR%B,KPOINTS%VKPT(:,NKP)),KPOINTS%WTKPT(NKP)
      ENDDO
      IF (CLINE=='L') THEN
      WRITE(IO,'(" LINE MODE: INPUT SPECIAL KPTS IN DIR:")')
      WRITE(IO,'(I4,3F8.3,F10.5)')1,KPOINTS%VKPT(:,1)
      DO NKP=NINTER,KPOINTS%NKPTS,NINTER
        WRITE(IO,'(I4,3F8.3,F10.5)')NKP,KPOINTS%VKPT(:,NKP)
      ENDDO
      WRITE(IO,'(" LINE MODE: INPUT SPECIAL KPTS IN KAR:")')
      WRITE(IO,'(I4,3F8.3,F10.5)')1,MATMUL(LATT_CUR%B,KPOINTS%VKPT(:,1))
      DO NKP=NINTER,KPOINTS%NKPTS,NINTER
        WRITE(IO,'(I4,3F8.3,F10.5)')NKP,MATMUL(LATT_CUR%B,KPOINTS%VKPT(:,NKP))
      ENDDO
      ENDIF
      RETURN

      END SUBROUTINE RD_KPOINTS

!**************** SUBROUTINE SET_UP_TET  ***********************************
      SUBROUTINE SET_UP_TET(KPOINTS,LTRS,IO)
      USE prec
      IMPLICIT NONE
      TYPE (kpoints_struct) KPOINTS
      LOGICAL LTRS ! Time-reversal symmetry
      INTEGER IO
! LOCAL
      REAL(gq) V(3)
      INTEGER IKPTD,I1,I2,I3,NKX,NKY,NKZ,IK
      INTEGER,parameter :: NTETD=9000000
      REAL(gq),PARAMETER :: TINY=1.E-6_gq

      IF(.NOT.KPOINTS%LTET)RETURN
      IKPTD=MAXVAL(KPOINTS%NKPT)
      NKX=KPOINTS%NKPT(1); NKY=KPOINTS%NKPT(2); NKZ=KPOINTS%NKPT(3)
      ALLOCATE(KPOINTS%IKPT(IKPTD,IKPTD,IKPTD),KPOINTS%IDTET(0:4,NTETD))
      KPOINTS%IKPT=0; KPOINTS%IDTET=0
! Go through the mesh and get the 'connection table'
      DO I3=0,NKZ-1; V(3)=REAL(I3-INT(NKZ/2),KIND=gq)/NKZ
      DO I2=0,NKY-1; V(2)=REAL(I2-INT(NKY/2),KIND=gq)/NKY
      DO I1=0,NKX-1; V(1)=REAL(I1-INT(NKX/2),KIND=gq)/NKX
      DO IK=1,KPOINTS%NKPTS
        IF(ABS(V(1)-KPOINTS%VKPT(1,IK)).GT.TINY)CYCLE
        IF(ABS(V(2)-KPOINTS%VKPT(2,IK)).GT.TINY)CYCLE
        IF(ABS(V(3)-KPOINTS%VKPT(3,IK)).GT.TINY)CYCLE
        KPOINTS%IKPT(I1+1,I2+1,I3+1)=IK
        EXIT
      ENDDO
      IF(KPOINTS%IKPT(I1+1,I2+1,I3+1).EQ.0.AND.LTRS)THEN
      DO IK=1,KPOINTS%NKPTS
        IF(ABS(V(1)+KPOINTS%VKPT(1,IK)).GT.TINY)CYCLE
        IF(ABS(V(2)+KPOINTS%VKPT(2,IK)).GT.TINY)CYCLE
        IF(ABS(V(3)+KPOINTS%VKPT(3,IK)).GT.TINY)CYCLE
        KPOINTS%IKPT(I1+1,I2+1,I3+1)=IK
        EXIT
      ENDDO
      ENDIF
      IF(KPOINTS%IKPT(I1+1,I2+1,I3+1).EQ.0)THEN
        WRITE(0,'(" CAN NOT SET CONNECTION FOR I:",3I4)')I1,I2,I3; STOP
      ENDIF
      ENDDO; ENDDO; ENDDO
      CALL TETIRR(KPOINTS%NTET,KPOINTS%IDTET,NTETD,LATT_CUR%B,NKX,NKY,NKZ,KPOINTS%IKPT,IKPTD,IO)
      KPOINTS%VOLWGT=1._gq/REAL(6*NKX*NKY*NKZ,KIND=gq)
      DEALLOCATE(KPOINTS%IKPT)
      RETURN

      END SUBROUTINE SET_UP_TET

!**************** SUBROUTINE LATTIC  ***********************************
!  subroutine for calculating the reciprocal lattice from the direct
!  lattice, in addition the norm of the lattice-vectors and the volume of
!  the basis-cell is calculated
!***********************************************************************
      SUBROUTINE LATTIC(Mylatt)
      USE prec
      IMPLICIT NONE

      TYPE(LATT) Mylatt
      REAL(gq) Omega
      INTEGER I,J
      INTRINSIC SUM

      CALL EXPRO(Mylatt%B(1:3,1),Mylatt%A(1:3,2),Mylatt%A(1:3,3))
      CALL EXPRO(Mylatt%B(1:3,2),Mylatt%A(1:3,3),Mylatt%A(1:3,1))
      CALL EXPRO(Mylatt%B(1:3,3),Mylatt%A(1:3,1),Mylatt%A(1:3,2))

      Omega =Mylatt%B(1,1)*Mylatt%A(1,1)+Mylatt%B(2,1)*Mylatt%A(2,1) &
     &      +Mylatt%B(3,1)*Mylatt%A(3,1)

      DO I=1,3; DO J=1,3
        Mylatt%B(I,J)=Mylatt%B(I,J)/Omega
      ENDDO; ENDDO

      DO I=1,3
        Mylatt%ANORM(I)=SQRT(SUM(Mylatt%A(:,I)*Mylatt%A(:,I)))
        Mylatt%BNORM(I)=SQRT(SUM(Mylatt%B(:,I)*Mylatt%B(:,I)))
      ENDDO
      Mylatt%Omega=Omega
      RETURN
      END SUBROUTINE LATTIC

!***********************************************************************
      SUBROUTINE READ_TBHS(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER IC,IQ1,IQ2,ISP,ILINE
      INTEGER NQ,FNQ,NC,ISPIN,NSP,NLINE,ITMP ! FNQ: frozen orbitalis
      COMPLEX(gq) ZES
      INTEGER IOS
      CHARACTER(3) STMP
      LOGICAL LORTHIN,LORTHOUT

      IOS=0
      OPEN(IU,FILE='TBHS.INP',STATUS='OLD')
      READ(IU,*,ERR=999)ISPIN,FNQ,LORTHIN,LORTHOUT
      WRITE(IO,'(" ISPIN,FNQ,LORTHIN,LORTHOUT:")')
      WRITE(IO,'(I2,I3,2L2,F8.4,I3)')ISPIN,FNQ,LORTHIN,LORTHOUT
      READ(IU,*); READ(IU,*)W%NET
      WRITE(IO,'(" NET=",F8.2)')W%NET
      IF(LORTHIN.AND.(.NOT.LORTHOUT))STOP'CONFLIT SETTING OF LORTHIN/LORTHOUT!'
      READ(IU,*); READ(IU,*)NC
      CELL%DIM=ABS(NC)
      WRITE(IO,'(" CELL%DIM=",I8)')NC
      ALLOCATE(CELL%NDEG(CELL%DIM),CELL%R(3,CELL%DIM))
      CELL%NDEG=1
      DO IC=1,CELL%DIM
        IF(NC.GT.0)THEN
          READ(IU,*)CELL%R(:,IC) ! RCUT type
        ELSE
          READ(IU,*)CELL%R(:,IC),CELL%NDEG(IC) ! Wigner_Seitz type
        ENDIF
      ENDDO
      READ(IU,*); READ(IU,*)NQ,NSP
      WRITE(IO,'(" NQ,NSP=",I5,I2)')NQ,NSP
      READ(IU,*); READ(IU,*)(ITMP,IC=1,NQ); READ(IU,*)(STMP,IC=1,NQ)
      ALLOCATE(QB%HR(NQ,NQ,CELL%DIM,NSP)); QB%HR=0
      READ(IU,*)
      ILINE=0
      DO
        READ(IU,*,ERR=102,END=102)ISP,IQ1,IQ2,IC,ZES
        IF(QB%HR(IQ1,IQ2,IC,ISP).NE.0)THEN
        WRITE(0,'(" ISP,IQ1,IQ2,IC",3I3,I6," ENTRY REPEATED IN TBHC.INP!")')ISP,IQ1,IQ2,IC
        STOP
        ENDIF
        QB%HR(IQ1,IQ2,IC,ISP)=ZES
        IF(IQ1/=IQ2)QB%HR(IQ2,IQ1,IC,ISP)=CONJG(ZES)
        ILINE=ILINE+1
      ENDDO
102   CONTINUE; BACKSPACE(IU)
      WRITE(IO,'(" TBH%LINES=",I7)')ILINE
      ILINE=0
      IF(.NOT.LORTHIN)THEN
      ALLOCATE(QB%SR(NQ,NQ,CELL%DIM,NSP)); QB%SR=0
      READ(IU,*,IOSTAT=IOS,ERR=101,END=101)
      DO
        READ(IU,*,ERR=101,END=101)ISP,IQ1,IQ2,IC,ZES
        IF(QB%SR(IQ1,IQ2,IC,ISP).NE.0)THEN
        WRITE(0,'(" ISP,IQ1,IQ2,IC",3I3,I6," ENTRY REPEATED IN TBHS.INP!")')ISP,IQ1,IQ2,IC
        STOP
        ENDIF
        QB%SR(IQ1,IQ2,IC,ISP)=ZES
        ILINE=ILINE+1
      ENDDO
      ENDIF
101   CONTINUE
      WRITE(IO,'(" TBS%LINES=",I7)')ILINE
      CLOSE(IU)

      DO IC=1,CELL%DIM
        IF(CELL%R(1,IC).EQ.0.AND.CELL%R(2,IC).EQ.0.AND.CELL%R(3,IC).EQ.0)THEN
          CELL%ORIG=IC
          GOTO 100
        ENDIF
      ENDDO
      STOP 'NO ORIGIN IN CELL!'
100   CONTINUE
      IF(IOS.NE.0)THEN
        DO IQ1=1,NQ; QB%SR(IQ1,IQ1,CELL%ORIG,:)=1._gq; ENDDO ! Set SR unit.
      ENDIF
      QB%DIM=NQ; QB%RDIM=NQ-FNQ; QB%NSP=NSP
      W%ISPIN=ISPIN; W%RSPIN=3-W%ISPIN
      QB%LORTHIN=LORTHIN; QB%LORTHOUT=LORTHOUT
      RETURN
999   WRITE(0,'(" NEED TO ADD ADDITIONAL CTRL PARAMETERS: ISPIN,FNQ,LORTHIN,LORTHOUT!")'); STOP

      END SUBROUTINE READ_TBHS

!***********************************************************************
      SUBROUTINE ALLOC_gqB()
      USE prec
      IMPLICIT NONE
! LOCAL
      INTEGER NQ,IQ

      NQ=QB%DIM
      ALLOCATE(QB%HK(NQ,NQ,KPOINTS%NKPTS,QB%NSP),QB%HAO(NQ,NQ,QB%NSP), &
              &QB%IDX(NQ),QB%PREROT(NQ,NQ))
      QB%HK=0; QB%HAO=0; QB%IDX=0; QB%PREROT=0
      DO IQ=1,NQ
        QB%HAO(IQ,IQ,:)=1._gq; QB%PREROT(IQ,IQ)=1._gq
        QB%IDX(IQ)=IQ
      ENDDO
      IF(.NOT.QB%LORTHOUT)THEN
      ALLOCATE(QB%SK(NQ,NQ,KPOINTS%NKPTS,QB%NSP)); QB%SK=0
      ENDIF
      RETURN

      END SUBROUTINE ALLOC_gqB

!***********************************************************************
      SUBROUTINE READ_HAO(IU)
      USE prec
      IMPLICIT NONE
      INTEGER IU
! LOCAL
      INTEGER NI
      INTEGER CRNG(2),DIM,NQBASE

      OPEN(IU,FILE='HAO.INP',STATUS='OLD')
      NQBASE=0
      T_INFO%NIONSL=0
      DO NI=1,T_INFO%NIONS
      READ(IU,*)
      READ(IU,*)DIM,CRNG
      IF(CRNG(2).GE.CRNG(1))THEN
        T_INFO%NIONSL=T_INFO%NIONSL+1
        TYZ(NI)%QO%HAO=>QB%HAO(NQBASE+CRNG(1):NQBASE+CRNG(2),NQBASE+CRNG(1):NQBASE+CRNG(2),:)
        READ(IU,*)TYZ(NI)%QO%HAO
      ENDIF
      TYZ(NI)%QO%DIM =DIM; TYZ(NI)%QO%CRNG=CRNG; TYZ(NI)%QO%CDIM=CRNG(2)-CRNG(1)+1
      NQBASE=NQBASE+DIM
      ENDDO ! NI
      CLOSE(IU)
      RETURN

      END SUBROUTINE READ_HAO

!***********************************************************************
      SUBROUTINE READ_PREROT(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER NI
      INTEGER CRNG(2),DIM,NQBASE

      OPEN(IU,FILE='PREROT.INP',STATUS='OLD',ERR=100)
      NQBASE=0
      DO NI=1,T_INFO%NIONS
      READ(IU,*)
      READ(IU,*)DIM
      IF(DIM.EQ.TYZ(NI)%QO%DIM)THEN
        TYZ(NI)%QO%PREROT=>QB%PREROT(NQBASE+1:NQBASE+DIM,NQBASE+1:NQBASE+DIM)
        READ(IU,*)TYZ(NI)%QO%PREROT
      ENDIF
      NQBASE=NQBASE+TYZ(NI)%QO%DIM
      ENDDO ! NI
      CLOSE(IU)
      WRITE(IO,'(" PREROT.INP READ IN.")')
100   RETURN

      END SUBROUTINE READ_PREROT

!***********************************************************************
      SUBROUTINE ROT_TBHS()
      USE prec
      IMPLICIT NONE
! LOCAL
      INTEGER IC,ISP

      DO ISP=1,QB%NSP; DO IC=1,CELL%DIM
        QB%HR(:,:,IC,ISP)=MATMUL(CONJG(TRANSPOSE(QB%HAO(:,:,ISP))), &
                         &MATMUL(CONJG(TRANSPOSE(QB%PREROT(:,:))), &
                         &MATMUL(QB%HR(:,:,IC,ISP), &
                         &MATMUL(QB%PREROT(:,:),QB%HAO(:,:,ISP)))))
        IF(.NOT.QB%LORTHIN)&
       &QB%SR(:,:,IC,ISP)=MATMUL(CONJG(TRANSPOSE(QB%HAO(:,:,ISP))), &
                         &MATMUL(CONJG(TRANSPOSE(QB%PREROT(:,:))), &
                         &MATMUL(QB%SR(:,:,IC,ISP), &
                         &MATMUL(QB%PREROT(:,:),QB%HAO(:,:,ISP)))))
      ENDDO; ENDDO
      RETURN

      END SUBROUTINE ROT_TBHS

!***********************************************************************
      SUBROUTINE SET_BAND()
      USE prec
      IMPLICIT NONE

      W%NBANDS=QB%RDIM
! WAVE
      ALLOCATE(W%CPTWFP(QB%RDIM,W%NBANDS,KPOINTS%NKPTS,W%ISPIN), &
               W%CELEN (        W%NBANDS,KPOINTS%NKPTS,W%ISPIN), &
               W%FERWE (        W%NBANDS,KPOINTS%NKPTS,W%ISPIN))
      W%CPTWFP=0; W%CELEN=0; W%FERWE=0
! DOS
      KPOINTS%NEDOS=501
      KPOINTS%EMIN= 100._gq
      KPOINTS%EMAX=-100._gq
      ALLOCATE(KPOINTS%DOS(KPOINTS%NEDOS,W%ISPIN),KPOINTS%DOSI(KPOINTS%NEDOS,W%ISPIN))
      KPOINTS%DOS=0; KPOINTS%DOSI=0
      RETURN

      END SUBROUTINE SET_BAND

!***********************************************************************
      SUBROUTINE SET_TBHSK()
      USE prec; USE gconstant
      IMPLICIT NONE
! LOCAL
      INTEGER ISP,IK
      COMPLEX(gq) SK(QB%DIM,QB%DIM)

      DO IK=1,KPOINTS%NKPTS
      DO ISP=1,QB%NSP
        CALL SET_HS1K(KPOINTS%VKPT(:,IK),QB%HK(:,:,IK,ISP),QB%HR(:,:,:,ISP))
        IF(.NOT.QB%LORTHIN)THEN
        CALL SET_HS1K(KPOINTS%VKPT(:,IK),   SK            ,QB%SR(:,:,:,ISP))
        IF(QB%LORTHOUT)THEN
          CALL ORTH_HS1K(QB%HK(:,:,IK,ISP),SK)
        ELSE
          QB%SK(:,:,IK,ISP)=SK
        ENDIF
        ENDIF
      ENDDO; ENDDO
      RETURN

      END SUBROUTINE SET_TBHSK

!***********************************************************************
      SUBROUTINE OUT_TBHSK(IU)
      USE prec
      IMPLICIT NONE
      INTEGER IU
! LOCAL
      INTEGER ISP,IK,I1,I2

      OPEN(IU,FILE='TBHK.dat',STATUS='REPLACE')
      WRITE(IU,*)KPOINTS%NKPTS
      WRITE(IU,*)QB%DIM
      DO ISP=1,QB%NSP; DO IK=1,KPOINTS%NKPTS
      DO I1=1,QB%DIM; DO I2=1,QB%DIM
        WRITE(IU,*)QB%HK(I1,I2,IK,ISP)
      ENDDO; ENDDO
      ENDDO; ENDDO
      CLOSE(IU)
      RETURN

      END SUBROUTINE OUT_TBHSK

!***********************************************************************
      SUBROUTINE SET_HS1K(VKPT,AK,AR)
      USE prec; USE gconstant
      IMPLICIT NONE
      REAL(gq) VKPT(3)
      COMPLEX(gq) AR(QB%DIM,QB%DIM,CELL%DIM),AK(QB%DIM,QB%DIM)
! LOCAL
      INTEGER IC,IQ1,IQ2
      COMPLEX(gq) PHASEK

      AK=0
      DO IC=1,CELL%DIM
      PHASEK=EXP( CITPI*SUM(CELL%R(:,IC)*VKPT) ) ! EXP(-iK.R)
      DO IQ1=1,QB%DIM; DO IQ2=IQ1,QB%DIM
        AK(IQ1,IQ2)=AK(IQ1,IQ2)+AR(IQ1,IQ2,IC)*PHASEK/CELL%NDEG(IC)
      ENDDO; ENDDO; ENDDO
      DO IQ1=1,QB%DIM; DO IQ2=1,IQ1-1
        AK(IQ1,IQ2)=CONJG(AK(IQ2,IQ1))
      ENDDO; ENDDO
      RETURN

      END SUBROUTINE SET_HS1K

!***********************************************************************
      SUBROUTINE CALC_FULL_BAND()
      USE prec
      IMPLICIT NONE
! LOCAL
      INTEGER ISP,K
      REAL(gq) KPT(3)

      DO ISP=1,W%ISPIN
      DO K=1,KPOINTS%NKPTS
        KPT=KPOINTS%VKPT(:,K)
        CALL EIGEN_SOL(KPT,K,W%CELEN(:,K,ISP),W%CPTWFP(:,:,K,ISP),QB%RDIM,'V',ISP)
      ENDDO; ENDDO ! ISP,K
      RETURN

      END SUBROUTINE CALC_FULL_BAND

!***********************************************************************
      SUBROUTINE EIGEN_SOL(KPT,K,W,A,N,JOBZ,ISP,PJ)
      USE prec
      IMPLICIT NONE
      INTEGER K,N,ISP
      REAL(gq) W(N),KPT(3)
      COMPLEX(gq) A(N,N)
      CHARACTER*1 JOBZ
      REAL(gq),OPTIONAL,INTENT(OUT)::PJ(N,N)
! Local
      INTEGER LWORK,INFO,IV,IVP
      REAL(gq) RWORK(3*N-2)
      COMPLEX(gq),ALLOCATABLE :: B(:,:),WORK(:),BP(:,:),SV(:)

      CALL  GET_EFF_HAMILT(A,N,KPT,K,ISP)
      LWORK=32*N
      ALLOCATE(WORK(LWORK)); WORK=0
      IF(.NOT.QB%LORTHOUT)THEN
      ALLOCATE(B(N,N)); B=0
      CALL  GET_EFF_OVERLAP(B,N,KPT,K,ISP)
      IF(PRESENT(PJ))THEN
      ALLOCATE(BP(N,N),SV(N)); BP=B
      ENDIF
      CALL ZHEGV(1,JOBZ,'L',N,A,N,B,N,W,WORK,LWORK,RWORK,INFO)
      IF(PRESENT(PJ))THEN
      DO IV=1,N
        SV=MATMUL(BP,A(:,IV))
        DO IVP=1,N
        PJ(IVP,IV)=CONJG(A(IVP,IV))*SV(IVP)
      ENDDO; ENDDO
      ENDIF

      ELSE
      CALL ZHEEV(    JOBZ,'L',N,A,N,    W,WORK,LWORK,RWORK,INFO)
      IF(PRESENT(PJ))THEN
        DO IV=1,N; DO IVP=1,N
        PJ(IVP,IV)=CONJG(A(IVP,IV))*A(IVP,IV)
        ENDDO; ENDDO
      ENDIF
      ENDIF
      DEALLOCATE(WORK)
      IF(.NOT.QB%LORTHOUT)THEN
      DEALLOCATE(B)
      IF(PRESENT(PJ))DEALLOCATE(BP)
      ENDIF
      IF(INFO.NE.0)THEN
        WRITE(0,'(" EIGEN_SOL: INFO=",I4)')INFO
        STOP
      ENDIF
      RETURN

      END SUBROUTINE EIGEN_SOL

!***********************************************************************
      SUBROUTINE GET_EFF_HAMILT(HMAT,NDIM,KPT,K,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(gq) KPT(3)
      COMPLEX(gq) HMAT(NDIM,NDIM)
! Local
      INTEGER N1,N2

      CALL GET_UCTBH(HMAT,NDIM,KPT,K,ISP)
      RETURN

      END SUBROUTINE GET_EFF_HAMILT

!***********************************************************************
      SUBROUTINE GET_UCTBH(HMAT,NDIM,KPT,K,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(gq) KPT(3)
      COMPLEX(gq) HMAT(NDIM,NDIM)
! LOCAL
      INTEGER IQ1,IQ2
      INTEGER ISP_
      COMPLEX(gq),ALLOCATABLE :: HK(:,:),SK(:,:)

      ISP_=MIN(ISP,QB%NSP)
      IF(K.GT.0)THEN
        DO IQ1=1,QB%RDIM; DO IQ2=1,QB%RDIM
          HMAT(IQ1,IQ2)=QB%HK(QB%IDX(IQ1),QB%IDX(IQ2),K,ISP_)
        ENDDO; ENDDO
      ELSE
        ALLOCATE(HK(QB%DIM,QB%DIM)); HK=0
        CALL SET_HS1K(KPT,HK,QB%HR(:,:,:,ISP_))
        IF(.NOT.QB%LORTHIN.AND.QB%LORTHOUT)THEN
        ALLOCATE(SK(QB%DIM,QB%DIM)); SK=0
        CALL SET_HS1K(KPT,SK,QB%SR(:,:,:,ISP_))
        CALL ORTH_HS1K(HK,SK)
        ENDIF
        DO IQ1=1,QB%RDIM; DO IQ2=1,QB%RDIM
          HMAT(IQ1,IQ2)=HK(QB%IDX(IQ1),QB%IDX(IQ2))
        ENDDO; ENDDO
        DEALLOCATE(HK)
      ENDIF
      RETURN

      END SUBROUTINE GET_UCTBH

!***********************************************************************
      SUBROUTINE GET_EFF_OVERLAP(SMAT,NDIM,KPT,K,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(gq) KPT(3)
      COMPLEX(gq) SMAT(NDIM,NDIM)
! Local
      INTEGER N1,N2

      CALL GET_UCTBS(SMAT,NDIM,KPT,K,ISP)
      RETURN

      END SUBROUTINE GET_EFF_OVERLAP

!***********************************************************************
      SUBROUTINE GET_UCTBS(SMAT,NDIM,KPT,K,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(gq) KPT(3)
      COMPLEX(gq) SMAT(NDIM,NDIM)
! LOCAL
      INTEGER IQ1,IQ2
      INTEGER ISP_
      COMPLEX(gq),ALLOCATABLE :: SK(:,:)

      ISP_=MIN(ISP,QB%NSP)
      IF(K.GT.0)THEN
        DO IQ1=1,QB%RDIM; DO IQ2=1,QB%RDIM
          SMAT(IQ1,IQ2)=QB%SK(QB%IDX(IQ1),QB%IDX(IQ2),K,ISP_)
        ENDDO; ENDDO
      ELSE
        ALLOCATE(SK(QB%DIM,QB%DIM)); SK=0
        CALL SET_HS1K(KPT,SK,QB%SR(:,:,:,ISP_))
        DO IQ1=1,QB%RDIM; DO IQ2=1,QB%RDIM
          SMAT(IQ1,IQ2)=SK(QB%IDX(IQ1),QB%IDX(IQ2))
        ENDDO; ENDDO
        DEALLOCATE(SK)
      ENDIF
      RETURN

      END SUBROUTINE GET_UCTBS

!***********************************************************************
      SUBROUTINE SET_FERMI_WT(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER ISP,K,N
      REAL(gq) EADD,E

      IF(KPOINTS%EMIN.GT.KPOINTS%EMAX)THEN
        DO ISP=1,W%ISPIN
        DO K=1,KPOINTS%NKPTS
        DO N=1,W%NBANDS
          E=W%CELEN(N,K,ISP)
          KPOINTS%EMAX=MAX(KPOINTS%EMAX,E)
          KPOINTS%EMIN=MIN(KPOINTS%EMIN,E)
        ENDDO; ENDDO; ENDDO
        EADD=(KPOINTS%EMAX-KPOINTS%EMIN)*0.05_gq
        EADD=MAX(EADD,10*ABS(KPOINTS%SIGMA))
        KPOINTS%EMIN=KPOINTS%EMIN-EADD
        KPOINTS%EMAX=KPOINTS%EMAX+EADD
      ENDIF
! Determine fermi level
      IF(KPOINTS%ISMEAR.GE.0)THEN
! Gaussian ...
        CALL DENMP(W%ISPIN,W%RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,W%NET,W%EFERMI, &
             & KPOINTS%ISMEAR,KPOINTS%SIGMA,W%NBANDS,KPOINTS%NKPTS,W%FERWE,W%CELEN, &
             & KPOINTS%NEDOS,KPOINTS%DOS,KPOINTS%DOSI,KPOINTS%WTKPT)
      ELSEIF(KPOINTS%ISMEAR.GE.-5)THEN
        CALL DENTET(IO,W%CELEN,KPOINTS%WTKPT,W%NBANDS,KPOINTS%NKPTS, &
           KPOINTS%DOS,KPOINTS%DOSI, &
           KPOINTS%NEDOS,W%ISPIN,W%RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET,  &
           KPOINTS%VOLWGT,W%NET,W%EFERMI,W%FERWE, &
           2,IO,KPOINTS%NKPTS)
      ELSE
        WRITE(0,'(" NOT SUPPORTED: ISMEAR=",I2)')KPOINTS%ISMEAR; STOP
      ENDIF
      WRITE(IO,'(" EFERMI=",F8.3)')W%EFERMI
      RETURN

      END SUBROUTINE SET_FERMI_WT

!********************************************************************
! PRINT OUT TDOS
!********************************************************************
      SUBROUTINE OUT_TDOS(IT)
      USE prec
      IMPLICIT NONE
      INTEGER IT,I
      REAL(gq) DELTAE

      DELTAE=(KPOINTS%EMAX-KPOINTS%EMIN)/(KPOINTS%NEDOS-1)
      WRITE(71,'("#TDOS: ITER=",I)')IT
      DO I=1,KPOINTS%NEDOS
        WRITE(71,'(3F)')KPOINTS%EMIN+(I-1)*DELTAE-W%EFERMI,KPOINTS%DOS(I,1),KPOINTS%DOSI(I,1)
      ENDDO
      WRITE(71,*)

      END SUBROUTINE OUT_TDOS

!********************************************************************
      SUBROUTINE TB_PAR_ANA_TRANS(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO

      CALL TB_PAR_DECAY_ANA(QB%HR,QB%DIM,CELL%DIM,QB%NSP,CELL%R,IO)
      RETURN

      END SUBROUTINE TB_PAR_ANA_TRANS

!********************************************************************
      SUBROUTINE ORTH_HS1K(HK,SK)
      USE prec
      IMPLICIT NONE
      COMPLEX(gq) HK(QB%DIM,QB%DIM),SK(QB%DIM,QB%DIM)
! LOCAL
      INTEGER LWORK,IERR,IQ,IQ1,IQ2
      COMPLEX(gq) SNHK(QB%DIM,QB%DIM),SK_(QB%DIM,QB%DIM)
      REAL(gq) RWORK(3*QB%DIM),W(QB%DIM)
      COMPLEX(gq) WORK(32*QB%DIM)

      LWORK=32*QB%DIM; SK_=SK
      CALL ZHEEV('V','U',QB%DIM,SK,QB%DIM,W,WORK,LWORK,RWORK,IERR)
      IF(IERR.NE.0)STOP 'ZHEEV ERROR-I IN RTH_HS1K!'
      IF(MINVAL(W)/MAXVAL(W).LT.1.E-5_gq)THEN ! Condition number
      WRITE(0,*)
      WRITE(0,*)'EIGENVALUES IN ORTH_HS1K:'
      WRITE(0,'(10f8.3)')W
      STOP 'ZHEEV ERROR-II IN RTH_HS1K'
      ENDIF
! Set up S^(-1/2)
      SNHK=0
      DO IQ1=1,QB%DIM; DO IQ2=1,QB%DIM; DO IQ=1,QB%DIM
        SNHK(IQ2,IQ1)=SNHK(IQ2,IQ1)+SK(IQ2,IQ)*CONJG(SK(IQ1,IQ))/SQRT(W(IQ))
      ENDDO; ENDDO; ENDDO
      HK=MATMUL(CONJG(TRANSPOSE(SNHK)),MATMUL(HK ,SNHK))
      SK=MATMUL(CONJG(TRANSPOSE(SNHK)),MATMUL(SK_,SNHK))
      RETURN

      END SUBROUTINE ORTH_HS1K

!********************************************************************
      SUBROUTINE TB_PAR_DECAY_ANA(HR,NQ,NC,NSP,R,IO)
      USE prec
      IMPLICIT NONE
      INTEGER NQ,NC,NSP,IO,R(3,NC)
      COMPLEX(gq) HR(NQ,NQ,NC,NSP)
! LOCAL
      INTEGER I,ISP,IQ1,IQ2,NI1,NI2,L1,L2,IC
      INTEGER,PARAMETER :: NDIV=23
      REAL(gq) RC(NDIV),EC(NDIV),DIST,HR1

      EC=0
      DO I=1,NDIV; RC(I)=I+1; ENDDO
      DO ISP=1,NSP
      IQ1=0
      DO NI1=1,T_INFO%NIONS; DO L1=1,TYZ(NI1)%QO%DIM
      IQ1=IQ1+1
      IQ2=0
      DO NI2=1,T_INFO%NIONS; DO L2=1,TYZ(NI2)%QO%DIM
      IQ2=IQ2+1
      DO IC=1,NC
      DIST=SQRT(SUM((MATMUL(LATT_CUR%A,(R(:,IC)+T_INFO%POSION(:,NI2)-T_INFO%POSION(:,NI1))))**2))
      HR1=ABS(HR(IQ1,IQ2,IC,ISP))
      DO I=1,NDIV
      IF(DIST.GT.RC(I))THEN
        IF(ABS(HR1).GT.ABS(EC(I)))EC(I)=HR1
      ENDIF
      ENDDO; ENDDO ! IC
      ENDDO; ENDDO ! L2,NI2
      ENDDO; ENDDO ! L1,NI1
      ENDDO  ! ISP
      WRITE(IO,'(" TB HOPPING PARAMETERS DECAY EXAMINATION:")')
      WRITE(IO,'(" RCUT          ECUT")')
      DO I=1,NDIV
        WRITE(IO,'(2F12.4)')RC(I),EC(I)
      ENDDO
      RETURN

      END SUBROUTINE TB_PAR_DECAY_ANA


      END MODULE TB_MODULE

! OPEN SUBROUTINES
!**************** SUBROUTINE EXPRO   ***********************************
! EXPRO: caclulates the x-product of two vectors
!***********************************************************************
      SUBROUTINE EXPRO(H,U1,U2)
      USE prec
      IMPLICIT NONE
      REAL(gq) H(3),U1(3),U2(3)

      H(1)=U1(2)*U2(3)-U1(3)*U2(2)
      H(2)=U1(3)*U2(1)-U1(1)*U2(3)
      H(3)=U1(1)*U2(2)-U1(2)*U2(1)

      RETURN
      END SUBROUTINE

!***********************SUBROUTINE DENNP *******************************
! if ISMEAR=0
! subroutine DENSTA calculates a continuous density of states in the
! interval (EMIN,EMAX) by applying a gaussian broadening to the discrete
! eigenvalue spectrum contained in CELEN(NBANDS,NKPTS). The width of the
! gaussians is SIGMA. The fermi energy EFERMI is calculated from the
! integrated dos
! correction to the variational energy is calculated (-SIGMA S)
! according to A.de Vita
! if ISMEAR>0 the generalized form of Methfessel and Paxton of order
!        N=ISMEAR will be used instead of Gaussians to get the dos ...
! routine is parallelized to get full speed ...
! initially it took 2 seconds to calculate occupancies (gK)
!***********************************************************************
      SUBROUTINE DENMP(ISPIN,RSPIN,EMIN,EMAX,NELECT,EFERMI, &
             & ISMEAR,SIGMA,NBANDS,NKPTS,FERWE,CELEN, &
             & NEDOS,DOS,DOSI,WTKPT)
      USE prec; USE gconstant
      IMPLICIT NONE
      INTEGER ISPIN,ISMEAR,NBANDS,NKPTS,NEDOS
      REAL(gq) DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN)
      REAL(gq) NELECT,EMIN,EMAX,EFERMI,SIGMA,RSPIN
      REAL(gq) FERWE(NBANDS,NKPTS,ISPIN),CELEN(NBANDS,NKPTS,ISPIN)
      REAL(gq) WTKPT(NKPTS)
! Local
      INTEGER ISP,K,N,NELOW,NEHIG,I,NITER
      REAL(gq) DELTAE,EPS,WEIGHT,SFUN_DONE,E,EPSDOS,DFUN,SFUN,DOSTOT
      REAL(gq) EF1,EF2,ELECT,X1
      LOGICAL LOWB,HIGHB

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
      DOS =0
      DOSI=0
!=======================================================================
! accumulate dos and integrated dos
!=======================================================================
      DO ISP=1,ISPIN
      DO K=1,NKPTS
      DO N=1,NBANDS
        EPS=CELEN(N,K,ISP)
        WEIGHT= RSPIN*WTKPT(K)
        NELOW=(EPS-8._gq*SIGMA-EMIN)/DELTAE+1
        NEHIG=(EPS+8._gq*SIGMA-EMIN)/DELTAE+1
        IF (NELOW<1)     NELOW=1
        IF (NELOW>NEDOS) NELOW=NEDOS
        IF (NEHIG<1)     NEHIG=1
        IF (NEHIG>NEDOS) NEHIG=NEDOS

        SFUN_DONE=0
        DO I=NELOW,NEHIG
          E=EMIN+DELTAE*(I-1)-EPS
          CALL DELSTP(ISMEAR,(E/SIGMA),DFUN,SFUN)
          EPSDOS=DFUN/SIGMA
!gK fix the DOS so that the integrated DOS yields accurate results
          EPSDOS=(SFUN-SFUN_DONE)/DELTAE
          SFUN_DONE=SFUN

          DOS(I,ISP) =DOS(I,ISP) +(WEIGHT*EPSDOS)
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT*SFUN
        ENDDO
        DO I=NEHIG+1,NEDOS
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT
        ENDDO
      ENDDO; ENDDO; ENDDO
!=======================================================================
! calculate approximated fermi energy
!=======================================================================
      DO I=1,NEDOS
        DOSTOT=SUM(DOSI(I,:))
        IF (ABS(DOSTOT-NELECT)<0.01_gq .OR.DOSTOT>NELECT) EXIT
      ENDDO
      EFERMI= EMIN+(I-1)*DELTAE
      IF (SIGMA<1E-5_gq) RETURN
!=======================================================================
! search now for exact Fermi-level using bisectioning
!=======================================================================
      EF1= EMIN+(I-2)*DELTAE
      LOWB =.FALSE.
      EF2= EMIN+(I-1)*DELTAE
      HIGHB=.FALSE.
      NITER=0

      setfermi: DO

      EFERMI=(EF1+EF2)/2
      NITER=NITER+1

      ELECT=0
      DO ISP=1,ISPIN
      DO K=1,NKPTS
      DO N=1,NBANDS
        EPS=CELEN(N,K,ISP)
        X1=(EFERMI-EPS)/SIGMA
        CALL DELSTP(ISMEAR,X1,DFUN,SFUN)
        FERWE(N,K,ISP)=SFUN
        ELECT=ELECT+FERWE(N,K,ISP)*WTKPT(K)
      ENDDO; ENDDO; ENDDO
      ELECT=ELECT*RSPIN

      ! compare with number of electrons

      IF ( ABS(ELECT-NELECT)<1E-10_gq) RETURN
      IF ( (ABS(EF1-EF2)/(ABS(EFERMI)+1.E-10_gq))<1E-14_gq) EXIT
      IF ( ELECT>NELECT) THEN
        IF (.NOT.LOWB)  EF1=EF1-DELTAE
        HIGHB=.TRUE.
        EF2  =EFERMI
      ELSE
        IF (.NOT.HIGHB) EF2=EF2+DELTAE
        LOWB =.TRUE.
        EF1  =EFERMI
      ENDIF
      ENDDO setfermi
      WRITE(0,*)'WARNING: DENMP: can''t reach specified precision'
      WRITE(0,*)'Number of Electrons is NELECT =',ELECT

      END SUBROUTINE DENMP

!******************** DELSTP    ****************************************
! Returns generalised delta and step functions (Methfessel & Paxton)
!
!  Input:
!      n > -1 : order of approximant; x : argument
!  Output:
!      D_n (x) ,  S_n (x)
!  Remarks:
!      D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
!      S_n (x) = (1 - erf x)/2 + exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
!      where H is a Hermite polynomial and
!      A_i = (-1)^i / ( i! 4^i sqrt(pi) )
!***********************************************************************
      SUBROUTINE DELSTP(N,X,D,S)
      USE prec
      USE gconstant
      IMPLICIT REAL(gq) (A-H,O-Z)

      IF (X<-1.E5_gq) THEN
         D=0._gq
         S=0._gq
         RETURN
      END IF
      IF (X>1.E5_gq) THEN
         D=0._gq
         S=1._gq
         RETURN
      END IF
!=======================================================================
!  If n < 0 : assume Gaussian type smearing
!  (must return  same as N=0 or ... )
!=======================================================================
      IF (N<0) THEN
         D=EXP(-(X*X))/SQRT(PI)
         S=0.5_gq+0.5_gq*ERRF(X)
         RETURN
      END IF
!=======================================================================
! Methfessel & Paxton
!=======================================================================
      EX2=EXP(-(X*X))
      S0=0.5_gq*ERRF(X)
      A=1._gq/SQRT(PI)
      K=0
      H1=1._gq
      H2=2._gq*X
      S=0._gq
      D=A
      DO I=1,N
         A=A/((-4._gq)*I)
         K=K+1
         H3=H1
         H1=H2
         H2=2._gq*X*H2-2*K*H3
         S=S+A*H1
         K=K+1
         H3=H1
         H1=H2
         H2=2._gq*X*H2-2*K*H3
         D=D+A*H1
      ENDDO
      D=D*EX2
      S=0.5_gq+S0-S*EX2
      RETURN
      END SUBROUTINE DELSTP

!**************** SUBROUTINE KARDIR ************************************
! transform a set of vectors from cartesian coordinates to
! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************
      SUBROUTINE KARDIR(NMAX,V,BASIS)
      USE prec
      IMPLICIT NONE
      INTEGER N,NMAX
      REAL(gq) V(3,NMAX),BASIS(3,3),V1,V2,V3

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(2,1)+V(3,N)*BASIS(3,1)
        V2=V(1,N)*BASIS(1,2)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(3,2)
        V3=V(1,N)*BASIS(1,3)+V(2,N)*BASIS(2,3)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE

!**************** SUBROUTINE TOPRIM ************************************
! bring all ions into the primitive cell
!***********************************************************************
      SUBROUTINE TOPRIM(NIONS,POSION)
      USE prec
      IMPLICIT NONE
      INTEGER NIONS,I
      REAL(gq) POSION(3,NIONS)

      DO I=1,NIONS
      POSION(1,I)=MOD(POSION(1,I)+60,1._gq)
      POSION(2,I)=MOD(POSION(2,I)+60,1._gq)
      POSION(3,I)=MOD(POSION(3,I)+60,1._gq)
      ENDDO
      RETURN
      END SUBROUTINE

!++++++++++++++++++++++++++++++++++++++++++++++++++
! TEST IF THE MATRIX IS DIAGONAL!
!++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CHK_MAT_DGL(A,LDA,LDGL)
      USE prec
      IMPLICIT NONE
      LOGICAL LDGL
      INTEGER LDA
      REAL(gq) A(LDA,LDA)
! LOCAL
      INTEGER I,J
      REAL(gq),PARAMETER :: TINY0=1.E-2_gq

      LDGL=.TRUE.
      DO I=1,LDA; DO J=1,LDA
      IF(I.EQ.J)CYCLE
      IF(ABS(A(I,J)).GT.TINY0)THEN
      LDGL=.FALSE.; RETURN
      ENDIF
      ENDDO; ENDDO

      END SUBROUTINE CHK_MAT_DGL

!******************** SUBROUTINE TETIRR ********************************
!
! subroutine tetirr finds inequivalent tetrahedra in an equally spaced
! k-mesh (setting also the correct weights) ...
!
! this routine needs the basis vectors of the k-lattice (BK) [cartesian]
!                    the size of the k-mesh (NKX,NKY,NKZ)
!                    the k-point 'connection list' IKPT
!                    the dimension parameter IKPTD for array IKPT
!                    the dimension parameter NTETD for array IDTET
!                    and the unit where to write informations (IU6)
! it returns:        the number of irreducible tetrahedra (NTET)
!                    the k-index of the four corners [IDTET(1-4,...)]
!                    and the weighting factors [IDTET(0,...)]
!
!***********************************************************************

      SUBROUTINE TETIRR(NTET,IDTET,NTETD,BK,NKX,NKY,NKZ,IKPT,IKPTD,IU6)
      USE prec

      IMPLICIT REAL(gq) (A-H,O-Z)
      DIMENSION IDTET(0:4,NTETD),BK(3,3),IKPT(IKPTD,IKPTD,IKPTD)
      DIMENSION KCUT0(3,4,6),KCUT(3,4,6),P(3,4),IQ(4),IMC(0:1,0:1,0:1)

      SAVE KCUT0
      DATA KCUT0/ &
     &         0,0,0, 0,1,0, 1,1,0, 1,1,1,  0,0,0, 1,0,0, 1,1,0, 1,1,1, &
     &         0,0,0, 1,0,0, 1,0,1, 1,1,1,  0,0,0, 0,1,0, 0,1,1, 1,1,1, &
     &         0,0,0, 0,0,1, 0,1,1, 1,1,1,  0,0,0, 0,0,1, 1,0,1, 1,1,1 /
      ANRM2(X,Y,Z)=X*X*1.00001E0_gq+Y*Y*1.00002E0_gq+Z*Z*1.00003E0_gq &
     &                           -X*0.000004E0_gq-Y*0.000003E0_gq-Z*0.000002E0_gq

! Setting up the tetrahedra will be done cutting a microcell with eight
! corners into six tetrahedra. The edges of the tetrahedra are given by
! three edges, two face diagonals and one space diagonal of the cell.
! Giving the space diagonal, the way how to choose the rest is uniquely
! determined ... . But now there are four different possibilities how to
! choose the space diagonal! Prefer the choice which gives the shortest
! edges for all tetrahedra ('the most compact tetrahedra') - just to
! avoid very long 'interpolation distances' ... !
      LZ=0
      LXX=0
      LYY=0
      EDGMAX=1.E30_gq
      EDGMIN=0._gq
! For the four choices ...
      DO 5 LX=0,1
       DO 5 LY=0,1
! ... we set up the 'trial division' of a given cell into 6 tetrahedra:
         DO 1 ITET=1,6
          DO 1 IC=1,4
            KCUT(1,IC,ITET)=KCUT0(1,IC,ITET)
            KCUT(2,IC,ITET)=KCUT0(2,IC,ITET)
            KCUT(3,IC,ITET)=KCUT0(3,IC,ITET)
            IF (LX==1) KCUT(1,IC,ITET)=1-KCUT0(1,IC,ITET)
            IF (LY==1) KCUT(2,IC,ITET)=1-KCUT0(2,IC,ITET)
    1    CONTINUE
         EDMIN=1.E30_gq
         EDMAX=0._gq
! For this trial setting, loop over all tetrahedra ...,
         DO 4 ITET=1,6
! ... set up the cartesian coordinates of the four corner points ...,
            DO 2 IC=1,4
               P(1,IC)=KCUT(1,IC,ITET)*BK(1,1)+ &
     &                 KCUT(2,IC,ITET)*BK(1,2)+ &
     &                 KCUT(3,IC,ITET)*BK(1,3)
               P(2,IC)=KCUT(1,IC,ITET)*BK(2,1)+ &
     &                 KCUT(2,IC,ITET)*BK(2,2)+ &
     &                 KCUT(3,IC,ITET)*BK(2,3)
               P(3,IC)=KCUT(1,IC,ITET)*BK(3,1)+ &
     &                 KCUT(2,IC,ITET)*BK(3,2)+ &
     &                 KCUT(3,IC,ITET)*BK(3,3)
    2       CONTINUE
! ... and get the shortest and longest distance between two points in
! each tetrahedron (minimum/maximum taken over all tetrahedra ...):
            DO 3 I=1,3
             DO 3 J=I+1,4
               XX=ANRM2(P(1,I)-P(1,J),P(2,I)-P(2,J),P(3,I)-P(3,J))
               EDMAX=MAX(EDMAX,XX)
               EDMIN=MIN(EDMIN,XX)
    3       CONTINUE
    4    CONTINUE
! Now look at the global maximum: Have we found a cut with smaller
! maximum distance between two points within one tetrahedron than
! before? If yes: store it  (until we find something better ...)!
         IF (EDMAX<EDGMAX) THEN
            LXX=LX
            LYY=LY
            EDGMAX=EDMAX
            EDGMIN=EDMIN
         END IF
    5 CONTINUE
! Now set up the 'correct' cutup giving the most compact tetrahdra ... :
      DO 6 ITET=1,6
       DO 6 IC=1,4
         KCUT(1,IC,ITET)=KCUT0(1,IC,ITET)
         KCUT(2,IC,ITET)=KCUT0(2,IC,ITET)
         KCUT(3,IC,ITET)=KCUT0(3,IC,ITET)
         IF (LXX==1) KCUT(1,IC,ITET)=1-KCUT0(1,IC,ITET)
         IF (LYY==1) KCUT(2,IC,ITET)=1-KCUT0(2,IC,ITET)
    6 CONTINUE
! Now start searching the tetrahedra ... :
      NTET=0
! For all k-points ...
      DO 14 I3=1,NKZ
       DO 14 I2=1,NKY
        DO 14 I1=1,NKX
! Set up microcell of 8 k-points (= 8 corners of unit cell of k-mesh):
         DO 7 K1=0,1
          J1=MOD(I1+K1-1,NKX)+1
          DO 7 K2=0,1
           J2=MOD(I2+K2-1,NKY)+1
           DO 7 K3=0,1
            J3=MOD(I3+K3-1,NKZ)+1
! Get the identifiers (the irreducible k-point connected to I1,I2,I3):
            IMC(K1,K2,K3)=IKPT(J1,J2,J3)
    7    CONTINUE
! From this microcell we can cut out six tetrahedra:
         DO 13 ITET=1,6
! Set the 4 corners (identifiers) of the actual tetrahedron:
            DO 8 IC=1,4
               K1=KCUT(1,IC,ITET)
               K2=KCUT(2,IC,ITET)
               K3=KCUT(3,IC,ITET)
               IQ(IC)=IMC(K1,K2,K3)
    8       CONTINUE
! Order the identifiers of the corners ...
            DO 9 J=1,3
             DO 9 I=1,4-J
               IF (IQ(I)>IQ(I+1)) THEN
                  II=IQ(I)
                  IQ(I)=IQ(I+1)
                  IQ(I+1)=II
               END IF
    9       CONTINUE
! First tetrahedron ...
            IF (NTET==0) GOTO 11
! Now test all tetrahedra found previously:
            DO 10 N=1,NTET
               IF ((IDTET(1,N)==IQ(1)) &
     &            .AND.(IDTET(2,N)==IQ(2)) &
     &               .AND.(IDTET(3,N)==IQ(3)) &
     &                  .AND.(IDTET(4,N)==IQ(4))) THEN
! We have found the same combination previously, so increment the
! counter for this type of tetrahedron ...
                  IDTET(0,N)=IDTET(0,N)+1
! ... and go to the next tetrahedron:
                  GOTO 13
               END IF
   10       CONTINUE
! New tetrahedron found if arriving here:
   11       CONTINUE
! Count it, ...
            NTET=NTET+1
! ... store the corner coordiantes (identifier) ...
            DO 12 I=1,4
   12       IDTET(I,NTET)=IQ(I)
! ... and initialize the counter for this new type of tetrahedron:
            IDTET(0,NTET)=1
   13    CONTINUE
   14 CONTINUE
      IF (NTET>NTETD) STOP ' TETIRR: NTET > NTETD!'
! Now tell us the result ... :
      IF (IU6>=0) WRITE(IU6,15) NTET,6*NKX*NKY*NKZ
   15 FORMAT(1X,'TETIRR: Found ',I6,' inequivalent tetrahedra from ',I8)
      RETURN

      END SUBROUTINE TETIRR

!***********************SUBROUTINE DENTET******************************
!
! subroutine DENTET calculates a continuous density of states in the
! interval (EMIN,EMAX) applying the tetrahedron method to the discrete
! eigenvalue spectrum in CELEN(NBANDS,NKPTS). EFERMI is calculated so
! that sum over occupation-numbers is equal to number of electrons
! executed on all nodes
!**********************************************************************

      SUBROUTINE DENTET(IU0,CELEN,WTKPT,NBANDS,NKPTS,DOS,DOSI,NEDOS, &
     &           ISPIN,RSPIN,EMIN,EMAX,IDTET,NTET,VOLWGT,NELECT,EFERMI,FERWE, &
     &           JOB,IU6,NKDIM)
      USE prec
      IMPLICIT COMPLEX(gq) (C)
      IMPLICIT REAL(gq) (A-B,D-H,O-Z)
      LOGICAL LOWB,HIGHB
      REAL(gq)   NELECT
      REAL(gq) CELEN(NBANDS,NKDIM,ISPIN),FERWE(NBANDS,NKDIM,ISPIN)
      DIMENSION DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN),IDTET(0:4,NTET)

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
!=======================================================================
! initialize arrays for dos and integr. dos
!=======================================================================
      DOS =0
      DOSI=0

!=======================================================================
! calculate dos and integrated dos
!=======================================================================
      CALL BZINTS(0,FERWE,CELEN,WTKPT,NBANDS,NBANDS,NKPTS,IDTET,NTET, &
     &     ISPIN,RSPIN,VOLWGT,EMIN,EMAX,DOS,DOSI,NEDOS,EFERMI,SUMWEI, &
     &     SUME,100,NKDIM)
!=======================================================================
! calculate approximated fermi energy
!=======================================================================

      DO I=1,NEDOS
        DOSTOT=DOSI(I,1)
        IF (ISPIN==2) DOSTOT=DOSTOT+DOSI(I,2)
        IF (ABS(DOSTOT-NELECT)<0.01_gq .OR.DOSTOT>NELECT) EXIT
      ENDDO

      EFERMI= EMIN+(I-1)*DELTAE
      IF (JOB==0) RETURN

!=======================================================================
! search now for exact Fermi-level
!=======================================================================
      EF1= EMIN+(I-2)*DELTAE
      LOWB =.FALSE.
      EF2= EMIN+(I-1)*DELTAE
      HIGHB=.FALSE.

!=======================================================================
! calculate fermi-weighting function, and their sum
!=======================================================================
   calcfermi: DO
      EFERMI=(EF1+EF2)/2
      CALL BZINTS(JOB,FERWE,CELEN,WTKPT,NBANDS,NBANDS,NKPTS,IDTET,NTET, &
     &       ISPIN,RSPIN,VOLWGT,EMIN,EMAX,DOS,DOSI,NEDOS,EFERMI,SUMWEI, &
     &       SUME,100,NKDIM)
      ELECT=SUMWEI*RSPIN

!=======================================================================
! compare now with Number of Electrons
!=======================================================================
      IF ( ABS(ELECT-NELECT)<1E-10_gq) GOTO 110
      IF ( (ABS(EF1-EF2)/(ABS(EFERMI)+1.E-10_gq))<1E-14_gq) GOTO 120
      IF ( ELECT>NELECT) THEN
        IF (.NOT.LOWB)  EF1=EF1-DELTAE
        HIGHB=.TRUE.
        EF2  =EFERMI
      ELSE
        IF (.NOT.HIGHB) EF2=EF2+DELTAE
        LOWB =.TRUE.
        EF1  =EFERMI
      ENDIF
      ENDDO calcfermi

  120 CONTINUE
      IF (IU0>=0) THEN
         WRITE(*,*)' WARNING: DENTET: can''t reach specified precision'
         WRITE(*,*)' Number of Electrons is NELECT =',ELECT
      ENDIF

  110 CONTINUE
! Final call to BZINTS (not really absolutely necessary ...) - can be
! commented out if no informational output/debugging desired ... !
      CALL BZINTS(JOB,FERWE,CELEN,WTKPT,NBANDS,NBANDS,NKPTS,IDTET,NTET, &
     &      ISPIN,RSPIN,VOLWGT,EMIN,EMAX,DOS,DOSI,NEDOS,EFERMI,SUMWEI, &
     &      SUME,IU6,NKDIM)

!=======================================================================
! How to calculate the correction term to the total Energy???? Is there
! a correction term at all ('no smearing', analytical interpolation and
! then integration with 'delta-function sampling'!!!) ?????????????????
! People say generally: NO! THERE IS NO ENTROPY! --> believe it or not!
!=======================================================================
      RETURN
      END SUBROUTINE DENTET

!****************** SUBROUTINE BZINTS **********************************
!
      SUBROUTINE BZINTS(JOB,FERWE,CELEN,WTKPT,NBAND,NBANDD,NKPTS,IDTET, &
     &           NTET,ISPIN,RSPIN,VOLWGT,EMIN,EMAX,DOS,DOSI,NEDOS,EFERMI,SUMWEI, &
     &           SUME,IU6,NKDIM)
      USE prec
      IMPLICIT REAL(gq) (A-B,D-H,O-Z)
      IMPLICIT COMPLEX(gq) (C)
!
!***********************************************************************
!
! This routine performs BZ-integrations by the tetrahedron method. It
! has two basic job modes (controlled by flag 'JOB'):
!    JOB=0:         Calculate the DOS and the integrated DOS
!    JOB=-2,-1,1,2: Calculate the Fermi-weights etc., here you have
!                   for submodes: JOB<0 = without Bloechl-correction,
!                   JOB>0 = with Bloechl-correction, for JOB=+/-2 some
!                   additional output is provided (band energy, number
!                   of electrons, Fermi-energy, Bloechl-correction ...)
!
! Input-parameters are CELEN: the band structure data (epsilon_i,ik)
!                      WTKPT: the weights of the k-points
!                      NBAND: the number of bands
!                      NKPTS: the number of irreducible k-points
!                      IDTET: weights/'coordinates' of all tetrahedra
!                      NTET:  number of tetrahedra
!                      EMAX,EMIN: energy window for DOS/DOSI (JOB=0)
!                      NEDOS: number of energy points for DOS/DOSI
!                      EFERMI: approximate/exact Fermi energy (JOB/=0)
!                      IU6: I/O-unit where to write data ...
! Output quantities FERWE:  the Fermi weights for each state (JOB/=0)
!                   DOS:    the density of states (JOB=0)
!                   DOSI:   the integrated density of states (JOB=0)
!                   SUMWEI: the sum of all Fermi weights (JOB/=0)
!                   SUME:   eigenvalue sum ['band energy'] (JOB/=0)
!
!***********************************************************************


      DIMENSION FERWE(NBANDD,NKDIM,ISPIN)
      REAL(gq) CELEN(NBANDD,NKDIM,ISPIN)
      DIMENSION IDTET(0:4,NTET),EC(4),WC(4,2),WTKPT(NKPTS),IQ(4)
      DIMENSION DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN)

! Fatal ERROR!
      IF (NKPTS<4) STOP ' BZINTS: Tetrahedron method fails (number of k-points < 4)!'
      IF ((JOB<(-2)).OR.(JOB>2)) STOP ' BZINTS: JOB must be +/-1 or +/-2 (make weights) or 0 (make DOS)!'
! Initialize arrays for DOS/integrated DOS (if JOB=0) ...
      IF (JOB==0) THEN
         DO 1 ISP=1,ISPIN
         DO 1 I=1,NEDOS
            DOS(I,ISP)=0._gq
    1    DOSI(I,ISP)=0._gq
      ELSE
! ... Fermi weights ...
         DO 2 ISP=1,ISPIN
          DO 2 IK=1,NKPTS
           DO 2 I=1,NBAND
    2    FERWE(I,IK,ISP)=0._gq
      END IF
! ... and eigenvalue sums:
      SEV1=0._gq
      SEV2=0._gq
! Start looping over tetrahedra:
      DO 20 ITET=1,NTET
! Get the four corner points:
       IQ(1)=IDTET(1,ITET)
       IQ(2)=IDTET(2,ITET)
       IQ(3)=IDTET(3,ITET)
       IQ(4)=IDTET(4,ITET)
       DO 20 ISP=1,ISPIN
       DO 20 I=1,NBAND
! Get the eigenvalues at each corner:
         EC(1)=CELEN(I,IQ(1),ISP)
         EC(2)=CELEN(I,IQ(2),ISP)
         EC(3)=CELEN(I,IQ(3),ISP)
         EC(4)=CELEN(I,IQ(4),ISP)
         IF (JOB==0) THEN
! Make the DOS:
            CALL SLINZ(VOLWGT*IDTET(0,ITET)*RSPIN,EC,EMIN,EMAX, &
     &            DOS(1,ISP),DOSI(1,ISP),NEDOS,IQ, &
     &            NKDIM,NBANDD,I)
         ELSE
! Make the weights:
            CALL FSWGTS(VOLWGT*IDTET(0,ITET),EC,EFERMI,WC,(JOB>0))
! Band occupations, band energy and number of electrons ... :
            DO 10 IC=1,4
               SEV1=SEV1+WC(IC,1)*EC(IC)*RSPIN
               SEV2=SEV2+WC(IC,2)*EC(IC)*RSPIN
               SUMWEI=(WC(IC,1)+WC(IC,2))/WTKPT(IQ(IC))
               FERWE(I,IQ(IC),ISP)=FERWE(I,IQ(IC),ISP)+SUMWEI
   10       CONTINUE
         END IF
   20 CONTINUE
! If JOB=+/-2 make additional checks and give some output (if desired):
      IF (ABS(JOB)==2) THEN
         SUME=SEV1+SEV2
! Check the sum of occupation numbers:
         SUMWEI=0._gq
         DO 30 ISP=1,ISPIN
          DO 30 IK=1,NKPTS
           DO 30 I=1,NBAND
   30    SUMWEI=SUMWEI+FERWE(I,IK,ISP)*WTKPT(IK)
         IF ((IU6>=0).AND.(IU6<=99)) &
     &                       WRITE(IU6,40) EFERMI,RSPIN*SUMWEI,SUME,SEV2
      END IF
   40 FORMAT(1X,'BZINTS: Fermi energy:',F10.6,';',F10.6,' electrons'/ &
     &       9X,'Band energy:',F11.6,';  BLOECHL correction:',F10.6)
      RETURN
      END


!****************** SUBROUTINE SLINZ ***********************************
!
      SUBROUTINE SLINZ(VOLWGT,EC,EMIN,EMAX,DOS,DOSI,NEDOS,IQ, &
     &                NKDIM,NBANDD,ISTATE)
      USE prec
      IMPLICIT REAL(gq) (A-H,O-Z)
!
!***********************************************************************
!
! This subroutine adds up the contributions to the density of
! states and the number of states for one single tetrahedron
!
!  Input-parameters are VOLWGT: weight on tetrahedron
!                       EC: energies at corners of tetrahedron
!                       EMIN,EMAX: energy window
!                       NEDOS: number of energy points for DOS/DOSI
!  Output quantities DOS(K): DOS at E(K)=EMIN+(K-1)(EMAX-EMIN)/(NEDOS-1)
!                    DOSI(K): Integrated DOS at E(K)
!
!***********************************************************************



      DIMENSION EC(4),DOS(NEDOS),DOSI(NEDOS),ES(4),EC1(4),IQ(4)

      DE=(EMAX-EMIN)/REAL(NEDOS-1,KIND=gq)
! Sort the energies at the four corners (array EC) into array ES
      DO 1 I=1,4
    1 EC1(I)=EC(I)
      DO 3 I=1,4
         I00=1
         DO 2 J=2,4
    2    IF (EC1(J)<EC1(I00)) I00=J
         ES(I)=EC1(I00)
         EC1(I00)=1.E30_gq
    3 CONTINUE
! Lowest energy still above EMAX ---> no contributions to DOS/DOSI ... :
      IF (ES(1)>=(EMAX+0.00000001_gq*DE)) RETURN
! Highest energy still below EMIN ---> no contribution to DOS and
! contribution of complete tetrahedron to DOSI (1*VOLWGT) ... :
      IF (ES(4)<=(EMIN-0.00000001_gq*DE)) THEN
         DO 4 I=1,NEDOS
            DOSI(I)=DOSI(I)+VOLWGT
    4    CONTINUE
         RETURN
      END IF
! Now the rest ...
      E1=ES(1)
      E2=ES(2)
      E3=ES(3)
      E4=ES(4)
! Now get the minimum and maximum index for the range we have to update
! DOS(I) and DOSI(I) [so that EMIN>E(ISTART) and EMAX<E(ISTOP)] ... :
      ISTART=MAX((INT((E1-EMIN)/DE-0.00000001_gq)),1)
      ISTART=MIN(ISTART,NEDOS)
      ISTOP=MIN((INT((E4-EMIN)/DE+0.00000001_gq))+2,NEDOS)
      ISTOP=MAX(ISTOP,1)
! Some constants occuring in the integration formulas ... :
      IF ((E3-E2)>0._gq) THEN
         C3=VOLWGT*(E1+E2-E3-E4)/((E3-E1)*(E4-E1)*(E3-E2)*(E4-E2))
         C2=VOLWGT*3._gq/((E3-E1)*(E4-E1))
      ELSE
         C3=0._gq
         C2=0._gq
      ENDIF
      C1=C2*(E2-E1)
      C0=C1*(E2-E1)/3._gq
      IF ((E2-E1)>0._gq) THEN
         CC12=VOLWGT/((E2-E1)*(E3-E1)*(E4-E1))
      ELSE
         CC12=0._gq
      ENDIF
      IF ((E4-E3)>0._gq) THEN
         CC34=VOLWGT/((E3-E4)*(E2-E4)*(E1-E4))
      ELSE
         CC34=0._gq
      ENDIF
      DO 7 I=ISTART,ISTOP
         EACT=EMIN+(I-1)*DE
         ADDDOS=0._gq
! Case EACT between E2,E3:
         IF ((E2<EACT).AND.(EACT<=E3)) THEN
            X=EACT-E2
            ADDDOS=C1+X*(2._gq*C2+3._gq*X*C3)
            DOSI(I)=DOSI(I)+C0+X*(C1+X*(C2+X*C3))
! Case EACT between E1,E2:
         ELSE IF ((E1<EACT).AND.(EACT<=E2)) THEN
            X=EACT-E1
            ADDDOS=3._gq*CC12*X*X
            DOSI(I)=DOSI(I)+CC12*X*X*X
! Case EACT between E3,E4:
         ELSE IF ((E3<EACT).AND.(EACT<=E4)) THEN
            X=EACT-E4
            ADDDOS=-3._gq*CC34*X*X
            DOSI(I)=DOSI(I)+VOLWGT-CC34*X*X*X
! Case EACT greater than E4 (might probably happen for I=ISTOP):
         ELSE IF (E4<=EACT) THEN
            DOSI(I)=DOSI(I)+VOLWGT
         END IF
         DOS(I)=DOS(I)+ADDDOS
    7 CONTINUE
! All energies higer than E(ISTOP) give same contribution to DOSI as
! in the case EACT greater than E4 above ... :
      IF (ISTOP<NEDOS) THEN
         DO 10 I=ISTOP+1,NEDOS
            DOSI(I)=DOSI(I)+VOLWGT
   10    CONTINUE
      END IF
      RETURN
      END


!****************** SUBROUTINE FSWGTS **********************************
!
      SUBROUTINE FSWGTS(VOLWGT,EC,EFERMI,W,BLOECH)
      USE prec
      IMPLICIT REAL(gq) (A-H,O-Z)
!
!***********************************************************************
!
! This routine makes the weights for integration up to EFERMI for
! one single tetrahedron ('contributions to occupation numbers').
!
!  Input-parameters are VOLWGT: weight on tetrahedron
!                       EC: energies at corners of tetrahedron
!                       EFERMI: Fermi energy
!  Output quantities W(I,1): 'normal' weights
!                    W(I,2): Bloechl-corrections to the weights
!
!***********************************************************************



      LOGICAL BLOECH
      DIMENSION EC(4),EC1(4),ES(4),W(4,2),W1(4),W2(4),ISORT(4)

! Sort the energies at the four corners (array EC) into array ES
      DO 1 I=1,4
    1 EC1(I)=EC(I)
      DO 3 I=1,4
         I00=1
         DO 2 J=2,4
    2    IF (EC1(J)<EC1(I00)) I00=J
         ES(I)=EC1(I00)
         ISORT(I)=I00
         EC1(I00)=1.E30_gq
    3 CONTINUE
! Initialise weights and some other things:
      DF=0._gq
      DO 4 I=1,4
         W1(I)=0._gq
         W2(I)=0._gq
         W(I,1)=0._gq
    4 W(I,2)=0._gq
! Each k-point belongs to four tetrahedra (and will be touched 4 times
! when looping over all tetrahedra -> 'redefine weight factor' -> VW4):
      VW4=VOLWGT/4._gq
! Lowest energy still >=EFERMI ---> no contributions to weights ... :
      IF (ES(1)>=EFERMI) RETURN
! Highest energy still <=EFERMI ---> just add up full weight (VW4):
      IF (ES(4)<=EFERMI) THEN
         DO 5 I=1,4
    5    W(I,1)=VW4
         RETURN
      END IF
! Now the rest ... :
      E1=ES(1)
      E2=ES(2)
      E3=ES(3)
      E4=ES(4)
! Case EFERMI between E2,E3:
      IF ((E2<EFERMI).AND.(EFERMI<=E3)) THEN
         A31=(EFERMI-E1)/(E3-E1)
         A41=(EFERMI-E1)/(E4-E1)
         A32=(EFERMI-E2)/(E3-E2)
         A42=(EFERMI-E2)/(E4-E2)
         V1=A31*A41
         V2=A31*A42*(1._gq-A41)
         V3=A42*A32*(1._gq-A31)
         W1(1)=(V1*(3._gq-A31-A41)+V2*(2._gq-A31-A41)+V3*(1._gq-A31))*VW4
         W1(2)=(V1+V2*(2._gq-A42)+V3*(3._gq-A32-A42))*VW4
         W1(3)=(V1*A31+V2*A31+V3*(A31+A32))*VW4
         W1(4)=(V1*A41+V2*(A41+A42)+V3*A42)*VW4
         DF=((E1+E2-E3-E4)*A32*A42+2*EFERMI-E1-E2)/((E3-E1)*(E4-E1))
         DF=3._gq*VOLWGT*DF
! Case EFERMI between E1,E2:
      ELSE IF ((E1<EFERMI).AND.(EFERMI<=E2)) THEN
         A21=(EFERMI-E1)/(E2-E1)
         A31=(EFERMI-E1)/(E3-E1)
         A41=(EFERMI-E1)/(E4-E1)
         XXX=A21*A31*A41*VW4
         W1(1)=XXX*(4._gq-A21-A31-A41)
         W1(2)=XXX*A21
         W1(3)=XXX*A31
         W1(4)=XXX*A41
         DF=3._gq*VOLWGT*A31*A41/(E2-E1)
! Case EFERMI between E3,E4:
      ELSE IF ((E3<EFERMI).AND.(EFERMI<=E4)) THEN
         A14=(EFERMI-E4)/(E1-E4)
         A24=(EFERMI-E4)/(E2-E4)
         A34=(EFERMI-E4)/(E3-E4)
         XXX=A14*A24*A34*VW4
         W1(1)=VW4-XXX*A14
         W1(2)=VW4-XXX*A24
         W1(3)=VW4-XXX*A34
         W1(4)=VW4-XXX*(4._gq-A14-A24-A34)
         DF=3._gq*VOLWGT*A14*A24/(E4-E3)
      END IF
! Here the BLOECHL corrections to the weights (if desired) ... :
      IF (BLOECH) THEN
         DO 15 I=1,4
            W2(I)=0._gq
            DO 10 J=1,4
   10       W2(I)=W2(I)+(ES(J)-ES(I))*DF*0.025_gq
   15    CONTINUE
      END IF
! Now store the weights into W with the correct ordering ... :
      DO 20 I=1,4
         J=ISORT(I)
         W(J,1)=W1(I)
         W(J,2)=W2(I)
   20 CONTINUE
      RETURN
      END

