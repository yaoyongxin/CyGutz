SUBROUTINE L2MAIN(Qcomplex,nsymop,Qident,mode,projector,Qrenormalize,fUdmft)
  !------------------------------------------------------------------------------
  !-- Main routine for setting up DMFT transformation to local basis set       --
  !-- The key object, DMFTransform can transform self-energy to Kohn-Sham base --
  !-- as well as project the Kohn-Sham green's function to localized base      --
  !------------------------------------------------------------------------------  
  !--  mode: 'g' -> computes the green's function
  !--        'e' -> eigenvalues are printed
  !--        'u' -> printing transformation
  !--        'x' -> computes xqtl2, just for debugging
  USE param,   ONLY: LMAX2, LOMAX, NRAD, nloat, nrf, nmat, nume, IBLOCK, fastFilesystem
  USE struct,  ONLY: VOL, RMT, iord, iz, nat, mult, pos, tau, rotij, tauij, nsym, iatnr, aname
  USE case,    ONLY: cf, csize, maxsize, nl, cix, ll, Sigind, iatom, legend, maxdim, natom, ncix, shft, isort, crotloc, ifirst
  USE sym2,    ONLY: tmat, phase, idet
  USE com,     ONLY: MINWAV, MAXWAV, iso, iso2, emin, emax, EF, WEIGH, &
                    &elecn,BND_NE,EFMOD,DELTA
  USE abc,     ONLY: KX, KY, KZ, BKX, BKY, BKZ, E, A, ALM, ALML
  USE kpts,    ONLY: mweight, numkpt,MSX,MSY,MSZ,MKNAME
  USE w_atpar, ONLY: w_alo, w_nlo, w_nlov, w_nlon, w_ilo, w_loor, w_lapw, w_ri_mat, w_P, w_DP, w_jatom, w_FJ, w_DFJ, w_allocate, w_deallocate
  USE com_mpi, ONLY: nprocs, myrank, master, FilenameMPI, Gather_MPI, Reduce_MPI, FindMax_MPI, &
                    &cpuID, vector_para,vectors,nvector,fvectors,VECFN
  IMPLICIT NONE
  !------- Input variables
  LOGICAL, intent(in)      :: Qcomplex, Qident,Qrenormalize
  INTEGER, intent(in)      :: nsymop, projector
  CHARACTER*1, intent(in)  :: mode
  CHARACTER*200, intent(in):: fUdmft
  !------ External Functions
  INTEGER :: CountSelfenergy
  !------ Common Blocks
  COMMON /GENER/   BR1(3,3), BR2(3,3)
  REAL*8        :: BR1, BR2
  COMMON /XA/      R(NRAD), BK(3), BKROT(3), BKROT2(3), BKROT3(3), BKRLOC(3)
  REAL*8        :: R, BK, BKROT, BKROT2, BKROT3, BKRLOC
  common /loabc/   alo(0:lomax,2,nloat,nrf)                        ! changed by atpar
  REAL*8        :: alo
  common /lolog/   nlo, nlov, nlon, loor(0:lomax), ilo(0:lomax), lapw(0:lmax2) ! changed by atpar
  INTEGER       :: nlo, nlov, nlon, ilo
  logical       :: loor, lapw
  COMMON /RINTEG/  RI_MAT(0:lmax2,0:lmax2,nrf,nrf,2,2)             ! changed by atpar
  REAL*8        :: RI_MAT
  COMMON /ATSPDT/  P(0:LMAX2,2,nrf), DP(0:LMAX2,2,nrf)             ! changed by atpar
  REAL*8        :: P, DP
  !------- allocatable local varaibles
  COMPLEX(8),ALLOCATABLE :: DMFTU(:,:,:)
  real*8,     allocatable :: a_real(:)
  complex*16, allocatable :: Olapm0(:,:,:)
  integer,    allocatable :: noccur(:,:), nbandsk(:), nemink(:)
  integer,    allocatable :: cixdim(:), iSx(:,:), nindo(:), csort(:), cix_orb(:)
  complex*16, allocatable :: cfX(:,:,:,:),Rspin(:,:,:)

  !------- other local arrays
  complex*16  :: h_yl(2*LMAX2+1,iblock)
  complex*16  :: h_alyl(2*LMAX2+1,iblock,2)
  complex*16  :: h_blyl(2*LMAX2+1,iblock,2)
  complex*16  :: YL((LMAX2+1)*(LMAX2+1))
  integer     :: iorbital(natom,lmax2+1)
  !------ other local variables
  logical     :: QoffDiag, Tcompute
  character*10:: BNAME
  complex*16  :: PHSHEL, CFAC, PH_SPIN(2), csum, cvalue, Sigmaij, xomega
  real*8      :: dummy, PI, TWOPI, EMIST, S, T, Z, exxx, ARG1, ARG2, ARG3, ARGT, ARGT2, FAC
  integer     :: norbitals, cnemin, cnbands
  integer     :: lcase1, l1, N1, nidn, ndands, idt, num1, num2, num, iscf, maxdim2
  integer     :: m1, lms1, ind1, is2, m2, lms2, ind2, irf1, irf2
  integer     :: nomega, icix, LATEQ, JC
  integer     :: i, ii, jj, j, iom, icase, lcase, is, is1, itape, jtape, jatom, i3, lda, ldb, ldc, irf, JR, INDEX
  integer     :: iorb, L, nind, ip, iq, it, ioccur, ic, ikp, iikp, ivector, &
               & iks,nkp,N, NE, NEMIN, NEMAX, nbands, isym, M, IND_YL, &
               & ibb, LATOM
  integer     :: fh_sig, fh_dos, fh_gc, fh_dt, fh_eig 
  integer     :: wat, iat, im, isrt, iucase, wfirst, wicase, pr_proc, max_bands, natm, maxucase, iorb1, nind1
  integer     :: fh_p, numk, ind, lfirst, inum
  character*100 :: filename
  logical     :: pform
  !
  CHARACTER*200 :: FNAME
  complex*16  :: IMAG, gtc
  real*8, PARAMETER       :: Ry2eV = 13.60569193D0
  CHARACTER*3 STR
  INTEGER,PARAMETER :: GIU=90
  !------------------------------------------------------------------
  DATA IMAG/(0.0D0,1.0D0)/
  !------------------------------------------------------------------     
  QoffDiag = .FALSE.   !--- if off-diagonal ri_mat overlap matrix is needed (between different l-s) ---!
  fh_sig = 80
  fh_dos = 100
  fh_gc  = 120
  fh_dt  = 140
  fh_eig = 180
  !------ Some parameters which should be read from the input! Not yet implemented! -----!
  !gamma  = gamma/Ry2eV  !--- broadening for all orbitals to Ry ---------!
  
  !----------  Some old stuff ---------------------
  PI=ACOS(-1.0D0)
  TWOPI=2.D0*PI
  MAXWAV=0
  MINWAV=100000

  tmat(:,:,:iord) =  iz(:,:,:iord)

  natm = sum(mult) 
  ALLOCATE( csort(nat) )
  CALL Create_Atom_Arrays(csort, maxucase, isort, mult, iatom, nat, natm, natom)

  !----------- Find index to atom/L named iorbital(icase,lcase) -------------------------------!
  CALL Create_Orbital_Arrays(iorbital, norbitals, maxdim2, nl, ll, cix, natom, lmax2, ncix, iso)
  
  allocate( noccur(maxsize,ncix) )
  ALLOCATE( cixdim(ncix), iSx(maxdim2, norbitals), nindo(norbitals), cix_orb(norbitals), cfX(maxdim2,maxdim2,norbitals,norbitals) )
  cfX=0
  
  CALL Create_Other_Arrays(cixdim, iSx, noccur, nindo, cix_orb, cfX, CF, nl, ll, cix, iorbital, csize, Sigind, iso, natom, maxdim, lmax2, ncix, maxsize, norbitals, maxdim2)

  ALLOCATE( Rspin(2,2,norbitals) )
  CALL GetSpinRotation(Rspin,rotij,crotloc,norbitals,natom,natm,iso,lmax2,iatom,nl,cix,ll,iorbital)


  filename = 'BasicArrays.dat'
  if (mode.EQ.'u') CALL PrintSomeArrays(filename, nat, iso, norbitals, ncix, natom, numkpt, nmat, nume, Qcomplex, lmax2, maxdim2, maxdim, maxsize, nindo, cixdim, nl, ll, cix, iorbital, csize, iSx, Sigind, EF, VOL)

  IF(myrank.EQ.master) CALL GUTZ4_WRT(CFX,MAXDIM2,NORBITALS,NINDO)
  !-------------------------------------------------------------------------------------!
  !--------------- Reading potential parameters for all atoms.   -----------------------!
  !-- It stores  necessary parameters and uses them below, when looping over atoms -----!
  !-- The step is necessary because green's function of each orbital depends on the ----!
  !--  self-energy of all l's and all atoms. If more than one atom is correlated, we ---!
  !--  can not avoid looping over all correlated atoms in the inside loop --------------!
  !-------------------------------------------------------------------------------------!
  do is=1,iso2 !------ preparing input files case.vsp and case.vspn containing potential parameters -----!
     jtape=17+is
     READ(jtape,2032) ISCF
  end do
  
  call w_allocate(maxucase) !------ contains all variables which depend on atom / which are changed by atpar -----------!


  DO jatom=1,nat    !------  over all atoms (all sorts) even if we do not require qtl for them ------------------!
     iucase = csort(jatom) !--- this particular atom is requested in the input under index iucase --!
     if(iucase.gt.0) then                 !--- requested in the inut ----!
        if (myrank.EQ.master .OR. fastFilesystem) write(6,'(A,1x,I2,1x,A,1x,I2)')' Calculation for atom', jatom, 'which has iucase=', iucase
        w_jatom(iucase) = jatom
     else
        if (myrank.EQ.master .OR. fastFilesystem) write(6,'(A,1x,I2,1x,A)')' ATOM=', jatom, ' left out ' !-- not requested in the input ----!
     endif
     
     !----  mult(jatom) stands for multiplicity of each atom -> number of equivalent atoms of this sort --------!
     !--- Reading radial functions U(R), UE(R), for all atoms including those that are left out ----------------!
     do is=1,iso
        itape=8+is
        jtape=17+is     
        if (vector_para) then
           if (nvector.ge.1) then
              FNAME = fvectors(1,is)
           else   ! No k-point to compute, but still need linearization energies
              FNAME = TRIM(VECFN(1))//'_1'
           endif
           open(itape,FILE=FNAME,STATUS='old',FORM='unformatted')
        else
           !------------ vector files need to be read from the beginning -------------------!
           rewind(itape)
        endif
        if(is.eq.2.and.iso2.eq.1) then  ! fix for SO-coupling and non-spinpolarized
           jtape=18
           rewind(jtape)
           READ(jtape,2032) ISCF 
        endif
        CALL ATPAR(jatom,itape,jtape,is,iso,QoffDiag)
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     end do
          
     if (iucase.gt.0) then ! save values of potential parameters
        if (myrank.EQ.master .OR. fastFilesystem) write(6,'(A,I2,A,I2)')'Saving potential parameters for atom',jatom, ' into ', iucase
        w_alo(:,:,:,:,iucase) = alo(:,:,:,:)
        w_nlo(iucase)  = nlo
        w_nlov(iucase) = nlov
        w_nlon(iucase) = nlon
        w_loor(:,iucase) = loor(:)
        w_ilo(:,iucase)  = ilo(:)
        w_lapw(:,iucase) = lapw(:)
        w_ri_mat(:,:,:,:,:,:,iucase) = ri_mat(:,:,:,:,:,:)
        w_p(:,:,:,iucase) = P(:,:,:)
        w_dp(:,:,:,iucase) = DP(:,:,:)
     endif
  ENDDO

!!! ---------- Preparation of arrays for paralel executaion --------------
  if (vector_para) then
     pr_proc = sum(vectors(:,2))
  else
     pr_proc  = floor(numkpt/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor
  endif
  if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'tot-k=', numkpt

  if (Qrenormalize) then
     allocate(Olapm0(maxdim,maxdim,ncix) )
     CALL cmp_overlap(projector,Olapm0, Qcomplex, nsymop, csort, iorbital, cix_orb, nindo, cixdim, &
                     &iSx, noccur, cfX, Rspin, maxucase, maxdim2, norbitals, pr_proc)
  endif
  
  if (myrank.EQ.master .OR. fastFilesystem) then
     WRITE(6,*) '------------------------------------------------------------------------------'
  endif

  !------------ allocating some important arrays --------------!
  ALLOCATE( a_real(nmat) ) !---  for eigenvectors -------------!

  IF(myrank==master) CALL GUTZ5_WRT(numkpt,MWEIGHT,ELECN,BND_NE,MSX,MSY,MSZ,MKNAME,EFMOD,DELTA)

  !!! start-kpoints
  iikp=0
  DO ivector=1,nvector
     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           FNAME = fvectors(ivector,is)
           open(itape,FILE=FNAME,STATUS='old',FORM='unformatted')
        else
           !------------ vector files need to be read from the beginning -------------------!
           rewind(itape) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        DO I=1,NAT
           READ(itape) EMIST !---- At the beginninge we have linearization energies --!
           READ(itape) EMIST !---- We just skip them ---------------------------------!
        ENDDO
     END DO

     if (vector_para) then
        nkp = vectors(ivector,2)
        WRITE(STR,'(I3)')vectors(ivector,1)
     else
        nkp= numkpt
        WRITE(STR,'(I3)')myrank
     endif
     OPEN(GIU,FILE='BNDU_'//TRIM(ADJUSTL(STR))//'.INP',STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")

     DO iks=1,nkp   !------ Over all irreducible k-points ------!
        Tcompute=.TRUE.
        if (vector_para) then
           ikp = vectors(ivector,3)+iks  ! successive index in k-point table from case.klist
           iikp = iikp+1                 ! successive index in k-point table on this processor
        else
           ikp = iks                     ! successive index in k-point table from case.klist
           !--- We need to go over all k-points even though we will compute only some of them on this processor.
           !--- This is because we need to read vector file sequentially.
           iikp = ikp-myrank*pr_proc            ! The index of the point to compute. If negative, do not compute!
           if (iikp.GT.pr_proc) EXIT            ! Processor finished. The rest of the points will be taken care of by the other pro
           if (iikp.LE.0) Tcompute=.FALSE. 
        endif

     !------- reading from vector for both spins -----------!
        DO is=1,iso    !------ over up/dn ---------------------!
           itape=8+is
           READ(itape,END=998) S,T,Z,BNAME,N,NE  !--- Here we can jupm out of loop 4 and stop at 998 -----------------------------------!
           IF(N.GT.MAXWAV) MAXWAV=N                                         
           IF(N.LT.MINWAV) MINWAV=N                                         
           READ(itape) (KX(I),KY(I),KZ(I),I=1,N) !--- Reads all reciprocal vectors -----------------------------------------------------!

           NUM=0
           DO WHILE (NUM.NE.NE)
              READ(itape) NUM,exxx                !----- eigenvalues read -----!
              E(NUM)=exxx
              if (.not.Qcomplex) then
                 READ(itape) (A_real(I),I=1,N)    !----- eigenvector read -----!
                 a(:,num,is) = a_real(:)
              else
                 READ(itape) (A(I,NUM,is),I=1,N)  !--- eigenvector complex of size n (where n is number of reciprocal vectors) --!
              end if
           ENDDO
           NEMIN = BND_NE(2,IKP)
           NEMAX = BND_NE(3,IKP)

           IF (Tcompute) THEN
              DO I=1,N                              !----- Over reciprocal vectors --------------------------------------------------------!
                 BKX(I)=(S+KX(I))                   !----  Creates K+k where K is reciprocal vector and k is irreducible k-point ----------!
                 BKY(I)=(T+KY(I))
                 BKZ(I)=(Z+KZ(I))
              ENDDO

              DO iucase=1,maxucase
                 !--- Computes Bessel functions for all needed l's and atoms. --------------------------------!
                 !---  For optimization purposes, done only for irreducible k-points -------------------------!
                 !---  output: FJ, DFJ / Spherical Bessel functions and derivative, uses common block BR1/general ---!
                 CALL HARMON(N,BKX,BKY,BKZ,LMAX2,w_FJ(:,:,iucase),w_DFJ(:,:,iucase),RMT(w_jatom(iucase))) 
              ENDDO
           ENDIF
        ENDDO  ! is
        WRITE(GIU)E(1:NE)
        if (.not.Tcompute) CYCLE  ! This k-points was read, but will not be computed on this processor
        !---  Finished reading eigenvectors and eigenvalues for both spins -----!
        !---   E(iband) contains eigenvalues   ---------------------------------!
        !---   A(:,iband,is) contains eigenvectors -----------------------------!
        
        nbands = nemax-nemin+1
        allocate( DMFTU(nbands,maxdim2,norbitals) )
        !----------- sum over all k-points in the star of the irreducible k-point ---------------------------------!
        DO isym=1,nsymop  !-- Over all symmetry operarations -> all k-points in the start of irreducible k point --!
           DO is=1,iso
              PH_SPIN(is)=EXP((2*is-3)*IMAG*PHASE(isym)/2)  !-- simple phase factor for this k-point ----!
           END DO
           DMFTU=0
           DO icase=1,natom  !--------------- over all atoms requested in the input ------------------------!
              latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
              jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
              lfirst = ifirst(latom)
              iucase = csort(jatom)  ! The renumbert sorts, such that the required atoms from the input give continuous index
              !----  Setting all values of common-blocks which depend on atom and were computed by atpar ---!
              alo(:,:,:,:) = w_alo(:,:,:,:,iucase)
              nlo = w_nlo(iucase)                 
              nlov = w_nlov(iucase)               
              nlon = w_nlon(iucase)               
              loor(:) = w_loor(:,iucase)          
              ilo(:) = w_ilo(:,iucase)            
              lapw(:) = w_lapw(:,iucase)          
              ri_mat(:,:,:,:,:,:) = w_ri_mat(:,:,:,:,:,:,iucase)
              P(:,:,:) = w_p(:,:,:,iucase)        
              DP(:,:,:) = w_dp(:,:,:,iucase)      

              FAC=4.0D0*PI*RMT(jatom)**2/SQRT(VOL)
              
              do lcase=1,nl(icase) !----------- loop over L(jatom) requested in the ionput ---------------!
                 l=ll(icase,lcase) !------ current L --!
                 CFAC=IMAG**L      !------  (i)^l -----!
                 ALM = 0.0         !------  ALM(m,band,nrf,is) will hold product of eigenvectors and a/b expansion coefficients --!
                 
                 !--------- blocks are for efficiency. Matrix is multiplied in block form. This must be important in the past, while modern BLAS should do that better. I think it is obsolete.
                 DO ii=1,N-(nlo+nlon+nlov),iblock !------ iblock is 128 for 32-bit system -------!
                    !-------- nlo-number of local orbitals -----!
                    i3=0                   
                    do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))  ! 121
                       !---------  rotates ylm(k+K) to ylm(k'+K) where k' is in star of irreducible k. ------------!
                       i3=i3+1
                       BK(1)=BKX(I) !-----  reciprocal vector and irreducible vector: G=K+k ----!
                       BK(2)=BKY(I)
                       BK(3)=BKZ(I)
                       ! BKROT = R_a.(k+K) transforms to the reducible k-point
                       CALL ROTATE (BK,TMAT(1,1,isym),BKROT)  
                       ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom 
                       CALL ROTATE (BKROT, rotij(1,1,latom), BKROT2)
                       !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
                       ! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
                       BKROT3(1)=BKROT2(1)*BR1(1,1)+BKROT2(2)*BR1(1,2)+BKROT2(3)*BR1(1,3) 
                       BKROT3(2)=BKROT2(1)*BR1(2,1)+BKROT2(2)*BR1(2,2)+BKROT2(3)*BR1(2,3)
                       BKROT3(3)=BKROT2(1)*BR1(3,1)+BKROT2(2)*BR1(3,2)+BKROT2(3)*BR1(3,3)
                       !---- BKRLOC = crotloc.R_n.R_a.(k+K),  rotates according to the user specified local coordinate system.
                       CALL ROTATE (BKROT3,crotloc(1,1,icase),BKRLOC)
                       !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
                       CALL YLM (BKRLOC,LMAX2,YL)  ! 
                       ! (R_n.R_a.(k+K)) *  R(first) * 2pi
                       ARG1=BKROT2(1)*(POS(1,lfirst))*TWOPI
                       ARG2=BKROT2(2)*(POS(2,lfirst))*TWOPI
                       ARG3=BKROT2(3)*(POS(3,lfirst))*TWOPI
                       ! ARGT = (k+K)*tau(isym) * 2pi
                       ARGT =(BKX(I)*TAU(1,isym)+BKY(I)*TAU(2,isym)+BKZ(I)*TAU(3,isym))*TWOPI
                       ! ARGT2 = (R_a.(k+K)).tau_n * 2pi
                       ARGT2=(BKROT(1)*(tauij(1,latom)+shft(latom,1))+BKROT(2)*(tauij(2,latom)+shft(latom,2))+BKROT(3)*(tauij(3,latom)+shft(latom,3)))*TWOPI                    
                       ! PHSEHL = e^{I*2pi*( (R_a.(k+K))*tau_n + (K+k)*tau(isym) + (R_n.R_a.(k+K)*R(first)))}
                       PHSHEL=EXP(IMAG*(ARG1+ARG2+ARG3+ARGT+ARGT2))
                       DO  M=1,2*L+1
                          IND_YL=M+L*L
                          h_yl(M,i3)=conjg(yl(ind_yl))*phshel !----- h_yl is rotated yl when k is rotated to k' -----!
                          !WRITE(*,'(7I2,1x,2f20.15)') ikp, icase, is, ii, i, lcase, M, h_yl(M,i3)
                       END DO
                    enddo
                    
                    DO is=1,iso  !--- over both spins
                       i3=0
                       do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                          i3=i3+1
                          DO M=1,2*L+1
                             if (lapw(l)) then
                                h_ALYL(m,i3,is)=(w_DFJ(L,I,iucase)*P(l,is,2)-w_FJ(L,I,iucase)*DP(l,is,2))* h_yl(M,i3)*ph_spin(is) ! derivatives of bessel functions and spheric harmonics
                                h_BLYL(m,i3,is)=(w_FJ(L,I,iucase)*DP(l,is,1)-w_DFJ(L,I,iucase)*P(l,is,1))* h_yl(M,i3)*ph_spin(is) 
                             else
                                h_ALYL(m,i3,is)=w_FJ(L,I,iucase)/P(l,is,1)/RMT(jatom)**2*h_yl(M,i3)*ph_spin(is)
                                h_BLYL(m,i3,is)=(0.d0,0.d0)
                             end if
                             !WRITE(*,'(7I2,1x,2f20.15,1x,2f20.15,1x,2f20.15)') ikp, icase, is, ii, lcase, i, M, h_yl(M,i3), h_alyl(m,i3,is), PHASE(isym)
                          END DO
                       enddo
                       ibb=min(iblock,N-(nlo+nlon+nlov)-ii+1)
                       lda=2*LMAX2+1
                       ldc=lda
                       ldb=nmat                    
                       !---- h_alyl(2*lmax+1,iblock,is)  contains rotated Apw's, such that chi(r) = (Apw*u(r) + Bpw*udot(r))*Ylm(r)
                       !---- h_blyl(2*lmax+1,iblock,is)  contains rotated Bpw's
                       !---- A(:N,iband,is)              contains eigenvectors, also named C(k+G,inad,is)
                       !---- alm[2*l+1,iband][irf=1,is] += sum_{iK\in block} h_alyl[2*l+1,iK][is]*A[iK,iband][is]
                       !---- alm[2*l+1,iband][irf=2,is] += sum_{iK\in block} h_blyl[2*l+1,iK][is]*A[iK,iband][is]
                       !---- 
                       !---- The results is:
                       !---- alm[lm,iband][1,is] = sum_G Apw(lm,is,K+G) * C(k+G,iband,is)
                       !---- alm[lm,iband][2,is] = sum_G Bpw(lm,is,K+G) * C(k+G,iband,is)
                       !---- Where C(k+G,iband,is) are eigenvectors, and Apw and Bpw are expansion coefficients defined in Shick et.al., PRB 60, 10763 (1999).
                       call zgemm('N','N',2*l+1,nemax-nemin+1,ibb,(1.d0,0.d0), h_alyl(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), alm(1,nemin,1,is),ldc)
                       call zgemm('N','N',2*l+1,nemax-nemin+1,ibb,(1.d0,0.d0), h_blyl(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), alm(1,nemin,2,is),ldc)
                    ENDDO !----------- over both spins ---------!
                 ENDDO    !----------- over iblock -------------!
                 
                 !-------------- Adds localized orbitals to alm. -------------------!
                 if (nlo.ne.0) then
                    call lomain (nemin,nemax,lfirst,latom,n,jatom,isym,L,iso,icase)
                 end if
                 
                 DO M=1,2*L+1  
                    DO iNUM=NEMIN,NEMAX 
                       DO is=1, iso    
                          DO irf=1,nrf 
                             ALML(L,M,iNUM,irf,is)=ALM(M,iNUM,irf,is)*FAC*CFAC !--------- FACCFAC = (-1)^l*4.*PI*RMT(JATOM)**2/SQRT(VOL) --------!
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              enddo  !--------- end of atom L loop (alm)  ----------------!
              !----------------  ALML(l,m,iband,ifr,ispin) contains (A*C,B*C) for all l,m,iband,ifr=(1,2),ispin --------!
              
              !---------- Computing the DMFT transformation --------------!
              do lcase=1,nl(icase)
                 
                 l1=ll(icase,lcase)
                 nind=(2*l1+1)*iso
                 idt = idet(isym)
                 icix = cix(icase,lcase)
                 iorb = iorbital(icase,lcase)

                 CALL CmpDMFTRANS(DMFTU, alml, ri_mat, projector, cfX, Rspin, l1, iorb, cix_orb, nindo, iso, nemin, nbands, nind, maxdim2, maxdim, icix, idt, lmax2, nume, nrf, natom, ncix, nrf, norbitals)
              enddo
           ENDDO  ! over atom-cases
           !-------- computing correlated green's function --------!
           if (Qrenormalize) CALL RenormalizeTrans(DMFTU, Olapm0, cix_orb, nindo, iSx, nbands, maxdim2, norbitals, maxdim, ncix)
           DO I=1,NORBITALS; DO II=1,NINDO(I)
             WRITE(GIU)DMFTU(:,II,I)
           ENDDO; ENDDO
        ENDDO !------  isym: over star of the irreducible k-point ----!
        
        WRITE(*,'(I3,A,1x,I4,1x,3(A,I4))') myrank, ') Finished k-point number', ikp, 'with #bands=', nbands,' nemin=',nemin,' ne=',ne
        DEALLOCATE(DMFTU)
     ENDDO !---- over reducible k-point: ikp
998  CONTINUE !---- irreducible k-points end (jump) from reading
  !--- end k-points
     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     ENDDO
  CLOSE(GIU)
  END DO ! over different vector files

  DEALLOCATE( Rspin )

  !---- Deallocating some arrays before MPI_Gather
  CALL w_deallocate()  
  DEALLOCATE(a_real)

  deallocate( noccur )
  deallocate( cixdim, iSx, nindo, cix_orb, cfX )
  deallocate( csort )
  
  if (Qident) close(210)
  
  if (Qrenormalize) deallocate( Olapm0 )
  !------ Deallocation of memory ---------------
  
205 FORMAT(/,1X,' K-POINT:',3F8.4,1X,I5,I4,2X,A10)                 
200 FORMAT(i5,' atom population matrix')
201 FORMAT(' L1=',i2,', L2=',i2)
202 FORMAT(' L=',i2)
654 FORMAT(6i3,2f14.8)
655 FORMAT(2f14.8,6i3,' Ms1,Ms2,L1,L2,ML1,ML2')
971 FORMAT(1X,f14.8,2x,10f14.8)                               
2032 FORMAT(50X,I2,//)          
1060 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'=',F8.3,//   ,':FER  :',1X,'F E R M I - ENERGY',11X,'= ',F9.5)                                              
!********************************************************************* *PH*
      
  RETURN                    
950 CALL OUTERR('l2main','error reading parallel vectors')
999 CONTINUE
  STOP 'L2main - Error'
!
END SUBROUTINE L2MAIN


SUBROUTINE Create_Atom_Arrays(csort, maxucase, isort, mult, iatom, nat, natm, natom)
  USE param,   ONLY: fastFilesystem
  USE com_mpi, ONLY: myrank, master
  IMPLICIT NONE
  INTEGER, intent(out) :: csort(nat), maxucase
  INTEGER, intent(in)  :: isort(natm), mult(nat), iatom(natom), nat, natm, natom
  !local variables
  INTEGER :: im, iat, wfirst, wat, icase, wicase, iucase
  csort=0
  maxucase=1
  do icase=1,natom
     if (csort(isort(iatom(icase))).EQ.0) then
        csort(isort(iatom(icase))) = maxucase
        maxucase = maxucase + 1
     endif
  enddo
  maxucase = maxucase-1
  
  if (myrank.EQ.master .OR. fastFilesystem) then
     WRITE(6,*)
     WRITE(6,*) '********** Information about the atomic index arrays ************'
     WRITE(6,*) 'iatom=', iatom
     WRITE(6,*) 'isort=', isort
     WRITE(6,*) 'csort=', csort
     WRITE(6,*) 'maxucase=', maxucase
  
     DO icase=1,natom
        iucase = csort(isort(iatom(icase)))
        WRITE(6, '(A,6I6,1x)') 'icase, iucase, jatom, latom=', icase, iucase, isort(iatom(icase)), iatom(icase)
     ENDDO
  endif
  
END SUBROUTINE Create_Atom_Arrays

SUBROUTINE Create_Orbital_Arrays(iorbital, norbitals, maxdim2, nl, ll, cix, natom, lmax2, ncix, iso)
  USE param,   ONLY: fastFilesystem
  USE com_mpi, ONLY: myrank, master
  IMPLICIT NONE
  INTEGER, intent(out) :: iorbital(natom,lmax2+1), maxdim2, norbitals
  INTEGER, intent(in)  :: nl(natom), ll(natom,4), cix(natom,4), natom, lmax2, ncix, iso
  !----- locals
  INTEGER :: icase, lcase, l1, nind, icix
  !----------- Find index to atom/L named iorbital(icase,lcase) -------------------------------!
  iorbital=0
  maxdim2=0   !-- maximum size of the matrix over all atoms and l's requested in the input ----!
              !---- do not confuse with maxdim, which is maximum size for correlated orbitals -!
  norbitals=0 !-- number of atom/l cases, called orbitals -------------------------------------!
  do icase=1,natom
     do lcase=1,nl(icase)
        norbitals = norbitals+1
        iorbital(icase,lcase)=norbitals  !-- index to the orbital number ----!
        icix = cix(icase,lcase)
        
        l1=ll(icase,lcase)
        nind=(2*l1+1)*iso
        maxdim2 = max(maxdim2, nind )

     enddo
  enddo  

  if (myrank.EQ.master .OR. fastFilesystem) then
     WRITE(6,*)
     WRITE(6,*) '*****  Arrays for correlated blocks ******'
     WRITE(6,*) 'norbitals=', norbitals
     WRITE(6,*) 'maxdim2=', maxdim2
     DO icase=1,natom
        DO lcase=1,nl(icase)
           WRITE(6,'(A,I2,A,I2,A,I4)') 'iorbital(',icase,',',lcase,')=', iorbital(icase,lcase)
        ENDDO
     ENDDO
     DO icase=1,natom
        DO lcase=1,nl(icase)
           WRITE(6,'(A,I2,A,I2,A,I4)') 'cix(',icase,',',lcase,')=', cix(icase,lcase)
        ENDDO
     ENDDO
  endif
  !----------- Find index to atom/L named iorbital(icase,lcase) --------
END SUBROUTINE Create_Orbital_Arrays
  

SUBROUTINE Create_Other_Arrays(cixdim, iSx, noccur, nindo, cix_orb, cfX, CF, nl, ll, cix, iorbital, csize, Sigind, iso, natom, maxdim, lmax2, ncix, maxsize, norbitals, maxdim2)
  USE param,   ONLY: fastFilesystem
  USE com_mpi, ONLY: myrank, master
  IMPLICIT NONE
  INTEGER, intent(out)    :: cixdim(ncix), iSx(maxdim2, norbitals), noccur(maxsize,ncix), nindo(norbitals), cix_orb(norbitals)
  COMPLEX*16, intent(out) :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  COMPLEX*16, intent(in)  :: CF(maxdim,maxdim,ncix)
  INTEGER, intent(in)     :: nl(natom), ll(natom,4), cix(natom,4), iorbital(natom,lmax2+1), csize(ncix)
  INTEGER, intent(in)     :: Sigind(maxdim,maxdim,ncix)
  INTEGER, intent(in)     :: natom, maxdim, lmax2, ncix, iso, maxsize, norbitals, maxdim2
  !----- locals
  INTEGER :: icase, lcase, iorb, icix, L, nind, ip, iq, it, ioccur, i, ic, l1, nind1, nind2, ip1, iorb1, iorb2, ind1, ind2

  cixdim=0
  iSx=0
  DO icase=1,natom            
     do lcase=1,nl(icase)     
        icix = cix(icase,lcase)
        l1 = ll(icase,lcase)
        iorb = iorbital(icase,lcase)
        nind1 = (2*l1+1)*iso
        nindo(iorb) = nind1

        cix_orb(iorb) = icix
        
        if ( icix.EQ.0 ) CYCLE
        do ip1=1,nind1
           cixdim(icix) = cixdim(icix) + 1
           iSx(ip1,iorb) = cixdim(icix)
           if (myrank.EQ.master .OR. fastFilesystem) WRITE(6,'(A,7I5)') 'icase,lcase,icix,iorb,nind1,ip1,iSx=', icase, lcase, icix, iorb, nind1,ip1,cixdim(icix)
        enddo
     enddo
  ENDDO
  noccur=0
  DO icix=1,ncix
     DO ip=1,cixdim(icix)
        DO iq=1,cixdim(icix)
           it = Sigind(ip,iq,icix)
           if (it.gt.0)  noccur(it,icix) = noccur(it,icix) + 1
        ENDDO
     ENDDO
  ENDDO

  DO iorb1=1,norbitals
     nind1 = nindo(iorb1)
     if ( cix_orb(iorb1).EQ.0 ) CYCLE
     DO iorb2=1,norbitals
        nind2 = nindo(iorb2)
        if ( cix_orb(iorb1).NE.cix_orb(iorb2) ) CYCLE
        icix = cix_orb(iorb1)
        do ind1=1,nind1
           do ind2=1,nind2
              !WRITE(*,'(I2,I2,I3,I3,I3,I3,2x,I3,2f14.5)') iorb1, iorb2, ind1, ind2, iSx(ind1,iorb1), iSx(ind2,iorb2), icix, CF( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
              cfX(ind1,ind2,iorb1,iorb2) = CF( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
           enddo
        enddo
     ENDDO
  ENDDO
  
  
  if (myrank.EQ.master .OR. fastFilesystem) then
     do icix=1,ncix
        WRITE(6,'(A,I2,A,I2)') 'cixdim(', icix, ')=', cixdim(icix)
        DO it=1,csize(icix)
           WRITE(6,'(A,I2,A,I2,A,I4)') 'noccur(',it,',',icix,')=', noccur(it,icix)
        ENDDO
     enddo
     do iorb=1,norbitals
        do ip1=1,nindo(iorb)
           WRITE(6,'(A,I2,A,I2,A,I3)') 'iSx(', ip1, ',', iorb, ')=', iSx(ip1,iorb)
        enddo
     enddo
     do iorb=1,norbitals
        WRITE(6,'(A,I2,A,I3)') 'nindo(', iorb, ')=', nindo(iorb)
     enddo
     do iorb1=1,norbitals
        do iorb2=1,norbitals
           if ( cix_orb(iorb1).EQ.0 ) CYCLE
           if ( cix_orb(iorb1).NE.cix_orb(iorb2) ) CYCLE
           WRITE(6,'(A,I2,A,I2)') 'CF for iorb1=', iorb1, ' iorb2=', iorb2
           do ind1=1,nind1
              do ind2=1,nind2
                 WRITE(6,'(f14.8,1x,f14.8,3x)', advance='no') real(cfX(ind1,ind2,iorb1,iorb2)), aimag(cfX(ind1,ind2,iorb1,iorb2))
              enddo
              WRITE(6,*)
           enddo
        enddo
     enddo
  endif
END SUBROUTINE Create_Other_Arrays


SUBROUTINE PrintSomeArrays(filename, nat, iso, norbitals, ncix, natom, nkpt, nmat, nume, Qcomplex, lmax2, maxdim2, maxdim, maxsize, nindo, cixdim, nl, ll, cix, iorbital, csize, iSx, Sigind, EF, VOL)
  IMPLICIT NONE
  CHARACTER*100, intent(in) :: filename
  INTEGER, intent(in) :: nat, iso, norbitals, ncix, natom, nkpt, nmat, nume, lmax2, maxdim2, maxdim, maxsize
  LOGICAL, intent(in) :: Qcomplex
  INTEGER, intent(in) :: nindo(norbitals), cixdim(ncix), nl(natom), ll(natom,4), cix(natom,4), iorbital(natom,lmax2+1)
  INTEGER, intent(in) :: csize(ncix), iSx(maxdim2,norbitals), Sigind(maxdim,maxdim,ncix)
  REAL*8, intent(in)  :: EF, VOL
  ! locals
  INTEGER :: fh, icase, lcase, iorb, ip1, icix, ip, iq
  fh=996
  open(fh,file=TRIM(filename),status='unknown', form='formatted')
  WRITE(fh, *) 'nat, iso, norbitals, ncix, natom'
  WRITE(fh, *) nat, iso, norbitals, ncix, natom
  WRITE(fh, *) 'nkpt, nmat, nume'
  WRITE(fh, *) nkpt, nmat, nume
  WRITE(fh, *) 'Qcomplex'
  WRITE(fh, *) Qcomplex
  WRITE(fh, *) 'lmax2, maxdim2, maxdim, maxsize'
  WRITE(fh, *) lmax2, maxdim2, maxdim, maxsize
  WRITE(fh, *) 'nindo'
  WRITE(fh, *) (nindo(iorb),iorb=1,norbitals)
  WRITE(fh, *) 'cidim'
  WRITE(fh, *) (cixdim(icix),icix=1,ncix)
  WRITE(fh, *) 'nl'
  WRITE(fh, *) (nl(icase), icase=1,natom)
  WRITE(fh, *) 'll, iorbital, cix, nind'
  do icase=1,natom
     do lcase=1,nl(icase)
        WRITE(fh, *) ll(icase,lcase)
        WRITE(fh, *) iorbital(icase,lcase)
        WRITE(fh, *) cix(icase,lcase)
        WRITE(fh, *) (2*ll(icase,lcase)+1)*iso          ! nind
     enddo
  enddo
  
  WRITE(fh, *) 'csize'
  WRITE(fh, *) (csize(icix), icix=1,ncix)
  
  WRITE(fh, *) 'iSx'
  do iorb=1,norbitals
     do ip1=1,nindo(iorb)
        WRITE(fh, *) iSx(ip1,iorb)
     enddo
  enddo
  WRITE(fh, *) 'Sigind'
  DO icix=1,ncix
     DO ip=1,cixdim(icix)
        DO iq=1,cixdim(icix)
           WRITE(fh, *) Sigind(ip,iq,icix)
        ENDDO
     ENDDO
  ENDDO
  WRITE(fh, *) 'EF, VOL'
  WRITE(fh, *) EF, VOL
END SUBROUTINE PrintSomeArrays

SUBROUTINE RenormalizeTrans(DMFTU, Olapm0, cix_orb, nindo, iSx, nbands, maxdim2, norbitals, maxdim, ncix)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: DMFTU(nbands,maxdim2,norbitals)
  COMPLEX*16, intent(in)    :: Olapm0(maxdim,maxdim,ncix)
  INTEGER, intent(in)       :: cix_orb(norbitals), nindo(norbitals), iSx(maxdim2,norbitals)
  INTEGER, intent(in)       :: nbands, maxdim2, norbitals, maxdim, ncix
  ! locals
  INTEGER :: iorb1, icix, nind1, ind1
  REAL*8  :: olocef

  DO iorb1=1,norbitals
     icix = cix_orb(iorb1)
     if ( icix.EQ.0 ) CYCLE
     nind1 = nindo(iorb1)
     
     do ind1=1,nind1
        olocef = 1/sqrt(real(Olapm0( iSx(ind1,iorb1), iSx(ind1,iorb1), icix )))
        DMFTU(:,ind1,iorb1) = DMFTU(:,ind1,iorb1) * olocef
     enddo
     
  ENDDO
END SUBROUTINE RenormalizeTrans
