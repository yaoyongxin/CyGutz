SUBROUTINE cmp_overlap(projector,Olapm, Qcomplex, nsymop, csort, iorbital, cix_orb, nindo, cixdim, &
                      &iSx, noccur, cfX, Rspin, maxucase, maxdim2, norbitals,pr_proc)
  USE com_mpi, ONLY: nprocs, myrank, master, AllReduce_MPI,vectors, nvector, vector_para, fvectors
  USE com,     ONLY: MINWAV, MAXWAV, iso, emin, emax
  USE abc,     ONLY: KX, KY, KZ, BKX, BKY, BKZ, E, A, ALM, ALML
  USE param,   ONLY: nmat , LMAX2, iblock, nume, nrf, NRAD, LOMAX, nloat, fastFilesystem
  USE w_atpar, ONLY: w_FJ, w_DFJ, w_jatom, w_alo, w_nlo, w_nlov, w_nlon, w_ilo, w_loor, w_lapw, w_ri_mat, w_P, w_DP
  USE struct,  ONLY: RMT, VOL, pos, tau, nat, rotij, tauij
  USE case,    ONLY: iatom, isort, nl, ll, crotloc, shft, maxdim, ncix, natom, Sigind, csize, maxsize, cix, ifirst
  USE sym2,    ONLY: tmat, phase, idet
  USE kpts,    ONLY: mweight, numkpt
  IMPLICIT NONE
  INTEGER, intent(in)     :: projector
  COMPLEX*16, intent(out) :: Olapm(maxdim,maxdim,ncix)
  LOGICAL, intent(in)     :: Qcomplex
  INTEGER, intent(in)     :: nsymop, csort(nat)
  INTEGER, intent(in)     :: iorbital(natom,lmax2+1), cix_orb(norbitals), nindo(norbitals)
  INTEGER, intent(in)     :: cixdim(ncix), iSx(maxdim2, norbitals), noccur(maxsize,ncix)
  COMPLEX*16, intent(in)  :: cfX(maxdim2,maxdim2,norbitals,norbitals), Rspin(2,2,norbitals)
  INTEGER, intent(in)     :: maxucase, maxdim2, norbitals
  ! common blocks
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
  ! locals
  INTEGER :: ikp, iikp, iks, nkp, ivector, N, NE, NEMIN, NEMAX, nbands, NUM, isym, icase, lcase, lfirst
  INTEGER :: pr_proc, latom, jatom, iucase, i, IND_YL, ibb, lda, ldb, ldc, l1
  INTEGER :: nind, idt, icix, iorb, m1, lms1, ind1, ind2, is1, iorb1, iorb2, nind1, nind2
  INTEGER :: cixdm, ip, iq, it, is, itape, L, M, irf, ii, i3, N1, num1, irf1, num2, is2, m2, lms2, irf2
  logical :: Tcompute
  character*10:: KNAME
  real*8      :: S, T, Z, exxx, FAC, PI, TWOPI, ARG1, ARG2, ARG3, ARGT, EMIST, cor2, fct
  complex*16  :: h_yl(2*LMAX2+1,iblock)
  complex*16  :: h_alyl(2*LMAX2+1,iblock,2)
  complex*16  :: h_blyl(2*LMAX2+1,iblock,2)
  complex*16  :: YL((LMAX2+1)*(LMAX2+1))
  complex*16  :: PHSHEL, CFAC, PH_SPIN(2), IMAG, cc, csum, corc, sm1, sm2, ff
  real*8      :: Olapc(maxsize,ncix), ARGT2
  real*8, allocatable     :: a_real(:)
  complex*16, allocatable :: DMFTrans(:,:,:,:,:), DMFTU(:,:,:), olp(:,:)!, cft(:,:), DMFTrans(:,:,:,:,:)
  complex*16, allocatable :: Olapmk(:,:,:) !, Olapm(:,:,:)
  !complex*16, allocatable :: Olapm2(:,:), Olapmk2(:,:)
  COMPLEX*16, allocatable :: URx(:,:), tmp(:,:) !, Rx(:,:)
  COMPLEX*16 :: tuu, tud, tdu, tdd
  COMPLEX*16, allocatable :: Uu(:,:), Ud(:,:)
  INTEGER :: no1
  CHARACTER*200 :: FNAME
  !------------------------------------------------------------------
  DATA IMAG/(0.0D0,1.0D0)/
  
  PI=ACOS(-1.0D0)
  TWOPI=2.D0*PI

  ALLOCATE( a_real(nmat) ) !---  for eigenvectors -------------!
  allocate( olp(maxdim2,maxdim2) )
  allocate( Olapmk(maxdim,maxdim,ncix) )
  Olapm=0

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
     else
        nkp= numkpt
     endif

     DO iks=1,nkp   !------ Over all irreducible k-points ------!
        Tcompute=.TRUE.
        if (vector_para) then
           ikp = vectors(ivector,3)+iks  ! successive index in k-point table from case.klist
        else
           ikp = iks                     ! successive index in k-point table from case.klist
           !--- We need to go over all k-points even though we will compute only some of them on this processor.
           !--- This is because we need to read vector file sequentially.
           iikp = ikp-myrank*pr_proc            ! The index of the point to compute. If negative, do not compute!
           if (iikp.GT.pr_proc) EXIT            ! Processor finished. The rest of the points will be taken care of by the other processors.
           if (iikp.LE.0)Tcompute=.FALSE.
        endif
   
        !------- reading from vector for both spins -----------!
        DO is=1,iso    !------ over up/dn ---------------------!
           itape=8+is
           READ(itape,END=998) S,T,Z,KNAME,N,NE  !--- Here we can jupm out of loop 4 and stop at 998 -----------------------------------!
           IF(N.GT.MAXWAV) MAXWAV=N                                         
           IF(N.LT.MINWAV) MINWAV=N                                         
           READ(itape) (KX(I),KY(I),KZ(I),I=1,N) !--- Reads all reciprocal vectors -----------------------------------------------------!
           NEMIN=1
           NEMAX=0
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
              IF(E(NUM).LT.EMIN) NEMIN=NEMIN+1
              IF(E(NUM).LT.EMAX) NEMAX=NEMAX+1
           ENDDO
   
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
   
        if (.not.Tcompute) CYCLE  ! This k-points was read, but will not be computed on this processor
        !---  Finished reading eigenvectors and eigenvalues for both spins -----!
        !---   E(iband) contains eigenvalues   ---------------------------------!
        !---   A(:,iband,is) contains eigenvectors -----------------------------!
        
        nbands = nemax-nemin+1
        allocate( DMFTU(nbands,maxdim2,norbitals) ) 
        allocate( URx(nbands,maxdim2), tmp(nbands,maxdim2) )
   
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
                 L=ll(icase,lcase) !------ current L --!
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
                       CALL ROTATE (BK,TMAT(1,1,isym),BKROT)  ! apply one of the symmetry operations to BKROT=R.(k+K)
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
                       ! ARG1 + ARG2 + ARG3 = (R_g.(k+K)) *  R(iatom) * 2pi
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
                       END DO
                    enddo
                    
                    DO is=1,iso  !--- over both spins
                       i3=0
                       do  i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                          i3=i3+1
                          DO M=1,2*L+1
                             if (lapw(l)) then
                                h_ALYL(m,i3,is)=(w_DFJ(L,I,iucase)*P(l,is,2)-w_FJ(L,I,iucase)*DP(l,is,2))* h_yl(M,i3)*ph_spin(is) ! derivatives of bessel functions and spheric harmonics
                                h_BLYL(m,i3,is)=(w_FJ(L,I,iucase)*DP(l,is,1)-w_DFJ(L,I,iucase)*P(l,is,1))* h_yl(M,i3)*ph_spin(is) 
                             else
                                h_ALYL(m,i3,is)=w_FJ(L,I,iucase)/P(l,is,1)/RMT(jatom)**2*h_yl(M,i3)*ph_spin(is)
                                h_BLYL(m,i3,is)=(0.d0,0.d0)
                             end if
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
                    DO NUM=NEMIN,NEMAX 
                       DO is=1, iso    
                          DO irf=1,nrf 
                             ALML(L,M,NUM,irf,is)=ALM(M,NUM,irf,is)*FAC*CFAC !--------- FACCFAC = (-1)^l*4.*PI*RMT(JATOM)**2/SQRT(VOL) --------!
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
                 
                 if (icix.gt.0) then
                    URx=0
                    N1=2*l1+1
                    if (projector.eq.1) then
                       ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U
                       do num1=1,nbands      !nemin,nemax     ! over bands
                          do is1=1,iso            ! over spin-1
                             do m1=-l1,l1           ! over m-1
                                lms1=l1+1+m1
                                ind1=l1+1+idt*m1+N1*(is1-1)
                                cc=0
                                do irf1=1,nrf
                                   cc = cc + alml(l1,lms1,num1+nemin-1,irf1,is1)*ri_mat(l1,l1,irf1,1,is1,is1)*(idt)**m1
                                enddo
                                
                                if (idt.gt.0) then
                                   URx(num1,ind1) = cc
                                else
                                   URx(num1,ind1) = conjg(cc)
                                endif
                             enddo
                          enddo
                       enddo
                    elseif (projector.eq.2) then
                       ! L -> l1
                       ! nrfmax -> nrf
                       
                       ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U
                       do num1=1,nbands      !nemin,nemax     ! over bands
                          do is1=1,iso            ! over spin-1
                             do m1=-l1,l1           ! over m-1
                                lms1=l1+1+m1
                                ind1=l1+1+idt*m1+N1*(is1-1)
                                cc=0
                                sm1=0
                                sm2=0
                                do irf1=1,nrf
                                   do irf2=1,nrf
                                      ff = alml(l1,lms1,num1+nemin-1,irf1,is1)*conjg(alml(l1,lms1,num1+nemin-1,irf2,is1))
                                      sm1 = sm1 + ff * ri_mat(l1,l1,irf1,irf2,is1,is1)
                                      sm2 = sm2 + ff * ri_mat(l1,l1,irf1,1,is1,is1) * ri_mat(l1,l1,1,irf2,is1,is1)
                                   enddo
                                   cc = cc + alml(l1,lms1,num1+nemin-1,irf1,is1)*ri_mat(l1,l1,irf1,1,is1,is1)*(idt)**m1
                                enddo
                                cc = cc * sqrt(abs(sm1/sm2))
                                if (idt.gt.0) then
                                   URx(num1,ind1) = cc
                                else
                                   URx(num1,ind1) = conjg(cc)
                                endif
                             enddo
                          enddo
                       enddo
                    elseif (projector.eq.3) then
                       ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U
                       do num1=1,nbands      !nemin,nemax     ! over bands
                          do is1=1,iso            ! over spin-1
                             do m1=-l1,l1           ! over m-1
                                lms1=l1+1+m1
                                ind1=l1+1+idt*m1+N1*(is1-1)
                                
                                irf1 = 1
                                cc = alml(l1,lms1,num1+nemin-1,irf1,is1)*sqrt(ri_mat(l1,l1,irf1,irf1,is1,is1))*(idt)**m1
                                
                                ! Adding correction due to \cdot{u} and local orbitals
                                ! The results is  A(i,L,irf=1)*sqrt(o(L)) * sqrt( 1 + \sum_{irf>1} |A(i,L,irf)|^2*o(L,irf))
                                if ( abs(cc)>1e-6 ) then
                                   corc = 0
                                   do irf1=1,nrf
                                      do irf2=1,nrf
                                         if (irf1.eq.1 .and. irf2.eq.1) cycle
                                         corc = corc + alml(l1,lms1,num1+nemin-1,irf1,is1)*conjg(alml(l1,lms1,num1+nemin-1,irf2,is1))*ri_mat(l1,l1,irf1,irf2,is1,is1)
                                      enddo
                                   enddo
                                   cc = cc*sqrt(1 + abs(corc)/abs(cc)**2)
                                endif
                                if (idt.gt.0) then
                                   URx(num1,ind1) = cc
                                else
                                   URx(num1,ind1) = conjg(cc)
                                endif
                             enddo
                          enddo
                       enddo
                    else
                       print *, 'Only projector=[1,2,3] is allowed!'
                       stop
                    endif
   
                    nind1 = nind
                    iorb1=iorb
   
                    if (iso.eq.2) then
                       no1 = nind1/iso
                       ALLOCATE( Uu(nbands,no1), Ud(nbands,no1))
                    endif
   
                    DO iorb2=1,norbitals
                       if ( cix_orb(iorb2).NE. icix ) CYCLE
                       nind2 = nindo(iorb2)
   
                       if (iso.eq.2) then
                          tuu = Rspin(1,1,iorb1)
                          tud = Rspin(1,2,iorb1)
                          tdu = Rspin(2,1,iorb1)
                          tdd = Rspin(2,2,iorb1)
                          Uu(:,:) = URx(:,1:no1)
                          Ud(:,:) = URx(:,no1+1:2*no1)
                          URx(:,1:no1)       = Uu*tuu + Ud*tdu
                          URx(:,no1+1:2*no1) = Uu*tud + Ud*tdd
                       end if
   
                       call zgemm('N','C', nbands, nind2, nind1, (1.d0,0.d0), URx,nbands, cfX(:,:,iorb2,iorb1),maxdim2, (0.d0,0.d0), tmp,nbands)
                       DMFTU(:,:nind2,iorb2) = DMFTU(:,:nind2,iorb2) + conjg(tmp(:,:nind2))
                    ENDDO
   
                    if (iso.eq.2) then
                       DEALLOCATE( Uu, Ud )
                    endif
                    
                 endif
              enddo  ! over L-case
           ENDDO     ! over atom-cases
   
           Olapmk=0
           DO iorb1=1,norbitals
              icix = cix_orb(iorb1)
              if ( icix.EQ.0 ) CYCLE
              nind1 = nindo(iorb1)
              DO iorb2=1,norbitals
                 if ( icix.NE.cix_orb(iorb2) ) CYCLE
                 nind2 = nindo(iorb2)
                 call zgemm('C','N', nind1, nind2, nbands, (1.d0,0.d0), DMFTU(:,:,iorb1),nbands, DMFTU(:,:,iorb2),nbands, (0.d0,0.d0), olp,maxdim2)
                 do ind1=1,nind1
                    do ind2=1,nind2
                       Olapmk( iSx(ind1,iorb1), iSx(ind2,iorb2), icix ) = olp(ind1,ind2)
                    enddo
                 enddo
                       
              enddo
           ENDDO
           Olapm(:,:,:) = Olapm(:,:,:) + Olapmk(:,:,:) * (mweight(ikp)/nsymop) ! proper weight for the reducible k-point
        ENDDO !------  isym: over star of the irreducible k-point ----!
        
        !deallocate( DMFTrans )
        deallocate( DMFTU)
        deallocate( URx, tmp)
        

     WRITE(*,'(I3,A,1x,I3,1x,A,I4)') myrank, ') Finished k-point number', ikp, 'with #bands=', nbands
  ENDDO; ENDDO !---- over reducible k-point: ikp

998 CONTINUE !---- irreducible k-points end (jump) from reading
  !--- end k-points

     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     ENDDO

  CALL AllReduce_MPI(Olapm, maxdim, ncix)

  write(*,100)Olapm
100 format("Olapm0",100F8.4)

  Olapc=0
  DO icix=1,ncix
     cixdm = cixdim(icix)
     !---- Olap to vector form, and s_oo to matrix form
     DO ip=1,cixdm
        do iq=1,cixdm
           it = Sigind(ip,iq,icix)
           if (it.gt.0) then
              Olapc(it, icix) =  Olapc(it, icix) + real(Olapm(ip,iq,icix))
           endif
        enddo
     ENDDO
     DO it=1,csize(icix)
        Olapc(it, icix) = Olapc(it, icix)/noccur(it,icix)
     ENDDO
  ENDDO
  
  if (myrank.eq.master) then
     WRITE(*,*) 'Renormalizing Gloc to account for the interstitials'
     WRITE(*,*) 'Z due to missing-interstitials=', Olapc
  endif

  DEALLOCATE( a_real ) !---  for eigenvectors -------------!
  deallocate( olp )
  deallocate( Olapmk )
  
END SUBROUTINE cmp_overlap
