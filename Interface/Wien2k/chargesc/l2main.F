SUBROUTINE L2MAIN(cform,nnlo,coord,zz,so,NSPIN1,vresp_write,sumw,nwave,recomputeEF,WL,TCLM,TCLM_w,TFOUR,TFOUR_w)
  USE defs,  ONLY: PI, ZERO, TWO, ZEROC, IMAG
  USE param, ONLY: IBLOCK, LMAX2, LOMAX, NCOM, NLOAT, NRAD, NGAU, NMAT, NUME, NSYM, NKPT, fh_vec, fastFilesystem
  USE atspdt,ONLY: EL, P, DP, PE, DPE, PEI, e_store
  USE char,  ONLY: modus
  USE com,   ONLY: rel, NSPIN, NAT, NBAND, NK, MINWAV, MAXWAV, NB, EMIN, EF, ELECN, XWT, weigh
  USE charp, ONLY: init_charp, zero_charp, fini_charp
  USE chard, ONLY: init_chard, zero_chard, fini_chard
  USE charf, ONLY: init_charf, zero_charf, fini_charf
  USE lo,    ONLY: loor, lapw, nlo, nlov, nlon, ilo, a1lo, b1lo, elo_store, alo, blo, clo, pi12lo, pe12lo, pr12lo
  USE lohelp,ONLY: u21, ue21, u12, ue12, u22, init_lohelp, zerotc_lohelp, fini_lohelp, sum12, sum21, sum22, sume21, sume12
  USE struk, ONLY: iatnr, mult, isplit, jri, pos, rotij, tauij, rmt, v, rotloc, r0, dx, vol, br1, br2
  USE xa,    ONLY: LM, R, AVEC, BK, BKROT, BKRLOC, fj, dfj, E, WEIGHT, TC100, TCA100, TCB100, init_xa, fini_xa, SUMA, SUMB, SUMAB, SUMBA
  USE xa3,   ONLY: bkx, bky, bkz, bkxlo, bkylo, bkzlo, init_xa3, fini_xa3, kx, ky, kz, kxlo, kylo, kzlo, As, As_lo!, A, A_lo
  USE reclat,ONLY: inst, tauk, kzz
  USE com_mpi,ONLY: nprocs,myrank,master,FilenameMPI,Reduce_MPI,Gather_procs,FindMax_MPI,cpuID,pr_proc,pr_procr,pr_procs,vector_para,vectors,nvector,fvectors,VECFN,FindMaxK_MPI
  USE w_atpar
  USE dmf
  USE GMPI
  use ghdf5_base
  IMPLICIT NONE
  CHARACTER*4, intent(in)  :: cform
  INTEGER, intent(out)     :: nnlo
  CHARACTER*5, intent(in)  :: coord
  REAL*8, intent(in)       :: ZZ(*), sumw, WL
  LOGICAL, intent(in)      :: SO
  INTEGER, intent(in)      :: NSPIN1, nwave
  LOGICAL, intent(inout)   :: vresp_write
  INTEGER, intent(in)      :: recomputeEF
  REAL*8, intent(out)      :: TCLM,TCLM_w,TFOUR,TFOUR_w
  ! Functions
  INTEGER :: NOTRI
  REAL*8  :: GAUNT
  INTEGER :: CountSelfenergy
  REAL*8  :: ferm
  ! Common blocks
  COMMON /FACT/   FCT(100)
  REAL*8     :: FCT
  COMMON /RADFU/  RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2),RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)
  REAL*8     :: RRAD1, RADE1, RRAD2, RADE2
  COMMON /UHELP/  UA(NRAD),UB(NRAD),UBA(NRAD),UAB(NRAD)
  REAL*8     :: UA, UB, UBA, UAB
  common /xmean/  xwteh(0:3),xwtel(0:3),xwt1h(0:3),xwt1l(0:3) 
  REAL*8     :: xwteh, xwtel, xwt1h, xwt1l
  common /mkef/delef,ts2
  REAL*8     :: delef, ts2
  ! Locals
  real*8, PARAMETER       :: Ry2eV = 13.60569193
  !
  LOGICAL       :: debug, testLDA
  CHARACTER*10  :: KNAME
  CHARACTER*150 :: parseline
  LOGICAL       :: more_kpoints
  COMPLEX*16    :: YL((LMAX2+1)*(LMAX2+1))
  COMPLEX*16    :: PHSHEL,CFAC,IMAG1
  REAL*8        :: tc_buf(nrad),vtc_buf(nrad)
  CHARACTER*100 :: CDUMMY
  COMPLEX*16    :: cc, gtc, CCC
  INTEGER       :: I, J, N, JK, ISCF, LFIRST, LMX, jatom, ilm1, ilm3, n0, ne, nm
  INTEGER       :: jlm, lmmax, imax, nemin, nemax, DM_nemin, DM_nemax, DM_nemaxx, label, i1
  INTEGER       :: ii, ibb, i3, isize, latom, mu, index, l, m, maxx, lda, ldb, ldc, num, num0, ilm2
  INTEGER       :: jlo, ilm, li, mi, L1, MMAX, LP, LP1, mpmax, jlop, MP, MS, MPS, MTEST, LY, LPY
  REAL*8        :: arg, arg1, arg2, arg3, argt, argt2, emist, fac, rmt2, GNT
  REAL*8        :: TWOPI, TEST1, ETOT, ETOT2, EDL,SQRT2, SQRT3, SQFP, Y, CIN, jatombad, lbad
  !
  REAL*8     :: time_bl, time_bl_w, time_reduc, time_reduc_w, time_write, time_writeclm, time_writescf, time4_w, time3, time4
  REAL*8     :: time_write_w, time_writeclm_w, time_ilm, time_ilm_w, time_radprod, time_m, time_rad, time3_w
  REAL*8     :: time_rd_w, time_rd_c, time_atpar_w, time_atpar_c, time_writescf_w, time_r_w, time_r_c, trw, trc
  REAL*8     :: time_writeefg, time_writeeefg_w, time_writeefg_w, t1c, t1w, t2c, t2w, time1, time1_w, time2, time2_w
  REAL*8     :: time_dmf, time_dmfw, time_int, time_intw, time_ed, time_edw, time_dmf0, time_dmf0w, time_dmf1, time_dmf1w
  !
  COMPLEX*16 :: gint, wkp, w0, p_w0
  !COMPLEX*16 :: dh_alm(nume,(LMAX2+1)*(LMAX2+1)), dh_blm(nume,(LMAX2+1)*(LMAX2+1)), dh_clm(nume,(LMAX2+1)*(LMAX2+1),nloat)
  COMPLEX*16 :: tsa, tsb, tsab, tsba, ts12(nloat), tse12(nloat), ts21(nloat), tse21(nloat), ts22(nloat,nloat)
  COMPLEX*16 :: tsuma,tsumb,tsumab,tsumba,tsum21(nloat),tsume21(nloat),tsum12(nloat),tsume12(nloat),tsum22(nloat,nloat)
  COMPLEX*16 :: SA, SB, SAB, SBA, s12(nloat), se12(nloat), s21(nloat), se21(nloat), s22(nloat,nloat)
  COMPLEX*16 :: vsa, vsb, vsab, vsba, vs12(nloat), vse12(nloat), vs21(nloat), vse21(nloat), vs22(nloat,nloat)
  REAL*8     :: facv, EQBAD, OVRDIV, TC, EFGFACT, EFG20, EFG22, EFG2M, QXX, QYY, QZZ, QXY, VXX, VYY, VZZ
  INTEGER    :: LM_MAX, iks,ikp, is, itape,itape1,itape2, ios, icase, lcase, iw, nind, icix, num1, m1, ind, iom, max_nbands
  REAL*8     :: S,T,Z, vnorm1, BKROT2(3), BKROT3(3), zero3(3), ff, twgh, emax, wgamma
  integer    :: iorbital(natom,lmax2+1), norbitals, maxdim2, iorb, nbands, nbandsx, it, nomega, npomega, fh_sig, im
  real*8     :: h_al(iblock),h_bl(iblock)
  complex*16 :: h_yl((LMAX2+1)*(LMAX2+1),iblock), h_alyl((LMAX2+1)*(LMAX2+1),iblock), h_blyl((LMAX2+1)*(LMAX2+1),iblock)
  REAL*8     :: zw0(nume), LDA_weight(nume)
  REAL*8,     ALLOCATABLE :: zw2(:)
  COMPLEX*16, ALLOCATABLE :: alm(:,:),blm(:,:),clm(:,:,:), dh_alm(:,:,:), dh_blm(:,:,:), dh_clm(:,:,:,:)
  COMPLEX*16, ALLOCATABLE :: tmp(:,:), p_tmp(:,:), cfX(:,:,:,:), DMFTU(:,:,:), STrans(:,:,:,:), zw1(:), w1(:,:), Aweight(:,:), Olapm0(:,:,:), Wsum(:,:), Winf(:,:), Wferm(:,:)
  INTEGER,    ALLOCATABLE :: nindo(:), cix_orb(:), cixdim(:), iSx(:,:), noccur(:,:)
  complex*16, allocatable :: sigma(:,:,:), s_oo(:,:), Eij(:,:), zek(:), p_zek(:), gloc(:), gij(:,:),  zar(:,:), p_zar(:,:), zal(:,:), p_zal(:,:), zar_inf(:,:), zal_inf(:,:)
  real*8,     allocatable :: omega(:), zek_inf(:)
  ! interstitials
  INTEGER    :: keigen(3,nmat), kmax(3), iff1, iff2, iff3, iff1t, iff2t, iff3t, isig, ierr, ia1, ia2, jj, jx, ij, k
  REAL*8     :: volin, Tx, beta, ek_inf, dsum, drest, Om
  COMPLEX*16, ALLOCATABLE :: fft(:,:,:),tp(:,:,:),sumfft(:,:,:),vsumfft(:,:,:),tloc(:,:,:)
  DOUBLE PRECISION,ALLOCATABLE :: dwork(:)
  COMPLEX*16, ALLOCATABLE :: cwork(:)
  COMPLEX*16, ALLOCATABLE :: rho1(:),vrho1(:)
  COMPLEX*16, ALLOCATABLE :: rhok_buf(:),vrhok_buf(:)
  COMPLEX*16, ALLOCATABLE :: rhok(:),vrhok(:), Asc(:,:), w0inf(:), s_inf(:,:)
  COMPLEX*16 :: cek, omn, e0, e1, csum(nume)
  INTEGER :: qmax, j0, j1, j2, iomt, iikp, ivector
  LOGICAL :: Rho_Renormalize, Tcompute
  REAL*8,ALLOCATABLE  :: LowE(:,:,:)
  REAL(8) DM_EF,DM_EORB,DM_ESUM
  character str*27,fname*256
  
  Rho_Renormalize = .FALSE.
  testLDA=.FALSE.
  debug = .FALSE.
  wgamma = 1e-10
  
  fh_sig = 80
  emax = ef
  TWOPI=TWO*PI

  TEST1=0.0D0                                                        
  ETOT=0.0D0                                                        
  ETOT2=0.0; EDL=0
  XWT=0.0                                                           
  SQRT2=SQRT(2.0D0)                                                 
  SQRT3=SQRT(3.D0)                                                  
  SQFP=SQRT(4.D0*PI)                                                
  zero3(:) = 0.0

  CALL init_lohelp
  CALL init_xa(nat)
  
  vnorm1=1.d0
  IF((iso.EQ.2).AND.(nspin1.EQ.1)) vnorm1=0.5d0
  
  Y=1.0D0                                                           
  DO I=1,49
     J=2*I-1
     FCT(J)=Y  ! FCT(i) = ((i+1)/2)!
     Y=Y*I
  END DO
  
  CIN=1.d0/137.0359895d0**2                       
  IF (.NOT.REL) CIN=4.0*1.0D-22                                     
  
  READ(18,2032) ISCF                                                     

  if (myrank.eq.master) then
     WRITE(8,787) ISCF                        
     WRITE(8,*) '   NORM OF CLM(R) =          '                        
     WRITE(8,*)                                                 
     if (vresp_write)  then       
        WRITE(28,787) ISCF                        
        WRITE(28,*) '   NORM OF CLM(R) =          '                        
        WRITE(28,*)  
     endif
     WRITE(6,800)
  endif
  
  LFIRST=1
  nnlo=0
  time_bl=zero; time_bl_w=zero; time_reduc=zero; time_reduc_w=zero
  time_write=zero; time_writeclm=zero; time_writescf=zero; time_write_w=zero
  time_writeclm_w=zero; time_writescf_w=zero; time_ilm=zero; time_ilm_w=zero
  time_radprod=zero; time_m=zero; time_rad=zero; time_rd_w=zero; time_rd_c=zero
  time_r_w=zero; time_r_c=zero; time_atpar_w=zero; time_atpar_c=zero; time_writeefg=zero; time_writeefg_w=zero
  time_dmf=zero; time_dmfw=zero; time_int=zero; time_intw=zero; time_ed=zero; time_edw=zero
  
  ALLOCATE(alm((lmax2+1)*(lmax2+1),nume),blm((lmax2+1)*(lmax2+1),nume),clm((lmax2+1)*(lmax2+1),nume,nloat))
  ALLOCATE(dh_alm(nume,(lmax2+1)*(lmax2+1),iso),dh_blm(nume,(lmax2+1)*(lmax2+1),iso),dh_clm(nume,(lmax2+1)*(lmax2+1),nloat,iso)) ! Mar 27
  CALL init_charp(nume)
  CALL init_chard(nume)
  CALL init_charf(nume)
  
  LMX=LMAX2                                                    
  
  allocate (e_store(0:lmax2,nat),elo_store(0:lomax,1:nloat,nat))

  itape=30
  if (vector_para) then
     if (nvector.ge.1) then
        FNAME = fvectors(1,3)
     else   ! No k-point to compute, but still need linearization energies
        FNAME = TRIM(VECFN(4))//'_1'
     endif
     open(itape,FILE=FNAME,STATUS='old',FORM='formatted')
  else
     !------------ vector files need to be read from the beginning -------------------!
     rewind(itape)
  endif
  DO i=1,nat
     READ(itape,'(a)') parseline
     READ(parseline,'(100(f9.5))') e_store(0:lmax2,i)
     READ(itape,'(100(f9.5))') elo_store(0:lomax,1:nloat,i)
  enddo
  if (vector_para) then
     close(itape)
  else
     rewind(itape)
  endif
  
  CALL w_allocate0(nat)
  ! Reads case.in2 to determin LM array
  ! and maximum LM -> LM_MAX=max(LM(1,:))
  LM_MAX=0
  do jatom=1,nat
     ! neg L MEANS NEGATIVE SPHERICAL HARMONIC COMB (SEE KURKI-SUONIO)   
     read(5,1003) ( (lm(j,jlm),j=1,2), jlm=1,ncom )
     ! Find maximum LMMAX
     DO JLM=2,NCOM
        IF(LM(1,JLM).EQ.0) EXIT
     ENDDO
     LMMAX=JLM-1
     if (LMMAX.GT.LM_MAX) LM_MAX = LMMAX
     
     if (myrank.EQ.master .OR. fastFilesystem/=0) then
        write(6,*) ' atom',jatom,' ncomu',0,' lmmax',LMMAX
        WRITE(6, 1005) LMMAX,((LM(J,JLM),J=1,2), JLM=1,LMMAX)             
     endif
     if (myrank.EQ.master) then
        WRITE(21,13)   JATOM,IATNR(JATOM),(POS(I,lfirst),I=1,3),MULT(JATOM)
        WRITE(21,1005) LMMAX,((LM(J,JLM),J=1,2), JLM=1,lmmax)             
     endif
     w_lm(1:2,1:NCOM,jatom) = lm(1:2,1:NCOM)
     w_lmmax(jatom) = lmmax
  enddo
  
  ! allocates memory to store radial wave functions and 
  ! all results of subroutine "atpar", which is called outside
  ! the k - loop.
  CALL w_allocate(LM_MAX,nat)
  
  do jatom=1,nat
     ! Here we calculate all radial wave functions and other quantities
     ! computed by atpar, which do not depend on k vector.
     
     lfirst=1
     do i=1,jatom-1
        lfirst=lfirst + mult(i)
     enddo

     ! calculate radial functions U(r), UE(R), ...                       
     call atpar (rel,nat,jatom,lfirst,cform,zz(jatom))               
     
     IMAX=JRI(JATOM)
     DO I=1,IMAX
        R(I)=R0(JATOM)*EXP((I-1)*DX(JATOM)) ! Radial mesh
     ENDDO
     
     nnlo=nlo+nlon+nlov
     
     w_lfirst(jatom) = lfirst
     w_nlo(jatom) = nlo
     w_nlov(jatom) = nlov
     w_nlon(jatom) = nlon

     w_ilo(0:lmax2,jatom) = ilo(0:lmax2)
     w_loor(1:nloat,0:lomax,jatom) = loor(1:nloat,0:lomax)
     w_lapw(0:lmax2,jatom) = lapw(0:lmax2)
     w_alo(0:lomax,1:nloat,jatom) = alo(0:lomax,1:nloat)
     w_blo(0:lomax,1:nloat,jatom) = blo(0:lomax,1:nloat)
     w_clo(0:lomax,1:nloat,jatom) = clo(0:lomax,1:nloat)
     
     w_a1lo(1:nrad,1:nloat,0:lomax,jatom) = a1lo(1:nrad,1:nloat,0:lomax)
     w_b1lo(1:nrad,1:nloat,0:lomax,jatom) = b1lo(1:nrad,1:nloat,0:lomax)
     w_RRAD1(1:nrad,0:lmax2,jatom) = RRAD1(1:nrad,0:lmax2)
     w_RRAD2(1:nrad,0:lmax2,jatom) = RRAD2(1:nrad,0:lmax2)
     w_RADE1(1:nrad,0:lmax2,jatom) = RADE1(1:nrad,0:lmax2)
     w_RADE2(1:nrad,0:lmax2,jatom) = RADE2(1:nrad,0:lmax2) 
     
     w_P(0:LMAX2,jatom) = P(0:LMAX2)
     w_DP(0:LMAX2,jatom) = DP(0:LMAX2)
     w_PE(0:LMAX2,jatom) = PE(0:LMAX2)
     w_DPE(0:LMAX2,jatom) = DPE(0:LMAX2)
     w_PEI(0:LMAX2,jatom) = PEI(0:LMAX2)
     w_pi12lo(1:nloat,0:lomax,jatom) = pi12lo(1:nloat,0:lomax)
     w_pe12lo(1:nloat,0:lomax,jatom) = pe12lo(1:nloat,0:lomax)
     w_pr12lo(1:nloat,1:nloat,0:lomax,jatom) = pr12lo(1:nloat,1:nloat,0:lomax)
     
     w_R(1:nrad,jatom) = R(1:nrad)
  enddo

  !---------------------------------
  ! START LOOP OVER ALL k-POINTS
  !---------------------------------
  CALL init_xa3(nmat,nnlo,nume,iso)
  if (vector_para) then
     pr_proc = sum(vectors(:,2))
     pr_procr = pr_proc
     if (myrank.EQ.master .or. FastFilesystem/=0) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'tot-k=', nkpt
  else
     pr_proc  = floor(nkpt/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor
     ! correcting pr_proc is the number of k-points is not dividable
     pr_procr = pr_proc
     if ((myrank+1)*pr_proc .gt. nkpt) then
        if (nkpt-myrank*pr_proc.gt.0) then
           pr_procr = nkpt-myrank*pr_proc
        else
           pr_procr = 0
        endif
     endif
     if (myrank.EQ.master .OR. fastFilesystem/=0) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'pr_procr=', pr_procr, 'tot-k=', nkpt
  endif

  ALLOCATE( pr_procs(nprocs) )
  CALL Gather_procs(pr_procr, pr_procs, nprocs)
  if (myrank.EQ.master .and. SUM(pr_procs).NE.nkpt) then
     WRITE(6,*) 'ERROR: sum(pr_procs) should be nkpt, but is not', sum(pr_procs), nkpt
  endif

  max_nbands=0
  
  ! ---- for interstitial ----
  kmax(:)=0
  ! ..... FIND MAX RECIPROCAL LATTICE VECTORS

  DO ivector=1,nvector
     if (vector_para) then
        open(fh_vec,FILE=fvectors(ivector,1),STATUS='old',FORM='unformatted')
     else
        rewind(fh_vec) !--- both vector files: 9 and 10 rewind -------------------------!
     endif
     DO I=1,NAT
        READ(fh_vec) EMIST !---- At the beginninge we have linearization energies --!
        READ(fh_vec) EMIST !---- We just skip them ---------------------------------!
     ENDDO
     READ(fh_vec) s,t,z,kname,n,ne
     READ(fh_vec) (keigen(1,i),keigen(2,i),keigen(3,i),i=1,n)
     if (vector_para) then
        close(fh_vec)
     else
        rewind(fh_vec)
     endif
     DO i=n-nnlo,1,-1
        kmax(1)=MAX(keigen(1,i),kmax(1))
        kmax(2)=MAX(keigen(2,i),kmax(2))
        kmax(3)=MAX(keigen(3,i),kmax(3))
     ENDDO
  ENDDO
  CALL FindMaxK_MPI(kmax)

  IF (myrank.EQ.master .OR. fastFilesystem/=0) WRITE(6,*) 'max k indices:',kmax(1),kmax(2),kmax(3)
  kmax(1)=kmax(1)+1
  kmax(2)=kmax(2)+1
  kmax(3)=kmax(3)+1
  ! -----  set iff values for FFT calculation
  iff1=2*kmax(1)
  iff2=2*kmax(2)
  iff3=2*kmax(3)
  iff1=(iff1+1)*2
  iff2=(iff2+1)*2
  iff3=(iff3+1)*2
  iff1t=iff1
  iff2t=iff2
  iff3t=iff3
  CALL IFFLIM(iff1t,iff1)
  CALL IFFLIM(iff2t,iff2)
  CALL IFFLIM(iff3t,iff3)
  IF (myrank.EQ.master .OR. fastFilesystem/=0) WRITE(6,*) 'n,iff1,iff2,iff3',n-nnlo,iff1,iff2,iff3
  !
  ALLOCATE (fft(iff1,iff2,iff3),tp(iff1,iff2,iff3),sumfft(iff1,iff2,iff3),vsumfft(iff1,iff2,iff3),tloc(iff1,iff2,iff3))
  ALLOCATE (dwork(4*(iff1+iff2+iff3)+15),cwork(iff1+iff2+iff3))
  fft=0.0
  tp=0.0
  sumfft=0.0
  vsumfft=0.0
  tloc=0.0
  ! --- finished for interstitial

  w_vRHOLM(:,:,:)=0.0
  w_RHOLM(:,:,:)=0.0
  w_xwt1(:,:)=0.0
  w_xwt1l(:,:)=0.0
  w_xwt1h(:,:)=0.0
  w_xwteh(:,:)=0.0
  w_xwtel(:,:)=0.0
     
  REWIND 5
  READ(5,fmt='(//,1A)') CDUMMY
  
  rewind(30)
  rewind(18)
  READ(18,2032) ISCF                                                     
  
!!! ---------- Preparation of arrays for paralel executaion --------------                                                     
  IF (myrank.EQ.master .OR. fastFilesystem/=0) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'tot-k=', nkpt-1
  IF (myrank.EQ.master .OR. fastFilesystem/=0) print *, 'mpr_proc=', pr_proc, 'myrank=', myrank, 'nprocs=', nprocs, 'nkpt=', nkpt
  NB(:)=0

  call gh5_init()
  
  iikp=0
  DO ivector=1,nvector
     IF(GP%LKPVEC)THEN
       WRITE(STR,'(I3)')GP%KVEC(IVECTOR,1)
     ELSE
       WRITE(STR,'(I3)')GP%MYRANK
     ENDIF
     call gh5_open_r('KSWT_'//TRIM(ADJUSTL(STR))//'.h5', f_id)
     call gh5_read(DM_EF,'/e_fermi',f_id)
     call gh5_read(DM_EORB,'/e_gamma_dc',f_id)
     call gh5_read(DM_ESUM,'/e_band',f_id)

     IF(IVECTOR.EQ.1)THEN
       WRITE(21,'(":FER DMFT Fermi Energy",F10.5)')DM_EF
       WRITE(21,'(":EORB ",F20.9)')DM_EORB
       WRITE(21,'(":SUM SUM OF DMFT BAND ENERGIES ",F20.9)')DM_ESUM
     ENDIF

     DO is=1,iso    !------ over up/dn ---------------------!
        itape1=8+is
        itape2=29+is
        if (vector_para) then
           open(itape1,FILE=fvectors(ivector,is),STATUS='old',FORM='unformatted')
           open(itape2,FILE=fvectors(ivector,is+2),STATUS='old',FORM='formatted')
        else
           rewind(itape1) !--- both vector files: 9 and 10 rewind -------------------------!
           rewind(itape2) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        DO I=1,NAT
           READ(itape1) EMIST !---- At the beginninge we have linearization energies --!
           READ(itape1) EMIST !---- We just skip them ---------------------------------!
        ENDDO
     END DO
     
     DO iks=1,vectors(ivector,2) ! kpoint loop begin
        if (vector_para) then
           ikp = vectors(ivector,3)+iks  ! successive index in k-point table from case.klist
           iikp = iikp+1                 ! successive index in k-point table on this processor
        else
           ikp = iks                     ! successive index in k-point table from case.klist
           !--- We need to go over all k-points even though we will compute only some of them on this processor.
           !--- This is because we need to read vector file sequentially.
           iikp = ikp-myrank*pr_proc            ! The index of the point to compute. If negative, do not compute!
        endif
        Tcompute=.FALSE.                                                                                                                                      
        if (iikp.gt.0) Tcompute=.TRUE.       ! If Tcompute is true, the point needs to be computed.                                                           
        if (iikp.GT.pr_proc) EXIT            ! Processor finished. The rest of the points will be taken care of by the other processors.                      

        call cputim(t1c)
        call walltim(t1w)
        !!! Need to read both spins!
        CALL Read_Vec_Spin(ikp, E, twgh, As, As_lo, kx, ky, kz, kxlo, kylo, kzlo, bkx, bky, bkz, bkxlo, bkylo, bkzlo, nemin, nemax, more_kpoints, n0, emin, emax, iso, nmat, nume, nnlo)
   
        IF (.NOT.Tcompute) CYCLE
        IF (.NOT.more_kpoints) EXIT
   
        LDA_weight=0
        DO num=nemin,nemax
           IF(E(num).GE.emin.AND.E(num).LE.emax) THEN
              LDA_weight(num)=weigh(ikp,num)
              IF(ABS(LDA_weight(num)).GT.1.d-8) THEN
                 NB(ikp)=NB(ikp)+1
                 etot=etot+E(num)*LDA_weight(num)*vnorm1
              ENDIF
              IF(ABS(LDA_weight(num)).LE.1.d-8) nemax=nemax-1
           ENDIF
        ENDDO
        
        call cputim(t2c)
        call walltim(t2w)
        time_rd_c=time_rd_c+t2c-t1c
        time_rd_w=time_rd_w+t2w-t1w
   
        do jatom=1,nat
           isize=n0-nnlo
           CALL HARMON(isize,bkx,bky,bkz,lmax2,fj(:,:,jatom),dfj(:,:,jatom),rmt(jatom))
        enddo
   
        call cputim(t1c)
        call walltim(t1w)
 
        WRITE(STR,'(I7)')ikp
        call gh5_read(DM_nemin, '/IKP_'//trim(adjustl(str))//'/nemin',f_id)
        call gh5_read(DM_nemax, '/IKP_'//trim(adjustl(str))//'/nemax',f_id)
        nbands = DM_nemax-DM_nemin+1
   
        call cputim(t2c)
        call walltim(t2w)
        
        time_dmf=time_dmf+t2c-t1c
        time_dmfw=time_dmfw+t2w-t1w
        
        call cputim(t1c)
        call walltim(t1w)
        allocate( Aweight(nbands,nbands), zw2(nbands) )
        call gh5_read(Aweight,nbands,nbands,'/IKP_'//trim(adjustl(str))// &
                &'/KSWT',f_id)
        wkp = (twgh/sumw*2.0)
        zw0(:)=0
        DO num=nemin,DM_nemin-1
           zw0(num) = wkp ! These bands are fully filled
        ENDDO
        
        ! We use the following formula to implement the total energy:
        ! dE = \sum_{k,i} E_{k,i}^{LDA} * n_{i,i}^{DMFT}
        ! where E_{k,i} are Kohn-Sham eigenvalues
        ! and n_{i,i}^{DMFT}= T*sum_{iom}g_{k,ii}^{DMFT} is the DMFT density matrix in KS base
   
        DO i=nemin,DM_nemin-1
           ETOT2 = ETOT2 + E(i)*dble(wkp)*VNORM1
           EDL=EDL+E(i)*dble(wkp)*VNORM1
        ENDDO
        DO i=1,nbands
           ETOT2 = ETOT2 + E(i+DM_nemin-1)*dble(Aweight(i,i))*VNORM1
        ENDDO
        CALL DMFT_WEIGHTS(zw2, Aweight, nbands, debug)
        
        call cputim(t2c)
        call walltim(t2w)
        time_ed=time_ed+t2c-t1c
        time_edw=time_edw+t2w-t1w
        
        ! Weights for all bands, which are needed below
        weight=0
        weight(:DM_nemin-1) = zw0(:DM_nemin-1)
        weight(DM_nemin:DM_nemax) = zw2(:)
        
        
        DM_nemaxx=DM_nemax
        DO num=DM_nemin,DM_nemax
           IF(abs(zw2(num-DM_nemin+1)).LE.1.d-8) DM_nemaxx=DM_nemaxx-1
        ENDDO
        
        if (testLDA) DM_nemaxx = max(DM_nemaxx,nemax)
        ! we need to have alm up to DM_nemax
        nemax     = max(DM_nemax,nemax)
   
        do num=nemin,DM_nemaxx
           xwt=xwt+weight(num)*vnorm1
        enddo
   
        if (debug) then
           print *, 'xwt=', xwt, vnorm1
           print *, 'DM=', nemin, nemax, DM_nemin, DM_nemax, DM_nemaxx
           do num=nemin,nemax
              write(*,'(I3,1x,f20.10,1x)',advance='no') num, E(num)-EF
              write(*,'(f20.10,1x)', advance='no') LDA_weight(num)  !-(cc-conjg(cc))/(2*pi*IMAG)
              write(*,'(2f20.10,1x)') weight(num)
           enddo
        endif
        WRITE(*,'(I3,A,1x,I4,1x,I4,1x,A,I3,1x,A,I3,1x,A,I3,1x,A,I3,1x,A,I3,1x)') myrank, ') Finished k-point', ikp, iikp, 'with nemin=', nemin, 'nemax=', nemax, 'DM_nemin=', DM_nemin, 'DM_nemax=', DM_nemax, 'DM_nemaxx', DM_nemaxx
   
        deallocate( zw2 )
        
        if (TestLDA) then
           weight(:) = LDA_weight(:)
        endif
        
        nbandsx = DM_nemaxx-DM_nemin+1
        do jatom=1,nat
           ! Quantities computed by "atpar" which do not depend on k-point
           lfirst = w_lfirst(jatom)
           lm(1:2,1:NCOM) = w_lm(1:2,1:NCOM,jatom)
           lmmax = w_lmmax(jatom)
           nlo = w_nlo(jatom)
           nlov = w_nlov(jatom)
           nlon = w_nlon(jatom)
           ilo(0:lmax2) = w_ilo(0:lmax2,jatom)
           loor(1:nloat,0:lomax) = w_loor(1:nloat,0:lomax,jatom)
           lapw(0:lmax2) = w_lapw(0:lmax2,jatom)
           alo(0:lomax,1:nloat) = w_alo(0:lomax,1:nloat,jatom)
           blo(0:lomax,1:nloat) = w_blo(0:lomax,1:nloat,jatom)
           clo(0:lomax,1:nloat) = w_clo(0:lomax,1:nloat,jatom)
           a1lo(1:nrad,1:nloat,0:lomax) = w_a1lo(1:nrad,1:nloat,0:lomax,jatom)
           b1lo(1:nrad,1:nloat,0:lomax) = w_b1lo(1:nrad,1:nloat,0:lomax,jatom)
           RRAD1(1:nrad,0:lmax2) = w_RRAD1(1:nrad,0:lmax2,jatom)
           RRAD2(1:nrad,0:lmax2) = w_RRAD2(1:nrad,0:lmax2,jatom)
           RADE1(1:nrad,0:lmax2) = w_RADE1(1:nrad,0:lmax2,jatom)
           RADE2(1:nrad,0:lmax2) = w_RADE2(1:nrad,0:lmax2,jatom)
           P(0:LMAX2) = w_P(0:LMAX2,jatom)
           DP(0:LMAX2) = w_DP(0:LMAX2,jatom)
           PE(0:LMAX2) = w_PE(0:LMAX2,jatom)
           DPE(0:LMAX2) = w_DPE(0:LMAX2,jatom)
           PEI(0:LMAX2) = w_PEI(0:LMAX2,jatom)
           pi12lo(1:nloat,0:lomax) = w_pi12lo(1:nloat,0:lomax,jatom)
           pe12lo(1:nloat,0:lomax) = w_pe12lo(1:nloat,0:lomax,jatom)
           pr12lo(1:nloat,1:nloat,0:lomax) = w_pr12lo(1:nloat,1:nloat,0:lomax,jatom)
           R(1:nrad) = w_R(1:nrad,jatom)
              
           !.....CALCULATE  BESSELFUNCTIONS                                        
           !     YLM ARE CALCULATED FOR EACH ENERGY SEPARATLY TO CONSERVE CM-SPACE
           !     CALL HARMON(N-(nlo+nlon+nlov),bkx,bky,bkz,lmax2,fj,dfj,rmt(jatom))
           !isize=n0-nnlo
           !CALL HARMON(isize,bkx,bky,bkz,lmax2,fj(:,:,jatom),dfj(:,:,jatom),rmt(jatom))
           
           ! zero koeff. for charge analysis
           CALL zero_charp(nemin,nemax)
           CALL zero_chard(nemin,nemax)
           CALL zero_charf(nemin,nemax)
           CALL zerotc_lohelp(nemin,nemax)
           TC100(0:lmax2,nemin:nemax)=0.0
           TCA100(0:lmax2,nemin:nemax)=0.0
           TCB100(0:lmax2,nemin:nemax)=0.0
           
           ! CALCULATE ALM, BLM                                                
           fac=4.0d0*pi*rmt(jatom)**2/SQRT(vol)                              
           latom=lfirst-1                                                    
   
           equiv_atom_loop: DO mu=1,mult(jatom)
              
              latom=latom+1
   
              DO is=1,iso !!! over spin: Mar 25, 2010           
                 
                 alm(:,:)=0.0
                 blm(:,:)=0.0
                 clm(:,:,:)=0.0
   
                 CALL cputim(time1)
                 CALL walltim(time1_w)
              
                 blocked_loop: DO ii=1,n0-nnlo,iblock !iblock*ibpp,iblock
                    ibb=MIN(iblock,n0-(nlo+nlon+nlov)-ii+1)
                    if(ibb.le.0) EXIT
                    i3=0
                    DO i=ii,min(ii+iblock-1,n0-nnlo) !ii+iblock-1
                       !---------  rotates ylm(k+K) to ylm(k'+K) where k' is in star of irreducible k. ------------!
                       i3=i3+1
                       bk(1)=bkx(i) !-----  reciprocal vector and irreducible vector: G=K+k ----!
                       bk(2)=bky(i)
                       bk(3)=bkz(i)
                       ! BKROT = R_g.(k+K)
                       CALL ROTATE (bk,rotij(1,1,latom),bkrot)
                       !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
                       bk(1)=bkrot(1)*br1(1,1)+bkrot(2)*br1(1,2)+bkrot(3)*br1(1,3)   
                       bk(2)=bkrot(1)*br1(2,1)+bkrot(2)*br1(2,2)+bkrot(3)*br1(2,3)   
                       bk(3)=bkrot(1)*br1(3,1)+bkrot(2)*br1(3,2)+bkrot(3)*br1(3,3)   
                       !---- BKRLOC = Rotloc.R_g.(k+K),  is rotates Rotloc.R.(k+K) . Rotation Rotloc entered by user
                       CALL ROTATE (bk,rotloc(1,1,jatom),bkrloc)
                       !!!! CALL ROTATE (BK,crotloc(1,1,icase),BKRLOC) might beed to do that
                       !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
                       CALL YLM (bkrloc,lmax2,yl)
                       arg1=bkrot(1)*(pos(1,lfirst))*twopi
                       arg2=bkrot(2)*(pos(2,lfirst))*twopi
                       arg3=bkrot(3)*(pos(3,lfirst))*twopi
                       ! ARGT = (K+k)*tau(isym) 
                       argt=(bkx(i)*tauij(1,latom)+bky(i)*tauij(2,latom)+bkz(i)*tauij(3,latom))*twopi
                       ! PHSEHL = e^{I*2pi*( (R_g.(k+K)) *  R(iatom) + (K+k)*tau(isym))}
                       phshel=EXP( imag*(arg1+arg2+arg3+argt) )
                       index=0
                       DO  l=0,lmx
                          maxx=2*l+1
                          DO  m=1,maxx
                             index=index+1
                             h_yl(index,i3)=CONJG(yl(index))*phshel
                          ENDDO
                       ENDDO
                    ENDDO
                    
                    index=0
                    rmt2=1.D0/(rmt(jatom)**2)
                    DO l=0,lmx
                       i3=0
                       DO i=ii,ii+iblock-1
                          IF(i.GT.n0-(nlo+nlon+nlov)) EXIT
                          i3=i3+1
                          if(lapw(l)) then
                             h_al(i3)=dfj(l,i,jatom)*pe(l)-fj(l,i,jatom)*dpe(l) 
                             h_bl(i3)=fj(l,i,jatom)*dp(l)-dfj(l,i,jatom)*p(l)
                          ELSE
                             h_al(i3) = rmt2*fj(l,i,jatom)/p(l)
                             h_bl(i3) = 0.d0
                          ENDIF
                       ENDDO
                       maxx=2*l+1
                       DO m=1,maxx
                          index=index+1
                          i3=0
                          DO i=ii,ii+iblock-1
                             IF(i.GT.n0-(nlo+nlon+nlov)) EXIT
                             i3=i3+1
                             h_alyl(index,i3)=h_al(i3)*h_yl(index,i3)
                             h_blyl(index,i3)=h_bl(i3)*h_yl(index,i3)
                          ENDDO
                       ENDDO
                    ENDDO
                    !
                    lda=(LMAX2+1)*(LMAX2+1)
                    ldc=lda
                    ldb=nmat 
                    CALL zgemm('N','N',index,nemax-nemin+1,ibb,(1.d0,0.d0),h_alyl,lda,As(ii,nemin,is),ldb,(1.d0,0.d0),alm(1,nemin),ldc)
                    CALL zgemm('N','N',index,nemax-nemin+1,ibb,(1.d0,0.d0),h_blyl,lda,As(ii,nemin,is),ldb,(1.d0,0.d0),blm(1,nemin),ldc)
                    
                 ENDDO blocked_loop
                 
              
                 CALL cputim(time2)
                 CALL walltim(time2_w)
                 time_bl=time_bl+time2-time1
                 time_bl_w=time_bl_w+time2_w-time1_w
                 
                 CALL cputim(time1)
                 CALL walltim(time1_w)
                 
                 ldc=(LMAX2+1)*(LMAX2+1)
   
                 if (nlo.ne.0) call lomain (rotloc(:,:,jatom),is,nemin,nemax,lfirst,latom,n0,jatom,alm,blm,clm)
              
                 index=0
                 DO l=0,lmax2
                    maxx=2*l+1
                    cfac=imag**l
                    cfac=cfac*fac
                    DO m=1,maxx
                       index=index+1
                       DO num=nemin,nemax
                          alm(index,num)=alm(index,num)*cfac
                          blm(index,num)=blm(index,num)*cfac
                       ENDDO
                       DO jlo=1,ilo(l)
                          DO num=nemin,nemax
                             clm(index,num,jlo)=clm(index,num,jlo)*cfac
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
                 
                 ! Chenging alm,blm,clm to make them DMFT-like
                 alm(:,DM_nemin:DM_nemaxx)= matmul( alm(:,DM_nemin:DM_nemax), Aweight(:nbands,:nbandsx) )
                 blm(:,DM_nemin:DM_nemaxx)= matmul( blm(:,DM_nemin:DM_nemax), Aweight(:nbands,:nbandsx) )
                 DO jlo=1,nloat
                    clm(:,DM_nemin:DM_nemaxx,jlo)= matmul( clm(:,DM_nemin:DM_nemax,jlo), Aweight(:nbands,:nbandsx) )
                 ENDDO
                 
                 dh_alm(:,:,is) = transpose( alm(:,:) )
                 dh_blm(:,:,is) = transpose( blm(:,:) )
                 DO jlo=1,nloat
                    dh_clm(:,:,jlo,is) = transpose( clm(:,:,jlo) )
                 ENDDO
                 
              enddo  !!! Mar 27, 2010: loop over spin
   
              !!! Mar 27, 2010: we need alm for csplit below. Take the average between up and down
              if (iso.eq.2) then
                 alm(:,:) = transpose( (dh_alm(:,:,1)+dh_alm(:,:,2))/sqrt(2.) )
                 blm(:,:) = transpose( (dh_blm(:,:,1)+dh_blm(:,:,2))/sqrt(2.) )
                 DO jlo=1,nloat
                    clm(:,:,jlo) = transpose( (dh_clm(:,:,jlo,1)+dh_blm(:,:,2))/sqrt(2.) )
                 ENDDO
              endif
   
   
              CALL cputim(time2)
              CALL walltim(time2_w)
              time_reduc=time_reduc+time2-time1
              time_reduc_w=time_reduc_w+time2_w-time1_w
              
              !.....C(L,M) TO BE CALCULATED                                           
              CALL cputim(time1)
              CALL walltim(time1_w)
              
              ilm_loop: DO ILM=1,LMMAX
                 Li=IABS(LM(1,ILM))                                                
                 Mi=LM(2,ILM)                                                      
                 !     SEE KURKI-SUONI  FACTOR FOR  LM+ (1,0), LM-(0,-1), M.NE.2N *(-1)  
                 IMAG1=(1.0,0.0)                                                   
   
                 IF(LM(1,ILM).LT.0) IMAG1= IMAG
                 IF(MOD(Mi,2).EQ.1) IMAG1=-IMAG1
                 
                 l_sum: DO L=0,lmx
                    L1=L+1
                    MMAX=2*L+1
                    lp_sum: DO LP=0,lmx
                       CALL cputim(time3)
                       LP1=LP+1
                       IF(NOTRI(Li,L,LP).LT.0) CYCLE
                       MPMAX=2*LP+1
                       IMAX=JRI(JATOM)
                       
                       ! Radial functions are constructed
                       DO I=1,IMAX
                          UA(I) =RRAD1(I,L)*RRAD1(I,LP)+CIN*RRAD2(I,L)*RRAD2(I,LP)    
                          UB(I) =RADE1(I,L)*RADE1(I,LP)+CIN*RADE2(I,L)*RADE2(I,LP)    
                          UAB(I)=RRAD1(I,L)*RADE1(I,LP)+CIN*RRAD2(I,L)*RADE2(I,LP)   
                          UBA(I)=RRAD1(I,LP)*RADE1(I,L)+CIN*RRAD2(I,LP)*RADE2(I,L)   
                       ENDDO
                       ! Radial functions as product of local orbitals and u,dot{u}
                       DO jlo=1,ilo(l)
                          IF(loor(jlo,l)) THEN
                             DO I=1,IMAX
                                U21(I,jlo) =a1lo(I,jlo,L)*RRAD1(I,LP)+CIN*b1lo(I,jlo,L)*RRAD2(I,LP)    
                                Ue21(I,jlo)=a1lo(I,jlo,L)*RADE1(I,LP)+CIN*b1lo(I,jlo,L)*RADE2(I,LP)    
                             ENDDO
                          ENDIF
                       ENDDO
                       DO jlop=1,ilo(lp)
                          IF(loor(jlop,lp)) THEN
                             DO i=1,imax
                                u12(i,jlop) =RRAD1(i,l)*a1lo(i,jlop,lp)+CIN*RRAD2(i,l)*b1lo(i,jlop,lp)    
                                ue12(i,jlop)=RADE1(i,l)*a1lo(i,jlop,lp)+CIN*RADE2(i,l)*b1lo(i,jlop,lp) 
                             ENDDO
                          ENDIF
                       ENDDO
                       ! Radial functions for local orbitals
                       DO jlo=1,ilo(l)
                          DO jlop=1,ilo(lp) 
                             IF(loor(jlo,l).AND.loor(jlop,lp)) THEN
                                DO i=1,imax
                                   u22(i,jlop,jlo)=a1lo(i,jlo,l)*a1lo(i,jlop,lp)+cin*b1lo(i,jlo,l)*b1lo(i,jlop,lp) 
                                ENDDO
                             ENDIF
                          ENDDO
                       ENDDO
                       CALL cputim(time4)
                       time_radprod=time_radprod+time4-time3
                       !.....M SUM                           
   
                       sa=0.0d0; sb=0.0d0; sab=0.0d0; sba=0.0d0; s12=0.0d0; 
                       se12=0.0d0; s21=0.0d0; se21=0.0d0; s22=0.0d0
                       vsa=0.0d0; vsb=0.0d0; vsab=0.0d0; vsba=0.0d0; vs12=0.0d0
                       vse12=0.0d0; vs21=0.0d0; vse21=0.0d0; vs22=0.0d0
                       
                       
                       CALL cputim(time3)
                       DO num=nemin,DM_nemaxx
                          tsuma = 0.0; tsumb = 0.0; tsumab= 0.0; tsumba= 0.0
                          tsum21 (:ilo(l))  = 0.0;  tsume21(:ilo(l))  = 0.0
                          tsum12 (:ilo(lp)) = 0.0;  tsume12(:ilo(lp)) = 0.0
                          tsum22 (:ilo(lp),:ilo(l))=0.0
                          
                          M=-(L+1) 
                          MP=-(LP+1)
                          DO MS=1,2*L+1
                             M=M+1
                             MP=-LP1
                             DO MPS=1,2*LP+1
                                MP=MP+1
                                MTEST=-M+Mi+MP
                                IF(MTEST.NE.0) CYCLE
                                LY  = L*(L+1)+M+1
                                LPY = LP*(LP+1)+MP+1
                                GNT=GAUNT(L,Li,LP,M,Mi,MP)
                                gint = IMAG1*GNT
                                
                                do is=1,iso
                                   tSA = dh_ALM(num,LY,is)*dconjg(dh_ALM(num,LPY,is))*gint
                                   tSB = dh_BLM(num,LY,is)*dconjg(dh_BLM(num,LPY,is))*gint
                                   tSAB= dh_ALM(num,LY,is)*dconjg(dh_BLM(num,LPY,is))*gint
                                   tSBA= dh_BLM(num,LY,is)*dconjg(dh_ALM(num,LPY,is))*gint
                                   DO jlo=1,ilo(l)
                                      tS21(jlo)  = dh_cLM(num,LY,jlo,is)*dconjg(dh_ALM(num,LPY,is))*gint
                                      tSe21(jlo) = dh_cLM(num,LY,jlo,is)*dconjg(dh_BLM(num,LPY,is))*gint
                                   ENDDO
                                   DO jlop=1,ilo(lp)
                                      tS12(jlop)  = dh_ALM(num,LY,is)*dconjg(dh_cLM(num,LPY,jlop,is))*gint
                                      tSe12(jlop) = dh_BLM(num,LY,is)*dconjg(dh_cLM(num,LPY,jlop,is))*gint
                                   ENDDO
                                   DO jlo=1,ilo(l)
                                      DO jlop=1,ilo(lp)
                                         tS22(jlop,jlo) = dh_CLM(num,LY,jlo,is)*dconjg(dh_CLM(num,LPY,jlop,is))*gint
                                      ENDDO
                                   ENDDO
                                   
                                   tSUMA  = tSUMA  + tSA/iso
                                   tSUMB  = tSUMB  + tSB/iso
                                   tSUMAB = tSUMAB + tSAB/iso
                                   tSUMBA = tSUMBA + tSBA/iso
                                   DO jlo=1,ilo(l)
                                      tSUM21 (jlo) = tSUM21 (jlo) + tS21 (jlo)/iso
                                      tSUMe21(jlo) = tSUMe21(jlo) + tSe21(jlo)/iso
                                   ENDDO
                                   DO jlop=1,ilo(lp)
                                      tSUM12 (jlop) = tSUM12 (jlop) + tS12(jlop)/iso
                                      tSUMe12(jlop) = tSUMe12(jlop) + tSe12(jlop)/iso
                                   ENDDO
                                   DO jlo=1,ilo(l)
                                      DO jlop=1,ilo(lp)
                                         tSUM22(jlop,jlo) = tSUM22(jlop,jlo) + tS22(jlop,jlo)/iso
                                      ENDDO
                                   ENDDO
                                enddo
                                
                             ENDDO
                          ENDDO
                          
                          wkp = weight(num)
                          facv = (-e(num))*wkp
                          sa  = sa  + tsuma  * wkp
                          sb  = sb  + tsumb  * wkp
                          sab = sab + tsumab * wkp
                          sba = sba + tsumba * wkp
                          vsa = vsa + tsuma  * facv
                          vsb = vsb + tsumb  * facv
                          vsab= vsab+ tsumab * facv
                          vsba= vsba+ tsumba * facv
                          DO jlo=1,ilo(l)
                             s21(jlo)   = s21(jlo)  + tsum21 (jlo) * wkp
                             se21(jlo)  = se21(jlo) + tsume21(jlo) * wkp
                             vs21(jlo)  = vs21(jlo) + tsum21 (jlo) * facv
                             vse21(jlo) = vse21(jlo)+ tsume21(jlo) * facv
                          ENDDO
                          DO jlop=1,ilo(lp)
                             s12  (jlop) = s12  (jlop) + tsum12 (jlop) * wkp
                             se12 (jlop) = se12 (jlop) + tsume12(jlop) * wkp
                             vs12 (jlop) = vs12 (jlop) + tsum12 (jlop) * facv
                             vse12(jlop) = vse12(jlop) + tsume12(jlop) * facv
                          ENDDO
                          DO jlo=1,ilo(l)
                             DO jlop=1,ilo(lp)
                                s22 (jlop,jlo) = s22 (jlop,jlo) + tsum22(jlop,jlo) * wkp
                                vs22(jlop,jlo) = vs22(jlop,jlo) + tsum22(jlop,jlo) * facv
                             ENDDO
                          ENDDO
   
                          suma  (num) = tsuma
                          sumb  (num) = tsumb
                          sumab (num) = tsumab
                          sumba (num) = tsumba
                          sum21 (num,:ilo(l))  = tsum21 (:ilo(l))
                          sume21(num,:ilo(l))  = tsume21(:ilo(l))
                          sum12 (num,:ilo(lp)) = tsum12 (:ilo(lp))
                          sume12(num,:ilo(lp)) = tsume12(:ilo(lp))
                          sum22 (num,:ilo(lp),:ilo(l)) = tsum22 (:ilo(lp),:ilo(l))
                       ENDDO
                       CALL cputim(time4)
                       time_m=time_m+time4-time3
                       
                       IF(Li.eq.0) then
                          call csplit(nemin,DM_nemaxx,l,jatom,mu,alm,blm,clm,coord) 
                       ENDIF
                       CALL cputim(time3)
                       IMAX=JRI(JATOM)
                       DO I=1,IMAX
                          TC_buf(i)=SA*UA(I)+SB*UB(I)+SAB*UAB(I)+SBA*UBA(I)
                          vTC_buf(i)=vSA*UA(I)+vSB*UB(I)+vSAB*UAB(I)+vSBA*UBA(I)
                       ENDDO
                       DO jlo=1,ilo(l)
                          DO i=1,imax
                             tc_buf(i)=tc_buf(i)+S21(jlo)*U21(I,jlo)+Se21(jlo)*Ue21(I,jlo)
                             vtc_buf(i)=vtc_buf(i)+vS21(jlo)*U21(I,jlo)+vSe21(jlo)*Ue21(I,jlo)
                          ENDDO
                       ENDDO
                       DO jlop=1,ilo(lp)
                          DO i=1,imax
                             tc_buf(i)=tc_buf(i)+S12(jlop)*U12(I,jlop)+Se12(jlop)*Ue12(I,jlop)
                             vtc_buf(i)=vtc_buf(i)+vS12(jlop)*U12(I,jlop)+vSe12(jlop)*Ue12(I,jlop)
                          ENDDO
                       ENDDO
                       DO jlo=1,ilo(l)
                          DO jlop=1,ilo(lp)
                             DO i=1,imax
                                tc_buf(i)=tc_buf(i)+s22(jlop,jlo)*u22(i,jlop,jlo)
                                vtc_buf(i)=vtc_buf(i)+vs22(jlop,jlo)*u22(i,jlop,jlo)
                             ENDDO
                          ENDDO
                       ENDDO
                       DO i=1,imax
                          w_RHOLM(I,ILM,jatom)=w_RHOLM(I,ILM,jatom)+TC_buf(i)/MULT(JATOM)
                          w_vRHOLM(I,ILM,jatom)=w_vRHOLM(I,ILM,jatom)+vTC_buf(i)/MULT(JATOM)
                       ENDDO
                       CALL cputim(time4)
                       time_rad=time_rad+time4-time3
                    ENDDO lp_sum
                 ENDDO l_sum
              ENDDO ilm_loop
              
              CALL cputim(time2)
              CALL walltim(time2_w)
              time_ilm=time_ilm+time2-time1
              time_ilm_w=time_ilm_w+time2_w-time1_w
           ENDDO equiv_atom_loop
   
           call psplit(w_xwt1(:,jatom),w_xwteh(:,jatom),w_xwtel(:,jatom),w_xwt1h(:,jatom),w_xwt1l(:,jatom),jatom,nemin,DM_nemaxx,test1,EQBAD,jatombad,lbad)
        ENDDO
   
        CALL cputim(time1)
        CALL walltim(time1_w)
        ! INTERSTITIAL CHARGE
        ! FFT eigenvectors and accumulate weighted square in sumfft
        keigen(1,:)=kx(:)
        keigen(2,:)=ky(:)
        keigen(3,:)=kz(:)
        isig=-1
   
        ALLOCATE( Asc(n0-nnlo,1:DM_nemaxx) )
   
        !!!! Mar 27, 2010: Over spins in interstitals
        do is=1,iso
           Asc(:,:) = As(:n0-nnlo,:DM_nemaxx,is)
           Asc(:,DM_nemin:DM_nemaxx) = matmul( As(:n0-nnlo,DM_nemin:DM_nemax,is), Aweight(:nbands,:nbandsx) )
   
           DO num=nemin,DM_nemaxx
              !  puts eigenvector into fft array for FFT
              CALL setfft1c(n0-nnlo,iff1,iff2,iff3,Asc(:,num),fft,keigen)
              !  computes a(k,r) = FFT(a(k,K))
              CALL c3fft(iff1,iff2,iff3,fft,iff1,iff2,isig,cwork,dwork,ierr)  
              !  |a(k,r)|^2 w(k)
              CALL sumupfft(iff1,iff2,iff3, weight(num)/iso,       fft,sumfft)
              ! -|a(k,r)|^2 w(k)*e(k)
              CALL sumupfft(iff1,iff2,iff3,-weight(num)*E(num)/iso,fft,vsumfft)   
           ENDDO
        enddo
   
        DEALLOCATE( Asc, Aweight )
   
        CALL cputim(time2)
        CALL walltim(time2_w)
        time_int =time_int +time2-time1
        time_intw =time_intw +time2_w-time1_w
     ENDDO     !  kpoint loop end
     call gh5_close(f_id) 
     DO is=1,iso    !------ over up/dn ---------------------!
        itape1=8+is
        itape2=29+is
        if (vector_para) then
           close(itape1)
           close(itape2)
        else
           rewind(itape1)
           rewind(itape2)
        endif
     ENDDO
  ENDDO

  call gh5_end()

  CALL Reduce_MPI(xwt, etot2, w_RHOLM, w_vRHOLM, w_xwt1, w_xwteh, w_xwtel, w_xwt1h, w_xwt1l, sumfft, vsumfft, nomega, NRAD, LM_MAX, nat, iff1, iff2, iff3)
#ifdef MPI
  CALL DSUM1_ALL_MPI(etot)
#endif
  if (myrank.eq.master) then
     WRITE(*,*)'LDA Sum of eigenvalues=',ETOT, 'Ry.'
     WRITE(*,*)'    Sum of eigenvalues=', ETOT2, 'Ry.'
     WRITE(*,*)'  Ratio to renormalize=', elecn/xwt
     IF(ABS(elecn/xwt-1._8).GT.1.E-4_8)THEN
       WRITE(*,'(" SEVERE WARNING: Ratio to renormalize deviate too much from 1!")')
     ENDIF
     if (Rho_Renormalize) then
        w_RHOLM = w_RHOLM * (elecn/xwt)
        w_vRHOLM = w_vRHOLM * (elecn/xwt)
     endif
  endif

  if (myrank.eq.master) then

     DO jatom=1,nat
        lm(1:2,1:NCOM) = w_lm(1:2,1:NCOM,jatom)
        LMMAX = w_lmmax(jatom)
        R(:) = w_R(:,jatom)
        
        CALL cputim(time1)
        CALL walltim(time1_w)
        
        DO ILM=1,LMMAX
           OVRDIV=1.0
           Li=LM(1,ILM)
           Mi=LM(2,ILM)
           IF(Li.NE.0) THEN
              IF(Mi.NE.0) THEN
                 IMAX=JRI(JATOM)
                 DO I=1,IMAX
                    w_vRHOLM(I,ILM,jatom)=w_vRHOLM(I,ILM,jatom)*SQRT2
                    w_RHOLM(I,ILM,jatom)=w_RHOLM(I,ILM,jatom)*SQRT2
                 enddo
              ENDIF
              CYCLE
           ENDIF
           CALL CHARGE(R,OVRDIV,SQFP,w_RHOLM(1,ILM,jatom),DX(JATOM),JRI(JATOM),TC)
           WRITE(6,207)  JATOM,TC
           WRITE(21,720) jatom,jatom,TC
           IF(ISPLIT(JATOM).EQ.1) THEN
              WRITE(6,250) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,250) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.2) THEN                               
              WRITE(6,251) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,251) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.3) THEN                               
              WRITE(6,252) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,252)jatom,JATOM, JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.-2) THEN                              
              WRITE(6,253) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,253) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.4) THEN                               
              WRITE(6,254) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,254) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.5) THEN                               
              WRITE(6,255) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,255) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.6) THEN                               
              WRITE(6,256) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,256) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.8) THEN                               
              WRITE(6,257) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,257) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.15) THEN                               
              WRITE(6,258) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,21)
              WRITE(21,258) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,21)
           ENDIF
           do i=0,3
              if(w_xwt1l(i,jatom).lt.0.00005d0) then
                 w_XWT1l(i,jatom)=0.00001d0
                 w_XWTel(i,jatom)=w_XWT1l(i,jatom)*10.00001d0
              endif
              if(w_xwt1h(i,jatom).lt.0.00005d0) then
                 w_XWT1h(i,jatom)=0.00001d0
                 w_XWTeh(i,jatom)=w_XWT1h(i,jatom)*10.00001d0
              endif
           enddo
           write(6,260)
           WRITE(6,261) jatom,(w_XWT1l(I,jatom),w_XWTel(I,jatom)/w_XWT1l(I,jatom),I=0,3)
           write(6,262)
           WRITE(6,263) jatom,(w_XWT1h(I,jatom),w_XWTeh(I,jatom)/w_XWT1h(I,jatom),I=0,3)
           write(21,260)
           WRITE(21,261) jatom,(w_XWT1l(I,jatom),w_XWTel(I,jatom)/w_XWT1l(I,jatom),I=0,3)
           write(21,262)
           WRITE(21,263) jatom,(w_XWT1h(I,jatom),w_XWTeh(I,jatom)/w_XWT1h(I,jatom),I=0,3)
           !
           IF(Mi.NE.0) THEN
              IMAX=JRI(JATOM)                                                   
              DO I=1,IMAX                                                    
                 w_vRHOLM(I,ILM,jatom)=w_vRHOLM(I,ILM,jatom)*SQRT2
                 w_RHOLM(I,ILM,jatom)=w_RHOLM(I,ILM,jatom)*SQRT2
              enddo
           ENDIF
        enddo
     
        call cputim(time3)
        call walltim(time3_w)
        
        ! Writting case.clmval
        WRITE(8,1990) JATOM
        WRITE(8,2001) LMMAX
        DO ILM=1,LMMAX
           Li=LM(1,ILM)
           Mi=LM(2,ILM)
           WRITE(8,2011) Li,Mi
           WRITE(8,2022) ( w_RHOLM(I,ILM,jatom), I=1,IMAX )
           WRITE(8,2031)
        end DO
        WRITE(8,2030)
        
        if(vresp_write) then
           WRITE(28,1990) JATOM
           WRITE(28,2001) LMMAX
           DO ILM=1,LMMAX
              Li=LM(1,ILM)
              Mi=LM(2,ILM)
              WRITE(28,2011) Li,Mi
              WRITE(28,2022) ( w_vRHOLM(I,ILM,jatom), I=1,IMAX )
              WRITE(28,2031)
           end DO
           WRITE(28,2030)
        endif
        
        call cputim(time4)
        call walltim(time4_w)
        time_writeclm=time_writeclm+time4-time3
        time_writeclm_w=time_writeclm_w+time4_w-time3_w
        
        DO ILM1=1,LMMAX
           IF (LM(1,ILM1).EQ.2.AND.LM(2,ILM1).EQ.0) THEN
              DO I=1,IMAX                                                
                 w_RHOLM(I,ILM1,jatom)=w_RHOLM(I,ILM1,jatom)/R(I)**3
              enddo
              LABEL=1
              DO ILM2=ILM1+1,LMMAX
                 IF (LM(1,ILM2).EQ.2.AND.LM(2,ILM2).EQ.2) THEN
                    DO I=1,IMAX
                       w_RHOLM(I,ILM2,jatom)=w_RHOLM(I,ILM2,jatom)/R(I)**3
                    enddo
                    LABEL=2
                    DO ILM3=ILM2+1,LMMAX
                       IF (LM(1,ILM3).EQ.-2.AND.LM(2,ILM3).EQ.2) THEN
                          DO I=1,IMAX
                             w_RHOLM(I,ILM3,jatom)=w_RHOLM(I,ILM3,jatom)/R(I)**3
                          enddo
                          LABEL=3
                          GOTO 7372
                       ENDIF
                    ENDDO
                    GOTO 7372
                 ENDIF
              ENDDO
              GOTO 7372
           ENDIF
        ENDDO
     
        GOTO 7371
     
7372 CONTINUE 
     
        EFGFACT=0.02997925
        FAC=0.8D0*PI*324.14D0
        DO I1=0,jri(jatom)-5,1000   ! 7374
           ! integrates: TC = integrate(R*RHOLM)
           CALL CHARGE(R,1.D0,1.D0,w_RHOLM(1,ILM1,jatom),DX(JATOM),JRI(JATOM)-I1,TC) 
           TC=TC*FAC
           EFG20=TC
           IF (LABEL.GT.1) THEN
              CALL CHARGE(R,1.D0,1.D0,w_RHOLM(1,ILM2,jatom),DX(JATOM),JRI(JATOM)-I1,TC) 
              TC=TC*FAC
              EFG22=TC
              IF (LABEL.EQ.3) THEN
                 CALL CHARGE(R,1.D0,1.D0,w_RHOLM(1,ILM3,jatom),DX(JATOM),JRI(JATOM)-I1,TC)
                 TC=TC*FAC
                 EFG2M=TC
                 QXX=(-EFG20/SQRT(3.D0)+EFG22)*SQRT(15.D0/4.D0/PI)
                 QYY=(-EFG20/SQRT(3.D0)-EFG22)*SQRT(15.D0/4.D0/PI)
                 QZZ=(EFG20*2.D0/SQRT(3.D0))*SQRT(15.D0/4.D0/PI)
                 QXY=EFG2M*SQRT(15.D0/4.D0/PI)
                 QXX=-QXX*EFGFACT
                 QXY=-QXY*EFGFACT
                 QYY=-QYY*EFGFACT
                 QZZ=-QZZ*EFGFACT
                 IF (I1.EQ.0) WRITE(6,3019)
                 IF (I1.EQ.0) WRITE(21,3019)
                 WRITE(6,3121) jatom,QXX,QXY,QYY,QZZ,R(JRI(JATOM)-I1)
                 IF (I1.EQ.0) WRITE(21,3121) jatom,QXX,QXY,QYY,QZZ,R(JRI(JATOM)-I1)      
              ELSE
                 VXX=(-EFG20/SQRT(3.D0)+EFG22)*SQRT(15.D0/4.D0/PI)
                 VYY=(-EFG20/SQRT(3.D0)-EFG22)*SQRT(15.D0/4.D0/PI)
                 VZZ=(EFG20*2.D0/SQRT(3.D0))*SQRT(15.D0/4.D0/PI)
                 IF (I1.EQ.0) WRITE(6,3017)
                 VXX=-VXX*EFGFACT
                 VYY=-VYY*EFGFACT
                 VZZ=-VZZ*EFGFACT
                 IF (I1.EQ.0) WRITE(21,3017)
                 WRITE(6,3021) jatom,VXX,VYY,VZZ,R(JRI(JATOM)-I1)
                 IF (I1.EQ.0) WRITE(21,3021) jatom,VXX,VYY,VZZ,R(JRI(JATOM)-I1)
              ENDIF
           ELSE
              EFG20=EFG20*2.D0/SQRT(3.D0)*SQRT(15.D0/4.D0/PI)
              EFG20=-EFG20*EFGFACT
              WRITE(6,7207)  JATOM,EFG20,R(JRI(JATOM)-I1)
              IF (I1.EQ.0) WRITE(21,7720) jatom,JATOM,EFG20,R(JRI(JATOM)-I1)
           ENDIF
        ENDDO !7374       CONTINUE                                                          
     
     
7371 CONTINUE                                     
     
        REWIND fh_vec
        
        CALL cputim(time3) 
        call walltim(time3_w)
        time_writeefg=time_writeefg+time3-time4
        time_writeefg_w=time_writeefg_w+time3_w-time4_w
        
        
        CALL cputim(time2)
        CALL walltim(time2_w)
        time_writescf=time_writescf+time2-time3
        time_write=time_write+time2-time1
        time_writescf_w=time_writescf_w+time2_w-time3_w
        time_write_w=time_write_w+time2_w-time1_w
     ENDDO
     
     write(6,*)
     write(6,'(a,2f8.1)') '   atpar              (cpu,wall):',time_atpar_c,time_atpar_w
     write(6,'(a,2f8.1)') '   readvec            (cpu,wall):',time_rd_c,time_rd_w
     write(6,'(a,2f8.1)') '   readvec, read only (cpu,wall):',time_r_c,time_r_w
     WRITE(6,'(a,2f8.1)') '   blocked loop            (cpu):',time_bl,time_bl_w
     WRITE(6,'(a,2f8.1)') '   reduced                 (cpu):',time_reduc,time_reduc_w
     WRITE(6,'(a,2f8.1)') '   ilm-tot loop            (cpu):',time_ilm,time_ilm_w
     WRITE(6,'(a,f8.1)')  '   ilm-1 loop              (cpu):',time_radprod
     WRITE(6,'(a,f8.1)')  '   ilm-2 loop              (cpu):',time_m
     WRITE(6,'(a,f8.1)')  '   ilm-3 loop              (cpu):',time_rad
     WRITE(6,'(a,2f8.1)') '   write-clm          (cpu,wall):',time_writeclm,time_writeclm_w
     WRITE(6,'(a,2f8.1)') '   write-scf          (cpu,wall):',time_writescf,time_writescf_w
     WRITE(6,'(a,2f8.1)') '   write-efg          (cpu,wakk):',time_writeefg,time_writeefg_w
     WRITE(6,'(a,2f8.1)') '   write              (cpu,wall):',time_write,time_write_w
     WRITE(6,'(a,2f8.1)') '   dmft renorm-overlap(cpu,wall):',time_dmf0, time_dmf0w
     WRITE(6,'(a,2f8.1)') '   dmft chemical pot  (cpu,wall):',time_dmf1, time_dmf1w
     WRITE(6,'(a,2f8.1)') '   dmft transformation(cpu,wall):',time_dmf, time_dmfw
     WRITE(6,'(a,2f8.1)') '   dmft exact-diag    (cpu,wall):',time_ed, time_edw
     WRITE(6,'(a,2f8.1)') '   interstitial charge(cpu,wall):',time_int, time_intw
     
     WRITE(6,881)  XWT
     WRITE(21,881) XWT
     WRITE(21,725) ETOT2 !+ ts2
     WRITE(6,882) ETOT2  !+ ts2
     WRITE(6,'(" EN_DEEP LEVEL=",F12.7,"Ry.")')EDL
     !
  endif

  REWIND fh_vec


  deallocate (e_store,elo_store)
  DEALLOCATE(alm,blm,clm)
  DEALLOCATE(dh_alm,dh_blm,dh_clm)
  CALL fini_charp
  CALL fini_chard
  CALL fini_charf
  CALL fini_lohelp
  CALL fini_xa
  CALL fini_xa3
  CALL w_deallocate

  CALL cputim(time2)
  TCLM = time2
  CALL walltim(time2_w)
  TCLM_w = time2

  
  if (myrank.EQ.master) then
     
     ! ---- Interstitials ---
     !
     ! Back transform cumulative grids
     isig=1
     ! computes rho(K)=FFT(rho(r))
     CALL c3fft(iff1,iff2,iff3,sumfft,iff1,iff2,isig,cwork,dwork,ierr) 
     CALL c3fft(iff1,iff2,iff3,vsumfft,iff1,iff2,isig,cwork,dwork,ierr)
     sumfft(:,:,:)=sumfft(:,:,:)/dble(iff1*iff2*iff3)
     vsumfft(:,:,:)=vsumfft(:,:,:)/dble(iff1*iff2*iff3)
     
     ALLOCATE( rho1(nwave*nsym),vrho1(nwave*nsym) )
     rho1(:)=0.0 ; vrho1(:)=0.0
     
     call getfft(nwave,iff1,iff2,iff3,rho1,sumfft,kmax)
     call getfft(nwave,iff1,iff2,iff3,vrho1,vsumfft,kmax)
     
     ALLOCATE( rhok(nwave),vrhok(nwave) )
     vrhok(:)=0.0;  rhok(:)=0.0
     
     !.....SUM OVER ALL RECPR. LATTIC VECTORS OF ONE STAR
     ia1=1
     volin = 1.0d0/vol
     DO j=1,nwave
        ia2=ia1+inst(j)-1
        DO jj=ia1,ia2
           rhok(j)  = rhok(j)  + rho1(jj)*CONJG(tauk(jj))
           vrhok(j) = vrhok(j) + vrho1(jj)*CONJG(tauk(jj))
        ENDDO
        ia1=ia2+1
        vrhok(j) = vrhok(j) * volin
        rhok(j)  = rhok(j)  * volin
     ENDDO
     
     if (Rho_Renormalize) then
        rhok = rhok * (elecn/xwt)
        vrhok = vrhok * (elecn/xwt)
     endif
     
     WRITE(8,*) '   VALENCE CHARGE DENSITY IN INTERSTITIAL '
     WRITE(8,*)
     WRITE(8,2061)  NWAVE
     WRITE(8,2071)  ( (KZZ(JX,J),JX=1,3),RHOK(J), J=1,NWAVE)
     
     if(vresp_write) then
        WRITE(28,*) '   VALENCE CHARGE DENSITY IN INTERSTITIAL '
        WRITE(28,*)
        WRITE(28,2061)  NWAVE
        WRITE(28,2071)  ( (KZZ(JX,J),JX=1,3),VRHOK(J), J=1,NWAVE)
     endif
     WRITE(6,204) NWAVE,kmax
     
     close(8)
     
     DEALLOCATE(rhok,vrhok)
     DEALLOCATE(rho1,vrho1)
  endif

  DEALLOCATE (fft, tp, sumfft, vsumfft, tloc )
  DEALLOCATE (dwork, cwork)
  ! --- Interstitials ---

  CALL cputim(time2)
  TFOUR = time2
  CALL walltim(time2_w)
  TFOUR_w = time2

  RETURN
  

205 FORMAT(/,1X,' K-POINT:',3F8.4,1X,I5,I4,2X,A10)                      
207 FORMAT(1H0,' TOTAL CHARGE INSIDE SPHERE',I5,1H:,F12.6)            
250 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PZ,PXY',/,':QTL',i3.3,':',12F7.4)               
251 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,      ','D-EG,D-T2G ',/,':QTL',i3.3,':',12F7.4)  
252 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,      ','D-Z2,D-XY,X2Y2,D-XZ,YZ ',/,':QTL',i3.3,':',12F7.4)
253 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PZ,PXY,,','D-Z2,D-XY,D-X2Y2,D-XZ,YZ ',/,':QTL',i3.3,':',12F7.4)                          
254 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PZ,PXY, ','D-Z2,D-XY,X2Y2,D-XZ,YZ ',/,':QTL',i3.3,':',12F7.4)                            
255 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,      ','D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ ',/,':QTL',i3.3,':',12F7.4)                        
256 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PX,PY,PZ,',' ',/,':QTL',i3.3,':',12F7.4)
257 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PX,PY,PZ,','D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ ',/,':QTL',i3.3,':',12F7.4)                        
258 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PX,PY,PZ,D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ,F00,F11,F22',',F33,','F1M,F2M,F3M ',/,':QTL',i3.3,':',19F7.4)
260 FORMAT(8x,'Q-s-low E-s-low   Q-p-low E-p-low   Q-d-low E-d-low','   Q-f-low E-f-low')
261 FORMAT(':EPL',i3.3,':',4(2F8.4,2x))
262 FORMAT(8x,'Q-s-hi  E-s-hi    Q-p-hi  E-p-hi    Q-d-hi  E-d-hi ','   Q-f-hi  E-f-hi ')
263 FORMAT(':EPH',i3.3,':',4(2F8.4,2x))
501 FORMAT(2i5,' k-point, nband')
503 FORMAT(5e17.9)
720 FORMAT(/,':CHA',i3.3,':',1X,'TOTAL CHARGE INSIDE SPHERE ',I3,' = ',F12.6)         
725 FORMAT(/,':ESUM  :',1X,'SUM OF EIGENVALUES =  ',F20.9,/)            
787 FORMAT('    VALENCE CHARGE DENSITY   IN  MT SPHERES',5X,I3,' ITERATION')
800 FORMAT(/,10X,68(1H-),/,11X,'W A V E F U N C T I O N S','   AND   C H A R G E S   IN   S P H E R E S',/,10X,68(1H-))                                         
881 FORMAT(/,':CHA  :',' TOTAL CHARGE INSIDE UNIT CELL =',F15.6) 
882 FORMAT(//,3X,'SUM OF EIGENVALUES:'/24X,F15.6/)                    
1003 FORMAT(251(I3,I2))                                                 
1005 FORMAT(/,7X,'LMMAX',I3,/,7x,'LM= ',17(i3,I2),(/,7x,18(i3,i2)))          
1990 FORMAT(3X,'ATOMNUMBER   ',I3,5X,10A4)                             
2001 FORMAT(3X,'NUMBER OF LM',I3//)                                   
2011 FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)                         
2022 FORMAT(3X,4E19.12)                                                 
2030 FORMAT(///)                                                       
2031 FORMAT(/)                                                         
2032 FORMAT(49X,I3,//)  
2055 FORMAT(/,1X,' K-POINT:',3F14.10,1X,I5,I4,2X,A10)
3017 FORMAT(/,22X,'VXX',9X,'VYY',9X,'VZZ',7X,'UP TO R',/)
3018 FORMAT(/'  EFG COMPONENTS BEFORE TRANSFORMATION [10**21 V/(M*M)]:',/,15X,'QXX (-V20/2. + V22)     =',F10.3,/,15X,'QYY (-V20/2. - V22)     =',F10.3,/,15X,'QZZ ( V20)              =',F10.3,/,15X,'QXY ( V2M/2. * SQRT(3.))=',F10.3)  
3019 FORMAT(/,22X,'QXX',9X,'QXY',9X,'QYY',9X,'QZZ',7X,'UP TO R',/)
3020 FORMAT(/'  EFG COMPONENTS AFTER TRANSFORMATION [10**21 V/(M*M)]:',/,15X,'VXX                     =',F10.3,/,15X,'VYY                     =',F10.3,/,15X,'VZZ                     =',F10.3,/)
3121 FORMAT(':VZZ',i3.3,':',8x,4F12.5,F12.3)
3021 FORMAT(':VZZ',i3.3,':',8x,3F12.5,F12.3)
3022 FORMAT(  15X,'QXX=',F10.3,/,15X,'QYY=',F10.3,/,15X,'QZZ=',F10.3,15X,'QXY=',F10.3,' UP TO R= ',F10.5)
4645 FORMAT(2i4,3f15.10)
4646 FORMAT(8e20.10)
4893 FORMAT(3i4,4(2e16.8,2x))  
13   FORMAT(////,':POS',i3.3,':',1x,'AT.NR.',I4,1X,'POSITION =',3F8.5,2X,'MULTIPLICITY =',I3)
7207 FORMAT(1H ,' EFG INSIDE SPHERE',I5,1H:,F12.6,5X,' UP TO R = ',F10.7)                                                          
7720 FORMAT(':VZZ',i3.3,':',1X,'EFG INSIDE SPHERE ',I3,' = ',F12.6,5X,' UP TO R =',F10.5)


204 FORMAT(I10,' FOURIER COEFFICIENTS CALCULATED',/' LARGEST COMPONENTS:',3I6)
2061 FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')
2071 FORMAT(3X,3I5,2E19.12)

END SUBROUTINE L2MAIN
