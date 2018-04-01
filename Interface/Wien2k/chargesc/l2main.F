subroutine l2main(cform,nnlo,coord,zz,nwave,&
        &tclm,tclm_w,tfour,tfour_w)
  use defs,  only: pi, zero, two, zeroc, imag
  use param, only: iblock, lmax2, lomax, ncom, nloat, nrad, &
      &ngau, nmat, nume, nsym, nkpt, fh_vec, fastfilesystem
  use atspdt,only: el, p, dp, pe, dpe, pei, e_store
  use char,  only: modus
  use com,   only: rel, sspin, nat, nband, nk, minwav, maxwav, &
      &emin, elecn, xwt
  use charp, only: init_charp, zero_charp, fini_charp
  use chard, only: init_chard, zero_chard, fini_chard
  use charf, only: init_charf, zero_charf, fini_charf
  use lo,    only: loor, lapw, nlo, nlov, nlon, ilo, a1lo, b1lo, &
      &elo_store, alo, blo, clo, pi12lo, pe12lo, pr12lo
  use lohelp,only: u21, ue21, u12, ue12, u22, init_lohelp, &
      &zerotc_lohelp, fini_lohelp, sum12, sum21, sum22, sume21, sume12
  use struk, only: iatnr, mult, isplit, jri, pos, rotij, tauij, &
      &rmt, v, rotloc, r0, dx, vol, br1, br2
  use xa,    only: lm, r, avec, bk, bkrot, bkrloc, fj, dfj, e, &
      &weight, tc100, tca100, tcb100, init_xa, fini_xa, &
      &suma, sumb, sumab, sumba
  use xa3,   only: bkx, bky, bkz, bkxlo, bkylo, bkzlo, init_xa3, &
      &fini_xa3, kx, ky, kz, kxlo, kylo, kzlo, as, as_lo!, a, a_lo
  use reclat,only: inst, tauk, kzz
  use com_mpi,only: nprocs,myrank,master,filenamempi,reduce_mpi, &
      &gather_procs,cpuid,pr_proc,pr_procr,pr_procs,vector_para, &
      &vectors,nvector,fvectors,vecfn,findmaxk_mpi
  use w_atpar
  use dmf
  use gmpi,only:gp
  use ghdf5_base
  implicit none
  character*4, intent(in)  :: cform
  integer, intent(out)     :: nnlo
  character*5, intent(in)  :: coord
  real*8, intent(in)       :: zz(*)
  integer, intent(in)      :: nwave
  real*8, intent(out)      :: tclm,tclm_w,tfour,tfour_w
  ! functions
  integer :: notri
  real*8  :: gaunt
  real*8  :: ferm
  ! common blocks
  common /fact/   fct(100)
  real*8     :: fct
  common /radfu/  rrad1(nrad,0:lmax2),rade1(nrad,0:lmax2),rrad2(nrad,0:lmax2),rade2(nrad,0:lmax2)
  real*8     :: rrad1, rade1, rrad2, rade2
  common /uhelp/  ua(nrad),ub(nrad),uba(nrad),uab(nrad)
  real*8     :: ua, ub, uba, uab
  common /xmean/  xwteh(0:3),xwtel(0:3),xwt1h(0:3),xwt1l(0:3) 
  real*8     :: xwteh, xwtel, xwt1h, xwt1l
  ! locals
  real*8, parameter       :: ry2ev = 13.60569193
  !
  character*10  :: kname
  logical       :: more_kpoints
  complex*16    :: yl((lmax2+1)*(lmax2+1))
  complex*16    :: phshel,cfac,imag1
  real*8        :: tc_buf(nrad),vtc_buf(nrad)
  character*100 :: cdummy
  integer       :: i, j, n, jk, iscf, lfirst, lmx, jatom, ilm1, ilm3, n0, ne, nm
  integer       :: jlm, lmmax, imax, nemin, dm_nemin, dm_nemax, dm_nemaxx, label, i1
  integer       :: ii, ibb, i3, isize, latom, mu, index, l, m, maxx, lda, ldb, ldc, num, ilm2
  integer       :: jlo, ilm, li, mi, l1, mmax, lp, lp1, mpmax, jlop, mp, ms, mps, mtest, ly, lpy
  real*8        :: arg, arg1, arg2, arg3, argt, argt2, emist, fac, rmt2, gnt
  real*8        :: twopi, test1, sqrt2, sqrt3, sqfp, y, cin, jatombad, lbad
  !
  real*8     :: time_bl, time_bl_w, time_reduc, time_reduc_w, time_write, time_writeclm, time_writescf, time4_w, time3, time4
  real*8     :: time_write_w, time_writeclm_w, time_ilm, time_ilm_w, time_radprod, time_m, time_rad, time3_w
  real*8     :: time_rd_w, time_rd_c, time_atpar_w, time_atpar_c, time_writescf_w, time_r_w, time_r_c, trw, trc
  real*8     :: time_writeefg, time_writeeefg_w, time_writeefg_w, t1c, t1w, t2c, t2w, time1, time1_w, time2, time2_w
  real*8     :: time_dmf, time_dmfw, time_int, time_intw, time_ed, time_edw
  !
  complex*16 :: gint
  real(8) :: wkp
  complex*16 :: tsa, tsb, tsab, tsba, ts12(nloat), &
      &tse12(nloat), ts21(nloat), tse21(nloat), ts22(nloat,nloat)
  complex*16 :: tsuma,tsumb,tsumab,tsumba,tsum21(nloat),&
      &tsume21(nloat),tsum12(nloat),tsume12(nloat),tsum22(nloat,nloat)
  complex*16 :: sa, sb, sab, sba, s12(nloat), se12(nloat), s21(nloat), se21(nloat), s22(nloat,nloat)
  complex*16 :: vsa, vsb, vsab, vsba, vs12(nloat), vse12(nloat), vs21(nloat), vse21(nloat), vs22(nloat,nloat)
  real*8     :: facv, eqbad, ovrdiv, tc, efgfact, efg20, efg22, efg2m, qxx, qyy, qzz, qxy, vxx, vyy, vzz
  integer    :: lm_max, iks,ikp, is, itape,itape1,ios
  real*8     :: s,t,z, vnorm1, vnorm(nume), bkrot2(3), bkrot3(3), zero3(3), ff
  integer    :: nbands, nbandsx
  real*8     :: h_al(iblock),h_bl(iblock)
  complex*16 :: h_yl((lmax2+1)*(lmax2+1),iblock), h_alyl((lmax2+1)*(lmax2+1),iblock), h_blyl((lmax2+1)*(lmax2+1),iblock)
  real*8,     allocatable :: zw2(:)
  complex*16, allocatable :: alm(:,:),blm(:,:),clm(:,:,:), dh_alm(:,:,:), dh_blm(:,:,:), dh_clm(:,:,:,:)
  complex*16, allocatable :: aweight(:,:)
  ! interstitials
  integer    :: keigen(3,nmat), kmax(3), iff1, iff2, iff3, iff1t, iff2t, iff3t, isig, ierr, ia1, ia2, jj, jx, ij, k
  real*8     :: volin
  complex*16, allocatable :: fft(:,:,:),tp(:,:,:),sumfft(:,:,:),vsumfft(:,:,:),tloc(:,:,:)
  double precision,allocatable :: dwork(:)
  complex*16, allocatable :: cwork(:)
  complex*16, allocatable :: rho1(:),vrho1(:)
  complex*16, allocatable :: rhok_buf(:),vrhok_buf(:)
  complex*16, allocatable :: rhok(:),vrhok(:), asc(:,:)
  integer :: iikp, ivector
  logical :: tcompute
  real(8) dm_ef,dm_eorb,dm_esum
  character str*27,fname*256,sspin1*1
  
  twopi=two*pi

  test1=0.0d0                                                        
  xwt=0.0                                                           
  sqrt2=sqrt(2.0d0)                                                 
  sqrt3=sqrt(3.d0)                                                  
  sqfp=sqrt(4.d0*pi)                                                
  zero3(:) = 0.0


  call init_lohelp
  call init_xa(nat)
 
  if(iso.eq.2)then ! iso_orig=2 and ispin_pol=1
      vnorm1=0.5d0
  else
      vnorm1=1.d0
  endif

  if(iso_orig==2)then
      sspin1="1"
  else
      sspin1=sspin
  endif

  y=1.0d0                                                           
  do i=1,49
     j=2*i-1
     fct(j)=y  ! fct(i) = ((i+1)/2)!
     y=y*i
  end do
  
  cin=1.d0/137.0359895d0**2                       
  if (.not.rel) cin=4.0*1.0d-22                                     
  read(18,2032) iscf                                                     

  if (myrank.eq.master) then
     write(8,787) iscf                        
     write(8,*) '   norm of clm(r) =          '                        
     write(8,*)                                                 
     write(6,800)
  endif
  
  lfirst=1
  nnlo=0
  time_bl=zero; time_bl_w=zero; time_reduc=zero; time_reduc_w=zero
  time_write=zero; time_writeclm=zero; time_writescf=zero; time_write_w=zero
  time_writeclm_w=zero; time_writescf_w=zero; time_ilm=zero; time_ilm_w=zero
  time_radprod=zero; time_m=zero; time_rad=zero; time_rd_w=zero; time_rd_c=zero
  time_r_w=zero; time_r_c=zero; time_atpar_w=zero; time_atpar_c=zero; time_writeefg=zero; time_writeefg_w=zero
  time_dmf=zero; time_dmfw=zero; time_int=zero; time_intw=zero; time_ed=zero; time_edw=zero
 
  allocate(alm((lmax2+1)*(lmax2+1),nume),blm((lmax2+1)*(lmax2+1),nume), &
      &clm((lmax2+1)*(lmax2+1),nume,nloat))
  allocate(dh_alm(nume,(lmax2+1)*(lmax2+1),iso), &
      &dh_blm(nume,(lmax2+1)*(lmax2+1),iso), &
      &dh_clm(nume,(lmax2+1)*(lmax2+1),nloat,iso)) ! mar 27
  call init_charp(nume)
  call init_chard(nume)
  call init_charf(nume)
  
  lmx=lmax2                                                    
  allocate (e_store(0:lmax2,nat),elo_store(0:lomax,1:nloat,nat))
  itape=30
  if (vector_para) then
     if (nvector.ge.1) then
        fname = fvectors(1,3)
     else   ! no k-point to compute, but still need linearization energies
        fname = trim(vecfn(4))//'_1'
     endif
     open(itape,file=fname,status='old',form='formatted')
  else
     !------------ vector files need to be read from the beginning -------------------!
     rewind(itape)
  endif
  ! needed in atpar
  do i=1,nat
     read(itape,'(100(f9.5))') e_store(0:lmax2,i)
     read(itape,'(100(f9.5))') elo_store(0:lomax,1:nloat,i)
  enddo
  if (vector_para) then
     close(itape)
  else
     rewind(itape)
  endif
  
  call w_allocate0(nat)
  ! reads case.in2 to determin lm array
  ! and maximum lm -> lm_max=max(lm(1,:))
  lm_max=0
  do jatom=1,nat
     ! neg l means negative spherical harmonic comb (see kurki-suonio)   
     read(5,1003) ( (lm(j,jlm),j=1,2), jlm=1,ncom )
     ! find maximum lmmax
     do jlm=2,ncom
        if(lm(1,jlm).eq.0) exit
     enddo
     lmmax=jlm-1
     if (lmmax.gt.lm_max) lm_max = lmmax
     
     if (myrank.eq.master .or. fastfilesystem/=0) then
        write(6,*) ' atom',jatom,' ncomu',0,' lmmax',lmmax
        write(6, 1005) lmmax,((lm(j,jlm),j=1,2), jlm=1,lmmax)             
     endif
     if (myrank.eq.master) then
        write(21,13)   jatom,iatnr(jatom),(pos(i,lfirst),i=1,3),mult(jatom)
        write(21,1005) lmmax,((lm(j,jlm),j=1,2), jlm=1,lmmax)             
     endif
     w_lm(1:2,1:ncom,jatom) = lm(1:2,1:ncom)
     w_lmmax(jatom) = lmmax
  enddo
  
  ! allocates memory to store radial wave functions and 
  ! all results of subroutine "atpar", which is called outside
  ! the k - loop.
  call w_allocate(lm_max,nat)
  
  do jatom=1,nat
     ! here we calculate all radial wave functions and other quantities
     ! computed by atpar, which do not depend on k vector.
     
     lfirst=1
     do i=1,jatom-1
        lfirst=lfirst + mult(i)
     enddo

     ! calculate radial functions u(r), ue(r), ...                       
     call atpar (rel,nat,jatom,lfirst,cform,zz(jatom))               
     
     imax=jri(jatom)
     do i=1,imax
        r(i)=r0(jatom)*exp((i-1)*dx(jatom)) ! radial mesh
     enddo
     
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
     w_rrad1(1:nrad,0:lmax2,jatom) = rrad1(1:nrad,0:lmax2)
     w_rrad2(1:nrad,0:lmax2,jatom) = rrad2(1:nrad,0:lmax2)
     w_rade1(1:nrad,0:lmax2,jatom) = rade1(1:nrad,0:lmax2)
     w_rade2(1:nrad,0:lmax2,jatom) = rade2(1:nrad,0:lmax2) 
     
     w_p(0:lmax2,jatom) = p(0:lmax2)
     w_dp(0:lmax2,jatom) = dp(0:lmax2)
     w_pe(0:lmax2,jatom) = pe(0:lmax2)
     w_dpe(0:lmax2,jatom) = dpe(0:lmax2)
     w_pei(0:lmax2,jatom) = pei(0:lmax2)
     w_pi12lo(1:nloat,0:lomax,jatom) = pi12lo(1:nloat,0:lomax)
     w_pe12lo(1:nloat,0:lomax,jatom) = pe12lo(1:nloat,0:lomax)
     w_pr12lo(1:nloat,1:nloat,0:lomax,jatom) = pr12lo(1:nloat,1:nloat,0:lomax)
     
     w_r(1:nrad,jatom) = r(1:nrad)
  enddo

  !---------------------------------
  ! start loop over all k-points
  !---------------------------------
  call init_xa3(nmat,nnlo,nume,iso)
  if (vector_para) then
     pr_proc = sum(vectors(:,2))
     pr_procr = pr_proc
     if (myrank.eq.master .or. fastfilesystem/=0) write(6,'(a,i3,2x,a,i3)') 'pr_proc=', pr_proc, 'tot-k=', nkpt
  else
     pr_proc  = floor(nkpt/dble(nprocs)+0.999)  ! the maximum number of points calculated per processor
     ! correcting pr_proc is the number of k-points is not dividable
     pr_procr = pr_proc
     if ((myrank+1)*pr_proc .gt. nkpt) then
        if (nkpt-myrank*pr_proc.gt.0) then
           pr_procr = nkpt-myrank*pr_proc
        else
           pr_procr = 0
        endif
     endif
     if (myrank.eq.master .or. fastfilesystem/=0) write(6,'(a,i3,2x,a,i3)') 'pr_proc=', pr_proc, 'pr_procr=', pr_procr, 'tot-k=', nkpt
  endif

  allocate( pr_procs(nprocs) )
  call gather_procs(pr_procr, pr_procs, nprocs)
  if (myrank.eq.master .and. sum(pr_procs).ne.nkpt) then
     write(6,*) 'error: sum(pr_procs) should be nkpt, but is not', sum(pr_procs), nkpt
  endif

  ! ---- for interstitial ----
  call get_kmax(kmax,nat,nnlo)

  if (myrank.eq.master .or. fastfilesystem/=0) write(6,*) 'max k indices:',kmax(1),kmax(2),kmax(3)
  kmax(1)=kmax(1)+1
  kmax(2)=kmax(2)+1
  kmax(3)=kmax(3)+1
  ! -----  set iff values for fft calculation
  iff1=2*kmax(1)
  iff2=2*kmax(2)
  iff3=2*kmax(3)
  iff1=(iff1+1)*2
  iff2=(iff2+1)*2
  iff3=(iff3+1)*2
  iff1t=iff1
  iff2t=iff2
  iff3t=iff3
  call ifflim(iff1t,iff1)
  call ifflim(iff2t,iff2)
  call ifflim(iff3t,iff3)
  if (myrank.eq.master .or. fastfilesystem/=0) write(6,*) 'n,iff1,iff2,iff3',n-nnlo,iff1,iff2,iff3
  !
  allocate (fft(iff1,iff2,iff3),tp(iff1,iff2,iff3),&
      &sumfft(iff1,iff2,iff3),vsumfft(iff1,iff2,iff3),tloc(iff1,iff2,iff3))
  allocate (dwork(4*(iff1+iff2+iff3)+15),cwork(iff1+iff2+iff3))
  fft=0.0
  tp=0.0
  sumfft=0.0
  vsumfft=0.0
  tloc=0.0
  ! --- finished for interstitial

  w_vrholm(:,:,:)=0.0
  w_rholm(:,:,:)=0.0
  w_xwt1(:,:)=0.0
     
  rewind 5
  read(5,fmt='(//,1a)') cdummy
  
  rewind(30)
  rewind(18)
  read(18,2032) iscf                                                     
  
!!! ---------- preparation of arrays for paralel executaion --------------                                                     
  if (myrank.eq.master .or. fastfilesystem/=0) write(6,'(a,i3,2x,a,i3)') 'pr_proc=', pr_proc, 'tot-k=', nkpt-1
  if (myrank.eq.master .or. fastfilesystem/=0) print *, 'mpr_proc=', pr_proc, 'myrank=', myrank, 'nprocs=', nprocs, 'nkpt=', nkpt

  iikp=0
  do ivector=1,nvector
     if(gp%lkpvec)then
       write(str,'(i3)')gp%kvec(ivector,1)
     else
       write(str,'(i3)')gp%myrank
     endif
     call gh5_open_r('KSWT_'//trim(adjustl(str))//'.h5', f_id)
     if(ivector.eq.1)then
        if(sspin=="1")then
            call gh5_read(dm_ef,'/e_fermi',f_id)
            call gh5_read(dm_eorb,'/e_gamma_dc',f_id)
            write(21,'(":FER gutzwiller fermi energy",f10.5)')dm_ef
            write(21,'(":EORB ",f20.9)')dm_eorb
        endif
        if(ispin_pol==1)then
            call gh5_read(dm_esum,'/e_band',f_id)
        else
            call gh5_read(dm_esum,'/e_band_spin'//sspin,f_id)
        endif
        write(21,'(":SUM sum of gutz band energies ",f20.9)')dm_esum
     endif

     do is=1,iso    !------ over up/dn ---------------------!
        itape1=8+is
        if (vector_para) then
           open(itape1,file=fvectors(ivector,is),status='old',form='unformatted')
        else
           rewind(itape1) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        do i=1,nat
           read(itape1) emist !---- at the beginninge we have linearization energies --!
           read(itape1) emist !---- we just skip them ---------------------------------!
        enddo
     end do
   
     ! normso_x
     if (vector_para) then
        open(12,file=fvectors(ivector,5), &
                &status='old',form='formatted')
     endif

     do iks=1,vectors(ivector,2) ! kpoint loop begin
        if (vector_para) then
           ikp = vectors(ivector,3)+iks  ! successive index in k-point table from case.klist
           iikp = iikp+1                 ! successive index in k-point table on this processor
        else
           ikp = iks                     ! successive index in k-point table from case.klist
           !--- we need to go over all k-points even though we will compute only some of them on this processor.
           !--- this is because we need to read vector file sequentially.
           iikp = ikp-myrank*pr_proc            ! the index of the point to compute. if negative, do not compute!
        endif
        tcompute=.false.   

        if (iikp.gt.0) tcompute=.true.       
        ! if tcompute is true, the point needs to be computed.                                                           
        if (iikp.gt.pr_proc) exit            
        ! processor finished. the rest of the points will be taken care of by the other processors.                      

        call cputim(t1c)
        call walltim(t1w)
        !!! need to read both spins!
        call read_vec_spin(ikp, e, as, as_lo, kx, ky, kz, vnorm, &
            &kxlo, kylo, kzlo, &
            &bkx, bky, bkz, bkxlo, bkylo, bkzlo, more_kpoints, n0, &
            &emin, nemin, iso, nmat, nume, nnlo)
        if (.not.tcompute) cycle
        if (.not.more_kpoints) exit
        call cputim(t2c)
        call walltim(t2w)
        time_rd_c=time_rd_c+t2c-t1c
        time_rd_w=time_rd_w+t2w-t1w
   
        do jatom=1,nat
           isize=n0-nnlo
           call harmon(isize,bkx,bky,bkz,lmax2,fj(:,:,jatom),dfj(:,:,jatom),rmt(jatom))
        enddo
   
        call cputim(t1c)
        call walltim(t1w)

        write(str,'(i7)')ikp
        call gh5_read(dm_nemin,'/IKP_'//trim(adjustl(str))// &
                &'/nemin_spin'//sspin1,f_id)
        call gh5_read(dm_nemax,'/IKP_'//trim(adjustl(str))// &
                &'/nemax_spin'//sspin1,f_id)
        nbands = dm_nemax-dm_nemin+1
        call cputim(t2c)
        call walltim(t2w)
        
        time_dmf=time_dmf+t2c-t1c
        time_dmfw=time_dmfw+t2w-t1w
        
        call cputim(t1c)
        call walltim(t1w)
        allocate( aweight(nbands,nbands), zw2(nbands) )
        call gh5_read(wkp,'/IKP_'//trim(adjustl(str))// &
                &"/wt",f_id)
        call gh5_read(aweight,nbands,nbands,'/IKP_'//trim(adjustl(str))// &
                &'/KSWT_SPIN'//sspin1,f_id)
        if(ispin_pol==1)then
            wkp=wkp*2
        endif
        call dmft_weights(zw2, aweight, nbands)
       
        call cputim(t2c)
        call walltim(t2w)
        time_ed=time_ed+t2c-t1c
        time_edw=time_edw+t2w-t1w
        
        ! weights for all bands, which are needed below
        weight=0
        weight(:dm_nemin-1) = wkp
        weight(dm_nemin:dm_nemax) = zw2/vnorm1
        dm_nemaxx=dm_nemax
        do num=dm_nemin,dm_nemax
           ! zw2 already in descending order
           if(abs(zw2(num-dm_nemin+1)).le.1.d-8)then
               dm_nemaxx=num-1
               exit
           endif 
        enddo
       
        do num=nemin,dm_nemaxx
           xwt=xwt+weight(num)*vnorm1*vnorm(num)
        enddo
 
        write(*,'(i3,a,1x,i4,1x,i4,1x,a,i3,1x,a,i3,1x,a,i3)') &
            &myrank, ') finished k-point', ikp, iikp, &
            &'with dm_nemin=', dm_nemin, 'dm_nemax=', dm_nemax, &
            &'dm_nemaxx', dm_nemaxx
   
        deallocate(zw2)
        
        nbandsx = dm_nemaxx-dm_nemin+1
        do jatom=1,nat
           ! quantities computed by "atpar" which do not depend on k-point
           lfirst = w_lfirst(jatom)
           lm(1:2,1:ncom) = w_lm(1:2,1:ncom,jatom)
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
           rrad1(1:nrad,0:lmax2) = w_rrad1(1:nrad,0:lmax2,jatom)
           rrad2(1:nrad,0:lmax2) = w_rrad2(1:nrad,0:lmax2,jatom)
           rade1(1:nrad,0:lmax2) = w_rade1(1:nrad,0:lmax2,jatom)
           rade2(1:nrad,0:lmax2) = w_rade2(1:nrad,0:lmax2,jatom)
           p(0:lmax2) = w_p(0:lmax2,jatom)
           dp(0:lmax2) = w_dp(0:lmax2,jatom)
           pe(0:lmax2) = w_pe(0:lmax2,jatom)
           dpe(0:lmax2) = w_dpe(0:lmax2,jatom)
           pei(0:lmax2) = w_pei(0:lmax2,jatom)
           pi12lo(1:nloat,0:lomax) = w_pi12lo(1:nloat,0:lomax,jatom)
           pe12lo(1:nloat,0:lomax) = w_pe12lo(1:nloat,0:lomax,jatom)
           pr12lo(1:nloat,1:nloat,0:lomax) = w_pr12lo(1:nloat,1:nloat,0:lomax,jatom)
           r(1:nrad) = w_r(1:nrad,jatom)
              
           ! zero koeff. for charge analysis
           call zero_charp(nemin,dm_nemax)
           call zero_chard(nemin,dm_nemax)
           call zero_charf(nemin,dm_nemax)
           call zerotc_lohelp(nemin,dm_nemax)
           tc100(0:lmax2,nemin:dm_nemax)=0.0
           tca100(0:lmax2,nemin:dm_nemax)=0.0
           tcb100(0:lmax2,nemin:dm_nemax)=0.0
           
           ! calculate alm, blm                                                
           fac=4.0d0*pi*rmt(jatom)**2/sqrt(vol)                              
           latom=lfirst-1                                                    
   
           equiv_atom_loop: do mu=1,mult(jatom)
              
              latom=latom+1
   
              do is=1,iso ! iso == 1           
                 
                 alm(:,:)=0.0
                 blm(:,:)=0.0
                 clm(:,:,:)=0.0
   
                 call cputim(time1)
                 call walltim(time1_w)
              
                 blocked_loop: do ii=1,n0-nnlo,iblock !iblock*ibpp,iblock
                    ibb=min(iblock,n0-(nlo+nlon+nlov)-ii+1)
                    if(ibb.le.0) exit
                    i3=0
                    do i=ii,min(ii+iblock-1,n0-nnlo) !ii+iblock-1
                       !---------  rotates ylm(k+k) to ylm(k'+k) where k' is in star of irreducible k. ------------!
                       i3=i3+1
                       bk(1)=bkx(i) !-----  reciprocal vector and irreducible vector: g=k+k ----!
                       bk(2)=bky(i)
                       bk(3)=bkz(i)
                       ! bkrot = r_g.(k+k)
                       call rotate (bk,rotij(1,1,latom),bkrot)
                       !---- br1 transforms integer reciprocal lattice vectors, as given in the vectorlist of lapw1, into cartesian system ----!
                       bk(1)=bkrot(1)*br1(1,1)+bkrot(2)*br1(1,2)+bkrot(3)*br1(1,3)   
                       bk(2)=bkrot(1)*br1(2,1)+bkrot(2)*br1(2,2)+bkrot(3)*br1(2,3)   
                       bk(3)=bkrot(1)*br1(3,1)+bkrot(2)*br1(3,2)+bkrot(3)*br1(3,3)   
                       !---- bkrloc = rotloc.r_g.(k+k),  is rotates rotloc.r.(k+k) . rotation rotloc entered by user
                       call rotate (bk,rotloc(1,1,jatom),bkrloc)
                       !---- ylm = y_{l}(rotloc.r_g.(k+k))
                       call ylm (bkrloc,lmax2,yl)
                       arg1=bkrot(1)*(pos(1,lfirst))*twopi
                       arg2=bkrot(2)*(pos(2,lfirst))*twopi
                       arg3=bkrot(3)*(pos(3,lfirst))*twopi
                       ! argt = (k+k)*tau(isym) 
                       argt=(bkx(i)*tauij(1,latom)+bky(i)*tauij(2,latom)+bkz(i)*tauij(3,latom))*twopi
                       ! phsehl = e^{i*2pi*( (r_g.(k+k)) *  r(iatom) + (k+k)*tau(isym))}
                       phshel=exp( imag*(arg1+arg2+arg3+argt) )
                       index=0
                       do  l=0,lmx
                          maxx=2*l+1
                          do  m=1,maxx
                             index=index+1
                             h_yl(index,i3)=conjg(yl(index))*phshel
                          enddo
                       enddo
                    enddo
                    
                    index=0
                    rmt2=1.d0/(rmt(jatom)**2)
                    do l=0,lmx
                       i3=0
                       do i=ii,ii+iblock-1
                          if(i.gt.n0-(nlo+nlon+nlov)) exit
                          i3=i3+1
                          if(lapw(l)) then
                             h_al(i3)=dfj(l,i,jatom)*pe(l)-fj(l,i,jatom)*dpe(l) 
                             h_bl(i3)=fj(l,i,jatom)*dp(l)-dfj(l,i,jatom)*p(l)
                          else
                             h_al(i3) = rmt2*fj(l,i,jatom)/p(l)
                             h_bl(i3) = 0.d0
                          endif
                       enddo
                       maxx=2*l+1
                       do m=1,maxx
                          index=index+1
                          i3=0
                          do i=ii,ii+iblock-1
                             if(i.gt.n0-(nlo+nlon+nlov)) exit
                             i3=i3+1
                             h_alyl(index,i3)=h_al(i3)*h_yl(index,i3)
                             h_blyl(index,i3)=h_bl(i3)*h_yl(index,i3)
                          enddo
                       enddo
                    enddo
                    !
                    lda=(lmax2+1)*(lmax2+1)
                    ldc=lda
                    ldb=nmat 
                    call zgemm('n','n',index,dm_nemax-nemin+1,ibb,&
                        &(1.d0,0.d0),h_alyl,lda,as(ii,nemin,is),ldb, &
                        &(1.d0,0.d0),alm(1,nemin),ldc)
                    call zgemm('n','n',index,dm_nemax-nemin+1,ibb,&
                        &(1.d0,0.d0),h_blyl,lda,as(ii,nemin,is),ldb, &
                        &(1.d0,0.d0),blm(1,nemin),ldc)
                    
                 enddo blocked_loop
                 
              
                 call cputim(time2)
                 call walltim(time2_w)
                 time_bl=time_bl+time2-time1
                 time_bl_w=time_bl_w+time2_w-time1_w
                 
                 call cputim(time1)
                 call walltim(time1_w)
                 
                 ldc=(lmax2+1)*(lmax2+1)
   
                 if (nlo.ne.0) call lomain (rotloc(:,:,jatom),is,nemin,dm_nemax,lfirst,latom,n0,jatom,alm,blm,clm)
              
                 index=0
                 do l=0,lmax2
                    maxx=2*l+1
                    cfac=imag**l
                    cfac=cfac*fac
                    do m=1,maxx
                       index=index+1
                       do num=nemin,dm_nemax
                          alm(index,num)=alm(index,num)*cfac
                          blm(index,num)=blm(index,num)*cfac
                       enddo
                       do jlo=1,ilo(l)
                          do num=nemin,dm_nemax
                             clm(index,num,jlo)=clm(index,num,jlo)*cfac
                          enddo
                       enddo
                    enddo
                 enddo
                 ! chenging alm,blm,clm to make them dmft-like
                 alm(:,dm_nemin:dm_nemaxx)= matmul( alm(:,dm_nemin:dm_nemax), aweight(:nbands,:nbandsx) )
                 blm(:,dm_nemin:dm_nemaxx)= matmul( blm(:,dm_nemin:dm_nemax), aweight(:nbands,:nbandsx) )
                 do jlo=1,nloat
                    clm(:,dm_nemin:dm_nemaxx,jlo)= matmul( clm(:,dm_nemin:dm_nemax,jlo), aweight(:nbands,:nbandsx) )
                 enddo
                 
                 dh_alm(:,:,is) = transpose( alm(:,:) )
                 dh_blm(:,:,is) = transpose( blm(:,:) )
                 do jlo=1,nloat
                    dh_clm(:,:,jlo,is) = transpose( clm(:,:,jlo) )
                 enddo
                 
              enddo  !!! mar 27, 2010: loop over spin
   
              call cputim(time2)
              call walltim(time2_w)
              time_reduc=time_reduc+time2-time1
              time_reduc_w=time_reduc_w+time2_w-time1_w
              
              !.....c(l,m) to be calculated                                           
              call cputim(time1)
              call walltim(time1_w)
              
              ilm_loop: do ilm=1,lmmax
                 li=iabs(lm(1,ilm))                                                
                 mi=lm(2,ilm)                                                      
                 !     see kurki-suoni  factor for  lm+ (1,0), lm-(0,-1), m.ne.2n *(-1)  
                 imag1=(1.0,0.0)                                                   
   
                 if(lm(1,ilm).lt.0) imag1= imag
                 if(mod(mi,2).eq.1) imag1=-imag1
                 
                 l_sum: do l=0,lmx
                    l1=l+1
                    mmax=2*l+1
                    lp_sum: do lp=0,lmx
                       call cputim(time3)
                       lp1=lp+1
                       if(notri(li,l,lp).lt.0) cycle
                       mpmax=2*lp+1
                       imax=jri(jatom)
                       
                       ! radial functions are constructed
                       do i=1,imax
                          ua(i) =rrad1(i,l)*rrad1(i,lp)+cin*rrad2(i,l)*rrad2(i,lp)    
                          ub(i) =rade1(i,l)*rade1(i,lp)+cin*rade2(i,l)*rade2(i,lp)    
                          uab(i)=rrad1(i,l)*rade1(i,lp)+cin*rrad2(i,l)*rade2(i,lp)   
                          uba(i)=rrad1(i,lp)*rade1(i,l)+cin*rrad2(i,lp)*rade2(i,l)   
                       enddo
                       ! radial functions as product of local orbitals and u,dot{u}
                       do jlo=1,ilo(l)
                          if(loor(jlo,l)) then
                             do i=1,imax
                                u21(i,jlo) =a1lo(i,jlo,l)*rrad1(i,lp)+cin*b1lo(i,jlo,l)*rrad2(i,lp)    
                                ue21(i,jlo)=a1lo(i,jlo,l)*rade1(i,lp)+cin*b1lo(i,jlo,l)*rade2(i,lp)    
                             enddo
                          endif
                       enddo
                       do jlop=1,ilo(lp)
                          if(loor(jlop,lp)) then
                             do i=1,imax
                                u12(i,jlop) =rrad1(i,l)*a1lo(i,jlop,lp)+cin*rrad2(i,l)*b1lo(i,jlop,lp)    
                                ue12(i,jlop)=rade1(i,l)*a1lo(i,jlop,lp)+cin*rade2(i,l)*b1lo(i,jlop,lp) 
                             enddo
                          endif
                       enddo
                       ! radial functions for local orbitals
                       do jlo=1,ilo(l)
                          do jlop=1,ilo(lp) 
                             if(loor(jlo,l).and.loor(jlop,lp)) then
                                do i=1,imax
                                   u22(i,jlop,jlo)=a1lo(i,jlo,l)*a1lo(i,jlop,lp)+cin*b1lo(i,jlo,l)*b1lo(i,jlop,lp) 
                                enddo
                             endif
                          enddo
                       enddo
                       call cputim(time4)
                       time_radprod=time_radprod+time4-time3
                       !.....m sum                           
   
                       sa=0.0d0; sb=0.0d0; sab=0.0d0; sba=0.0d0; s12=0.0d0; 
                       se12=0.0d0; s21=0.0d0; se21=0.0d0; s22=0.0d0
                       vsa=0.0d0; vsb=0.0d0; vsab=0.0d0; vsba=0.0d0; vs12=0.0d0
                       vse12=0.0d0; vs21=0.0d0; vse21=0.0d0; vs22=0.0d0
                       
                       
                       call cputim(time3)
                       do num=nemin,dm_nemaxx
                          tsuma = 0.0; tsumb = 0.0; tsumab= 0.0; tsumba= 0.0
                          tsum21 (:ilo(l))  = 0.0;  tsume21(:ilo(l))  = 0.0
                          tsum12 (:ilo(lp)) = 0.0;  tsume12(:ilo(lp)) = 0.0
                          tsum22 (:ilo(lp),:ilo(l))=0.0
                          
                          m=-(l+1) 
                          mp=-(lp+1)
                          do ms=1,2*l+1
                             m=m+1
                             mp=-lp1
                             do mps=1,2*lp+1
                                mp=mp+1
                                mtest=-m+mi+mp
                                if(mtest.ne.0) cycle
                                ly  = l*(l+1)+m+1
                                lpy = lp*(lp+1)+mp+1
                                gnt=gaunt(l,li,lp,m,mi,mp)
                                gint = imag1*gnt
                                
                                do is=1,iso
                                   tsa = dh_alm(num,ly,is)*dconjg(dh_alm(num,lpy,is))*gint
                                   tsb = dh_blm(num,ly,is)*dconjg(dh_blm(num,lpy,is))*gint
                                   tsab= dh_alm(num,ly,is)*dconjg(dh_blm(num,lpy,is))*gint
                                   tsba= dh_blm(num,ly,is)*dconjg(dh_alm(num,lpy,is))*gint
                                   do jlo=1,ilo(l)
                                      ts21(jlo)  = dh_clm(num,ly,jlo,is)*dconjg(dh_alm(num,lpy,is))*gint
                                      tse21(jlo) = dh_clm(num,ly,jlo,is)*dconjg(dh_blm(num,lpy,is))*gint
                                   enddo
                                   do jlop=1,ilo(lp)
                                      ts12(jlop)  = dh_alm(num,ly,is)*dconjg(dh_clm(num,lpy,jlop,is))*gint
                                      tse12(jlop) = dh_blm(num,ly,is)*dconjg(dh_clm(num,lpy,jlop,is))*gint
                                   enddo
                                   do jlo=1,ilo(l)
                                      do jlop=1,ilo(lp)
                                         ts22(jlop,jlo) = dh_clm(num,ly,jlo,is)*dconjg(dh_clm(num,lpy,jlop,is))*gint
                                      enddo
                                   enddo
                                   
                                   tsuma  = tsuma  + tsa/iso
                                   tsumb  = tsumb  + tsb/iso
                                   tsumab = tsumab + tsab/iso
                                   tsumba = tsumba + tsba/iso
                                   do jlo=1,ilo(l)
                                      tsum21 (jlo) = tsum21 (jlo) + ts21 (jlo)/iso
                                      tsume21(jlo) = tsume21(jlo) + tse21(jlo)/iso
                                   enddo
                                   do jlop=1,ilo(lp)
                                      tsum12 (jlop) = tsum12 (jlop) + ts12(jlop)/iso
                                      tsume12(jlop) = tsume12(jlop) + tse12(jlop)/iso
                                   enddo
                                   do jlo=1,ilo(l)
                                      do jlop=1,ilo(lp)
                                         tsum22(jlop,jlo) = tsum22(jlop,jlo) + ts22(jlop,jlo)/iso
                                      enddo
                                   enddo
                                enddo
                                
                             enddo
                          enddo
                          
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
                          do jlo=1,ilo(l)
                             s21(jlo)   = s21(jlo)  + tsum21 (jlo) * wkp
                             se21(jlo)  = se21(jlo) + tsume21(jlo) * wkp
                             vs21(jlo)  = vs21(jlo) + tsum21 (jlo) * facv
                             vse21(jlo) = vse21(jlo)+ tsume21(jlo) * facv
                          enddo
                          do jlop=1,ilo(lp)
                             s12  (jlop) = s12  (jlop) + tsum12 (jlop) * wkp
                             se12 (jlop) = se12 (jlop) + tsume12(jlop) * wkp
                             vs12 (jlop) = vs12 (jlop) + tsum12 (jlop) * facv
                             vse12(jlop) = vse12(jlop) + tsume12(jlop) * facv
                          enddo
                          do jlo=1,ilo(l)
                             do jlop=1,ilo(lp)
                                s22 (jlop,jlo) = s22 (jlop,jlo) + tsum22(jlop,jlo) * wkp
                                vs22(jlop,jlo) = vs22(jlop,jlo) + tsum22(jlop,jlo) * facv
                             enddo
                          enddo
   
                          suma  (num) = tsuma
                          sumb  (num) = tsumb
                          sumab (num) = tsumab
                          sumba (num) = tsumba
                          sum21 (num,:ilo(l))  = tsum21 (:ilo(l))
                          sume21(num,:ilo(l))  = tsume21(:ilo(l))
                          sum12 (num,:ilo(lp)) = tsum12 (:ilo(lp))
                          sume12(num,:ilo(lp)) = tsume12(:ilo(lp))
                          sum22 (num,:ilo(lp),:ilo(l)) = tsum22 (:ilo(lp),:ilo(l))
                       enddo
                       call cputim(time4)
                       time_m=time_m+time4-time3
                       
                       if(li.eq.0) then
                          call csplit(nemin,dm_nemaxx,l,jatom,mu,alm,blm,clm,coord) 
                       endif
                       call cputim(time3)
                       imax=jri(jatom)
                       do i=1,imax
                          tc_buf(i)=sa*ua(i)+sb*ub(i)+sab*uab(i)+sba*uba(i)
                          vtc_buf(i)=vsa*ua(i)+vsb*ub(i)+vsab*uab(i)+vsba*uba(i)
                       enddo
                       do jlo=1,ilo(l)
                          do i=1,imax
                             tc_buf(i)=tc_buf(i)+s21(jlo)*u21(i,jlo)+se21(jlo)*ue21(i,jlo)
                             vtc_buf(i)=vtc_buf(i)+vs21(jlo)*u21(i,jlo)+vse21(jlo)*ue21(i,jlo)
                          enddo
                       enddo
                       do jlop=1,ilo(lp)
                          do i=1,imax
                             tc_buf(i)=tc_buf(i)+s12(jlop)*u12(i,jlop)+se12(jlop)*ue12(i,jlop)
                             vtc_buf(i)=vtc_buf(i)+vs12(jlop)*u12(i,jlop)+vse12(jlop)*ue12(i,jlop)
                          enddo
                       enddo
                       do jlo=1,ilo(l)
                          do jlop=1,ilo(lp)
                             do i=1,imax
                                tc_buf(i)=tc_buf(i)+s22(jlop,jlo)*u22(i,jlop,jlo)
                                vtc_buf(i)=vtc_buf(i)+vs22(jlop,jlo)*u22(i,jlop,jlo)
                             enddo
                          enddo
                       enddo
                       do i=1,imax
                          w_rholm(i,ilm,jatom)=w_rholm(i,ilm,jatom)+tc_buf(i)/mult(jatom)
                          w_vrholm(i,ilm,jatom)=w_vrholm(i,ilm,jatom)+vtc_buf(i)/mult(jatom)
                       enddo
                       call cputim(time4)
                       time_rad=time_rad+time4-time3
                    enddo lp_sum
                 enddo l_sum
              enddo ilm_loop
              
              call cputim(time2)
              call walltim(time2_w)
              time_ilm=time_ilm+time2-time1
              time_ilm_w=time_ilm_w+time2_w-time1_w
           enddo equiv_atom_loop
   
           call psplit(w_xwt1(:,jatom),jatom,nemin,dm_nemaxx,test1,eqbad,jatombad,lbad)
        enddo
   
        call cputim(time1)
        call walltim(time1_w)
        ! interstitial charge
        ! fft eigenvectors and accumulate weighted square in sumfft
        keigen(1,:)=kx(:)
        keigen(2,:)=ky(:)
        keigen(3,:)=kz(:)
        isig=-1
   
        allocate( asc(n0-nnlo,1:dm_nemaxx) )
   
        !!!! mar 27, 2010: over spins in interstitals
        do is=1,iso
           asc(:,:) = as(:n0-nnlo,:dm_nemaxx,is)
           asc(:,dm_nemin:dm_nemaxx) = matmul( as(:n0-nnlo,dm_nemin:dm_nemax,is), aweight(:nbands,:nbandsx) )
   
           do num=nemin,dm_nemaxx
              !  puts eigenvector into fft array for fft
              call setfft1c(n0-nnlo,iff1,iff2,iff3,asc(:,num),fft,keigen)
              !  computes a(k,r) = fft(a(k,k))
              call c3fft(iff1,iff2,iff3,fft,iff1,iff2,isig,cwork,dwork,ierr)  
              !  |a(k,r)|^2 w(k)
              call sumupfft(iff1,iff2,iff3, weight(num)/iso,       fft,sumfft)
              ! -|a(k,r)|^2 w(k)*e(k)
              call sumupfft(iff1,iff2,iff3,-weight(num)*e(num)/iso,fft,vsumfft)   
           enddo
        enddo
   
        deallocate( asc, aweight )
   
        call cputim(time2)
        call walltim(time2_w)
        time_int =time_int +time2-time1
        time_intw =time_intw +time2_w-time1_w
     enddo     !  kpoint loop end
     call gh5_close(f_id) 
     do is=1,iso    !------ over up/dn ---------------------!
        itape1=8+is
        if (vector_para) then
           close(itape1)
        else
           rewind(itape1)
        endif
     enddo
     if (vector_para) then
        close(12)
     endif
  enddo

  call reduce_mpi(xwt, w_rholm, w_vrholm, w_xwt1, sumfft, vsumfft, &
      &nrad, lm_max, nat, iff1, iff2, iff3)

  if (myrank.eq.master) then

     do jatom=1,nat
        lm(1:2,1:ncom) = w_lm(1:2,1:ncom,jatom)
        lmmax = w_lmmax(jatom)
        r(:) = w_r(:,jatom)
        
        call cputim(time1)
        call walltim(time1_w)
        
        do ilm=1,lmmax
           ovrdiv=1.0
           li=lm(1,ilm)
           mi=lm(2,ilm)
           if(li.ne.0) then
              if(mi.ne.0) then
                 imax=jri(jatom)
                 do i=1,imax
                    w_vrholm(i,ilm,jatom)=w_vrholm(i,ilm,jatom)*sqrt2
                    w_rholm(i,ilm,jatom)=w_rholm(i,ilm,jatom)*sqrt2
                 enddo
              endif
              cycle
           endif
           call charge(r,ovrdiv,sqfp,w_rholm(1,ilm,jatom),dx(jatom),jri(jatom),tc)
           write(6,207)  jatom,tc
           write(21,720) jatom,jatom,tc
           if(isplit(jatom).eq.1) then
              write(6,250) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
              write(21,250) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
           else if(isplit(jatom).eq.2) then                               
              write(6,251) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
              write(21,251) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
           else if(isplit(jatom).eq.3) then                               
              write(6,252) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
              write(21,252)jatom,jatom, jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
           else if(isplit(jatom).eq.-2) then                              
              write(6,253) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
              write(21,253) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
           else if(isplit(jatom).eq.4) then                               
              write(6,254) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
              write(21,254) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
           else if(isplit(jatom).eq.5) then                               
              write(6,255) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
              write(21,255) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
           else if(isplit(jatom).eq.6) then                               
              write(6,256) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
              write(21,256) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
           else if(isplit(jatom).eq.8) then                               
              write(6,257) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
              write(21,257) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,14)
           else if(isplit(jatom).eq.15) then                               
              write(6,258) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,21)
              write(21,258) jatom,jatom,jatom,(w_xwt1(i,jatom),i=0,3),(w_xwt1(i,jatom),i=7,21)
           endif
           !
           if(mi.ne.0) then
              imax=jri(jatom)                                                   
              do i=1,imax                                                    
                 w_vrholm(i,ilm,jatom)=w_vrholm(i,ilm,jatom)*sqrt2
                 w_rholm(i,ilm,jatom)=w_rholm(i,ilm,jatom)*sqrt2
              enddo
           endif
        enddo
     
        call cputim(time3)
        call walltim(time3_w)
        
        ! writting case.clmval
        write(8,1990) jatom
        write(8,2001) lmmax
        do ilm=1,lmmax
           li=lm(1,ilm)
           mi=lm(2,ilm)
           write(8,2011) li,mi
           write(8,2022) ( w_rholm(i,ilm,jatom), i=1,imax )
           write(8,2031)
        end do
        write(8,2030)
        
        call cputim(time4)
        call walltim(time4_w)
        time_writeclm=time_writeclm+time4-time3
        time_writeclm_w=time_writeclm_w+time4_w-time3_w
        
        do ilm1=1,lmmax
           if (lm(1,ilm1).eq.2.and.lm(2,ilm1).eq.0) then
              do i=1,imax                                                
                 w_rholm(i,ilm1,jatom)=w_rholm(i,ilm1,jatom)/r(i)**3
              enddo
              label=1
              do ilm2=ilm1+1,lmmax
                 if (lm(1,ilm2).eq.2.and.lm(2,ilm2).eq.2) then
                    do i=1,imax
                       w_rholm(i,ilm2,jatom)=w_rholm(i,ilm2,jatom)/r(i)**3
                    enddo
                    label=2
                    do ilm3=ilm2+1,lmmax
                       if (lm(1,ilm3).eq.-2.and.lm(2,ilm3).eq.2) then
                          do i=1,imax
                             w_rholm(i,ilm3,jatom)=w_rholm(i,ilm3,jatom)/r(i)**3
                          enddo
                          label=3
                          goto 7372
                       endif
                    enddo
                    goto 7372
                 endif
              enddo
              goto 7372
           endif
        enddo
     
        goto 7371
     
7372 continue 
     
        efgfact=0.02997925
        fac=0.8d0*pi*324.14d0
        do i1=0,jri(jatom)-5,1000   ! 7374
           ! integrates: tc = integrate(r*rholm)
           call charge(r,1.d0,1.d0,w_rholm(1,ilm1,jatom),dx(jatom),jri(jatom)-i1,tc) 
           tc=tc*fac
           efg20=tc
           if (label.gt.1) then
              call charge(r,1.d0,1.d0,w_rholm(1,ilm2,jatom),dx(jatom),jri(jatom)-i1,tc) 
              tc=tc*fac
              efg22=tc
              if (label.eq.3) then
                 call charge(r,1.d0,1.d0,w_rholm(1,ilm3,jatom),dx(jatom),jri(jatom)-i1,tc)
                 tc=tc*fac
                 efg2m=tc
                 qxx=(-efg20/sqrt(3.d0)+efg22)*sqrt(15.d0/4.d0/pi)
                 qyy=(-efg20/sqrt(3.d0)-efg22)*sqrt(15.d0/4.d0/pi)
                 qzz=(efg20*2.d0/sqrt(3.d0))*sqrt(15.d0/4.d0/pi)
                 qxy=efg2m*sqrt(15.d0/4.d0/pi)
                 qxx=-qxx*efgfact
                 qxy=-qxy*efgfact
                 qyy=-qyy*efgfact
                 qzz=-qzz*efgfact
                 if (i1.eq.0) write(6,3019)
                 if (i1.eq.0) write(21,3019)
                 write(6,3121) jatom,qxx,qxy,qyy,qzz,r(jri(jatom)-i1)
                 if (i1.eq.0) write(21,3121) jatom,qxx,qxy,qyy,qzz,r(jri(jatom)-i1)      
              else
                 vxx=(-efg20/sqrt(3.d0)+efg22)*sqrt(15.d0/4.d0/pi)
                 vyy=(-efg20/sqrt(3.d0)-efg22)*sqrt(15.d0/4.d0/pi)
                 vzz=(efg20*2.d0/sqrt(3.d0))*sqrt(15.d0/4.d0/pi)
                 if (i1.eq.0) write(6,3017)
                 vxx=-vxx*efgfact
                 vyy=-vyy*efgfact
                 vzz=-vzz*efgfact
                 if (i1.eq.0) write(21,3017)
                 write(6,3021) jatom,vxx,vyy,vzz,r(jri(jatom)-i1)
                 if (i1.eq.0) write(21,3021) jatom,vxx,vyy,vzz,r(jri(jatom)-i1)
              endif
           else
              efg20=efg20*2.d0/sqrt(3.d0)*sqrt(15.d0/4.d0/pi)
              efg20=-efg20*efgfact
              write(6,7207)  jatom,efg20,r(jri(jatom)-i1)
              if (i1.eq.0) write(21,7720) jatom,jatom,efg20,r(jri(jatom)-i1)
           endif
        enddo !7374       continue                                                          
     
     
7371 continue                                     
     
        rewind fh_vec
        
        call cputim(time3) 
        call walltim(time3_w)
        time_writeefg=time_writeefg+time3-time4
        time_writeefg_w=time_writeefg_w+time3_w-time4_w
        
        call cputim(time2)
        call walltim(time2_w)
        time_writescf=time_writescf+time2-time3
        time_write=time_write+time2-time1
        time_writescf_w=time_writescf_w+time2_w-time3_w
        time_write_w=time_write_w+time2_w-time1_w
     enddo
     
     write(6,*)
     write(6,'(a,2f8.1)') '   atpar              (cpu,wall):',time_atpar_c,time_atpar_w
     write(6,'(a,2f8.1)') '   readvec            (cpu,wall):',time_rd_c,time_rd_w
     write(6,'(a,2f8.1)') '   readvec, read only (cpu,wall):',time_r_c,time_r_w
     write(6,'(a,2f8.1)') '   blocked loop            (cpu):',time_bl,time_bl_w
     write(6,'(a,2f8.1)') '   reduced                 (cpu):',time_reduc,time_reduc_w
     write(6,'(a,2f8.1)') '   ilm-tot loop            (cpu):',time_ilm,time_ilm_w
     write(6,'(a,f8.1)')  '   ilm-1 loop              (cpu):',time_radprod
     write(6,'(a,f8.1)')  '   ilm-2 loop              (cpu):',time_m
     write(6,'(a,f8.1)')  '   ilm-3 loop              (cpu):',time_rad
     write(6,'(a,2f8.1)') '   write-clm          (cpu,wall):',time_writeclm,time_writeclm_w
     write(6,'(a,2f8.1)') '   write-scf          (cpu,wall):',time_writescf,time_writescf_w
     write(6,'(a,2f8.1)') '   write-efg          (cpu,wakk):',time_writeefg,time_writeefg_w
     write(6,'(a,2f8.1)') '   write              (cpu,wall):',time_write,time_write_w
     write(6,'(a,2f8.1)') '   dmft transformation(cpu,wall):',time_dmf, time_dmfw
     write(6,'(a,2f8.1)') '   dmft exact-diag    (cpu,wall):',time_ed, time_edw
     write(6,'(a,2f8.1)') '   interstitial charge(cpu,wall):',time_int, time_intw
     
     write(6,881)  xwt
     write(21,881) xwt
  endif

  rewind fh_vec


  deallocate (e_store,elo_store)
  deallocate(alm,blm,clm)
  deallocate(dh_alm,dh_blm,dh_clm)
  call fini_charp
  call fini_chard
  call fini_charf
  call fini_lohelp
  call fini_xa
  call fini_xa3
  call w_deallocate

  call cputim(time2)
  tclm = time2
  call walltim(time2_w)
  tclm_w = time2
  
  if (myrank.eq.master) then
     
     ! ---- interstitials ---
     !
     ! back transform cumulative grids
     isig=1
     ! computes rho(k)=fft(rho(r))
     call c3fft(iff1,iff2,iff3,sumfft,iff1,iff2,isig,cwork,dwork,ierr) 
     call c3fft(iff1,iff2,iff3,vsumfft,iff1,iff2,isig,cwork,dwork,ierr)
     sumfft(:,:,:)=sumfft(:,:,:)/dble(iff1*iff2*iff3)
     vsumfft(:,:,:)=vsumfft(:,:,:)/dble(iff1*iff2*iff3)
     
     allocate( rho1(nwave*nsym),vrho1(nwave*nsym) )
     rho1(:)=0.0 ; vrho1(:)=0.0
     
     call getfft(nwave,iff1,iff2,iff3,rho1,sumfft,kmax)
     call getfft(nwave,iff1,iff2,iff3,vrho1,vsumfft,kmax)
     
     allocate( rhok(nwave),vrhok(nwave) )
     vrhok(:)=0.0;  rhok(:)=0.0
     
     !.....sum over all recpr. lattic vectors of one star
     ia1=1
     volin = 1.0d0/vol
     do j=1,nwave
        ia2=ia1+inst(j)-1
        do jj=ia1,ia2
           rhok(j)  = rhok(j)  + rho1(jj)*conjg(tauk(jj))
           vrhok(j) = vrhok(j) + vrho1(jj)*conjg(tauk(jj))
        enddo
        ia1=ia2+1
        vrhok(j) = vrhok(j) * volin
        rhok(j)  = rhok(j)  * volin
     enddo
     
     write(8,*) '   valence charge density in interstitial '
     write(8,*)
     write(8,2061)  nwave
     write(8,2071)  ( (kzz(jx,j),jx=1,3),rhok(j), j=1,nwave)
     close(8)

     write(6,204) nwave,kmax
     deallocate(rhok,vrhok)
     deallocate(rho1,vrho1)
  endif

  deallocate (fft, tp, sumfft, vsumfft, tloc )
  deallocate (dwork, cwork)
  ! --- interstitials ---

  call cputim(time2)
  tfour = time2
  call walltim(time2_w)
  tfour_w = time2

  return
  

205 format(/,1x,' k-point:',3f8.4,1x,i5,i4,2x,a10)                      
207 format(1h0,' total charge inside sphere',i5,1h:,f12.6)            
250 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,pz,pxy',/,':QTL',i3.3,':',12f7.4)               
251 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,      ','d-eg,d-t2g ',/,':QTL',i3.3,':',12f7.4)  
252 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,      ','d-z2,d-xy,x2y2,d-xz,yz ',/,':QTL',i3.3,':',12f7.4)
253 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,pz,pxy,,','d-z2,d-xy,d-x2y2,d-xz,yz ',/,':QTL',i3.3,':',12f7.4)                          
254 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,pz,pxy, ','d-z2,d-xy,x2y2,d-xz,yz ',/,':QTL',i3.3,':',12f7.4)                            
255 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,      ','d-z2,d-x2y2,d-xy,d-xz,d-yz ',/,':QTL',i3.3,':',12f7.4)                        
256 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,px,py,pz,',' ',/,':QTL',i3.3,':',12f7.4)
257 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,px,py,pz,','d-z2,d-x2y2,d-xy,d-xz,d-yz ',/,':QTL',i3.3,':',12f7.4)                        
258 format(':PCS',i3.3,':',1x,'partial charges sphere =',i3,' s,p,d,f,px,py,pz,d-z2,d-x2y2,d-xy,d-xz,d-yz,f00,f11,f22',',f33,','f1m,f2m,f3m ',/,':QTL',i3.3,':',19f7.4)
720 format(/,':CHA',i3.3,':',1x,'total charge inside sphere ',i3,' = ',f12.6)         
787 format('    valence charge density   in  mt spheres',5x,i3,' iteration')
800 format(/,10x,68(1h-),/,11x,'w a v e f u n c t i o n s','   and   c h a r g e s   in   s p h e r e s',/,10x,68(1h-))                                         
881 format(/,':CHA  :',' total charge inside unit cell =',f15.6) 
1003 format(251(i3,i2))                                                 
1005 format(/,7x,'lmmax',i3,/,7x,'lm= ',17(i3,i2),(/,7x,18(i3,i2)))          
1990 format(3x,'atomnumber   ',i3,5x,10a4)                             
2001 format(3x,'number of lm',i3//)                                   
2011 format(3x,'clm(r) for l',i3,3x,'m=',i2/)                         
2022 format(3x,4e19.12)                                                 
2030 format(///)                                                       
2031 format(/)                                                         
2032 format(49x,i3,//)  
2055 format(/,1x,' k-point:',3f14.10,1x,i5,i4,2x,a10)
3017 format(/,22x,'vxx',9x,'vyy',9x,'vzz',7x,'up to r',/)
3018 format(/'  efg components before transformation [10**21 v/(m*m)]:',/,15x,'qxx (-v20/2. + v22)     =',f10.3,/,15x,'qyy (-v20/2. - v22)     =',f10.3,/,15x,'qzz ( v20)              =',f10.3,/,15x,'qxy ( v2m/2. * sqrt(3.))=',f10.3)  
3019 format(/,22x,'qxx',9x,'qxy',9x,'qyy',9x,'qzz',7x,'up to r',/)
3020 format(/'  efg components after transformation [10**21 v/(m*m)]:',/,15x,'vxx                     =',f10.3,/,15x,'vyy                     =',f10.3,/,15x,'vzz                     =',f10.3,/)
3121 format(':VZZ',i3.3,':',8x,4f12.5,f12.3)
3021 format(':VZZ',i3.3,':',8x,3f12.5,f12.3)
3022 format(  15x,'qxx=',f10.3,/,15x,'qyy=',f10.3,/,15x,'qzz=',f10.3,15x,'qxy=',f10.3,' up to r= ',f10.5)
4645 format(2i4,3f15.10)
4646 format(8e20.10)
4893 format(3i4,4(2e16.8,2x))  
13   format(////,':POS',i3.3,':',1x,'at.nr.',i4,1x,'position =',3f8.5,2x,'multiplicity =',i3)
7207 format(1h ,' efg inside sphere',i5,1h:,f12.6,5x,' up to r = ',f10.7)                                                          
7720 format(':VZZ',i3.3,':',1x,'efg inside sphere ',i3,' = ',f12.6,5x,' up to r =',f10.5)


204 format(i10,' fourier coefficients calculated',/' largest components:',3i6)
2061 format(1x,'        ',i10,1x,'number of pw')
2071 format(3x,3i5,2e19.12)

end subroutine l2main
