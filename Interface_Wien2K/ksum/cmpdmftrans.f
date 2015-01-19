SUBROUTINE CmpDMFTRANS(DMFTU, alml, ri_mat, projector, cfX, Rspin, L, iorb1, cix_orb, nindo, iso, nemin, nbands, nind, maxdim2, maxdim, icix, idt, lmax2, nume, nrf, natom, ncix, nrfmax, norbitals)
  ! The DMFT transformation is created, which transforms between Kohn-Sham basis
  ! and DMFT basis.
  ! It is created by :
  ! DMFTransform(i,j,L,M) =
  ! <Y_{L}(r_1)|psi_{ki}(r_1)><psi_{kj}(r_2)|Y_{M}(r_2)>*delta(r_1-r_2)
  ! Where <|> stands for the integral over space angles of r_1 and r_2 only,
  ! while delta(r_1-r_2) is delta function in radial distance only.
  !
  ! The transformation has the following important properties:
  !
  ! \sum_i DMFTrans(i,i,L1,L2) = delta(L1,L2)
  !
  ! \sum_L DMFTrans(i,j,L,L) = delta(i,j) . Be careful! For this to be walid,
  ! one needs to sum over all L to high cutoff.
  !                                         Something we never really need to
  !                                         do.
  !
  ! sum_{ij} DMFTrans(i,j,L1,L2) * DMFTrans(j,i,L2p,L1p) =
  ! delta(L1,L1p)*delta(L2,L2p)
  ! which makes shure that a quantity G_{LL'} when transformed to G_{ij} and
  ! back to G_{LL'} takes the same form
  !
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: DMFTU(nbands,maxdim2,norbitals)
  COMPLEX*16, intent(in)  :: alml(0:lmax2,2*lmax2+1,nume,nrf,2)
  REAL*8,     intent(in)  :: ri_mat(0:lmax2,0:lmax2,nrf,nrf,2,2)
  COMPLEX*16, intent(in)  :: cfX(maxdim2,maxdim2,norbitals,norbitals), Rspin(2,2,norbitals)
  INTEGER,    intent(in)  :: L, iorb1, cix_orb(norbitals), nindo(norbitals), iso, nemin, nbands, nind, icix, idt
  INTEGER,    intent(in)  :: lmax2, nume, nrf, natom, maxdim, maxdim2, ncix, nrfmax, norbitals
  INTEGER,    intent(in)  :: projector
  ! locals
  INTEGER :: N1, num1, num2, is1, m1, lms1, ind1, is2, m2, lms2, ind2, irf1, irf2, it, nind1, nind2, iorb2
  COMPLEX*16 :: cvalue, csum, cc, corc, sm1, sm2, ff
  REAL*8     :: fct
  COMPLEX*16, allocatable :: Rx(:,:), cft(:,:), URx(:,:), tmp(:,:)
  REAL*8 :: small, cor2
  COMPLEX*16 :: tuu, tud, tdu, tdd
  COMPLEX*16, allocatable :: Uu(:,:), Ud(:,:)
  INTEGER :: no1

  small = 1e-6   ! cutoff for ri_mat
  N1=2*L+1

  ALLOCATE( cft(maxdim2,maxdim2) )
  ALLOCATE( tmp(nbands,maxdim2) )

  if (icix.gt.0) then
     allocate( URx(nbands,nind) )


     if (projector.eq.1) then
        
        ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U
        do num1=1,nbands      !nemin,nemax     ! over bands
           do is1=1,iso            ! over spin-1
              do m1=-L,L           ! over m-1
                 lms1=L+1+m1
                 ind1=L+1+idt*m1+N1*(is1-1)
                 cc=0
                 do irf1=1,nrfmax
                    cc = cc + alml(L,lms1,num1+nemin-1,irf1,is1)*ri_mat(L,L,irf1,1,is1,is1)*(idt)**m1
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
        
        ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U
        do num1=1,nbands      !nemin,nemax     ! over bands
           do is1=1,iso            ! over spin-1
              do m1=-L,L           ! over m-1
                 lms1=L+1+m1
                 ind1=L+1+idt*m1+N1*(is1-1)
                 cc=0
                 sm1=0
                 sm2=0
                 do irf1=1,nrfmax
                    do irf2=1,nrfmax
                       ff = alml(L,lms1,num1+nemin-1,irf1,is1)*conjg(alml(L,lms1,num1+nemin-1,irf2,is1))
                       sm1 = sm1 + ff * ri_mat(L,L,irf1,irf2,is1,is1)
                       sm2 = sm2 + ff * ri_mat(L,L,irf1,1,is1,is1) * ri_mat(L,L,1,irf2,is1,is1)
                    enddo
                    cc = cc + alml(L,lms1,num1+nemin-1,irf1,is1)*ri_mat(L,L,irf1,1,is1,is1)*(idt)**m1
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
              do m1=-L,L           ! over m-1
                 lms1=L+1+m1
                 ind1=L+1+idt*m1+N1*(is1-1)

                 irf1 = 1
                 cc = alml(L,lms1,num1+nemin-1,irf1,is1)*sqrt(ri_mat(L,L,irf1,irf1,is1,is1))*(idt)**m1


                 ! Adding correction due to \cdot{u} and local orbitals
                 ! The results is  A(i,L,irf=1)*sqrt(o(L)) * sqrt( 1 + \sum_{irf>1} |A(i,L,irf)|^2*o(L,irf))
                 if ( abs(cc)>1e-6 ) then
                    corc = 0
                    do irf1=1,nrfmax
                       do irf2=1,nrfmax
                          if (irf1.eq.1 .and. irf2.eq.1) cycle
                          corc = corc + alml(L,lms1,num1+nemin-1,irf1,is1)*conjg(alml(L,lms1,num1+nemin-1,irf2,is1))*ri_mat(L,L,irf1,irf2,is1,is1)
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

     nind1 = (2*L+1)*iso
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
     deallocate(URx)
  endif

  deallocate( cft )
  deallocate( tmp )
end SUBROUTINE CmpDMFTrans
