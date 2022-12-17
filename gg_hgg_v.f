!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_hgg_v(p,msq)
      implicit none
      include 'types.f'
c--- Virtual matrix element squared averaged over initial colors and spins

c     g(-p1)+g(-p2) -->  H(p3)+g(p_iglue1=5)+g(p_iglue2=6)

c    Calculation is fully analytic

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      include 'scheme.f'
      include 'nflav.f'
      include 'deltar.f'
      include 'hdecaymode.f'
      include 'blha.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),s34
      real(dp):: hdecay,Asq,fac,msqhgamgam
      real(dp):: qrqr,qarb,aqbr,abab,qbra,bqar
      real(dp):: qaqa,aqaq,qqqq,aaaa
      real(dp):: qagg,aqgg,qgqg,gqqg,agag,gaag,ggqa
      real(dp):: gggg
      real(dp):: Hqarbvsqanal
      real(dp):: Hqaqavsqanal
      real(dp):: HAQggvsqanal
      real(dp):: Hggggvsqanal
      logical:: qb1,qb2,ab1,ab2,gb1,gb2
      logical:: CheckEGZ
      common/CheckEGZ/CheckEGZ
!$omp threadprivate(/CheckEGZ/)
c***************************************************
      scheme='tH-V'
c***************************************************

      if     (scheme == 'dred') then
        deltar=0._dp
      elseif (scheme == 'tH-V') then
        deltar=1._dp
      else
        write(6,*) 'Invalid scheme in gg_hgg_v.f'
        stop
      endif

c--- Set this to true to check squared matrix elements against
c--- hep-ph/0506196 using the point specified in Eq. (51)
      CheckEGZ=.false.

c--- Set up spinor products
      call spinoru(6,p,za,zb)

      Asq=(as/(3._dp*pi))**2/vevsq

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

      if (hdecaymode == 'none') then
         hdecay = 1._dp
      else
c   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqhgamgam(s34)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      endif

      fac=ason2pi*Asq*gsq**2*hdecay

c--- for checking EGZ
      if (CheckEGZ) then
        call CheckEGZres
      endif

c--- for checking scheme dependence of amplitudes
c      call CheckScheme(1,2,5,6)

      qb1=.true.
      qb2=.true.
      ab1=.true.
      ab2=.true.
      gb1=.true.
      gb2=.true.

c--- Note that Hqarbvsqanal(1,2,5,6)=Hqarbvsqanal(6,5,2,1)
c--- and the basic process is q(-k6)+r(-k2)-->q(-k5)+r(-k1)

c--- FOUR-QUARK PROCESSES WITH NON-IDENTICAL QUARKS

c---quark-quark
c     q(1)+r(2)->q(5)+r(6)
      if (useblha==0.and.(qb1 .and. qb2)) then
        qrqr=Hqarbvsqanal(6,2,5,1)
      else
        qrqr=zip
      endif

c----quark-antiquark annihilation (6-->5-->2-->6) wrt q(1)+r(2)->q(5)+r(6)
c     q(1)+a(2)->r(5)+b(6)
      if (useblha==0.and.((qb1 .and. ab2) .or. (ab1 .and. qb2))) then
        qarb=Hqarbvsqanal(5,6,2,1)
      else
        qarb=zip
      endif

c----antiquark-quark annihilation (1<-->2, 5<-->6) wrt to the above
c     a(1)+q(2)->b(5)+r(6)
c      aqbr=Hqarbvsqanal(6,5,1,2)
      aqbr=qarb

      if (useblha==0.and.((qb1 .and. ab2) .or. (ab1 .and. qb2))) then
c----quark-antiquark scattering (6<-->2) wrt q(1)+r(2)->q(5)+r(6)
c     q(1)+b(2)->r(5)+a(6)
        qbra=Hqarbvsqanal(2,6,5,1)
      else
        qbra=zip
      endif

c----antiquark-quark scattering
c     b(1)+q(2)->a(5)+r(6) (1<-->2, 5<-->6) wrt to the above
c      bqar=Hqarbvsqanal(1,5,6,2)
      bqar=qbra

      if ((useblha==0.and.(ab1 .and. ab2))
     &  .or.(blhatype==8).or.(blhatype==10)) then
c---antiquark-antiquark scattering (1<-->5,2<-->6) wrt q(1)+r(2)->q(5)+r(6)
c     a(1)+b(2)->a(5)+b(6)
        abab=Hqarbvsqanal(2,6,1,5)
      else
        abab=zip
      endif

c--- FOUR-QUARK PROCESSES WITH IDENTICAL QUARKS

      if (useblha==0.and.(qb1 .and. qb2)) then
c     q(1)+q(2)->q(5)+q(6)
        qqqq=qrqr+Hqarbvsqanal(5,2,6,1)+Hqaqavsqanal(6,2,5,1)
      else
        qqqq=zip
      endif

      if ((useblha==0.and.(ab1 .and. ab2))
     &  .or.(blhatype==8)) then
c     a(1)+a(2)->a(5)+a(6) (1<-->5,2<-->6) wrt q(1)+q(2)->q(5)+q(6)
        aaaa=abab+Hqarbvsqanal(2,5,1,6)+Hqaqavsqanal(2,6,1,5)
      else
        aaaa=zip
      endif

      if (useblha==0.and.((qb1 .and. ab2) .or. (ab1 .and. qb2))) then
c     q(1)+a(2)->q(5)+a(6) (2<-->6) wrt q(1)+q(2)->q(5)+q(6)
        qaqa=qbra+qarb+Hqaqavsqanal(2,6,5,1)
      else
        qaqa=zip
      endif

c     a(1)+q(2)->a(5)+q(6) (1<-->2, 5<-->6) wrt the above
c      aqqa=qbra+qarb+Hqaqavsqanal(1,5,6,2)
      aqaq=qaqa

c--- TWO-QUARK, TWO GLUON PROCESSES

      if (useblha==0.and.((qb1 .and. ab2) .or. (ab1 .and. qb2))) then
c     a(1)+q(2)->g(3)+g(4)
        aqgg=+HAQggvsqanal(2,1,5,6)
      else
        aqgg=zip
      endif

      if (useblha==0.and.((qb1 .and. gb2) .or. (ab1 .and. gb2))) then
c     q(1)+g(2)->q(5)+g(6)
        qgqg=+HAQggvsqanal(1,5,2,6)
      else
        qgqg=zip
      endif

      if (useblha==0.and.((gb1 .and. qb2) .or. (gb1 .and. ab2))) then
c     g(1)+q(2)->q(5)+g(6)
        gqqg=+HAQggvsqanal(2,5,1,6)
      else
        gqqg=zip
      endif

c     q(1)+a(2)->g(5)+g(6)
c      qagg=+HAQggvsqanal(1,2,5,6)
      qagg=aqgg

c     a(1)+g(2)->a(5)+g(6)
c      agag=+HAQggvsqanal(5,1,2,6)
      agag=qgqg

c     g(1)+a(2)->a(5)+g(6)
c      gaag=+HAQggvsqanal(5,2,1,6)
      gaag=gqqg

      if ((useblha==0.and.(gb1 .and. gb2))
     &  .or.(blhatype==2)) then
c     g(1)+g(2)->q(5)+a(6)
        ggqa=+HAQggvsqanal(6,5,1,2)
      else
        ggqa=zip
      endif

      if ((useblha==0.and.(gb1 .and. gb2))
     &  .or.(blhatype==1)) then
c--- FOUR GLUON PROCESS
        gggg=+Hggggvsqanal(1,2,5,6)
      else
        gggg=zip
      endif


c--- DEBUGGING OUTPUT
c      write(6,*) 'qrqr',qrqr
c      write(6,*) 'qarb',qarb
c      write(6,*) 'aqrb',aqrb
c      write(6,*) 'abab',abab
c      write(6,*) 'qbra',qbra
c      write(6,*) 'bqra',bqra

c      write(6,*) 'Identical'
c      write(6,*) 'qaqa',qaqa
c      write(6,*) 'aqqa',aqqa
c      write(6,*) 'qqqq',qqqq
c      write(6,*) 'aaaa',aaaa


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp

      if ((j==0).and.(k==0)) then
c---gg - all poles cancelled
         msq(j,k)=fac*avegg*(half*gggg+real(nflav,dp)*ggqa)
         if (useblha==1) then
            if ((blhafl(5)==0).and.(blhafl(6)==0)) then
               msq(j,k)=fac*avegg*half*gggg
            else
               msq(j,k)=fac*avegg*ggqa
            endif
         endif
      elseif ((j>0).and.(k>0)) then
c---qq - all poles cancelled
         if (j==k) then
         msq(j,k)=aveqq*fac*half*qqqq
         else
         msq(j,k)=aveqq*fac*qrqr
         endif

      elseif ((j<0).and.(k<0)) then
c---aa - all poles cancelled
         if (j==k) then
         msq(j,k)=aveqq*fac*half*aaaa
         else
         msq(j,k)=aveqq*fac*abab
         endif

      elseif ((j>0).and.(k<0)) then
c----qa scattering - all poles cancelled
         if (j==-k) then
         msq(j,k)=aveqq*fac*(real(nflav-1,dp)*qarb+qaqa+half*qagg)
         if (useblha==1) then
            if ((blhafl(5)==0).and.(blhafl(6)==0)) then
               msq(j,k)=aveqq*fac*half*qagg
            elseif ((blhafl(5)==blhafl(1)).and.
     &              (blhafl(6)==blhafl(2))) then
               msq(j,k)=aveqq*fac*qaqa
            else
               msq(j,k)=aveqq*fac*qarb
            endif
         endif
             else
         msq(j,k)=aveqq*fac*qbra
         endif

      elseif ((j<0).and.(k>0)) then
c----aq scattering - all poles cancelled
         if (j==-k) then
         msq(j,k)=aveqq*fac*(real(nflav-1,dp)*aqbr+aqaq+half*aqgg)
         if (useblha==1) then
            if ((blhafl(5)==0).and.(blhafl(6)==0)) then
               msq(j,k)=aveqq*fac*half*aqgg
            elseif ((blhafl(5)==blhafl(1)).and.
     &              (blhafl(6)==blhafl(2))) then
               msq(j,k)=aveqq*fac*aqaq
            else
               msq(j,k)=aveqq*fac*aqbr
            endif
         endif
             else
         msq(j,k)=aveqq*fac*bqar
         endif

      elseif ((j==0).and.(k>0)) then
c----gq scattering - all poles cancelled
         msq(j,k)=aveqg*fac*gqqg

      elseif ((j==0).and.(k<0)) then
c----ga scattering - all poles cancelled
         msq(j,k)=aveqg*fac*gaag

      elseif ((j>0).and.(k==0)) then
c----qg scattering - all poles cancelled
         msq(j,k)=aveqg*fac*qgqg

      elseif ((j<0).and.(k==0)) then
c----ag scattering - all poles cancelled
         msq(j,k)=aveqg*fac*agag
      endif

      enddo
      enddo

      return
      end




