!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine h2q3g(p1,p2,p3,p4,p5,Hqaggg)
      implicit none
      include 'types.f'

c-----calculates the matrix element squared
c-----for q(p1)+q~(p2)-->g(p3)+g(p4)+g(p5)+H
c----- using the results of S. Badger and
c----   V.~Del Duca, A.~Frizzo and F.~Maltoni,
c----   %``Higgs boson production in association with three jets,''
c----   JHEP {\bf 0405}, 064 (2004)
c----   [arXiv:hep-ph/0404013].
c-----has to be preceded by a call to spinoru to set up the za,zb
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      real(dp)::Hqaggg,sign,factor
      complex(dp)::qamp(6,2,2,2,2),temp(6,2,2,2,2)
      complex(dp)::
     & na2q3g_mmmpp,na2q3g_mmpmp,na2q3g_mpmmp,na2q3g_mpppp
c     & A2q3g_mpppp,
c     & A2q3g_mmmpp,
c     & A2q3g_mmpmp,
c     & A2q3g_mpmmp,
c      complex(dp)::a0hqbqggg_mp_ppp,
c     & a0hqbqggg_mpppm,a0hqbqggg_mp_pmp


      real(dp)::xa,xb,xc,xd
c--- Full result
      parameter(xa=xn*cf**2,xb=-cf/two,xc=0.25_dp/xn,xd=(xn**2+one)*xc)

c--- Leading in colour
c    parameter(xa=xn**3/four,xb=zip,xc=zip,xd=zip)
c--- 1/N^2 suppressed
c    parameter(xa=-xn/two,xb=-xn/four,xc=zip,xd=xn/four)
c--- 1/N^4 suppressed
c   parameter(xa=one/four/xn,xb=one/four/xn,xc=one/four/xn,xd=one/four/xn)
c--- Note that the definition of matrix DJK here differs from the one
c--- in (B.22) of the paper by an overall factor of Cf which is
c--- restored in the sum below
      integer::j,k,h1,h3,h4,h5,p1,p2,p3,p4,p5,n(5)
      integer, parameter, dimension(6) ::i3=(/3,3,4,4,5,5/)
      integer, parameter, dimension(6) ::i4=(/4,5,3,5,3,4/)
      integer, parameter, dimension(6) ::i5=(/5,4,5,3,4,3/)

      real(dp), parameter ::DJK(6,6)=reshape(
     & (/xa,xb,xb,xc,xc,xd,
     &   xb,xa,xc,xd,xb,xc,
     &   xb,xc,xa,xb,xd,xc,
     &   xc,xd,xb,xa,xc,xb,
     &   xc,xb,xd,xc,xa,xb,
     &   xd,xc,xc,xb,xb,xa/),(/6,6/))

      n(1)=p1
      n(2)=p2
      n(3)=p3
      n(4)=p4
      n(5)=p5

c--- To be used when relating amplitudes by symmetry:
c--- additional sign for the crossed amplitudes with initial gluon,
c--- but no sign for two gluons
      if ((p4 == 4) .and. (p3 < 3)) then
        sign=-one
      else
        sign=+one
      endif

c----definition of helicities is for outgoing lines
c----labelling is as follows
c     temp(j,h1,h3,h4,h5) since h2 can be obtained from h1

      do j=1,6

c      otmp(j,2,2,2,2)=A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),
c     & za,zb)
      temp(j,2,2,2,2)=-nA2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j))
     & ,za,zb)

c--------
c      write(6,*) '2q3g:o:+- +++',
c     & A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),za,zb)
c      write(6,*) '2q3g:n:+- +++',
c     & -nA2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),za,zb)
c--------

c      otmp(j,2,1,1,1)=A2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j)),
c     & zb,za)
      temp(j,2,1,1,1)=-nA2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j))
     & ,zb,za)

c      otmp(j,2,2,2,1) = A2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     & n(2),zb,za)
      temp(j,2,2,2,1) =-nA2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

c----------
c      write(6,*) '2q3g:o:+- ++-',
c     & A2q3g_mmmpp(n(2),n(i3(j)),n(i4(j)),n(i5(j)),n(1),zb,za)
c      write(6,*) '2q3g:n:+- ++-',
c     & -nA2q3g_mmmpp(n(2),n(i3(j)),n(i4(j)),n(i5(j)),n(1),zb,za)
c------------

c      otmp(j,2,2,1,1) = A2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     & n(1),za,zb)
      temp(j,2,2,1,1) =-nA2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

c      otmp(j,2,2,1,2) = A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j))
c     & ,n(2),zb,za)
      temp(j,2,2,1,2) = -nA2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

c--------
c      write(6,*) '2q3g:o:+- +-+',
c     & A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
c      write(6,*) '2q3g:n:+- +-+',
c     & -nA2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
c--------

c      otmp(j,2,1,2,1) = A2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     & n(1),za,zb)
      temp(j,2,1,2,1) = -nA2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

c      otmp(j,2,1,2,2) = A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     & n(2),zb,za)
      temp(j,2,1,2,2) = -nA2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     & n(2),zb,za)

c--------
c      write(6,*) '2q3g:o:+- -++',
c     & A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
c      write(6,*) '2q3g:n:+- -++',
c     & -nA2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),n(2),zb,za)
c---------

c      otmp(j,2,1,1,2) = A2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     &n(1),za,zb)
      temp(j,2,1,1,2) = -nA2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     & n(1),za,zb)

c       pause

c      temp(j,1,1,1,1)=A2q3g_mpppp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),
c     & zb,za)
c      temp(j,1,2,2,2)=A2q3g_mpppp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j)),
c     & za,zb)

c      temp(j,1,2,1,1) = A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     & n(2),za,zb)
c      temp(j,1,2,2,1) = A2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     & n(1),zb,za)

c      temp(j,1,1,1,2) = A2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     & n(2),za,zb)
c      temp(j,1,1,2,2) = A2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     & n(1),zb,za)

c      temp(j,1,1,2,1)= A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
c     & n(2),za,zb)
c      temp(j,1,2,1,2)= A2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
c     & n(1),zb,za)

c--- fastest to obtain remaining amplitudes by symmetry
c      otmp(j,1,1,1,1)=-sign*conjg(otmp(j,2,2,2,2))
c      otmp(j,1,2,2,2)=-sign*conjg(otmp(j,2,1,1,1))
c      otmp(j,1,2,1,1)=-sign*conjg(otmp(j,2,1,2,2))
c      otmp(j,1,2,2,1)=-sign*conjg(otmp(j,2,1,1,2))
c      otmp(j,1,1,1,2)=-sign*conjg(otmp(j,2,2,2,1))
c      otmp(j,1,1,2,2)=-sign*conjg(otmp(j,2,2,1,1))
c      otmp(j,1,1,2,1)=-sign*conjg(otmp(j,2,2,1,2))
c      otmp(j,1,2,1,2)=-sign*conjg(otmp(j,2,1,2,1))

      temp(j,1,1,1,1)=-sign*conjg(temp(j,2,2,2,2))
      temp(j,1,2,2,2)=-sign*conjg(temp(j,2,1,1,1))
      temp(j,1,2,1,1)=-sign*conjg(temp(j,2,1,2,2))
      temp(j,1,2,2,1)=-sign*conjg(temp(j,2,1,1,2))
      temp(j,1,1,1,2)=-sign*conjg(temp(j,2,2,2,1))
      temp(j,1,1,2,2)=-sign*conjg(temp(j,2,2,1,1))
      temp(j,1,1,2,1)=-sign*conjg(temp(j,2,2,1,2))
      temp(j,1,2,1,2)=-sign*conjg(temp(j,2,1,2,1))

      enddo

c      do j=1,6
c      do h1=1,2
c      do h3=1,2
c      do h4=1,2
c      do h5=1,2
c      write(6,*) j,h1,h3,h4,h5,temp(j,h1,h3,h4,h5)/
c     & otmp(j,h1,h3,h4,h5)
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      pause


c----At this stage we have setup the amplitudes but failed
c----to assign the helicities properly. So we now reshuffle
c----to get these right.
      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      do j=1,6
      if (j== 1) qamp(j,h1,h3,h4,h5)=temp(j,h1,h3,h4,h5)
      if (j== 2) qamp(j,h1,h3,h4,h5)=temp(j,h1,h3,h5,h4)
      if (j== 3) qamp(j,h1,h3,h4,h5)=temp(j,h1,h4,h3,h5)
      if (j== 4) qamp(j,h1,h3,h4,h5)=temp(j,h1,h4,h5,h3)
      if (j== 5) qamp(j,h1,h3,h4,h5)=temp(j,h1,h5,h3,h4)
      if (j== 6) qamp(j,h1,h3,h4,h5)=temp(j,h1,h5,h4,h3)
      enddo
      enddo
      enddo
      enddo
      enddo

c--- now perform the sum with the appropriate weights,
c--- c.f. Eq. (B.20); note that the factor of (gsq)**3
c--- is included in the wrapping routine

c--- NB: use symmetry to slightly improve speed, summing over
c---     diagonal and above in matrix
      Hqaggg=zip
      do j=1,6
      do k=j,6
      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      if (j == k) then
      factor=one
      else
      factor=two
      endif
      Hqaggg=Hqaggg+two*Cf*djk(j,k)
     & *real(qamp(j,h1,h3,h4,h5)*conjg(qamp(k,h1,h3,h4,h5)))*factor
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
