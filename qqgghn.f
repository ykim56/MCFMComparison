!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqgghn(p1,p2,p3,p4,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
      implicit none
      include 'types.f'

c---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> H((p5+p6)+g(p3)+g(p4)
c   with momentum 4 contracted with the vector n(mu),
c   separated into colour orderings AB, BA and the sum
c     calculated by the program qqgghn.frm
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4
      real(dp):: p(mxpart,4),n(4),nDn,nDp1,nDp2,nDp3
      real(dp):: is,s123,s124,qqgghn_ab,qqgghn_ba,qqgghn_sym
      is(p1,p2)=one/s(p1,p2)

      nDp1=n(4)*p(p1,4)-n(3)*p(p1,3)-n(2)*p(p1,2)-n(1)*p(p1,1)
      nDp2=n(4)*p(p2,4)-n(3)*p(p2,3)-n(2)*p(p2,2)-n(1)*p(p2,1)
      nDp3=n(4)*p(p3,4)-n(3)*p(p3,3)-n(2)*p(p3,2)-n(1)*p(p3,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2
      s123=s(p1,p2)+s(p2,p3)+s(p3,p1)
      s124=s(p1,p2)+s(p2,p4)+s(p4,p1)

      call checkndotp(p,n,p4)
c--- appropriate scale is approx 1.e-3_dp*energy(incoming)
c--- so of order(1) for the Tevatron
c      if (abs(nDp4)>1.e-3_dp*sqrt(abs(p(1,4)*n(4)))) then
c        write(*,*) 'Error in qqgghn for :',p1,p2,p3,p4
c        write(*,*) 'cutoff',1.e-3_dp*sqrt(abs(p(1,4)*n(4)))
c        write(6,*) 'nDp4',nDp4
c        call flush(6)
c        stop
c      endif

      qqgghn_ab=  + s123**(-2)*xn * ( 6.D0*s(p1,p4)**2*nDn + 12.D0*s(p1
     &    ,p4)*s(p2,p4)*nDn + 12.D0*s(p1,p4)*s(p3,p4)*nDn + 6.D0*s(p2,
     &    p4)**2*nDn + 12.D0*s(p2,p4)*s(p3,p4)*nDn + 6.D0*s(p3,p4)**2*
     &    nDn )
      qqgghn_ab = qqgghn_ab + s123**(-1)*s124**(-1)*xn * (  - 16.D0*s(
     &    p1,p2)*nDp1*nDp3 - 8.D0*s(p1,p2)*s(p1,p4)*nDn - 4.D0*s(p1,p2)
     &    *s(p3,p4)*nDn - 16.D0*s(p1,p3)*nDp2*nDp3 - 32.D0*s(p1,p3)*
     &    nDp1*nDp3 - 8.D0*s(p1,p3)*s(p1,p4)*nDn + 32.D0*s(p1,p4)*
     &    nDp3**2 - 16.D0*s(p1,p4)*nDp2*nDp3 + 16.D0*s(p1,p4)*s(p3,p4)*
     &    nDn - 32.D0*s(p3,p4)*nDp1*nDp3 - 16.D0*s(p3,p4)*nDp1*nDp2 -
     &    48.D0*s(p3,p4)*nDp1**2 + 8.D0*s(p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + s123**(-1)*xn * ( 48.D0*nDp2*nDp3 + 24.D0
     &    *nDp1*nDp3 + 16.D0*nDp1*nDp2 + 2.D0*s(p1,p2)*nDn + 2.D0*s(p1,
     &    p3)*nDn + 36.D0*s(p1,p4)*nDn + 16.D0*s(p2,p4)*nDn + 16.D0*s(
     &    p3,p4)*nDn )
      qqgghn_ab = qqgghn_ab + s124**(-2)*xn * ( 2.D0*s(p1,p3)**2*nDn +
     &    4.D0*s(p1,p3)*s(p2,p3)*nDn + 4.D0*s(p1,p3)*s(p3,p4)*nDn + 2.D0
     &    *s(p2,p3)**2*nDn + 4.D0*s(p2,p3)*s(p3,p4)*nDn + 2.D0*s(p3,p4)
     &    **2*nDn )
      qqgghn_ab = qqgghn_ab + s124**(-1)*xn * (  - 16.D0*nDp2*nDp3 - 16.
     &    D0*nDp2**2 - 16.D0*nDp1*nDp3 - 32.D0*nDp1**2 - 6.D0*s(p1,p2)*
     &    nDn + 16.D0*s(p1,p3)*nDn + 8.D0*s(p1,p4)*nDn + 16.D0*s(p2,p3)
     &    *nDn + 20.D0*s(p3,p4)*nDn )
      qqgghn_ab = qqgghn_ab + xn * ( 6.D0*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*s123**(-2)*xn * ( 8.D0*s(p1,p3)
     &    *s(p1,p4)**2*nDn + 16.D0*s(p1,p3)*s(p1,p4)*s(p2,p4)*nDn + 16.D
     &    0*s(p1,p3)*s(p1,p4)*s(p3,p4)*nDn + 8.D0*s(p1,p3)*s(p2,p4)**2*
     &    nDn + 16.D0*s(p1,p3)*s(p2,p4)*s(p3,p4)*nDn + 8.D0*s(p1,p3)*s(
     &    p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*s123**(-1)*s124**(-1)*xn * (
     &     - 16.D0*s(p1,p3)**2*nDp2*nDp3 - 16.D0*s(p1,p3)**2*nDp1*nDp3
     &     + 16.D0*s(p1,p3)*s(p1,p4)*nDp3**2 - 16.D0*s(p1,p3)*s(p1,p4)*
     &    nDp2*nDp3 - 16.D0*s(p1,p3)*s(p1,p4)*nDp1*nDp3 + 16.D0*s(p1,p3
     &    )*s(p1,p4)*s(p3,p4)*nDn + 16.D0*s(p1,p3)*s(p3,p4)*nDp2*nDp3
     &     - 16.D0*s(p1,p3)*s(p3,p4)*nDp1*nDp3 - 48.D0*s(p1,p3)*s(p3,p4
     &    )*nDp1*nDp2 - 48.D0*s(p1,p3)*s(p3,p4)*nDp1**2 + 32.D0*s(p1,p4
     &    )**2*nDp3**2 + 16.D0*s(p1,p4)**2*nDp2*nDp3 + 16.D0*s(p1,p4)**
     &    2*nDp1*nDp3 + 32.D0*s(p1,p4)*s(p3,p4)*nDp2*nDp3 - 16.D0*s(p1,
     &    p4)*s(p3,p4)*nDp2**2 - 32.D0*s(p1,p4)*s(p3,p4)*nDp1*nDp3 - 80.
     &    D0*s(p1,p4)*s(p3,p4)*nDp1*nDp2 - 64.D0*s(p1,p4)*s(p3,p4)*
     &    nDp1**2 - 8.D0*s(p1,p4)*s(p3,p4)**2*nDn - 32.D0*s(p3,p4)**2*
     &    nDp1*nDp2 - 4.D0*s(p3,p4)**3*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*s123**(-1)*xn * ( 64.D0*s(p1,p3
     &    )*nDp2*nDp3 + 48.D0*s(p1,p3)*nDp1*nDp3 + 16.D0*s(p1,p3)*nDp1*
     &    nDp2 + 36.D0*s(p1,p3)*s(p1,p4)*nDn + 18.D0*s(p1,p3)*s(p2,p4)*
     &    nDn + 8.D0*s(p1,p3)*s(p3,p4)*nDn - 32.D0*s(p1,p4)*nDp3**2 +
     &    56.D0*s(p1,p4)*nDp2*nDp3 - 16.D0*s(p1,p4)*nDp2**2 - 24.D0*s(
     &    p1,p4)*nDp1*nDp3 - 48.D0*s(p1,p4)*nDp1*nDp2 - 14.D0*s(p1,p4)
     &    **2*nDn - 20.D0*s(p1,p4)*s(p2,p4)*nDn - 26.D0*s(p1,p4)*s(p3,
     &    p4)*nDn + 8.D0*s(p2,p4)*nDp1*nDp3 + 16.D0*s(p2,p4)*nDp1*nDp2
     &     + 48.D0*s(p2,p4)*nDp1**2 - 6.D0*s(p2,p4)**2*nDn - 10.D0*s(p2
     &    ,p4)*s(p3,p4)*nDn + 32.D0*s(p3,p4)*nDp1*nDp3 - 8.D0*s(p3,p4)*
     &    nDp1*nDp2 + 56.D0*s(p3,p4)*nDp1**2 - 10.D0*s(p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*s124**(-1)*xn * ( 32.D0*s(p1,p3
     &    )*nDp2*nDp3 + 16.D0*s(p1,p3)*nDp2**2 - 16.D0*s(p1,p3)*nDp1*
     &    nDp2 - 4.D0*s(p1,p3)**2*nDn + 16.D0*s(p1,p3)*s(p1,p4)*nDn - 8.
     &    D0*s(p1,p3)*s(p2,p3)*nDn - 10.D0*s(p1,p3)*s(p3,p4)*nDn - 16.D0
     &    *s(p1,p4)*nDp3**2 + 32.D0*s(p1,p4)*nDp2*nDp3 + 16.D0*s(p1,p4)
     &    *nDp2**2 - 64.D0*s(p1,p4)*nDp1*nDp3 - 48.D0*s(p1,p4)*nDp1*
     &    nDp2 - 64.D0*s(p1,p4)*nDp1**2 - 8.D0*s(p1,p4)*s(p3,p4)*nDn -
     &    32.D0*s(p2,p3)*nDp2*nDp3 - 16.D0*s(p2,p3)*nDp1*nDp2 + 16.D0*
     &    s(p2,p3)*nDp1**2 - 4.D0*s(p2,p3)**2*nDn - 10.D0*s(p2,p3)*s(p3
     &    ,p4)*nDn - 24.D0*s(p3,p4)*nDp2*nDp3 + 8.D0*s(p3,p4)*nDp1*nDp3
     &     - 32.D0*s(p3,p4)*nDp1*nDp2 + 32.D0*s(p3,p4)*nDp1**2 - 10.D0*
     &    s(p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*xn * ( 8.D0*nDp3**2 - 32.D0*
     &    nDp2*nDp3 + 8.D0*nDp2**2 + 16.D0*nDp1*nDp3 - 8.D0*nDp1*nDp2
     &     + 32.D0*nDp1**2 - 20.D0*s(p1,p3)*nDn - 42.D0*s(p1,p4)*nDn -
     &    16.D0*s(p2,p3)*nDn - 24.D0*s(p2,p4)*nDn - 24.D0*s(p3,p4)*nDn
     &     )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*s123**(-2)*xn * ( 4.D0*s(p1,
     &    p3)**2*s(p1,p4)**2*nDn + 8.D0*s(p1,p3)**2*s(p1,p4)*s(p2,p4)*
     &    nDn + 8.D0*s(p1,p3)**2*s(p1,p4)*s(p3,p4)*nDn + 4.D0*s(p1,p3)
     &    **2*s(p2,p4)**2*nDn + 8.D0*s(p1,p3)**2*s(p2,p4)*s(p3,p4)*nDn
     &     + 4.D0*s(p1,p3)**2*s(p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*s123**(-1)*s124**(-1)*xn
     &  * (  - 16.D0*s(p1,p3)**2*s(p3,p4)*nDp2**2 - 32.D0*s(p1,p3)**2*
     &    s(p3,p4)*nDp1*nDp2 - 16.D0*s(p1,p3)**2*s(p3,p4)*nDp1**2 - 32.D
     &    0*s(p1,p3)*s(p1,p4)*s(p3,p4)*nDp2**2 - 64.D0*s(p1,p3)*s(p1,p4
     &    )*s(p3,p4)*nDp1*nDp2 - 32.D0*s(p1,p3)*s(p1,p4)*s(p3,p4)*
     &    nDp1**2 - 8.D0*s(p1,p3)*s(p1,p4)*s(p3,p4)**2*nDn - 16.D0*s(p1
     &    ,p4)**2*s(p3,p4)*nDp2**2 - 32.D0*s(p1,p4)**2*s(p3,p4)*nDp1*
     &    nDp2 - 16.D0*s(p1,p4)**2*s(p3,p4)*nDp1**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*s123**(-1)*xn * ( 32.D0*s(p1
     &    ,p3)**2*nDp2*nDp3 + 32.D0*s(p1,p3)**2*nDp1*nDp3 + 16.D0*s(p1,
     &    p3)**2*s(p1,p4)*nDn + 16.D0*s(p1,p3)**2*s(p2,p4)*nDn + 8.D0*
     &    s(p1,p3)**2*s(p3,p4)*nDn + 32.D0*s(p1,p3)*s(p1,p4)*nDp2*nDp3
     &     - 32.D0*s(p1,p3)*s(p1,p4)*nDp2**2 - 32.D0*s(p1,p3)*s(p1,p4)*
     &    nDp1*nDp2 - 8.D0*s(p1,p3)*s(p1,p4)**2*nDn - 8.D0*s(p1,p3)*s(
     &    p1,p4)*s(p2,p4)*nDn - 16.D0*s(p1,p3)*s(p1,p4)*s(p3,p4)*nDn -
     &    32.D0*s(p1,p3)*s(p2,p4)*nDp1*nDp3 + 32.D0*s(p1,p3)*s(p2,p4)*
     &    nDp1*nDp2 + 32.D0*s(p1,p3)*s(p2,p4)*nDp1**2 + 32.D0*s(p1,p3)*
     &    s(p3,p4)*nDp1*nDp2 + 32.D0*s(p1,p3)*s(p3,p4)*nDp1**2 - 32.D0*
     &    s(p1,p4)**2*nDp2**2 + 64.D0*s(p1,p4)*s(p2,p4)*nDp1*nDp2 + 32.D
     &    0*s(p1,p4)*s(p3,p4)*nDp1*nDp2 + 16.D0*s(p1,p4)*s(p3,p4)*
     &    nDp1**2 - 32.D0*s(p2,p4)**2*nDp1**2 - 16.D0*s(p2,p4)*s(p3,p4)
     &    *nDp1**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*s124**(-2)*xn * ( 4.D0*s(p1,
     &    p3)**2*s(p1,p4)**2*nDn + 8.D0*s(p1,p3)*s(p1,p4)**2*s(p2,p3)*
     &    nDn + 8.D0*s(p1,p3)*s(p1,p4)**2*s(p3,p4)*nDn + 4.D0*s(p1,p4)
     &    **2*s(p2,p3)**2*nDn + 8.D0*s(p1,p4)**2*s(p2,p3)*s(p3,p4)*nDn
     &     + 4.D0*s(p1,p4)**2*s(p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*s124**(-1)*xn * (  - 8.D0*s(
     &    p1,p3)**2*s(p1,p4)*nDn - 32.D0*s(p1,p3)*s(p1,p4)*nDp2*nDp3 -
     &    32.D0*s(p1,p3)*s(p1,p4)*nDp2**2 - 32.D0*s(p1,p3)*s(p1,p4)*
     &    nDp1*nDp2 + 16.D0*s(p1,p3)*s(p1,p4)**2*nDn - 8.D0*s(p1,p3)*s(
     &    p1,p4)*s(p2,p3)*nDn - 16.D0*s(p1,p3)*s(p1,p4)*s(p3,p4)*nDn +
     &    32.D0*s(p1,p3)*s(p3,p4)*nDp1*nDp2 + 16.D0*s(p1,p3)*s(p3,p4)*
     &    nDp1**2 - 32.D0*s(p1,p4)**2*nDp2*nDp3 - 32.D0*s(p1,p4)**2*
     &    nDp2**2 - 32.D0*s(p1,p4)**2*nDp1*nDp3 - 64.D0*s(p1,p4)**2*
     &    nDp1*nDp2 - 32.D0*s(p1,p4)**2*nDp1**2 + 16.D0*s(p1,p4)**2*s(
     &    p2,p3)*nDn + 8.D0*s(p1,p4)**2*s(p3,p4)*nDn + 32.D0*s(p1,p4)*
     &    s(p2,p3)*nDp1*nDp3 + 32.D0*s(p1,p4)*s(p2,p3)*nDp1*nDp2 + 32.D0
     &    *s(p1,p4)*s(p2,p3)*nDp1**2 + 32.D0*s(p1,p4)*s(p3,p4)*nDp1*
     &    nDp2 + 32.D0*s(p1,p4)*s(p3,p4)*nDp1**2 - 16.D0*s(p2,p3)*s(p3,
     &    p4)*nDp1**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*xn * (  - 32.D0*s(p1,p3)*
     &    nDp1*nDp3 + 32.D0*s(p1,p3)*nDp1*nDp2 + 4.D0*s(p1,p3)**2*nDn
     &     - 32.D0*s(p1,p3)*s(p1,p4)*nDn + 32.D0*s(p1,p4)*nDp1*nDp3 +
     &    96.D0*s(p1,p4)*nDp1*nDp2 + 32.D0*s(p1,p4)*nDp1**2 + 4.D0*s(p1
     &    ,p4)**2*nDn - 32.D0*s(p2,p3)*nDp1**2 - 64.D0*s(p2,p4)*nDp1**2
     &     - 32.D0*s(p3,p4)*nDp1**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*is(p3,p4)*s123**(-1)*xn * (
     &     - 16.D0*s(p1,p3)**2*s(p1,p4)*nDp3**2 + 8.D0*s(p1,p3)**2*s(p1
     &    ,p4)**2*nDn + 16.D0*s(p1,p3)**2*s(p1,p4)*s(p2,p4)*nDn - 16.D0
     &    *s(p1,p3)**2*s(p2,p4)*nDp3**2 + 8.D0*s(p1,p3)**2*s(p2,p4)**2*
     &    nDn + 32.D0*s(p1,p3)*s(p1,p4)**2*nDp2*nDp3 + 32.D0*s(p1,p3)*
     &    s(p1,p4)*s(p2,p4)*nDp2*nDp3 - 32.D0*s(p1,p3)*s(p1,p4)*s(p2,p4
     &    )*nDp1*nDp3 - 32.D0*s(p1,p3)*s(p2,p4)**2*nDp1*nDp3 - 16.D0*s(
     &    p1,p4)**3*nDp2**2 - 16.D0*s(p1,p4)**2*s(p2,p4)*nDp2**2 + 32.D0
     &    *s(p1,p4)**2*s(p2,p4)*nDp1*nDp2 + 32.D0*s(p1,p4)*s(p2,p4)**2*
     &    nDp1*nDp2 - 16.D0*s(p1,p4)*s(p2,p4)**2*nDp1**2 - 16.D0*s(p2,
     &    p4)**3*nDp1**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*is(p3,p4)*s124**(-1)*xn * (
     &    8.D0*s(p1,p3)**2*s(p1,p4)**2*nDn - 16.D0*s(p1,p3)*s(p1,p4)**2
     &    *nDp3**2 - 32.D0*s(p1,p3)*s(p1,p4)**2*nDp2*nDp3 - 16.D0*s(p1,
     &    p3)*s(p1,p4)**2*nDp2**2 - 32.D0*s(p1,p3)*s(p1,p4)**2*nDp1*
     &    nDp3 - 32.D0*s(p1,p3)*s(p1,p4)**2*nDp1*nDp2 - 16.D0*s(p1,p3)*
     &    s(p1,p4)**2*nDp1**2 + 16.D0*s(p1,p3)*s(p1,p4)**2*s(p2,p3)*nDn
     &     - 16.D0*s(p1,p4)**2*s(p2,p3)*nDp3**2 - 32.D0*s(p1,p4)**2*s(
     &    p2,p3)*nDp2*nDp3 - 16.D0*s(p1,p4)**2*s(p2,p3)*nDp2**2 - 32.D0
     &    *s(p1,p4)**2*s(p2,p3)*nDp1*nDp3 - 32.D0*s(p1,p4)**2*s(p2,p3)*
     &    nDp1*nDp2 - 16.D0*s(p1,p4)**2*s(p2,p3)*nDp1**2 + 8.D0*s(p1,p4
     &    )**2*s(p2,p3)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*is(p3,p4)*xn * (  - 8.D0*s(
     &    p1,p3)**2*s(p1,p4)*nDn + 8.D0*s(p1,p3)**2*s(p2,p4)*nDn + 32.D0
     &    *s(p1,p3)*s(p1,p4)*nDp3**2 + 32.D0*s(p1,p3)*s(p1,p4)*nDp2*
     &    nDp3 + 32.D0*s(p1,p3)*s(p1,p4)*nDp1*nDp3 + 32.D0*s(p1,p3)*s(
     &    p1,p4)*nDp1*nDp2 + 16.D0*s(p1,p3)*s(p1,p4)*nDp1**2 - 8.D0*s(
     &    p1,p3)*s(p1,p4)**2*nDn - 16.D0*s(p1,p3)*s(p1,p4)*s(p2,p3)*nDn
     &     - 16.D0*s(p1,p3)*s(p1,p4)*s(p2,p4)*nDn - 32.D0*s(p1,p3)*s(p2
     &    ,p4)*nDp1*nDp3 - 16.D0*s(p1,p3)*s(p2,p4)*nDp1**2 - 32.D0*s(p1
     &    ,p4)**2*nDp2*nDp3 - 32.D0*s(p1,p4)**2*nDp2**2 + 8.D0*s(p1,p4)
     &    **2*s(p2,p3)*nDn + 32.D0*s(p1,p4)*s(p2,p3)*nDp1*nDp3 + 32.D0*
     &    s(p1,p4)*s(p2,p3)*nDp1*nDp2 + 16.D0*s(p1,p4)*s(p2,p3)*nDp1**2
     &     + 32.D0*s(p1,p4)*s(p2,p4)*nDp1*nDp3 + 64.D0*s(p1,p4)*s(p2,p4
     &    )*nDp1*nDp2 - 16.D0*s(p2,p3)*s(p2,p4)*nDp1**2 - 32.D0*s(p2,p4
     &    )**2*nDp1**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)**2*is(p3,p4)**2*xn * ( 4.D0*s(
     &    p1,p3)**2*s(p2,p4)**2*nDn - 8.D0*s(p1,p3)*s(p1,p4)*s(p2,p3)*
     &    s(p2,p4)*nDn + 4.D0*s(p1,p4)**2*s(p2,p3)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p1,p3)*s123**(-1)*s124**(-1)
     & *xn * ( 16.D0*s(p1,p4)**3*nDp3**2 + 16.D0*s(p1,p4)**3*nDp2*nDp3
     &     + 16.D0*s(p1,p4)**3*nDp1*nDp3 + 16.D0*s(p1,p4)**2*s(p3,p4)*
     &    nDp2*nDp3 - 16.D0*s(p1,p4)**2*s(p3,p4)*nDp1*nDp3 - 16.D0*s(p1
     &    ,p4)**2*s(p3,p4)*nDp1*nDp2 - 16.D0*s(p1,p4)**2*s(p3,p4)*
     &    nDp1**2 - 16.D0*s(p1,p4)*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p1,p3)*s123**(-1)*xn * (  -
     &    8.D0*s(p1,p4)**2*nDp3**2 + 8.D0*s(p1,p4)**2*nDp2*nDp3 - 16.D0
     &    *s(p1,p4)**2*nDp2**2 - 16.D0*s(p1,p4)**2*nDp1*nDp3 + 8.D0*s(
     &    p1,p4)*s(p2,p4)*nDp3**2 + 8.D0*s(p1,p4)*s(p2,p4)*nDp2*nDp3 -
     &    8.D0*s(p1,p4)*s(p2,p4)*nDp1*nDp3 + 16.D0*s(p1,p4)*s(p2,p4)*
     &    nDp1*nDp2 - 8.D0*s(p1,p4)*s(p3,p4)*nDp2*nDp3 - 8.D0*s(p1,p4)*
     &    s(p3,p4)*nDp2**2 - 8.D0*s(p1,p4)*s(p3,p4)*nDp1*nDp2 + 16.D0*
     &    s(p1,p4)*s(p3,p4)*nDp1**2 - 8.D0*s(p2,p4)**2*nDp1*nDp3 - 8.D0
     &    *s(p2,p4)*s(p3,p4)*nDp1*nDp3 + 8.D0*s(p2,p4)*s(p3,p4)*nDp1*
     &    nDp2 + 8.D0*s(p2,p4)*s(p3,p4)*nDp1**2 + 8.D0*s(p3,p4)**2*nDp1
     &    *nDp2 + 8.D0*s(p3,p4)**2*nDp1**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p1,p3)*s124**(-1)*xn * ( 2.D0
     &    *s(p1,p4)*s(p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p1,p3)*xn * (  - 2.D0*s(p1,
     &    p4)**2*nDn - 2.D0*s(p1,p4)*s(p2,p4)*nDn + 2.D0*s(p1,p4)*s(p3,
     &    p4)*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p1,p3)*is(p3,p4)*s123**(-1)*
     & xn * (  - 16.D0*s(p1,p4)**3*nDp2**2 + 16.D0*s(p1,p4)**2*s(p2,p4)
     &    *nDp2*nDp3 + 32.D0*s(p1,p4)**2*s(p2,p4)*nDp1*nDp2 - 16.D0*s(
     &    p1,p4)*s(p2,p4)**2*nDp1*nDp3 - 16.D0*s(p1,p4)*s(p2,p4)**2*
     &    nDp1**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p1,p3)*is(p3,p4)*xn * (  - 2.
     &    D0*s(p1,p4)**3*nDn - 4.D0*s(p1,p4)**2*s(p2,p4)*nDn - 2.D0*s(
     &    p1,p4)*s(p2,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p2,p4)*s123**(-1)*s124**(-1)
     & *xn * (  - 16.D0*s(p1,p3)**2*s(p3,p4)*nDp2**2 - 16.D0*s(p1,p3)**
     &    2*s(p3,p4)*nDp1*nDp2 + 16.D0*s(p1,p3)*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p2,p4)*s123**(-1)*xn * ( 16.D
     &    0*s(p1,p3)**2*nDp2*nDp3 - 16.D0*s(p1,p3)*s(p1,p4)*nDp2**2 + 2.
     &    D0*s(p1,p3)*s(p1,p4)**2*nDn - 16.D0*s(p1,p3)*s(p3,p4)*nDp2*
     &    nDp3 - 32.D0*s(p1,p3)*s(p3,p4)*nDp2**2 - 2.D0*s(p1,p3)*s(p3,
     &    p4)**2*nDn + 8.D0*s(p1,p4)**2*nDp2*nDp3 - 8.D0*s(p1,p4)*s(p3,
     &    p4)*nDp2**2 - 8.D0*s(p1,p4)*s(p3,p4)*nDp1*nDp2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p2,p4)*s124**(-1)*xn * ( 8.D0
     &    *s(p1,p3)**2*nDp2**2 + 8.D0*s(p1,p3)**2*nDp1*nDp2 + 32.D0*s(
     &    p1,p3)*s(p3,p4)*nDp1*nDp2 + 8.D0*s(p2,p3)**2*nDp2**2 + 8.D0*
     &    s(p2,p3)**2*nDp1*nDp2 + 16.D0*s(p2,p3)*s(p3,p4)*nDp2**2 + 8.D0
     &    *s(p3,p4)**2*nDp2**2 - 8.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p2,p4)*xn * (  - 32.D0*s(p1,
     &    p3)*nDp2*nDp3 - 16.D0*s(p1,p3)*nDp2**2 + 16.D0*s(p1,p3)*nDp1*
     &    nDp2 + 4.D0*s(p1,p3)*s(p2,p3)*nDn - 16.D0*s(p1,p4)*nDp2*nDp3
     &     - 8.D0*s(p1,p4)*nDp2**2 + 8.D0*s(p1,p4)*nDp1*nDp2 - 2.D0*s(
     &    p1,p4)**2*nDn + 4.D0*s(p1,p4)*s(p2,p3)*nDn + 16.D0*s(p2,p3)*
     &    nDp2**2 + 16.D0*s(p2,p3)*nDp1*nDp2 + 4.D0*s(p2,p3)**2*nDn + 2.
     &    D0*s(p2,p3)*s(p3,p4)*nDn + 8.D0*s(p3,p4)*nDp2*nDp3 + 16.D0*s(
     &    p3,p4)*nDp2**2 + 16.D0*s(p3,p4)*nDp1*nDp2 + 2.D0*s(p3,p4)**2*
     &    nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p2,p4)*is(p3,p4)*xn * (  - 8.
     &    D0*s(p1,p3)**2*nDp2*nDp3 - 16.D0*s(p1,p3)*s(p1,p4)*nDp2*nDp3
     &     + 4.D0*s(p1,p3)*s(p1,p4)*s(p2,p3)*nDn - 8.D0*s(p1,p4)**2*
     &    nDp2*nDp3 + 2.D0*s(p1,p4)**2*s(p2,p3)*nDn - 16.D0*s(p1,p4)*s(
     &    p2,p3)*nDp2*nDp3 - 16.D0*s(p1,p4)*s(p2,p3)*nDp2**2 + 4.D0*s(
     &    p1,p4)*s(p2,p3)**2*nDn - 8.D0*s(p2,p3)**2*nDp2*nDp3 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p3,p4)*s123**(-1)*xn * (  -
     &    16.D0*s(p1,p3)*s(p1,p4)*nDp3**2 + 16.D0*s(p1,p3)*s(p1,p4)*
     &    nDp2*nDp3 + 16.D0*s(p1,p3)*s(p1,p4)**2*nDn + 24.D0*s(p1,p3)*
     &    s(p1,p4)*s(p2,p4)*nDn - 32.D0*s(p1,p3)*s(p2,p4)*nDp3**2 - 16.D
     &    0*s(p1,p3)*s(p2,p4)*nDp1*nDp3 + 8.D0*s(p1,p3)*s(p2,p4)**2*nDn
     &     + 48.D0*s(p1,p4)**2*nDp2*nDp3 - 2.D0*s(p1,p4)**3*nDn - 6.D0*
     &    s(p1,p4)**2*s(p2,p4)*nDn - 32.D0*s(p1,p4)*s(p2,p4)*nDp3**2 +
     &    16.D0*s(p1,p4)*s(p2,p4)*nDp2*nDp3 - 48.D0*s(p1,p4)*s(p2,p4)*
     &    nDp1*nDp3 - 6.D0*s(p1,p4)*s(p2,p4)**2*nDn - 16.D0*s(p2,p4)**2
     &    *nDp1*nDp3 - 2.D0*s(p2,p4)**3*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p3,p4)*s124**(-1)*xn * (  -
     &    8.D0*s(p1,p3)**2*nDp2*nDp3 - 8.D0*s(p1,p3)**2*nDp1*nDp3 + 8.D0
     &    *s(p1,p3)**2*s(p1,p4)*nDn - 16.D0*s(p1,p3)*s(p1,p4)*nDp2*nDp3
     &     - 48.D0*s(p1,p3)*s(p1,p4)*nDp1*nDp3 - 32.D0*s(p1,p3)*s(p1,p4
     &    )*nDp1*nDp2 - 32.D0*s(p1,p3)*s(p1,p4)*nDp1**2 + 8.D0*s(p1,p3)
     &    *s(p1,p4)*s(p2,p3)*nDn + 32.D0*s(p1,p4)**2*nDp3**2 + 16.D0*s(
     &    p1,p4)**2*nDp2*nDp3 + 16.D0*s(p1,p4)**2*nDp1*nDp3 - 16.D0*s(
     &    p1,p4)*s(p2,p3)*nDp3**2 + 16.D0*s(p1,p4)*s(p2,p3)*nDp2*nDp3
     &     + 16.D0*s(p1,p4)*s(p2,p3)*nDp2**2 - 16.D0*s(p1,p4)*s(p2,p3)*
     &    nDp1*nDp3 - 16.D0*s(p1,p4)*s(p2,p3)*nDp1**2 - 8.D0*s(p2,p3)**
     &    2*nDp2*nDp3 - 8.D0*s(p2,p3)**2*nDp1*nDp3 )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p3,p4)*xn * ( 16.D0*s(p1,p3)
     &    *nDp1**2 - 4.D0*s(p1,p3)**2*nDn - 10.D0*s(p1,p3)*s(p1,p4)*nDn
     &     - 4.D0*s(p1,p3)*s(p2,p3)*nDn - 10.D0*s(p1,p3)*s(p2,p4)*nDn
     &     - 16.D0*s(p1,p4)*nDp3**2 - 8.D0*s(p1,p4)*nDp2*nDp3 - 24.D0*
     &    s(p1,p4)*nDp1*nDp3 - 16.D0*s(p1,p4)**2*nDn - 6.D0*s(p1,p4)*s(
     &    p2,p3)*nDn - 24.D0*s(p1,p4)*s(p2,p4)*nDn + 16.D0*s(p2,p3)*
     &    nDp3**2 - 32.D0*s(p2,p3)*nDp2*nDp3 - 8.D0*s(p2,p3)*s(p2,p4)*
     &    nDn + 16.D0*s(p2,p4)*nDp3**2 - 16.D0*s(p2,p4)*nDp2*nDp3 + 8.D0
     &    *s(p2,p4)*nDp1*nDp3 - 8.D0*s(p2,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p2)*is(p3,p4)**2*xn * ( 8.D0*s(p1,
     &    p3)**2*nDp3**2 + 16.D0*s(p1,p3)*s(p1,p4)*nDp3**2 + 8.D0*s(p1,
     &    p4)**2*nDp3**2 + 8.D0*s(p2,p3)**2*nDp3**2 + 16.D0*s(p2,p3)*s(
     &    p2,p4)*nDp3**2 + 8.D0*s(p2,p4)**2*nDp3**2 )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*s123**(-2)*xn * ( 2.D0*s(p1,p2)
     &    *s(p1,p4)**2*nDn + 4.D0*s(p1,p2)*s(p1,p4)*s(p2,p4)*nDn + 4.D0
     &    *s(p1,p2)*s(p1,p4)*s(p3,p4)*nDn + 2.D0*s(p1,p2)*s(p2,p4)**2*
     &    nDn + 4.D0*s(p1,p2)*s(p2,p4)*s(p3,p4)*nDn + 2.D0*s(p1,p2)*s(
     &    p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*s123**(-1)*s124**(-1)*xn * (
     &     - 2.D0*s(p1,p2)**2*s(p1,p4)*nDn - 2.D0*s(p1,p2)**2*s(p3,p4)*
     &    nDn + 16.D0*s(p1,p2)*s(p1,p4)*nDp3**2 + 16.D0*s(p1,p2)*s(p1,
     &    p4)*nDp1*nDp3 + 4.D0*s(p1,p2)*s(p1,p4)*s(p3,p4)*nDn - 16.D0*
     &    s(p1,p2)*s(p3,p4)*nDp1*nDp3 - 16.D0*s(p1,p2)*s(p3,p4)*nDp1**2
     &     + 4.D0*s(p1,p2)*s(p3,p4)**2*nDn + 32.D0*s(p1,p4)**2*nDp3**2
     &     + 16.D0*s(p1,p4)**2*nDp2*nDp3 + 32.D0*s(p1,p4)**2*nDp1*nDp3
     &     + 16.D0*s(p1,p4)*s(p3,p4)*nDp2*nDp3 - 32.D0*s(p1,p4)*s(p3,p4
     &    )*nDp1*nDp3 - 16.D0*s(p1,p4)*s(p3,p4)*nDp1*nDp2 - 32.D0*s(p1,
     &    p4)*s(p3,p4)*nDp1**2 - 2.D0*s(p1,p4)*s(p3,p4)**2*nDn - 16.D0*
     &    s(p3,p4)**2*nDp1*nDp2 - 2.D0*s(p3,p4)**3*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*s123**(-1)*xn * ( 16.D0*s(p1,p2
     &    )*nDp2*nDp3 + 8.D0*s(p1,p2)*nDp1*nDp2 + 12.D0*s(p1,p2)*s(p1,
     &    p4)*nDn + 4.D0*s(p1,p2)*s(p2,p4)*nDn + 8.D0*s(p1,p2)*s(p3,p4)
     &    *nDn - 24.D0*s(p1,p4)*nDp3**2 + 16.D0*s(p1,p4)*nDp2*nDp3 - 24.
     &    D0*s(p1,p4)*nDp1*nDp3 - 16.D0*s(p1,p4)*nDp1*nDp2 - 6.D0*s(p1,
     &    p4)**2*nDn - 6.D0*s(p1,p4)*s(p2,p4)*nDn - 8.D0*s(p1,p4)*s(p3,
     &    p4)*nDn + 16.D0*s(p2,p4)*nDp1*nDp3 + 16.D0*s(p2,p4)*nDp1**2
     &     - 2.D0*s(p2,p4)*s(p3,p4)*nDn + 24.D0*s(p3,p4)*nDp1*nDp3 - 16.
     &    D0*s(p3,p4)*nDp1*nDp2 + 24.D0*s(p3,p4)*nDp1**2 - 4.D0*s(p3,p4
     &    )**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*s124**(-1)*xn * ( 2.D0*s(p1,p2)
     &    *s(p1,p4)*nDn - 4.D0*s(p1,p4)*s(p3,p4)*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*xn * ( 8.D0*nDp1*nDp3 - 8.D0*
     &    nDp1*nDp2 - 8.D0*s(p1,p4)*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*is(p2,p4)*s123**(-1)*xn * (  -
     &    8.D0*s(p1,p2)*s(p1,p4)*nDp2**2 + 2.D0*s(p1,p2)*s(p1,p4)**2*
     &    nDn + 4.D0*s(p1,p2)*s(p1,p4)*s(p3,p4)*nDn - 8.D0*s(p1,p2)*s(
     &    p3,p4)*nDp2**2 + 2.D0*s(p1,p2)*s(p3,p4)**2*nDn + 8.D0*s(p1,p4
     &    )**2*nDp2*nDp3 + 8.D0*s(p1,p4)*s(p3,p4)*nDp2*nDp3 - 8.D0*s(p1
     &    ,p4)*s(p3,p4)*nDp1*nDp2 - 8.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*is(p2,p4)*xn * ( 8.D0*s(p1,p4)*
     &    nDp2*nDp3 + 8.D0*s(p1,p4)*nDp2**2 - 2.D0*s(p1,p4)**2*nDn - 2.D
     &    0*s(p1,p4)*s(p3,p4)*nDn - 8.D0*s(p3,p4)*nDp1*nDp2 )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*is(p3,p4)*s123**(-1)*xn * ( 8.D0
     &    *s(p1,p2)*s(p1,p4)*nDp2*nDp3 + 4.D0*s(p1,p2)*s(p1,p4)**2*nDn
     &     + 4.D0*s(p1,p2)*s(p1,p4)*s(p2,p4)*nDn - 8.D0*s(p1,p2)*s(p2,
     &    p4)*nDp3**2 - 8.D0*s(p1,p2)*s(p2,p4)*nDp1*nDp3 + 16.D0*s(p1,
     &    p4)**2*nDp2*nDp3 - 2.D0*s(p1,p4)**3*nDn - 4.D0*s(p1,p4)**2*s(
     &    p2,p4)*nDn - 16.D0*s(p1,p4)*s(p2,p4)*nDp3**2 - 16.D0*s(p1,p4)
     &    *s(p2,p4)*nDp1*nDp3 - 2.D0*s(p1,p4)*s(p2,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p1,p3)*is(p3,p4)*xn * (  - 8.D0*s(p1,
     &    p4)*nDp3**2 - 8.D0*s(p1,p4)*nDp2*nDp3 - 4.D0*s(p1,p4)**2*nDn
     &     - 4.D0*s(p1,p4)*s(p2,p4)*nDn + 8.D0*s(p2,p4)*nDp1*nDp3 )
      qqgghn_ab = qqgghn_ab + is(p2,p4)*s123**(-1)*s124**(-1)*xn * ( 2.D
     &    0*s(p1,p2)**3*nDn + 2.D0*s(p1,p2)**2*s(p1,p3)*nDn - 6.D0*s(p1
     &    ,p2)**2*s(p3,p4)*nDn - 16.D0*s(p1,p2)*s(p1,p3)*nDp2*nDp3 - 4.D
     &    0*s(p1,p2)*s(p1,p3)*s(p3,p4)*nDn + 6.D0*s(p1,p2)*s(p3,p4)**2*
     &    nDn - 16.D0*s(p1,p3)**2*nDp2*nDp3 + 16.D0*s(p1,p3)*s(p3,p4)*
     &    nDp2*nDp3 - 16.D0*s(p1,p3)*s(p3,p4)*nDp1*nDp2 + 2.D0*s(p1,p3)
     &    *s(p3,p4)**2*nDn - 2.D0*s(p3,p4)**3*nDn )
      qqgghn_ab = qqgghn_ab + is(p2,p4)*s123**(-1)*xn * (  - 2.D0*s(p1,
     &    p2)**2*nDn - 2.D0*s(p1,p2)*s(p1,p3)*nDn + 2.D0*s(p1,p2)*s(p1,
     &    p4)*nDn + 6.D0*s(p1,p2)*s(p3,p4)*nDn + 16.D0*s(p1,p3)*nDp2*
     &    nDp3 + 2.D0*s(p1,p3)*s(p1,p4)*nDn + 4.D0*s(p1,p3)*s(p3,p4)*
     &    nDn + 8.D0*s(p1,p4)*nDp2*nDp3 - 16.D0*s(p1,p4)*nDp2**2 + 4.D0
     &    *s(p1,p4)**2*nDn + 4.D0*s(p1,p4)*s(p3,p4)*nDn - 24.D0*s(p3,p4
     &    )*nDp2**2 - 8.D0*s(p3,p4)*nDp1*nDp2 - 2.D0*s(p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p2,p4)*s124**(-2)*xn * ( 2.D0*s(p1,p2)
     &    *s(p1,p3)**2*nDn + 4.D0*s(p1,p2)*s(p1,p3)*s(p2,p3)*nDn + 4.D0
     &    *s(p1,p2)*s(p1,p3)*s(p3,p4)*nDn + 2.D0*s(p1,p2)*s(p2,p3)**2*
     &    nDn + 4.D0*s(p1,p2)*s(p2,p3)*s(p3,p4)*nDn + 2.D0*s(p1,p2)*s(
     &    p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p2,p4)*s124**(-1)*xn * (  - 16.D0*s(p1
     &    ,p2)*nDp2*nDp3 - 16.D0*s(p1,p2)*nDp2**2 - 4.D0*s(p1,p2)**2*
     &    nDn + 4.D0*s(p1,p2)*s(p1,p3)*nDn + 12.D0*s(p1,p2)*s(p2,p3)*
     &    nDn + 12.D0*s(p1,p2)*s(p3,p4)*nDn + 32.D0*s(p1,p3)*nDp2*nDp3
     &     + 16.D0*s(p1,p3)*nDp2**2 - 2.D0*s(p1,p3)**2*nDn - 8.D0*s(p1,
     &    p3)*s(p2,p3)*nDn - 4.D0*s(p1,p3)*s(p3,p4)*nDn - 32.D0*s(p2,p3
     &    )*nDp2*nDp3 - 16.D0*s(p2,p3)*nDp1*nDp2 - 6.D0*s(p2,p3)**2*nDn
     &     - 10.D0*s(p2,p3)*s(p3,p4)*nDn - 24.D0*s(p3,p4)*nDp2*nDp3 -
     &    16.D0*s(p3,p4)*nDp1*nDp2 - 6.D0*s(p3,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p2,p4)*xn * ( 16.D0*nDp2*nDp3 + 16.D0*
     &    nDp2**2 + 4.D0*s(p1,p2)*nDn - 4.D0*s(p1,p3)*nDn - 4.D0*s(p1,
     &    p4)*nDn - 12.D0*s(p2,p3)*nDn - 8.D0*s(p3,p4)*nDn )
      qqgghn_ab = qqgghn_ab + is(p2,p4)*is(p3,p4)*s124**(-1)*xn * (  -
     &    2.D0*s(p1,p2)**2*s(p2,p3)*nDn + 4.D0*s(p1,p2)*s(p1,p3)*s(p2,
     &    p3)*nDn - 16.D0*s(p1,p2)*s(p2,p3)*nDp2*nDp3 - 16.D0*s(p1,p2)*
     &    s(p2,p3)*nDp2**2 + 4.D0*s(p1,p2)*s(p2,p3)**2*nDn - 8.D0*s(p1,
     &    p3)**2*nDp2*nDp3 - 8.D0*s(p2,p3)**2*nDp2*nDp3 )
      qqgghn_ab = qqgghn_ab + is(p2,p4)*is(p3,p4)*xn * ( 2.D0*s(p1,p2)*
     &    s(p2,p3)*nDn - 4.D0*s(p1,p3)*s(p2,p3)*nDn - 2.D0*s(p1,p4)*s(
     &    p2,p3)*nDn + 16.D0*s(p2,p3)*nDp2*nDp3 + 16.D0*s(p2,p3)*
     &    nDp2**2 - 4.D0*s(p2,p3)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p3,p4)*s123**(-1)*xn * (  - 8.D0*s(p1,
     &    p4)*nDp3**2 + 16.D0*s(p1,p4)*nDp2*nDp3 + 12.D0*s(p1,p4)**2*
     &    nDn + 16.D0*s(p1,p4)*s(p2,p4)*nDn - 24.D0*s(p2,p4)*nDp3**2 -
     &    16.D0*s(p2,p4)*nDp1*nDp3 + 4.D0*s(p2,p4)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p3,p4)*s124**(-1)*xn * (  - 2.D0*s(p1,
     &    p2)*s(p1,p3)*nDn - 2.D0*s(p1,p2)*s(p2,p3)*nDn - 16.D0*s(p1,p3
     &    )*nDp1*nDp3 - 16.D0*s(p1,p3)*nDp1**2 + 4.D0*s(p1,p3)**2*nDn
     &     + 8.D0*s(p1,p3)*s(p2,p3)*nDn + 16.D0*s(p1,p4)*nDp3**2 + 16.D0
     &    *s(p1,p4)*nDp1*nDp3 - 16.D0*s(p2,p3)*nDp2*nDp3 - 16.D0*s(p2,
     &    p3)*nDp2**2 + 4.D0*s(p2,p3)**2*nDn )
      qqgghn_ab = qqgghn_ab + is(p3,p4)*xn * ( 2.D0*s(p1,p3)*nDn )
      qqgghn_ab = qqgghn_ab + (nDp2*is(p2,p4))**2*s124**(-1)*xn * ( 8.D0*s(p1,
     &    p3)**2 + 8.D0*s(p2,p3)**2 + 16.D0*s(p2,p3)*s(
     &    p3,p4) + 8.D0*s(p3,p4)**2 )

      qqgghn_ba=  + s123**(-2)*xn * ( 6.D0*s(p1,p4)**2*nDn + 12.D0*s(p1
     &    ,p4)*s(p2,p4)*nDn + 12.D0*s(p1,p4)*s(p3,p4)*nDn + 6.D0*s(p2,
     &    p4)**2*nDn + 12.D0*s(p2,p4)*s(p3,p4)*nDn + 6.D0*s(p3,p4)**2*
     &    nDn )
      qqgghn_ba = qqgghn_ba + s123**(-1)*s124**(-1)*xn * (  - 16.D0*s(
     &    p1,p2)*nDp2*nDp3 - 8.D0*s(p1,p2)*s(p2,p4)*nDn - 4.D0*s(p1,p2)
     &    *s(p3,p4)*nDn - 32.D0*s(p2,p3)*nDp2*nDp3 - 16.D0*s(p2,p3)*
     &    nDp1*nDp3 - 8.D0*s(p2,p3)*s(p2,p4)*nDn + 32.D0*s(p2,p4)*
     &    nDp3**2 - 16.D0*s(p2,p4)*nDp1*nDp3 + 16.D0*s(p2,p4)*s(p3,p4)*
     &    nDn - 32.D0*s(p3,p4)*nDp2*nDp3 - 48.D0*s(p3,p4)*nDp2**2 - 16.D
     &    0*s(p3,p4)*nDp1*nDp2 + 8.D0*s(p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn * ( 24.D0*nDp2*nDp3 + 48.D0
     &    *nDp1*nDp3 + 16.D0*nDp1*nDp2 + 2.D0*s(p1,p2)*nDn + 16.D0*s(p1
     &    ,p4)*nDn + 2.D0*s(p2,p3)*nDn + 36.D0*s(p2,p4)*nDn + 16.D0*s(
     &    p3,p4)*nDn )
      qqgghn_ba = qqgghn_ba + s124**(-2)*xn * ( 2.D0*s(p1,p3)**2*nDn +
     &    4.D0*s(p1,p3)*s(p2,p3)*nDn + 4.D0*s(p1,p3)*s(p3,p4)*nDn + 2.D0
     &    *s(p2,p3)**2*nDn + 4.D0*s(p2,p3)*s(p3,p4)*nDn + 2.D0*s(p3,p4)
     &    **2*nDn )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn * (  - 16.D0*nDp2*nDp3 - 32.
     &    D0*nDp2**2 - 16.D0*nDp1*nDp3 - 16.D0*nDp1**2 - 6.D0*s(p1,p2)*
     &    nDn + 16.D0*s(p1,p3)*nDn + 16.D0*s(p2,p3)*nDn + 8.D0*s(p2,p4)
     &    *nDn + 20.D0*s(p3,p4)*nDn )
      qqgghn_ba = qqgghn_ba + xn * ( 6.D0*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*s123**(-2)*xn * ( 8.D0*s(p1,p4)
     &    **2*s(p2,p3)*nDn + 16.D0*s(p1,p4)*s(p2,p3)*s(p2,p4)*nDn + 16.D
     &    0*s(p1,p4)*s(p2,p3)*s(p3,p4)*nDn + 8.D0*s(p2,p3)*s(p2,p4)**2*
     &    nDn + 16.D0*s(p2,p3)*s(p2,p4)*s(p3,p4)*nDn + 8.D0*s(p2,p3)*s(
     &    p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*s123**(-1)*s124**(-1)*xn * (
     &     - 16.D0*s(p2,p3)**2*nDp2*nDp3 - 16.D0*s(p2,p3)**2*nDp1*nDp3
     &     + 16.D0*s(p2,p3)*s(p2,p4)*nDp3**2 - 16.D0*s(p2,p3)*s(p2,p4)*
     &    nDp2*nDp3 - 16.D0*s(p2,p3)*s(p2,p4)*nDp1*nDp3 + 16.D0*s(p2,p3
     &    )*s(p2,p4)*s(p3,p4)*nDn - 16.D0*s(p2,p3)*s(p3,p4)*nDp2*nDp3
     &     - 48.D0*s(p2,p3)*s(p3,p4)*nDp2**2 + 16.D0*s(p2,p3)*s(p3,p4)*
     &    nDp1*nDp3 - 48.D0*s(p2,p3)*s(p3,p4)*nDp1*nDp2 + 32.D0*s(p2,p4
     &    )**2*nDp3**2 + 16.D0*s(p2,p4)**2*nDp2*nDp3 + 16.D0*s(p2,p4)**
     &    2*nDp1*nDp3 - 32.D0*s(p2,p4)*s(p3,p4)*nDp2*nDp3 - 64.D0*s(p2,
     &    p4)*s(p3,p4)*nDp2**2 + 32.D0*s(p2,p4)*s(p3,p4)*nDp1*nDp3 - 80.
     &    D0*s(p2,p4)*s(p3,p4)*nDp1*nDp2 - 16.D0*s(p2,p4)*s(p3,p4)*
     &    nDp1**2 - 8.D0*s(p2,p4)*s(p3,p4)**2*nDn - 32.D0*s(p3,p4)**2*
     &    nDp1*nDp2 - 4.D0*s(p3,p4)**3*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*s123**(-1)*xn * ( 8.D0*s(p1,p4)
     &    *nDp2*nDp3 + 48.D0*s(p1,p4)*nDp2**2 + 16.D0*s(p1,p4)*nDp1*
     &    nDp2 - 6.D0*s(p1,p4)**2*nDn + 18.D0*s(p1,p4)*s(p2,p3)*nDn -
     &    20.D0*s(p1,p4)*s(p2,p4)*nDn - 10.D0*s(p1,p4)*s(p3,p4)*nDn +
     &    48.D0*s(p2,p3)*nDp2*nDp3 + 64.D0*s(p2,p3)*nDp1*nDp3 + 16.D0*
     &    s(p2,p3)*nDp1*nDp2 + 36.D0*s(p2,p3)*s(p2,p4)*nDn + 8.D0*s(p2,
     &    p3)*s(p3,p4)*nDn - 32.D0*s(p2,p4)*nDp3**2 - 24.D0*s(p2,p4)*
     &    nDp2*nDp3 + 56.D0*s(p2,p4)*nDp1*nDp3 - 48.D0*s(p2,p4)*nDp1*
     &    nDp2 - 16.D0*s(p2,p4)*nDp1**2 - 14.D0*s(p2,p4)**2*nDn - 26.D0
     &    *s(p2,p4)*s(p3,p4)*nDn + 32.D0*s(p3,p4)*nDp2*nDp3 + 56.D0*s(
     &    p3,p4)*nDp2**2 - 8.D0*s(p3,p4)*nDp1*nDp2 - 10.D0*s(p3,p4)**2*
     &    nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*s124**(-1)*xn * ( 16.D0*s(p1,p3
     &    )*nDp2**2 - 32.D0*s(p1,p3)*nDp1*nDp3 - 16.D0*s(p1,p3)*nDp1*
     &    nDp2 - 4.D0*s(p1,p3)**2*nDn - 8.D0*s(p1,p3)*s(p2,p3)*nDn - 10.
     &    D0*s(p1,p3)*s(p3,p4)*nDn + 32.D0*s(p2,p3)*nDp1*nDp3 - 16.D0*
     &    s(p2,p3)*nDp1*nDp2 + 16.D0*s(p2,p3)*nDp1**2 - 4.D0*s(p2,p3)**
     &    2*nDn + 16.D0*s(p2,p3)*s(p2,p4)*nDn - 10.D0*s(p2,p3)*s(p3,p4)
     &    *nDn - 16.D0*s(p2,p4)*nDp3**2 - 64.D0*s(p2,p4)*nDp2*nDp3 - 64.
     &    D0*s(p2,p4)*nDp2**2 + 32.D0*s(p2,p4)*nDp1*nDp3 - 48.D0*s(p2,
     &    p4)*nDp1*nDp2 + 16.D0*s(p2,p4)*nDp1**2 - 8.D0*s(p2,p4)*s(p3,
     &    p4)*nDn + 8.D0*s(p3,p4)*nDp2*nDp3 + 32.D0*s(p3,p4)*nDp2**2 -
     &    24.D0*s(p3,p4)*nDp1*nDp3 - 32.D0*s(p3,p4)*nDp1*nDp2 - 10.D0*
     &    s(p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*xn * ( 8.D0*nDp3**2 + 16.D0*
     &    nDp2*nDp3 + 32.D0*nDp2**2 - 32.D0*nDp1*nDp3 - 8.D0*nDp1*nDp2
     &     + 8.D0*nDp1**2 - 16.D0*s(p1,p3)*nDn - 24.D0*s(p1,p4)*nDn -
     &    20.D0*s(p2,p3)*nDn - 42.D0*s(p2,p4)*nDn - 24.D0*s(p3,p4)*nDn
     &     )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*s123**(-2)*xn * ( 4.D0*s(p1,
     &    p4)**2*s(p2,p3)**2*nDn + 8.D0*s(p1,p4)*s(p2,p3)**2*s(p2,p4)*
     &    nDn + 8.D0*s(p1,p4)*s(p2,p3)**2*s(p3,p4)*nDn + 4.D0*s(p2,p3)
     &    **2*s(p2,p4)**2*nDn + 8.D0*s(p2,p3)**2*s(p2,p4)*s(p3,p4)*nDn
     &     + 4.D0*s(p2,p3)**2*s(p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*s123**(-1)*s124**(-1)*xn
     &  * (  - 16.D0*s(p2,p3)**2*s(p3,p4)*nDp2**2 - 32.D0*s(p2,p3)**2*
     &    s(p3,p4)*nDp1*nDp2 - 16.D0*s(p2,p3)**2*s(p3,p4)*nDp1**2 - 32.D
     &    0*s(p2,p3)*s(p2,p4)*s(p3,p4)*nDp2**2 - 64.D0*s(p2,p3)*s(p2,p4
     &    )*s(p3,p4)*nDp1*nDp2 - 32.D0*s(p2,p3)*s(p2,p4)*s(p3,p4)*
     &    nDp1**2 - 8.D0*s(p2,p3)*s(p2,p4)*s(p3,p4)**2*nDn - 16.D0*s(p2
     &    ,p4)**2*s(p3,p4)*nDp2**2 - 32.D0*s(p2,p4)**2*s(p3,p4)*nDp1*
     &    nDp2 - 16.D0*s(p2,p4)**2*s(p3,p4)*nDp1**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*s123**(-1)*xn * (  - 32.D0*
     &    s(p1,p4)**2*nDp2**2 - 32.D0*s(p1,p4)*s(p2,p3)*nDp2*nDp3 + 32.D
     &    0*s(p1,p4)*s(p2,p3)*nDp2**2 + 32.D0*s(p1,p4)*s(p2,p3)*nDp1*
     &    nDp2 + 16.D0*s(p1,p4)*s(p2,p3)**2*nDn - 8.D0*s(p1,p4)*s(p2,p3
     &    )*s(p2,p4)*nDn + 64.D0*s(p1,p4)*s(p2,p4)*nDp1*nDp2 - 16.D0*s(
     &    p1,p4)*s(p3,p4)*nDp2**2 + 32.D0*s(p2,p3)**2*nDp2*nDp3 + 32.D0
     &    *s(p2,p3)**2*nDp1*nDp3 + 16.D0*s(p2,p3)**2*s(p2,p4)*nDn + 8.D0
     &    *s(p2,p3)**2*s(p3,p4)*nDn + 32.D0*s(p2,p3)*s(p2,p4)*nDp1*nDp3
     &     - 32.D0*s(p2,p3)*s(p2,p4)*nDp1*nDp2 - 32.D0*s(p2,p3)*s(p2,p4
     &    )*nDp1**2 - 8.D0*s(p2,p3)*s(p2,p4)**2*nDn - 16.D0*s(p2,p3)*s(
     &    p2,p4)*s(p3,p4)*nDn + 32.D0*s(p2,p3)*s(p3,p4)*nDp2**2 + 32.D0
     &    *s(p2,p3)*s(p3,p4)*nDp1*nDp2 - 32.D0*s(p2,p4)**2*nDp1**2 + 16.
     &    D0*s(p2,p4)*s(p3,p4)*nDp2**2 + 32.D0*s(p2,p4)*s(p3,p4)*nDp1*
     &    nDp2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*s124**(-2)*xn * ( 4.D0*s(p1,
     &    p3)**2*s(p2,p4)**2*nDn + 8.D0*s(p1,p3)*s(p2,p3)*s(p2,p4)**2*
     &    nDn + 8.D0*s(p1,p3)*s(p2,p4)**2*s(p3,p4)*nDn + 4.D0*s(p2,p3)
     &    **2*s(p2,p4)**2*nDn + 8.D0*s(p2,p3)*s(p2,p4)**2*s(p3,p4)*nDn
     &     + 4.D0*s(p2,p4)**2*s(p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*s124**(-1)*xn * (  - 8.D0*s(
     &    p1,p3)*s(p2,p3)*s(p2,p4)*nDn + 32.D0*s(p1,p3)*s(p2,p4)*nDp2*
     &    nDp3 + 32.D0*s(p1,p3)*s(p2,p4)*nDp2**2 + 32.D0*s(p1,p3)*s(p2,
     &    p4)*nDp1*nDp2 + 16.D0*s(p1,p3)*s(p2,p4)**2*nDn - 16.D0*s(p1,
     &    p3)*s(p3,p4)*nDp2**2 - 8.D0*s(p2,p3)**2*s(p2,p4)*nDn - 32.D0*
     &    s(p2,p3)*s(p2,p4)*nDp1*nDp3 - 32.D0*s(p2,p3)*s(p2,p4)*nDp1*
     &    nDp2 - 32.D0*s(p2,p3)*s(p2,p4)*nDp1**2 + 16.D0*s(p2,p3)*s(p2,
     &    p4)**2*nDn - 16.D0*s(p2,p3)*s(p2,p4)*s(p3,p4)*nDn + 16.D0*s(
     &    p2,p3)*s(p3,p4)*nDp2**2 + 32.D0*s(p2,p3)*s(p3,p4)*nDp1*nDp2
     &     - 32.D0*s(p2,p4)**2*nDp2*nDp3 - 32.D0*s(p2,p4)**2*nDp2**2 -
     &    32.D0*s(p2,p4)**2*nDp1*nDp3 - 64.D0*s(p2,p4)**2*nDp1*nDp2 -
     &    32.D0*s(p2,p4)**2*nDp1**2 + 8.D0*s(p2,p4)**2*s(p3,p4)*nDn +
     &    32.D0*s(p2,p4)*s(p3,p4)*nDp2**2 + 32.D0*s(p2,p4)*s(p3,p4)*
     &    nDp1*nDp2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*xn * (  - 32.D0*s(p1,p3)*
     &    nDp2**2 - 64.D0*s(p1,p4)*nDp2**2 - 32.D0*s(p2,p3)*nDp2*nDp3
     &     + 32.D0*s(p2,p3)*nDp1*nDp2 + 4.D0*s(p2,p3)**2*nDn - 32.D0*s(
     &    p2,p3)*s(p2,p4)*nDn + 32.D0*s(p2,p4)*nDp2*nDp3 + 32.D0*s(p2,
     &    p4)*nDp2**2 + 96.D0*s(p2,p4)*nDp1*nDp2 + 4.D0*s(p2,p4)**2*nDn
     &     - 32.D0*s(p3,p4)*nDp2**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*is(p3,p4)*s123**(-1)*xn * (
     &     - 16.D0*s(p1,p4)**3*nDp2**2 - 32.D0*s(p1,p4)**2*s(p2,p3)*
     &    nDp2*nDp3 + 8.D0*s(p1,p4)**2*s(p2,p3)**2*nDn - 16.D0*s(p1,p4)
     &    **2*s(p2,p4)*nDp2**2 + 32.D0*s(p1,p4)**2*s(p2,p4)*nDp1*nDp2
     &     - 16.D0*s(p1,p4)*s(p2,p3)**2*nDp3**2 + 16.D0*s(p1,p4)*s(p2,
     &    p3)**2*s(p2,p4)*nDn - 32.D0*s(p1,p4)*s(p2,p3)*s(p2,p4)*nDp2*
     &    nDp3 + 32.D0*s(p1,p4)*s(p2,p3)*s(p2,p4)*nDp1*nDp3 + 32.D0*s(
     &    p1,p4)*s(p2,p4)**2*nDp1*nDp2 - 16.D0*s(p1,p4)*s(p2,p4)**2*
     &    nDp1**2 - 16.D0*s(p2,p3)**2*s(p2,p4)*nDp3**2 + 8.D0*s(p2,p3)
     &    **2*s(p2,p4)**2*nDn + 32.D0*s(p2,p3)*s(p2,p4)**2*nDp1*nDp3 -
     &    16.D0*s(p2,p4)**3*nDp1**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*is(p3,p4)*s124**(-1)*xn * (
     &    8.D0*s(p1,p3)**2*s(p2,p4)**2*nDn + 16.D0*s(p1,p3)*s(p2,p3)*s(
     &    p2,p4)**2*nDn - 16.D0*s(p1,p3)*s(p2,p4)**2*nDp3**2 - 32.D0*s(
     &    p1,p3)*s(p2,p4)**2*nDp2*nDp3 - 16.D0*s(p1,p3)*s(p2,p4)**2*
     &    nDp2**2 - 32.D0*s(p1,p3)*s(p2,p4)**2*nDp1*nDp3 - 32.D0*s(p1,
     &    p3)*s(p2,p4)**2*nDp1*nDp2 - 16.D0*s(p1,p3)*s(p2,p4)**2*
     &    nDp1**2 + 8.D0*s(p2,p3)**2*s(p2,p4)**2*nDn - 16.D0*s(p2,p3)*
     &    s(p2,p4)**2*nDp3**2 - 32.D0*s(p2,p3)*s(p2,p4)**2*nDp2*nDp3 -
     &    16.D0*s(p2,p3)*s(p2,p4)**2*nDp2**2 - 32.D0*s(p2,p3)*s(p2,p4)
     &    **2*nDp1*nDp3 - 32.D0*s(p2,p3)*s(p2,p4)**2*nDp1*nDp2 - 16.D0*
     &    s(p2,p3)*s(p2,p4)**2*nDp1**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*is(p3,p4)*xn * (  - 16.D0*s(
     &    p1,p3)*s(p1,p4)*nDp2**2 - 16.D0*s(p1,p3)*s(p2,p3)*s(p2,p4)*
     &    nDn + 32.D0*s(p1,p3)*s(p2,p4)*nDp2*nDp3 + 16.D0*s(p1,p3)*s(p2
     &    ,p4)*nDp2**2 + 32.D0*s(p1,p3)*s(p2,p4)*nDp1*nDp2 + 8.D0*s(p1,
     &    p3)*s(p2,p4)**2*nDn - 32.D0*s(p1,p4)**2*nDp2**2 - 32.D0*s(p1,
     &    p4)*s(p2,p3)*nDp2*nDp3 - 16.D0*s(p1,p4)*s(p2,p3)*nDp2**2 + 8.D
     &    0*s(p1,p4)*s(p2,p3)**2*nDn - 16.D0*s(p1,p4)*s(p2,p3)*s(p2,p4)
     &    *nDn + 32.D0*s(p1,p4)*s(p2,p4)*nDp2*nDp3 + 64.D0*s(p1,p4)*s(
     &    p2,p4)*nDp1*nDp2 - 8.D0*s(p2,p3)**2*s(p2,p4)*nDn + 32.D0*s(p2
     &    ,p3)*s(p2,p4)*nDp3**2 + 32.D0*s(p2,p3)*s(p2,p4)*nDp2*nDp3 +
     &    16.D0*s(p2,p3)*s(p2,p4)*nDp2**2 + 32.D0*s(p2,p3)*s(p2,p4)*
     &    nDp1*nDp3 + 32.D0*s(p2,p3)*s(p2,p4)*nDp1*nDp2 - 8.D0*s(p2,p3)
     &    *s(p2,p4)**2*nDn - 32.D0*s(p2,p4)**2*nDp1*nDp3 - 32.D0*s(p2,
     &    p4)**2*nDp1**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)**2*is(p3,p4)**2*xn * ( 4.D0*s(
     &    p1,p3)**2*s(p2,p4)**2*nDn - 8.D0*s(p1,p3)*s(p1,p4)*s(p2,p3)*
     &    s(p2,p4)*nDn + 4.D0*s(p1,p4)**2*s(p2,p3)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p1,p4)*s123**(-1)*s124**(-1)
     & *xn * (  - 16.D0*s(p2,p3)**2*s(p3,p4)*nDp1*nDp2 - 16.D0*s(p2,p3)
     &    **2*s(p3,p4)*nDp1**2 + 16.D0*s(p2,p3)*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p1,p4)*s123**(-1)*xn * ( 16.D
     &    0*s(p2,p3)**2*nDp1*nDp3 - 16.D0*s(p2,p3)*s(p2,p4)*nDp1**2 + 2.
     &    D0*s(p2,p3)*s(p2,p4)**2*nDn - 16.D0*s(p2,p3)*s(p3,p4)*nDp1*
     &    nDp3 - 32.D0*s(p2,p3)*s(p3,p4)*nDp1**2 - 2.D0*s(p2,p3)*s(p3,
     &    p4)**2*nDn + 8.D0*s(p2,p4)**2*nDp1*nDp3 - 8.D0*s(p2,p4)*s(p3,
     &    p4)*nDp1*nDp2 - 8.D0*s(p2,p4)*s(p3,p4)*nDp1**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p1,p4)*s124**(-1)*xn * ( 8.D0
     &    *s(p1,p3)**2*nDp1*nDp2 + 8.D0*s(p1,p3)**2*nDp1**2 + 16.D0*s(
     &    p1,p3)*s(p3,p4)*nDp1**2 + 8.D0*s(p2,p3)**2*nDp1*nDp2 + 8.D0*
     &    s(p2,p3)**2*nDp1**2 + 32.D0*s(p2,p3)*s(p3,p4)*nDp1*nDp2 - 8.D0
     &    *s(p3,p4)**2*nDp1*nDp2 + 8.D0*s(p3,p4)**2*nDp1**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p1,p4)*xn * ( 16.D0*s(p1,p3)
     &    *nDp1*nDp2 + 16.D0*s(p1,p3)*nDp1**2 + 4.D0*s(p1,p3)**2*nDn +
     &    4.D0*s(p1,p3)*s(p2,p3)*nDn + 4.D0*s(p1,p3)*s(p2,p4)*nDn + 2.D0
     &    *s(p1,p3)*s(p3,p4)*nDn - 32.D0*s(p2,p3)*nDp1*nDp3 + 16.D0*s(
     &    p2,p3)*nDp1*nDp2 - 16.D0*s(p2,p3)*nDp1**2 - 16.D0*s(p2,p4)*
     &    nDp1*nDp3 + 8.D0*s(p2,p4)*nDp1*nDp2 - 8.D0*s(p2,p4)*nDp1**2
     &     - 2.D0*s(p2,p4)**2*nDn + 8.D0*s(p3,p4)*nDp1*nDp3 + 16.D0*s(
     &    p3,p4)*nDp1*nDp2 + 16.D0*s(p3,p4)*nDp1**2 + 2.D0*s(p3,p4)**2*
     &    nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p1,p4)*is(p3,p4)*xn * (  - 8.
     &    D0*s(p1,p3)**2*nDp1*nDp3 + 4.D0*s(p1,p3)**2*s(p2,p4)*nDn + 4.D
     &    0*s(p1,p3)*s(p2,p3)*s(p2,p4)*nDn - 16.D0*s(p1,p3)*s(p2,p4)*
     &    nDp1*nDp3 - 16.D0*s(p1,p3)*s(p2,p4)*nDp1**2 + 2.D0*s(p1,p3)*
     &    s(p2,p4)**2*nDn - 8.D0*s(p2,p3)**2*nDp1*nDp3 - 16.D0*s(p2,p3)
     &    *s(p2,p4)*nDp1*nDp3 - 8.D0*s(p2,p4)**2*nDp1*nDp3 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p2,p3)*s123**(-1)*s124**(-1)
     & *xn * ( 16.D0*s(p2,p4)**3*nDp3**2 + 16.D0*s(p2,p4)**3*nDp2*nDp3
     &     + 16.D0*s(p2,p4)**3*nDp1*nDp3 - 16.D0*s(p2,p4)**2*s(p3,p4)*
     &    nDp2*nDp3 - 16.D0*s(p2,p4)**2*s(p3,p4)*nDp2**2 + 16.D0*s(p2,
     &    p4)**2*s(p3,p4)*nDp1*nDp3 - 16.D0*s(p2,p4)**2*s(p3,p4)*nDp1*
     &    nDp2 - 16.D0*s(p2,p4)*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p2,p3)*s123**(-1)*xn * (  -
     &    8.D0*s(p1,p4)**2*nDp2*nDp3 + 8.D0*s(p1,p4)*s(p2,p4)*nDp3**2
     &     - 8.D0*s(p1,p4)*s(p2,p4)*nDp2*nDp3 + 8.D0*s(p1,p4)*s(p2,p4)*
     &    nDp1*nDp3 + 16.D0*s(p1,p4)*s(p2,p4)*nDp1*nDp2 - 8.D0*s(p1,p4)
     &    *s(p3,p4)*nDp2*nDp3 + 8.D0*s(p1,p4)*s(p3,p4)*nDp2**2 + 8.D0*
     &    s(p1,p4)*s(p3,p4)*nDp1*nDp2 - 8.D0*s(p2,p4)**2*nDp3**2 - 16.D0
     &    *s(p2,p4)**2*nDp2*nDp3 + 8.D0*s(p2,p4)**2*nDp1*nDp3 - 16.D0*
     &    s(p2,p4)**2*nDp1**2 + 16.D0*s(p2,p4)*s(p3,p4)*nDp2**2 - 8.D0*
     &    s(p2,p4)*s(p3,p4)*nDp1*nDp3 - 8.D0*s(p2,p4)*s(p3,p4)*nDp1*
     &    nDp2 - 8.D0*s(p2,p4)*s(p3,p4)*nDp1**2 + 8.D0*s(p3,p4)**2*
     &    nDp2**2 + 8.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p2,p3)*s124**(-1)*xn * ( 2.D0
     &    *s(p2,p4)*s(p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p2,p3)*xn * (  - 2.D0*s(p1,
     &    p4)*s(p2,p4)*nDn - 2.D0*s(p2,p4)**2*nDn + 2.D0*s(p2,p4)*s(p3,
     &    p4)*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p2,p3)*is(p3,p4)*s123**(-1)*
     & xn * (  - 16.D0*s(p1,p4)**2*s(p2,p4)*nDp2*nDp3 - 16.D0*s(p1,p4)
     &    **2*s(p2,p4)*nDp2**2 + 16.D0*s(p1,p4)*s(p2,p4)**2*nDp1*nDp3
     &     + 32.D0*s(p1,p4)*s(p2,p4)**2*nDp1*nDp2 - 16.D0*s(p2,p4)**3*
     &    nDp1**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p2,p3)*is(p3,p4)*xn * (  - 2.
     &    D0*s(p1,p4)**2*s(p2,p4)*nDn - 4.D0*s(p1,p4)*s(p2,p4)**2*nDn
     &     - 2.D0*s(p2,p4)**3*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p3,p4)*s123**(-1)*xn * (  -
     &    16.D0*s(p1,p4)**2*nDp2*nDp3 - 2.D0*s(p1,p4)**3*nDn + 8.D0*s(
     &    p1,p4)**2*s(p2,p3)*nDn - 6.D0*s(p1,p4)**2*s(p2,p4)*nDn - 32.D0
     &    *s(p1,p4)*s(p2,p3)*nDp3**2 - 16.D0*s(p1,p4)*s(p2,p3)*nDp2*
     &    nDp3 + 24.D0*s(p1,p4)*s(p2,p3)*s(p2,p4)*nDn - 32.D0*s(p1,p4)*
     &    s(p2,p4)*nDp3**2 - 48.D0*s(p1,p4)*s(p2,p4)*nDp2*nDp3 + 16.D0*
     &    s(p1,p4)*s(p2,p4)*nDp1*nDp3 - 6.D0*s(p1,p4)*s(p2,p4)**2*nDn
     &     - 16.D0*s(p2,p3)*s(p2,p4)*nDp3**2 + 16.D0*s(p2,p3)*s(p2,p4)*
     &    nDp1*nDp3 + 16.D0*s(p2,p3)*s(p2,p4)**2*nDn + 48.D0*s(p2,p4)**
     &    2*nDp1*nDp3 - 2.D0*s(p2,p4)**3*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p3,p4)*s124**(-1)*xn * (  -
     &    8.D0*s(p1,p3)**2*nDp2*nDp3 - 8.D0*s(p1,p3)**2*nDp1*nDp3 + 8.D0
     &    *s(p1,p3)*s(p2,p3)*s(p2,p4)*nDn - 16.D0*s(p1,p3)*s(p2,p4)*
     &    nDp3**2 - 16.D0*s(p1,p3)*s(p2,p4)*nDp2*nDp3 - 16.D0*s(p1,p3)*
     &    s(p2,p4)*nDp2**2 + 16.D0*s(p1,p3)*s(p2,p4)*nDp1*nDp3 + 16.D0*
     &    s(p1,p3)*s(p2,p4)*nDp1**2 - 8.D0*s(p2,p3)**2*nDp2*nDp3 - 8.D0
     &    *s(p2,p3)**2*nDp1*nDp3 + 8.D0*s(p2,p3)**2*s(p2,p4)*nDn - 48.D0
     &    *s(p2,p3)*s(p2,p4)*nDp2*nDp3 - 32.D0*s(p2,p3)*s(p2,p4)*
     &    nDp2**2 - 16.D0*s(p2,p3)*s(p2,p4)*nDp1*nDp3 - 32.D0*s(p2,p3)*
     &    s(p2,p4)*nDp1*nDp2 + 32.D0*s(p2,p4)**2*nDp3**2 + 16.D0*s(p2,
     &    p4)**2*nDp2*nDp3 + 16.D0*s(p2,p4)**2*nDp1*nDp3 )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p3,p4)*xn * ( 16.D0*s(p1,p3)
     &    *nDp3**2 - 32.D0*s(p1,p3)*nDp1*nDp3 - 8.D0*s(p1,p3)*s(p1,p4)*
     &    nDn - 4.D0*s(p1,p3)*s(p2,p3)*nDn - 6.D0*s(p1,p3)*s(p2,p4)*nDn
     &     + 16.D0*s(p1,p4)*nDp3**2 + 8.D0*s(p1,p4)*nDp2*nDp3 - 16.D0*
     &    s(p1,p4)*nDp1*nDp3 - 8.D0*s(p1,p4)**2*nDn - 10.D0*s(p1,p4)*s(
     &    p2,p3)*nDn - 24.D0*s(p1,p4)*s(p2,p4)*nDn + 16.D0*s(p2,p3)*
     &    nDp2**2 - 4.D0*s(p2,p3)**2*nDn - 10.D0*s(p2,p3)*s(p2,p4)*nDn
     &     - 16.D0*s(p2,p4)*nDp3**2 - 24.D0*s(p2,p4)*nDp2*nDp3 - 8.D0*
     &    s(p2,p4)*nDp1*nDp3 - 16.D0*s(p2,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p2)*is(p3,p4)**2*xn * ( 8.D0*s(p1,
     &    p3)**2*nDp3**2 + 16.D0*s(p1,p3)*s(p1,p4)*nDp3**2 + 8.D0*s(p1,
     &    p4)**2*nDp3**2 + 8.D0*s(p2,p3)**2*nDp3**2 + 16.D0*s(p2,p3)*s(
     &    p2,p4)*nDp3**2 + 8.D0*s(p2,p4)**2*nDp3**2 )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*s123**(-1)*s124**(-1)*xn * ( 2.D
     &    0*s(p1,p2)**3*nDn + 2.D0*s(p1,p2)**2*s(p2,p3)*nDn - 6.D0*s(p1
     &    ,p2)**2*s(p3,p4)*nDn - 16.D0*s(p1,p2)*s(p2,p3)*nDp1*nDp3 - 4.D
     &    0*s(p1,p2)*s(p2,p3)*s(p3,p4)*nDn + 6.D0*s(p1,p2)*s(p3,p4)**2*
     &    nDn - 16.D0*s(p2,p3)**2*nDp1*nDp3 + 16.D0*s(p2,p3)*s(p3,p4)*
     &    nDp1*nDp3 - 16.D0*s(p2,p3)*s(p3,p4)*nDp1*nDp2 + 2.D0*s(p2,p3)
     &    *s(p3,p4)**2*nDn - 2.D0*s(p3,p4)**3*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*s123**(-1)*xn * (  - 2.D0*s(p1,
     &    p2)**2*nDn - 2.D0*s(p1,p2)*s(p2,p3)*nDn + 2.D0*s(p1,p2)*s(p2,
     &    p4)*nDn + 6.D0*s(p1,p2)*s(p3,p4)*nDn + 16.D0*s(p2,p3)*nDp1*
     &    nDp3 + 2.D0*s(p2,p3)*s(p2,p4)*nDn + 4.D0*s(p2,p3)*s(p3,p4)*
     &    nDn + 8.D0*s(p2,p4)*nDp1*nDp3 - 16.D0*s(p2,p4)*nDp1**2 + 4.D0
     &    *s(p2,p4)**2*nDn + 4.D0*s(p2,p4)*s(p3,p4)*nDn - 8.D0*s(p3,p4)
     &    *nDp1*nDp2 - 24.D0*s(p3,p4)*nDp1**2 - 2.D0*s(p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*s124**(-2)*xn * ( 2.D0*s(p1,p2)
     &    *s(p1,p3)**2*nDn + 4.D0*s(p1,p2)*s(p1,p3)*s(p2,p3)*nDn + 4.D0
     &    *s(p1,p2)*s(p1,p3)*s(p3,p4)*nDn + 2.D0*s(p1,p2)*s(p2,p3)**2*
     &    nDn + 4.D0*s(p1,p2)*s(p2,p3)*s(p3,p4)*nDn + 2.D0*s(p1,p2)*s(
     &    p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*s124**(-1)*xn * (  - 16.D0*s(p1
     &    ,p2)*nDp1*nDp3 - 16.D0*s(p1,p2)*nDp1**2 - 4.D0*s(p1,p2)**2*
     &    nDn + 12.D0*s(p1,p2)*s(p1,p3)*nDn + 4.D0*s(p1,p2)*s(p2,p3)*
     &    nDn + 12.D0*s(p1,p2)*s(p3,p4)*nDn - 32.D0*s(p1,p3)*nDp1*nDp3
     &     - 16.D0*s(p1,p3)*nDp1*nDp2 - 6.D0*s(p1,p3)**2*nDn - 8.D0*s(
     &    p1,p3)*s(p2,p3)*nDn - 10.D0*s(p1,p3)*s(p3,p4)*nDn + 32.D0*s(
     &    p2,p3)*nDp1*nDp3 + 16.D0*s(p2,p3)*nDp1**2 - 2.D0*s(p2,p3)**2*
     &    nDn - 4.D0*s(p2,p3)*s(p3,p4)*nDn - 24.D0*s(p3,p4)*nDp1*nDp3
     &     - 16.D0*s(p3,p4)*nDp1*nDp2 - 6.D0*s(p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*xn * ( 16.D0*nDp1*nDp3 + 16.D0*
     &    nDp1**2 + 4.D0*s(p1,p2)*nDn - 12.D0*s(p1,p3)*nDn - 4.D0*s(p2,
     &    p3)*nDn - 4.D0*s(p2,p4)*nDn - 8.D0*s(p3,p4)*nDn )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*is(p2,p3)*s123**(-1)*xn * (  -
     &    8.D0*s(p1,p2)*s(p2,p4)*nDp1**2 + 2.D0*s(p1,p2)*s(p2,p4)**2*
     &    nDn + 4.D0*s(p1,p2)*s(p2,p4)*s(p3,p4)*nDn - 8.D0*s(p1,p2)*s(
     &    p3,p4)*nDp1**2 + 2.D0*s(p1,p2)*s(p3,p4)**2*nDn + 8.D0*s(p2,p4
     &    )**2*nDp1*nDp3 + 8.D0*s(p2,p4)*s(p3,p4)*nDp1*nDp3 - 8.D0*s(p2
     &    ,p4)*s(p3,p4)*nDp1*nDp2 - 8.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*is(p2,p3)*xn * ( 8.D0*s(p2,p4)*
     &    nDp1*nDp3 + 8.D0*s(p2,p4)*nDp1**2 - 2.D0*s(p2,p4)**2*nDn - 2.D
     &    0*s(p2,p4)*s(p3,p4)*nDn - 8.D0*s(p3,p4)*nDp1*nDp2 )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*is(p3,p4)*s124**(-1)*xn * (  -
     &    2.D0*s(p1,p2)**2*s(p1,p3)*nDn - 16.D0*s(p1,p2)*s(p1,p3)*nDp1*
     &    nDp3 - 16.D0*s(p1,p2)*s(p1,p3)*nDp1**2 + 4.D0*s(p1,p2)*s(p1,
     &    p3)**2*nDn + 4.D0*s(p1,p2)*s(p1,p3)*s(p2,p3)*nDn - 8.D0*s(p1,
     &    p3)**2*nDp1*nDp3 - 8.D0*s(p2,p3)**2*nDp1*nDp3 )
      qqgghn_ba = qqgghn_ba + is(p1,p4)*is(p3,p4)*xn * ( 2.D0*s(p1,p2)*
     &    s(p1,p3)*nDn + 16.D0*s(p1,p3)*nDp1*nDp3 + 16.D0*s(p1,p3)*
     &    nDp1**2 - 4.D0*s(p1,p3)**2*nDn - 4.D0*s(p1,p3)*s(p2,p3)*nDn
     &     - 2.D0*s(p1,p3)*s(p2,p4)*nDn )
      qqgghn_ba = qqgghn_ba + is(p2,p3)*s123**(-2)*xn * ( 2.D0*s(p1,p2)
     &    *s(p1,p4)**2*nDn + 4.D0*s(p1,p2)*s(p1,p4)*s(p2,p4)*nDn + 4.D0
     &    *s(p1,p2)*s(p1,p4)*s(p3,p4)*nDn + 2.D0*s(p1,p2)*s(p2,p4)**2*
     &    nDn + 4.D0*s(p1,p2)*s(p2,p4)*s(p3,p4)*nDn + 2.D0*s(p1,p2)*s(
     &    p3,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p2,p3)*s123**(-1)*s124**(-1)*xn * (
     &     - 2.D0*s(p1,p2)**2*s(p2,p4)*nDn - 2.D0*s(p1,p2)**2*s(p3,p4)*
     &    nDn + 16.D0*s(p1,p2)*s(p2,p4)*nDp3**2 + 16.D0*s(p1,p2)*s(p2,
     &    p4)*nDp2*nDp3 + 4.D0*s(p1,p2)*s(p2,p4)*s(p3,p4)*nDn - 16.D0*
     &    s(p1,p2)*s(p3,p4)*nDp2*nDp3 - 16.D0*s(p1,p2)*s(p3,p4)*nDp2**2
     &     + 4.D0*s(p1,p2)*s(p3,p4)**2*nDn + 32.D0*s(p2,p4)**2*nDp3**2
     &     + 32.D0*s(p2,p4)**2*nDp2*nDp3 + 16.D0*s(p2,p4)**2*nDp1*nDp3
     &     - 32.D0*s(p2,p4)*s(p3,p4)*nDp2*nDp3 - 32.D0*s(p2,p4)*s(p3,p4
     &    )*nDp2**2 + 16.D0*s(p2,p4)*s(p3,p4)*nDp1*nDp3 - 16.D0*s(p2,p4
     &    )*s(p3,p4)*nDp1*nDp2 - 2.D0*s(p2,p4)*s(p3,p4)**2*nDn - 16.D0*
     &    s(p3,p4)**2*nDp1*nDp2 - 2.D0*s(p3,p4)**3*nDn )
      qqgghn_ba = qqgghn_ba + is(p2,p3)*s123**(-1)*xn * ( 16.D0*s(p1,p2
     &    )*nDp1*nDp3 + 8.D0*s(p1,p2)*nDp1*nDp2 + 4.D0*s(p1,p2)*s(p1,p4
     &    )*nDn + 12.D0*s(p1,p2)*s(p2,p4)*nDn + 8.D0*s(p1,p2)*s(p3,p4)*
     &    nDn + 16.D0*s(p1,p4)*nDp2*nDp3 + 16.D0*s(p1,p4)*nDp2**2 - 6.D0
     &    *s(p1,p4)*s(p2,p4)*nDn - 2.D0*s(p1,p4)*s(p3,p4)*nDn - 24.D0*
     &    s(p2,p4)*nDp3**2 - 24.D0*s(p2,p4)*nDp2*nDp3 + 16.D0*s(p2,p4)*
     &    nDp1*nDp3 - 16.D0*s(p2,p4)*nDp1*nDp2 - 6.D0*s(p2,p4)**2*nDn
     &     - 8.D0*s(p2,p4)*s(p3,p4)*nDn + 24.D0*s(p3,p4)*nDp2*nDp3 + 24.
     &    D0*s(p3,p4)*nDp2**2 - 16.D0*s(p3,p4)*nDp1*nDp2 - 4.D0*s(p3,p4
     &    )**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p2,p3)*s124**(-1)*xn * ( 2.D0*s(p1,p2)
     &    *s(p2,p4)*nDn - 4.D0*s(p2,p4)*s(p3,p4)*nDn )
      qqgghn_ba = qqgghn_ba + is(p2,p3)*xn * ( 8.D0*nDp2*nDp3 - 8.D0*
     &    nDp1*nDp2 - 8.D0*s(p2,p4)*nDn )
      qqgghn_ba = qqgghn_ba + is(p2,p3)*is(p3,p4)*s123**(-1)*xn * (  -
     &    8.D0*s(p1,p2)*s(p1,p4)*nDp3**2 - 8.D0*s(p1,p2)*s(p1,p4)*nDp2*
     &    nDp3 + 4.D0*s(p1,p2)*s(p1,p4)*s(p2,p4)*nDn + 8.D0*s(p1,p2)*s(
     &    p2,p4)*nDp1*nDp3 + 4.D0*s(p1,p2)*s(p2,p4)**2*nDn - 2.D0*s(p1,
     &    p4)**2*s(p2,p4)*nDn - 16.D0*s(p1,p4)*s(p2,p4)*nDp3**2 - 16.D0
     &    *s(p1,p4)*s(p2,p4)*nDp2*nDp3 - 4.D0*s(p1,p4)*s(p2,p4)**2*nDn
     &     + 16.D0*s(p2,p4)**2*nDp1*nDp3 - 2.D0*s(p2,p4)**3*nDn )
      qqgghn_ba = qqgghn_ba + is(p2,p3)*is(p3,p4)*xn * ( 8.D0*s(p1,p4)*
     &    nDp2*nDp3 - 4.D0*s(p1,p4)*s(p2,p4)*nDn - 8.D0*s(p2,p4)*
     &    nDp3**2 - 8.D0*s(p2,p4)*nDp1*nDp3 - 4.D0*s(p2,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p3,p4)*s123**(-1)*xn * (  - 24.D0*s(p1
     &    ,p4)*nDp3**2 - 16.D0*s(p1,p4)*nDp2*nDp3 + 4.D0*s(p1,p4)**2*
     &    nDn + 16.D0*s(p1,p4)*s(p2,p4)*nDn - 8.D0*s(p2,p4)*nDp3**2 +
     &    16.D0*s(p2,p4)*nDp1*nDp3 + 12.D0*s(p2,p4)**2*nDn )
      qqgghn_ba = qqgghn_ba + is(p3,p4)*s124**(-1)*xn * (  - 2.D0*s(p1,
     &    p2)*s(p1,p3)*nDn - 2.D0*s(p1,p2)*s(p2,p3)*nDn - 16.D0*s(p1,p3
     &    )*nDp1*nDp3 - 16.D0*s(p1,p3)*nDp1**2 + 4.D0*s(p1,p3)**2*nDn
     &     + 8.D0*s(p1,p3)*s(p2,p3)*nDn - 16.D0*s(p2,p3)*nDp2*nDp3 - 16.
     &    D0*s(p2,p3)*nDp2**2 + 4.D0*s(p2,p3)**2*nDn + 16.D0*s(p2,p4)*
     &    nDp3**2 + 16.D0*s(p2,p4)*nDp2*nDp3 )
      qqgghn_ba = qqgghn_ba + is(p3,p4)*xn * ( 2.D0*s(p2,p3)*nDn )
      qqgghn_ba = qqgghn_ba + (nDp1*is(p1,p4))**2*s124**(-1)*xn * ( 8.D0*s(p1,
     &    p3)**2 + 16.D0*s(p1,p3)*s(p3,p4) + 8.D0*s(p2,
     &    p3)**2 + 8.D0*s(p3,p4)**2 )

      qqgghn_sym=  + s123**(-2)*xn**(-1) * (  - 4.D0*s(p1,p4)**2*nDn -
     &    8.D0*s(p1,p4)*s(p2,p4)*nDn - 8.D0*s(p1,p4)*s(p3,p4)*nDn - 4.D0
     &    *s(p2,p4)**2*nDn - 8.D0*s(p2,p4)*s(p3,p4)*nDn - 4.D0*s(p3,p4)
     &    **2*nDn )
      qqgghn_sym = qqgghn_sym + s123**(-1)*s124**(-1)*xn**(-1) * ( 16.D0
     &    *s(p1,p2)*nDp3**2 - 16.D0*s(p1,p2)*nDp2*nDp3 - 16.D0*s(p1,p2)
     &    *nDp1*nDp3 - 8.D0*s(p1,p2)**2*nDn + 16.D0*s(p1,p2)*s(p3,p4)*
     &    nDn - 32.D0*s(p3,p4)*nDp2**2 - 64.D0*s(p3,p4)*nDp1*nDp2 - 32.D
     &    0*s(p3,p4)*nDp1**2 - 8.D0*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1) * (  - 16.D0*
     &    nDp3**2 - 32.D0*nDp1*nDp2 + 4.D0*s(p1,p2)*nDn - 6.D0*s(p1,p4)
     &    *nDn - 6.D0*s(p2,p4)*nDn - 8.D0*s(p3,p4)*nDn )
      qqgghn_sym = qqgghn_sym + s124**(-2)*xn**(-1) * (  - 4.D0*s(p1,p3
     &    )**2*nDn - 8.D0*s(p1,p3)*s(p2,p3)*nDn - 8.D0*s(p1,p3)*s(p3,p4
     &    )*nDn - 4.D0*s(p2,p3)**2*nDn - 8.D0*s(p2,p3)*s(p3,p4)*nDn - 4.
     &    D0*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + s124**(-1)*xn**(-1) * (  - 16.D0*
     &    nDp3**2 - 16.D0*nDp2*nDp3 - 16.D0*nDp2**2 - 16.D0*nDp1*nDp3
     &     - 32.D0*nDp1*nDp2 - 16.D0*nDp1**2 + 8.D0*s(p1,p2)*nDn - 8.D0
     &    *s(p3,p4)*nDn )
      qqgghn_sym = qqgghn_sym + xn**(-1) * (  - 8.D0*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*s123**(-2)*xn**(-1) * (  - 2.D
     &    0*s(p1,p2)*s(p1,p4)**2*nDn - 4.D0*s(p1,p2)*s(p1,p4)*s(p2,p4)*
     &    nDn - 4.D0*s(p1,p2)*s(p1,p4)*s(p3,p4)*nDn - 2.D0*s(p1,p2)*s(
     &    p2,p4)**2*nDn - 4.D0*s(p1,p2)*s(p2,p4)*s(p3,p4)*nDn - 2.D0*s(
     &    p1,p2)*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*s123**(-1)*s124**(-1)*
     & xn**(-1) * ( 16.D0*s(p1,p2)**2*nDp3**2 + 16.D0*s(p1,p2)**2*nDp1*
     &    nDp3 - 2.D0*s(p1,p2)**3*nDn + 4.D0*s(p1,p2)**2*s(p3,p4)*nDn
     &     + 32.D0*s(p1,p2)*s(p1,p4)*nDp3**2 + 16.D0*s(p1,p2)*s(p1,p4)*
     &    nDp2*nDp3 + 32.D0*s(p1,p2)*s(p1,p4)*nDp1*nDp3 + 16.D0*s(p1,p2
     &    )*s(p3,p4)*nDp2*nDp3 - 32.D0*s(p1,p2)*s(p3,p4)*nDp1*nDp3 - 16.
     &    D0*s(p1,p2)*s(p3,p4)*nDp1*nDp2 - 32.D0*s(p1,p2)*s(p3,p4)*
     &    nDp1**2 - 2.D0*s(p1,p2)*s(p3,p4)**2*nDn + 16.D0*s(p1,p4)**2*
     &    nDp3**2 + 16.D0*s(p1,p4)**2*nDp2*nDp3 + 16.D0*s(p1,p4)**2*
     &    nDp1*nDp3 + 16.D0*s(p1,p4)*s(p3,p4)*nDp2*nDp3 - 16.D0*s(p1,p4
     &    )*s(p3,p4)*nDp1*nDp3 - 16.D0*s(p1,p4)*s(p3,p4)*nDp1*nDp2 - 16.
     &    D0*s(p1,p4)*s(p3,p4)*nDp1**2 - 16.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*s123**(-1)*xn**(-1) * (  - 16.
     &    D0*s(p1,p2)*nDp3**2 - 8.D0*s(p1,p2)*nDp2*nDp3 - 16.D0*s(p1,p2
     &    )*nDp1*nDp3 - 16.D0*s(p1,p2)*nDp1*nDp2 + 2.D0*s(p1,p2)**2*nDn
     &     - 2.D0*s(p1,p2)*s(p1,p4)*nDn - 4.D0*s(p1,p2)*s(p3,p4)*nDn -
     &    8.D0*s(p1,p4)*nDp3**2 - 8.D0*s(p1,p4)*nDp2*nDp3 - 8.D0*s(p1,
     &    p4)*nDp1*nDp3 - 2.D0*s(p1,p4)**2*nDn - 4.D0*s(p1,p4)*s(p2,p4)
     &    *nDn - 4.D0*s(p1,p4)*s(p3,p4)*nDn + 8.D0*s(p2,p4)*nDp3**2 + 8.
     &    D0*s(p2,p4)*nDp1*nDp3 - 2.D0*s(p2,p4)**2*nDn - 4.D0*s(p2,p4)*
     &    s(p3,p4)*nDn - 8.D0*s(p3,p4)*nDp2*nDp3 + 8.D0*s(p3,p4)*nDp1*
     &    nDp3 + 8.D0*s(p3,p4)*nDp1*nDp2 + 8.D0*s(p3,p4)*nDp1**2 - 2.D0
     &    *s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*s124**(-1)*xn**(-1) * ( 2.D0*
     &    s(p1,p2)**2*nDn - 4.D0*s(p1,p2)*s(p3,p4)*nDn + 2.D0*s(p3,p4)
     &    **2*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*xn**(-1) * ( 8.D0*nDp1*nDp3
     &     + 16.D0*nDp1*nDp2 - 2.D0*s(p1,p2)*nDn + 2.D0*s(p1,p4)*nDn +
     &    2.D0*s(p3,p4)*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*is(p1,p4)*s123**(-1)*
     & s124**(-1)*xn**(-1) * (  - 2.D0*s(p1,p2)**3*s(p3,p4)*nDn - 16.D0
     &    *s(p1,p2)**2*s(p3,p4)*nDp1*nDp3 - 16.D0*s(p1,p2)**2*s(p3,p4)*
     &    nDp1**2 + 4.D0*s(p1,p2)**2*s(p3,p4)**2*nDn - 16.D0*s(p1,p2)*
     &    s(p3,p4)**2*nDp1*nDp2 - 2.D0*s(p1,p2)*s(p3,p4)**3*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*is(p1,p4)*s123**(-1)*xn**(-1)
     &  * ( 2.D0*s(p1,p2)**2*s(p3,p4)*nDn + 8.D0*s(p1,p2)*s(p2,p4)*nDp1
     &    *nDp3 + 8.D0*s(p1,p2)*s(p2,p4)*nDp1**2 + 16.D0*s(p1,p2)*s(p3,
     &    p4)*nDp1*nDp3 - 8.D0*s(p1,p2)*s(p3,p4)*nDp1*nDp2 + 16.D0*s(p1
     &    ,p2)*s(p3,p4)*nDp1**2 - 2.D0*s(p1,p2)*s(p3,p4)**2*nDn - 8.D0*
     &    s(p2,p4)*s(p3,p4)*nDp1*nDp3 - 8.D0*s(p2,p4)*s(p3,p4)*nDp1**2
     &     + 8.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*is(p1,p4)*xn**(-1) * (  - 8.D0
     &    *s(p2,p4)*nDp1**2 - 8.D0*s(p3,p4)*nDp1**2 )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*is(p2,p3)*s123**(-1)*xn**(-1)
     &  * ( 8.D0*s(p1,p4)**2*nDp3**2 + 16.D0*s(p1,p4)**2*nDp2*nDp3 + 16.
     &    D0*s(p1,p4)**2*nDp2**2 - 16.D0*s(p1,p4)*s(p2,p4)*nDp2*nDp3 -
     &    16.D0*s(p1,p4)*s(p2,p4)*nDp1*nDp3 - 32.D0*s(p1,p4)*s(p2,p4)*
     &    nDp1*nDp2 + 16.D0*s(p1,p4)*s(p3,p4)*nDp2**2 - 16.D0*s(p1,p4)*
     &    s(p3,p4)*nDp1*nDp3 - 16.D0*s(p1,p4)*s(p3,p4)*nDp1*nDp2 + 8.D0
     &    *s(p2,p4)**2*nDp3**2 + 16.D0*s(p2,p4)**2*nDp1*nDp3 + 16.D0*s(
     &    p2,p4)**2*nDp1**2 - 16.D0*s(p2,p4)*s(p3,p4)*nDp2*nDp3 - 16.D0
     &    *s(p2,p4)*s(p3,p4)*nDp1*nDp2 + 16.D0*s(p2,p4)*s(p3,p4)*
     &    nDp1**2 + 8.D0*s(p3,p4)**2*nDp2**2 + 8.D0*s(p3,p4)**2*nDp1**2
     &     )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*is(p2,p3)*xn**(-1) * ( 4.D0*
     &    s(p1,p4)**2*nDn + 8.D0*s(p1,p4)*s(p2,p4)*nDn + 8.D0*s(p1,p4)*
     &    s(p3,p4)*nDn + 4.D0*s(p2,p4)**2*nDn + 8.D0*s(p2,p4)*s(p3,p4)*
     &    nDn + 4.D0*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*is(p2,p4)*s123**(-1)*xn**(-1)
     &  * ( 8.D0*s(p1,p2)*s(p1,p4)*nDp2**2 - 2.D0*s(p1,p2)*s(p1,p4)**2*
     &    nDn - 4.D0*s(p1,p2)*s(p1,p4)*s(p3,p4)*nDn + 8.D0*s(p1,p2)*s(
     &    p3,p4)*nDp2**2 - 2.D0*s(p1,p2)*s(p3,p4)**2*nDn - 8.D0*s(p1,p4
     &    )**2*nDp2*nDp3 - 8.D0*s(p1,p4)*s(p3,p4)*nDp2*nDp3 + 8.D0*s(p1
     &    ,p4)*s(p3,p4)*nDp1*nDp2 + 8.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_sym = qqgghn_sym + is(p1,p3)*is(p2,p4)*xn**(-1) * (  - 8.D0
     &    *s(p1,p4)*nDp2*nDp3 - 8.D0*s(p1,p4)*nDp2**2 + 2.D0*s(p1,p4)**
     &    2*nDn + 2.D0*s(p1,p4)*s(p3,p4)*nDn + 8.D0*s(p3,p4)*nDp1*nDp2
     &     )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*s123**(-1)*s124**(-1)*
     & xn**(-1) * (  - 16.D0*s(p1,p2)**2*nDp1*nDp3 - 2.D0*s(p1,p2)**3*
     &    nDn + 4.D0*s(p1,p2)**2*s(p3,p4)*nDn - 16.D0*s(p1,p2)*s(p1,p3)
     &    *nDp1*nDp3 - 16.D0*s(p1,p2)*s(p3,p4)*nDp1*nDp3 - 16.D0*s(p1,
     &    p2)*s(p3,p4)*nDp1*nDp2 - 32.D0*s(p1,p2)*s(p3,p4)*nDp1**2 - 2.D
     &    0*s(p1,p2)*s(p3,p4)**2*nDn - 16.D0*s(p1,p3)*s(p3,p4)*nDp1*
     &    nDp2 - 16.D0*s(p1,p3)*s(p3,p4)*nDp1**2 - 16.D0*s(p3,p4)**2*
     &    nDp1*nDp2 )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*s123**(-1)*xn**(-1) * ( 16.D0
     &    *s(p1,p2)*nDp1*nDp3 + 2.D0*s(p1,p2)**2*nDn - 2.D0*s(p1,p2)*s(
     &    p2,p4)*nDn - 4.D0*s(p1,p2)*s(p3,p4)*nDn + 16.D0*s(p1,p3)*nDp1
     &    *nDp3 + 16.D0*s(p2,p4)*nDp1**2 - 2.D0*s(p2,p4)**2*nDn + 16.D0
     &    *s(p3,p4)*nDp1*nDp3 + 32.D0*s(p3,p4)*nDp1**2 + 2.D0*s(p3,p4)
     &    **2*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*s124**(-2)*xn**(-1) * (  - 2.D
     &    0*s(p1,p2)*s(p1,p3)**2*nDn - 4.D0*s(p1,p2)*s(p1,p3)*s(p2,p3)*
     &    nDn - 4.D0*s(p1,p2)*s(p1,p3)*s(p3,p4)*nDn - 2.D0*s(p1,p2)*s(
     &    p2,p3)**2*nDn - 4.D0*s(p1,p2)*s(p2,p3)*s(p3,p4)*nDn - 2.D0*s(
     &    p1,p2)*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*s124**(-1)*xn**(-1) * (  - 16.
     &    D0*s(p1,p2)*nDp1**2 + 16.D0*s(p1,p3)*nDp1*nDp3 + 2.D0*s(p1,p3
     &    )**2*nDn + 4.D0*s(p1,p3)*s(p2,p3)*nDn - 16.D0*s(p2,p3)*nDp1*
     &    nDp3 + 2.D0*s(p2,p3)**2*nDn + 16.D0*s(p3,p4)*nDp1*nDp3 + 16.D0
     &    *s(p3,p4)*nDp1**2 - 2.D0*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*xn**(-1) * ( 16.D0*nDp1**2 )
      qqgghn_sym = qqgghn_sym + is(p1,p4)**2*s124**(-1)*xn**(-1) * (
     &     - 8.D0*s(p1,p3)**2*nDp1**2 - 16.D0*s(p1,p3)*s(p3,p4)*nDp1**2
     &     - 8.D0*s(p2,p3)**2*nDp1**2 - 8.D0*s(p3,p4)**2*nDp1**2 )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*is(p2,p3)*s123**(-1)*xn**(-1)
     &  * ( 8.D0*s(p1,p2)*s(p2,p4)*nDp1**2 - 2.D0*s(p1,p2)*s(p2,p4)**2*
     &    nDn - 4.D0*s(p1,p2)*s(p2,p4)*s(p3,p4)*nDn + 8.D0*s(p1,p2)*s(
     &    p3,p4)*nDp1**2 - 2.D0*s(p1,p2)*s(p3,p4)**2*nDn - 8.D0*s(p2,p4
     &    )**2*nDp1*nDp3 - 8.D0*s(p2,p4)*s(p3,p4)*nDp1*nDp3 + 8.D0*s(p2
     &    ,p4)*s(p3,p4)*nDp1*nDp2 + 8.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*is(p2,p3)*xn**(-1) * (  - 8.D0
     &    *s(p2,p4)*nDp1*nDp3 - 8.D0*s(p2,p4)*nDp1**2 + 2.D0*s(p2,p4)**
     &    2*nDn + 2.D0*s(p2,p4)*s(p3,p4)*nDn + 8.D0*s(p3,p4)*nDp1*nDp2
     &     )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*is(p2,p4)*s124**(-1)*xn**(-1)
     &  * ( 16.D0*s(p1,p3)**2*nDp1*nDp2 + 16.D0*s(p1,p3)*s(p3,p4)*nDp1*
     &    nDp2 + 16.D0*s(p2,p3)**2*nDp1*nDp2 + 16.D0*s(p2,p3)*s(p3,p4)*
     &    nDp1*nDp2 )
      qqgghn_sym = qqgghn_sym + is(p1,p4)*is(p2,p4)*xn**(-1) * ( 4.D0*
     &    s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p3)*s123**(-2)*xn**(-1) * (  - 2.D
     &    0*s(p1,p2)*s(p1,p4)**2*nDn - 4.D0*s(p1,p2)*s(p1,p4)*s(p2,p4)*
     &    nDn - 4.D0*s(p1,p2)*s(p1,p4)*s(p3,p4)*nDn - 2.D0*s(p1,p2)*s(
     &    p2,p4)**2*nDn - 4.D0*s(p1,p2)*s(p2,p4)*s(p3,p4)*nDn - 2.D0*s(
     &    p1,p2)*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p3)*s123**(-1)*s124**(-1)*
     & xn**(-1) * (  - 2.D0*s(p1,p2)**3*nDn + 4.D0*s(p1,p2)**2*s(p3,p4)
     &    *nDn + 16.D0*s(p1,p2)*s(p1,p4)*nDp1*nDp3 - 16.D0*s(p1,p2)*s(
     &    p3,p4)*nDp2*nDp3 - 16.D0*s(p1,p2)*s(p3,p4)*nDp2**2 - 2.D0*s(
     &    p1,p2)*s(p3,p4)**2*nDn + 16.D0*s(p1,p4)**2*nDp3**2 + 16.D0*s(
     &    p1,p4)**2*nDp2*nDp3 + 16.D0*s(p1,p4)**2*nDp1*nDp3 + 16.D0*s(
     &    p1,p4)*s(p3,p4)*nDp2*nDp3 + 16.D0*s(p1,p4)*s(p3,p4)*nDp2**2
     &     - 16.D0*s(p1,p4)*s(p3,p4)*nDp1*nDp3 + 16.D0*s(p1,p4)*s(p3,p4
     &    )*nDp1*nDp2 - 16.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_sym = qqgghn_sym + is(p2,p3)*s123**(-1)*xn**(-1) * (  - 8.D
     &    0*s(p1,p2)*nDp1*nDp3 - 16.D0*s(p1,p2)*nDp1*nDp2 + 2.D0*s(p1,
     &    p2)**2*nDn - 2.D0*s(p1,p2)*s(p2,p4)*nDn - 4.D0*s(p1,p2)*s(p3,
     &    p4)*nDn - 8.D0*s(p1,p4)*nDp3**2 - 8.D0*s(p1,p4)*nDp2*nDp3 -
     &    16.D0*s(p1,p4)*nDp1*nDp3 - 2.D0*s(p1,p4)**2*nDn - 4.D0*s(p1,
     &    p4)*s(p2,p4)*nDn - 4.D0*s(p1,p4)*s(p3,p4)*nDn + 8.D0*s(p2,p4)
     &    *nDp3**2 + 8.D0*s(p2,p4)*nDp2*nDp3 + 8.D0*s(p2,p4)*nDp1*nDp3
     &     - 2.D0*s(p2,p4)**2*nDn - 4.D0*s(p2,p4)*s(p3,p4)*nDn - 8.D0*
     &    s(p3,p4)*nDp2*nDp3 - 8.D0*s(p3,p4)*nDp2**2 + 8.D0*s(p3,p4)*
     &    nDp1*nDp3 - 8.D0*s(p3,p4)*nDp1*nDp2 - 2.D0*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p3)*s124**(-1)*xn**(-1) * ( 2.D0*
     &    s(p1,p2)**2*nDn - 4.D0*s(p1,p2)*s(p3,p4)*nDn + 2.D0*s(p3,p4)
     &    **2*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p3)*xn**(-1) * ( 8.D0*nDp2*nDp3
     &     + 16.D0*nDp1*nDp2 - 2.D0*s(p1,p2)*nDn + 2.D0*s(p2,p4)*nDn +
     &    2.D0*s(p3,p4)*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p3)*is(p2,p4)*s123**(-1)*
     & s124**(-1)*xn**(-1) * (  - 2.D0*s(p1,p2)**3*s(p3,p4)*nDn - 16.D0
     &    *s(p1,p2)**2*s(p3,p4)*nDp2*nDp3 - 16.D0*s(p1,p2)**2*s(p3,p4)*
     &    nDp2**2 + 4.D0*s(p1,p2)**2*s(p3,p4)**2*nDn - 16.D0*s(p1,p2)*
     &    s(p3,p4)**2*nDp1*nDp2 - 2.D0*s(p1,p2)*s(p3,p4)**3*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p3)*is(p2,p4)*s123**(-1)*xn**(-1)
     &  * ( 2.D0*s(p1,p2)**2*s(p3,p4)*nDn + 8.D0*s(p1,p2)*s(p1,p4)*nDp2
     &    *nDp3 + 8.D0*s(p1,p2)*s(p1,p4)*nDp2**2 + 16.D0*s(p1,p2)*s(p3,
     &    p4)*nDp2*nDp3 + 16.D0*s(p1,p2)*s(p3,p4)*nDp2**2 - 8.D0*s(p1,
     &    p2)*s(p3,p4)*nDp1*nDp2 - 2.D0*s(p1,p2)*s(p3,p4)**2*nDn - 8.D0
     &    *s(p1,p4)*s(p3,p4)*nDp2*nDp3 - 8.D0*s(p1,p4)*s(p3,p4)*nDp2**2
     &     + 8.D0*s(p3,p4)**2*nDp1*nDp2 )
      qqgghn_sym = qqgghn_sym + is(p2,p3)*is(p2,p4)*xn**(-1) * (  - 8.D0
     &    *s(p1,p4)*nDp2**2 - 8.D0*s(p3,p4)*nDp2**2 )
      qqgghn_sym = qqgghn_sym + is(p2,p4)*s123**(-1)*s124**(-1)*
     & xn**(-1) * (  - 2.D0*s(p1,p2)**3*nDn + 4.D0*s(p1,p2)**2*s(p3,p4)
     &    *nDn + 16.D0*s(p1,p2)*s(p1,p3)*nDp2*nDp3 - 16.D0*s(p1,p2)*s(
     &    p3,p4)*nDp2*nDp3 - 16.D0*s(p1,p2)*s(p3,p4)*nDp2**2 - 2.D0*s(
     &    p1,p2)*s(p3,p4)**2*nDn + 16.D0*s(p1,p3)*s(p3,p4)*nDp2**2 + 16.
     &    D0*s(p1,p3)*s(p3,p4)*nDp1*nDp2 - 16.D0*s(p3,p4)**2*nDp1*nDp2
     &     )
      qqgghn_sym = qqgghn_sym + is(p2,p4)*s123**(-1)*xn**(-1) * ( 2.D0*
     &    s(p1,p2)**2*nDn - 2.D0*s(p1,p2)*s(p1,p4)*nDn - 4.D0*s(p1,p2)*
     &    s(p3,p4)*nDn - 16.D0*s(p1,p3)*nDp2*nDp3 + 16.D0*s(p1,p4)*
     &    nDp2**2 - 2.D0*s(p1,p4)**2*nDn + 16.D0*s(p3,p4)*nDp2*nDp3 +
     &    32.D0*s(p3,p4)*nDp2**2 + 2.D0*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p4)*s124**(-2)*xn**(-1) * (  - 2.D
     &    0*s(p1,p2)*s(p1,p3)**2*nDn - 4.D0*s(p1,p2)*s(p1,p3)*s(p2,p3)*
     &    nDn - 4.D0*s(p1,p2)*s(p1,p3)*s(p3,p4)*nDn - 2.D0*s(p1,p2)*s(
     &    p2,p3)**2*nDn - 4.D0*s(p1,p2)*s(p2,p3)*s(p3,p4)*nDn - 2.D0*s(
     &    p1,p2)*s(p3,p4)**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p4)*s124**(-1)*xn**(-1) * (  - 16.
     &    D0*s(p1,p2)*nDp2*nDp3 - 16.D0*s(p1,p2)*nDp2**2 - 16.D0*s(p1,
     &    p3)*nDp2*nDp3 + 2.D0*s(p1,p3)**2*nDn + 4.D0*s(p1,p3)*s(p2,p3)
     &    *nDn + 16.D0*s(p2,p3)*nDp2*nDp3 + 2.D0*s(p2,p3)**2*nDn + 16.D0
     &    *s(p3,p4)*nDp2*nDp3 - 16.D0*s(p3,p4)*nDp1*nDp2 - 2.D0*s(p3,p4
     &    )**2*nDn )
      qqgghn_sym = qqgghn_sym + is(p2,p4)*xn**(-1) * ( 16.D0*nDp2*nDp3
     &     + 16.D0*nDp2**2 )
      qqgghn_sym = qqgghn_sym + is(p2,p4)**2*s124**(-1)*xn**(-1) * (
     &     - 8.D0*s(p1,p3)**2*nDp2**2 - 8.D0*s(p2,p3)**2*nDp2**2 - 16.D0
     &    *s(p2,p3)*s(p3,p4)*nDp2**2 - 8.D0*s(p3,p4)**2*nDp2**2 )

      return
      end
