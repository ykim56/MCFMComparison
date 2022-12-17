!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine h4qn(p1,p2,p3,p4,ampsq)
      implicit none
      include 'types.f'

c     Taken from Kauffman,Desai,Risal
c     PRD 55 1997 (4009)
c     q(-p1)+qp(-p2)--> h --> q(p3)+qp(p4)
c     returns overall matrix element squared
c     summed over colors and spins with factor of g^4*A^2 removed
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      real(dp):: ampsq,ammsq
      integer:: p1,p2,p3,p4
c      include 'zprods_com.f'
c      integer:: j1,j2
c      complex(dp):: a(2,2),b(2,2),amm
c  Eq 29 (suitably modified for my momentum configuration)
c====statement function
c--left-left amplitude
c      amm(p1,p2,p3,p4)=
c     & +za(p3,p4)**2/(za(p1,p3)*za(p2,p4))
c     & +zb(p1,p2)**2/(zb(p1,p3)*zb(p2,p4))
c--The above amplitude squared + color factor
      ammsq(p1,p2,p3,p4)=V/4._dp
     & *((s(p1,p2)-s(p3,p4))**2+
     & (s(p1,p3)*s(p2,p4)+s(p3,p4)*s(p1,p2)-s(p1,p4)*s(p2,p3))**2
     &  /(s(p1,p3)*s(p2,p4)))/(s(p1,p3)*s(p2,p4))
      ampsq=ammsq(p1,p2,p3,p4)+ammsq(p3,p2,p1,p4)
     &     +ammsq(p1,p4,p3,p2)+ammsq(p3,p4,p1,p2)

      return
      end

      subroutine h4qi(p1,p2,p3,p4,ampsqid,ampsq_a,ampsq_b,ampsq_i)
      implicit none
      include 'types.f'

c     Taken from Kauffman,Desai,Risal
c     PRD 55 1997 (4009)
c     q(-p1)+qp(-p2)--> h --> q(p3)+qp(p4)
c     returns overall matrix element squared
c     summed over colors and spins with factor of g^4*A^2 removed
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      real(dp):: ampsqid,ammsq,ammsqi,ampsq_a,ampsq_b,ampsq_i
      integer:: p1,p2,p3,p4
c      include 'zprods_com.f'
c      integer:: j1,j2
c      complex(dp):: a(2,2),b(2,2),amm
c  Eq 29 (suitably modified for my momentum configuration)
c====statement function
c--left-left amplitude
c      amm(p1,p2,p3,p4)=
c     & +za(p3,p4)**2/(za(p1,p3)*za(p2,p4))
c     & +zb(p1,p2)**2/(zb(p1,p3)*zb(p2,p4))
c--The above amplitude squared + color factor
      ammsq(p1,p2,p3,p4)=V/4._dp
     & *((s(p1,p2)-s(p3,p4))**2+
     & (s(p1,p3)*s(p2,p4)+s(p3,p4)*s(p1,p2)-s(p1,p4)*s(p2,p3))**2
     &  /(s(p1,p3)*s(p2,p4)))/(s(p1,p3)*s(p2,p4))
c--The interference with color factor Eq.(A20)
c-- Note: sign of this term changed by JMC on 8/8/05 to agree with
c-- the other routine in Hqaqasq.f; one minus sign comes from the
c-- colour factor, the other from a single fermion loop
      ammsqi(p1,p2,p3,p4)=+cf/2._dp
     & *((s(p1,p2)-s(p3,p4))**2
     & *(s(p1,p3)*s(p2,p4)+s(p1,p4)*s(p2,p3)-s(p1,p2)*s(p3,p4))
     & -2._dp*(s(p1,p2)*s(p3,p4)+s(p2,p3)*s(p1,p4)-s(p1,p3)*s(p2,p4))
     &     *(s(p3,p4)*s(p1,p2)+s(p1,p3)*s(p2,p4)-s(p1,p4)*s(p2,p3))
     &  )/(s(p1,p3)*s(p2,p3)*s(p1,p4)*s(p2,p4))

c--- calculate the matrix element as the sum of two pieces,
c--- related by the interchange of p3 and p4, plus an interference
c--- term that is colour-suppressed

      ampsq_a=ammsq(p1,p2,p3,p4)+ammsq(p3,p2,p1,p4)
     &       +ammsq(p1,p4,p3,p2)+ammsq(p3,p4,p1,p2)
      ampsq_b=ammsq(p1,p2,p4,p3)+ammsq(p4,p2,p1,p3)
     &       +ammsq(p1,p3,p4,p2)+ammsq(p4,p3,p1,p2)
      ampsq_i=ammsqi(p1,p2,p3,p4)+ammsqi(p3,p4,p1,p2)

      ampsqid=ampsq_a+ampsq_b+ampsq_i

      return
      end
