!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function q4ghpmp3(p1,p2,p3,p4,p5,za,zb)

c     This is the reduced matrix element squared
c     for the process
c     q(p1)+qbar(p2) --> H((p5+p6)+Q(p3)+qbar(p4)+g(p5)
c     calculated by the program Hq4g.frm
      include 'types.f'
      complex(dp)::q4ghpmp3

      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5
      real(dp)::s5h,t134,t234,t125
      integer::i1,i2,i3
      complex(dp)::zab
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)

      t134=s(p1,p3)+s(p3,p4)+s(p4,p1)
      t234=s(p2,p3)+s(p3,p4)+s(p4,p2)
      t125=s(p1,p2)+s(p2,p5)+s(p5,p1)
      s5h=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &            +s(p2,p3)+s(p2,p4)
     &                     +s(p3,p4)

c %\cite{DelDuca:2004wt}
c \bibitem{DelDuca:2004wt}
c V.~Del Duca, A.~Frizzo and F.~Maltoni,
c %``Higgs boson production in association with three jets,''
c JHEP {\bf 0405}, 064 (2004)
c [arXiv:hep-ph/0404013].
c %%CITATION = HEP-PH 0404013;%%
c Eq. B15
      q4ghpmp3=
     & -za(p2,p3)*zb(p1,p5)/(s(p3,p4)*t234*s5h)
     & *(-zb(p1,p5)*(zab(p1,p2,p4)+zab(p1,p3,p4))
     & -zb(p4,p5)*t234)
      q4ghpmp3=q4ghpmp3
     & -one/(za(p1,p5)*s(p3,p4)*t125)
     & *(za(p2,p3)*zb(p3,p4)*(za(p3,p1)*zb(p1,p5)+za(p3,p2)*zb(p2,p5))
     & +za(p3,p4)*zb(p4,p5)*(zab(p2,p1,p4)+zab(p2,p5,p4)))
      q4ghpmp3=q4ghpmp3
     & +za(p1,p2)/(za(p1,p5)*za(p2,p5)*s(p3,p4)*t125)
     & *(za(p2,p3)*zb(p3,p4)*(zab(p3,p2,p1)+zab(p3,p5,p1))
     & -za(p3,p4)*zb(p1,p4)*(zab(p2,p1,p4)+zab(p2,p5,p4)))
     & +zb(p1,p4)/(s(p3,p4)*t134*s5h)
     & *(-zab(p2,p1,p5)-zab(p2,p3,p5)-zab(p2,p4,p5))
     & *(zab(p3,p1,p5)+zab(p3,p4,p5))

      return
      end
