!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      SUBROUTINE h5g(ANS)

c RETURNS AMPLITUDE SQUARED SUMMED OVER COLORS
c AND HELICITIES FOR THE POINT IN PHASE SPACE P(mxpart,4)
c FOR PROCESS : g g -> g g g h

c overall coupling factors have been removed


      include 'types.f'
      include 'constants.f'

c CONSTANTS

      integer,parameter:: NEXTERNAL=5, NCOMB= 16


c ARGUMENTS

      real(dp)::ANS

c LOCAL VARIABLES

      integer::NHEL(NEXTERNAL,NCOMB)
      real(dp)::T,GG_GGG
      integer:: IHEL
      include 'hels.f'

c ----------
c BEGIN CODE
c ----------


      ANS=zip

c--  sum over helicities
      DO IHEL=1,16
              T=GG_GGG(NHEL(1,IHEL))
              ANS=ANS+T
      ENDDO

c--  Multiply by two to account for the other 16 helicity configs
      ANS=ANS*two

      RETURN
      END



