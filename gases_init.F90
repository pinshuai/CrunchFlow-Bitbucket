!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:58:09
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE gases_init(ncomp,ngas,tempc)
USE crunchtype
USE params
USE concentration
USE temperature

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: ngas
REAL(DP), INTENT(IN)                                       :: tempc

!  Internal variables

REAL(DP)                                                   :: tempk
REAL(DP)                                                   :: denmol
REAL(DP)                                                   :: sum

INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: kk


tempk = tempc + 273.15
denmol = DLOG(1.e05/(8.314*tempk))           ! P/RT = n/V, with pressure converted from bars to Pascals

DO kk = 1,ngas

  IF (ikh2o /= 0) THEN

    sum = 0.0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
        sum = sum + mugas(kk,i)*(gamtmp(i))
      ELSE
        sum = sum + mugas(kk,i)*(sptmp(i) + gamtmp(i))
      END IF
    END DO

  ELSE

    sum = 0.0
    DO i = 1,ncomp
      sum = sum + mugas(kk,i)*(sptmp(i) + gamtmp(i))
    END DO

  END IF
  spgastmp(kk) = keqgas_tmp(kk) + sum 
  spgastmp10(kk) = DEXP(spgastmp(kk))       !  mole fraction of gases
END DO

RETURN
END SUBROUTINE gases_init
!  **************************************************************
