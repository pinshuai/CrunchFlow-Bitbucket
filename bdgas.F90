!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:53:05
 
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

SUBROUTINE bdgas(ncomp,nspec,nrct,ngas,nbnd,scorr)
USE crunchtype
USE concentration
USE params

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec
INTEGER(I4B), INTENT(IN)                      :: nrct
INTEGER(I4B), INTENT(IN)                      :: ngas
INTEGER(I4B), INTENT(IN)                      :: nbnd
REAL(DP), DIMENSION(:), INTENT(OUT)           :: scorr

! Internal arrays and variables

REAL(DP)                                      :: sum

INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: kk

!  This routine calculates the gas concentrations at the
!  ghost zone outside the domain.

DO i = 1,ncomp
  sum = 0.0
  DO kk = 1,ngas
    sum = sum + mugas(kk,i)*spbgas(kk,nbnd)
  END DO
  scorr(i) = sum     ! Assumes gas is not a basis species
END DO

RETURN
END SUBROUTINE bdgas
!**************************************************************
