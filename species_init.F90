!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:08:24
 
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

SUBROUTINE species_init(ncomp,nspec)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec

!  Internal variables and arrays

INTEGER(I4B)                                                :: ksp
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: nk

REAL(DP)                                                    :: sum

DO ksp = 1,nspec
  sum = 0.0D0
  DO i = 1,ncomp
    IF (ulab(i) == 'H2O') THEN
      sum = sum + muaq(ksp,i)*(gamtmp(i))
    ELSE
      sum = sum + muaq(ksp,i)*(sptmp(i) + gamtmp(i))
    END IF
  END DO
  nk = ksp + ncomp
  sptmp(nk) = keqaq_tmp(ksp) - gamtmp(nk) + sum
  sptmp10(nk) = EXP(sptmp(nk))
END DO

RETURN
END SUBROUTINE species_init
!  **************************************************************
