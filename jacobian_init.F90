!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:59:47
 
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

SUBROUTINE jacobian_init(ncomp,nspec)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nspec

!  Internal variables

REAL(DP)                                             :: sum

INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: ksp

DO i = 1,ncomp
  DO i2 = 1,ncomp
    sum = 0.0
    DO ksp = 1,nspec
      sum = sum + muaq(ksp,i)*muaq(ksp,i2)*sptmp10(ksp+ncomp)
    END DO
    dpsi(i,i2) = sum
  END DO
  dpsi(i,i) = dpsi(i,i) + sptmp10(i)
END DO

RETURN
END SUBROUTINE jacobian_init
!***********************************************************
