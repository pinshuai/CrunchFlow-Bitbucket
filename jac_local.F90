!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:59:43
 
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

SUBROUTINE jac_local(ncomp,nspec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: jx
INTEGER(I4B), INTENT(IN)                             :: jy
INTEGER(I4B), INTENT(IN)                             :: jz

!  Internal variables

REAL(DP)                                             :: sum
REAL(DP)                                             :: spec_conc
REAL(DP)                                             :: mutemp

INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: kk
INTEGER(I4B)                                         :: ksp

fjac_loc = 0.0

DO ksp = 1,nspec
  spec_conc = sp10(ksp+ncomp,jx,jy,jz)
  DO i = 1,ncomp
    IF (muaq(ksp,i) /= 0.0) THEN
      mutemp = muaq(ksp,i)
      DO i2 = 1,i-1
        fjac_loc(i2,i) = fjac_loc(i2,i) +       &
          mutemp*muaq(ksp,i2)*spec_conc
        fjac_loc(i,i2) = fjac_loc(i2,i)
      END DO
      fjac_loc(i,i) = fjac_loc(i,i) +       &
          mutemp*mutemp*spec_conc
    END IF
  END DO
END DO

DO I = 1,ncomp
  fjac_loc(i,i) = fjac_loc(i,i) + sp10(i,jx,jy,jz)
END DO


RETURN
END SUBROUTINE jac_local
!***********************************************************
