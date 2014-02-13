!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:52:52
 
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

SUBROUTINE bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
USE crunchtype
USE params
USE concentration
USE transport

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec
INTEGER(I4B), INTENT(IN)                      :: nbnd
REAL(DP), INTENT(OUT)                         :: sumchgbd

! Internal arrays and variables

REAL(DP)                                      :: sum1
REAL(DP)                                      :: sum2
REAL(DP)                                      :: sum3
REAL(DP)                                      :: correct
REAL(DP)                                      :: tmp

INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: ksp


IF (idiffus == 0) THEN
  d_25 = dzero
ELSE
  d_25 = dcoeff
END IF

sum3 = 0.0
DO i = 1,ncomp
  sum1 = 0.0
  sum2 = 0.0
  DO ksp = 1,nspec
    correct = d_sp(ksp+ncomp)/d_25
    tmp = correct*muaq(ksp,i)*spb(ksp+ncomp,nbnd)
    sum1 = sum1 + tmp
    sum2 = sum2 + tmp*chg(ksp+ncomp)
  END DO
  correct = d_sp(i)/d_25
  tmp = correct*spb(i,nbnd)
  sdsp(i) = sum1 + tmp
  schg(i) = sum2 + tmp*chg(i)
  sum3 = sum3 + chg(i)*schg(i)
END DO

sumchgbd = sum3

RETURN
END SUBROUTINE bd_diffuse
!**************************************************************
