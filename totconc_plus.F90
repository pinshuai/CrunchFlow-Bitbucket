!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:09:21
 
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

SUBROUTINE totconc_plus(ncomp,nspec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE transport

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz

!  Internal variables and arrays

INTEGER(I4B)                                                :: ksp
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: kk

REAL(DP)                                                    :: sum1
REAL(DP)                                                    :: sum2
REAL(DP)                                                    :: sum3
REAL(DP)                                                    :: sum4
REAL(DP)                                                    :: tmp

sum4 = 0.0
DO i = 1,ncomp
  sum1 = 0.0
  sum2 = 0.0
  sum3 = 0.0
  DO ksp = 1,nspec
    kk = ksp + ncomp
    tmp = muaq(ksp,i)*sp10(kk,jx,jy,jz)
    sum1 = sum1 + tmp
    sum2 = sum2 + tmp*d_correct(kk)
    sum3 = sum3 + tmp*d_correct(kk)*chg(kk)
  END DO
  tmp = sp10(i,jx,jy,jz)
  s(i,jx,jy,jz)     = sum1 + tmp
  s_dsp(i,jx,jy,jz) = sum2 + tmp*d_correct(i)
  s_chg(i,jx,jy,jz) = sum3 + tmp*d_correct(i)*chg(i)
  sum4 = sum4 + s_chg(i,jx,jy,jz)*chg(i)
END DO

sumwtchg(jx,jy,jz) = sum4

RETURN
END SUBROUTINE totconc_plus
!  **************************************************************
