!******************        GIMRT98      ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:00:15
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*********************** C.I. Steefel *******************
!                      All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE keqcalc2_init(ncomp,nrct,nspec,ngas,nsurf_sec,tempc)
USE crunchtype
USE params
USE concentration
USE mineral
USE temperature

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nrct
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: ngas
INTEGER(I4B), INTENT(IN)                             :: nsurf_sec

REAL(DP), INTENT(IN)                                 :: tempc

!  Internal variables

REAL(DP)                                             :: temp
REAL(DP)                                             :: temp2
REAL(DP)                                             :: x1
REAL(DP)                                             :: x2
REAL(DP)                                             :: x3
REAL(DP)                                             :: x4
REAL(DP)                                             :: x5

INTEGER(I4B)                                         :: ksp
INTEGER(I4B)                                         :: kg
INTEGER(I4B)                                         :: k
INTEGER(I4B)                                         :: msub
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: np

temp = tempc + 273.15
temp2 = temp*temp

DO ksp = 1,nspec
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqaq_tmp(ksp) = -clg*eqhom(ksp)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqaq_tmp(ksp) = -clg*(x1*DLOG(temp) + x2 +  &
        x3*temp + x4/temp + x5/(temp2))
  END IF
END DO

DO kg = 1,ngas
  ksp = kg + nspec
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqgas_tmp(kg) = -clg*eqgas(kg)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqgas_tmp(kg) = -clg*(x1*DLOG(temp) + x2 +  &
        x3*temp + x4/temp + x5/(temp2))
  END IF
END DO

msub = 0
DO k = 1,nrct
  DO np = 1,nreactmin(k)
    msub = msub + 1
    ksp = msub + ngas + nspec
    IF (ntemp == 1 .OR. RunIsothermal) THEN
      keqmin_tmp(np,k) = clg*alnk(msub)
    ELSE
      x1 = as1(ksp,1)
      x2 = as1(ksp,2)
      x3 = as1(ksp,3)
      x4 = as1(ksp,4)
      x5 = as1(ksp,5)
      keqmin_tmp(np,k) = clg*(x1*DLOG(temp) + x2 +  &
          x3*temp + x4/temp + x5/(temp2))
    END IF
  END DO
END DO

DO ns = 1,nsurf_sec
  ksp = msub + ngas + nspec + ns
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqsurf_tmp(ns) = -clg*eqsurf(ns)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqsurf_tmp(ns) = -clg*(x1*DLOG(temp) + x2 +  &
        x3*temp + x4/temp + x5/(temp2))
  END IF
END DO


RETURN
END SUBROUTINE keqcalc2_init
!*********************************************************************
