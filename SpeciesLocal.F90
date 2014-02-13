!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:08:21
 
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

SUBROUTINE SpeciesLocal(ncomp,nspec,jx,jy,jz)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz

!  Internal variables and arrays

INTEGER(I4B)                                                :: ksp
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: nk

REAL(DP)                                                    :: sum

DO ksp = 1,nspec

  IF (ikh2o /= 0) THEN

    sum = 0.0D0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
        sum = sum + muaq(ksp,i)*(gam(i,jx,jy,jz))
      ELSE
        sum = sum + muaq(ksp,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
      END IF
    END DO

  ELSE

    sum = 0.0D0
    DO i = 1,ncomp
      sum = sum + muaq(ksp,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
    END DO

  END IF


  nk = ncomp + ksp
  sp(nk,jx,jy,jz) = keqaq(ksp,jx,jy,jz) - gam(nk,jx,jy,jz) + sum
  sp10(nk,jx,jy,jz) = EXP(sp(nk,jx,jy,jz))
END DO

RETURN
END SUBROUTINE SpeciesLocal
!  **************************************************************
