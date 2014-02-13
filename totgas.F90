!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:09:35
 
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

SUBROUTINE totgas(ncomp,nspec,ngas,jx,jy,jz)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: ngas
INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz

!  Internal variables and arrays

INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: kk

REAL(DP)                                                    :: sum

DO i = 1,ncomp
  sum=0.0
  DO kk = 1,ngas
    sum = sum + mugas(kk,i)*spgas10(kk,jx,jy,jz)
  END DO
  sgas(i,jx,jy,jz) = sum   ! Assumes gas is not primary species
END DO

RETURN
END SUBROUTINE totgas
!  **************************************************************
