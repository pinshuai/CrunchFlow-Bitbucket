!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:59:07
 
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

SUBROUTINE GasInjection(ncomp,nspec,ngas,scorr,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE flow

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nspec
INTEGER(I4B), INTENT(IN)                           :: ngas
INTEGER(I4B), INTENT(IN)                           :: jz
INTEGER(I4B), INTENT(IN)                           :: jx
INTEGER(I4B), INTENT(IN)                           :: jy

REAL(DP), DIMENSION(:), INTENT(OUT)                :: scorr

!  Internal variables and arrays

INTEGER(I4B)                                       :: nco
INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: kk

REAL(DP)                                           :: sum

nco = intbndgas(jx,jy,jz)

DO  i = 1,ncomp
  sum=0.0
  DO kk = 1,ngas
    sum = sum + mugas(kk,i)*spcondgas10(kk,nco)
  END DO
  scorr(i) = sum     !! NOTE:  Assumes gas is not a basis species
END DO

RETURN
END SUBROUTINE GasInjection
!**************************************************************
