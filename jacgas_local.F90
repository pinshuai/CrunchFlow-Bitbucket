!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:59:38
 
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

SUBROUTINE jacgas_local(ncomp,ngas,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE
!fp! auto_par_loops = 0;
!  External variables

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: ngas
INTEGER(I4B), INTENT(IN)                           :: jx
INTEGER(I4B), INTENT(IN)                           :: jy
INTEGER(I4B), INTENT(IN)                           :: jz

!  Internal variables

INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: i2
INTEGER(I4B)                                       :: kk

REAL(DP)                                           :: mutemp
REAL(DP)                                           :: gastemp

fgas_local = 0.0d0

DO kk = 1,ngas
  gastemp = spgas10(kk,jx,jy,jz)
  DO i = 1,ncomp
    IF (mugas(kk,i) /= 0.0) THEN
      mutemp = mugas(kk,i)
      DO i2 = 1,ncomp
        fgas_local(i2,i) = fgas_local(i2,i) + mutemp*mugas(kk,i2)*gastemp
      END DO
    END IF
  END DO
END DO

RETURN
END SUBROUTINE jacgas_local
!***********************************************************
