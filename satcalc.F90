!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:07:50
 
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

SUBROUTINE satcalc(ncomp,nrct,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE mineral
USE medium
USE temperature


IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                        :: ncomp
INTEGER(I4B), INTENT(IN)                                        :: nrct
INTEGER(I4B), INTENT(IN)                                        :: jx
INTEGER(I4B), INTENT(IN)                                        :: jy
INTEGER(I4B), INTENT(IN)                                        :: jz

!  Internal variables and arrays

REAL(DP)                                                        :: tk
REAL(DP)                                                        :: tkinv
REAL(DP)                                                        :: reft
REAL(DP)                                                        :: sum
REAL(DP)                                                        :: sumiap 

INTEGER(I4B)                                                    :: k
INTEGER(I4B)                                                    :: np
INTEGER(I4B)                                                    :: i
INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: kd
INTEGER(I4B)                                                    :: isotope

tk = t(jx,jy,jz) + 273.15D0
tkinv = 1.0/tk
reft = 1.0/298.15

DO k = 1,nrct
  DO np = 1,nreactmin(k)

    IF (ikh2o /= 0) THEN

      sumiap = 0.0D0
      DO i = 1,ncomp
        IF (ulab(i) == 'H2O') THEN
          sumiap = sumiap + decay_correct(i,k)*mumin(np,k,i)*(gam(i,jx,jy,jz))
        ELSE
          sumiap = sumiap + decay_correct(i,k)*mumin(np,k,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
        END IF
      END DO

    ELSE

      sumiap = 0.0D0
      DO i = 1,ncomp
        sumiap = sumiap + decay_correct(i,k)*mumin(np,k,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
      END DO

    END IF

    silog(np,k) = (sumiap - keqmin(np,k,jx,jy,jz))/clg
  END DO
END DO

RETURN
END SUBROUTINE satcalc
!********************************************************
