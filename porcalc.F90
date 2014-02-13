!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:07
 
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

SUBROUTINE porcalc(nx,ny,nz,nkin,jpor)
USE crunchtype
USE params
USE concentration
USE mineral
USE medium

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nx
INTEGER(I4B), INTENT(IN)                                        :: ny
INTEGER(I4B), INTENT(IN)                                        :: nz
INTEGER(I4B), INTENT(IN)                                        :: nkin
INTEGER(I4B), INTENT(IN)                                        :: jpor

!  Internal variables and arrays

REAL(DP)                                                        :: sum
REAL(DP)                                                        :: vinit
REAL(DP)                                                        :: porfactor

INTEGER(I4B)                                                    :: jx
INTEGER(I4B)                                                    :: jy
INTEGER(I4B)                                                    :: jz
INTEGER(I4B)                                                    :: k

!  NOTE:  Initial volume fractions should be advected as well
!         Assumed immobile for the moment

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      sum = 0.0
      DO k = 1,nkin
        vinit = volin(k,jinit(jx,jy,jz))
        IF (vinit == 0.0) THEN
          area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))*(volfx(k,jx,jy,jz)/0.01)**0.66666
          sum = sum + volfx(k,jx,jy,jz)
        ELSE
          sum = sum + volfx(k,jx,jy,jz)
          area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))* (volfx(k,jx,jy,jz)/vinit)**0.6666
        END IF
      END DO
      IF (jpor /= 1) THEN
        por(jx,jy,jz) = porin(jx,jy,jz)
      ELSE
        por(jx,jy,jz) = 1.0 - sum
        IF (por(jx,jy,jz) <= 0.0) THEN
          por(jx,jy,jz) = 1.d-14
        END IF
      END IF
      porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**0.6666666666666
      DO k = 1,nkin
        area(k,jx,jy,jz) = area(k,jx,jy,jz)*porfactor
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE porcalc
!  *******************************************************
