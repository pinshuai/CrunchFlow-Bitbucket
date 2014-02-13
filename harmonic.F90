!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:58:43
 
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

SUBROUTINE harmonic(nx,ny,nz)
USE crunchtype
USE params
USE medium
USE flow

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz

!  Internal variables

INTEGER(I4B)                                       :: jx
INTEGER(I4B)                                       :: jy
INTEGER(I4B)                                       :: jz


!   generate average (harmonic mean) conductivities at grid cell faces

DO jz = 1,nz 
  DO jy = 1,ny
    DO jx = 0,nx
      IF (permx(jx,jy,jz) == 0.0 .AND. permx(jx+1,jy,jz) == 0.0) THEN
        harx(jx,jy,jz) = 0.0
      ELSE
        harx(jx,jy,jz) = permx(jx,jy,jz)*permx(jx+1,jy,jz) *(dxx(jx)+dxx(jx+1))/(dxx(jx)*permx(jx+1,jy,jz) + dxx(jx+1)*permx(jx,jy,jz))
        continue
      END IF
    END DO
  END DO
END DO

DO jz = 1,nz 
  DO jy = 0,ny
    DO jx = 1,nx
      IF (permy(jx,jy,jz) == 0.0 .AND. permy(jx,jy+1,jz) == 0.0) THEN
        hary(jx,jy,jz) = 0.0
      ELSE
        hary(jx,jy,jz) = permy(jx,jy,jz)*permy(jx,jy+1,jz) *(dyy(jy)+dyy(jy+1))  &
          /(dyy(jy)*permy(jx,jy+1,jz) + dyy(jy+1)*permy(jx,jy,jz))
      END IF
    END DO
  END DO
END DO

DO jz = 0,nz 
  DO jy = 1,ny
    DO jx = 1,nx
      IF (permz(jx,jy,jz) == 0.0 .AND. permz(jx,jy,jz+1) == 0.0) THEN
        harz(jx,jy,jz) = 0.0
      ELSE
        harz(jx,jy,jz) = permz(jx,jy,jz)*permz(jx,jy,jz+1) *(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))  &
          /(dzz(jx,jy,jz)*permz(jx,jy,jz+1) + dzz(jx,jy,jz+1)*permz(jx,jy,jz))
      END IF
    END DO
  END DO
END DO

RETURN
END SUBROUTINE harmonic


