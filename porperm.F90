!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:43:27
 
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

SUBROUTINE porperm(nx,ny,nz)
USE crunchtype
USE params
USE medium
USE flow
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nx
INTEGER(I4B), INTENT(IN)                                        :: ny
INTEGER(I4B), INTENT(IN)                                        :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                    :: jx
INTEGER(I4B)                                                    :: jy
INTEGER(I4B)                                                    :: jz

permxOld = permx
permyOld = permy
permzOld = permz

do jz = 1,nz
DO jy = 1,ny
  DO jx = 1,nx
!!    permx(jx,jy,jz) = permxOld(jx,jy,jz) * ( por(jx,jy,jz)          )**3.0/ ( ( 1.0d0-por(jx,jy,jz) )**2.0 ) * &
!!                                           ( 1.0d0-porOld(jx,jy,jz) )**2.0/ ( (porOld(jx,jy,jz)     )**3.0 ) 
    permx(jx,jy,jz) = permxOld(jx,jy,jz)*( por(jx,jy,jz)/porOld(jx,jy,jz) )**3  * ( (1.0d0-porOld(jx,jy,jz))/(1.0d0-por(jx,jy,jz)) )**2.0d0
    permy(jx,jy,jz) = permyOld(jx,jy,jz) * ( por(jx,jy,jz)          )**3.0/ ( ( 1.0d0-por(jx,jy,jz) )**2.0 ) * &
                                           ( 1.0d0-porOld(jx,jy,jz) )**2.0/ ( (porOld(jx,jy,jz)     )**3.0 )  
    permz(jx,jy,jz) = permzOld(jx,jy,jz) * ( por(jx,jy,jz)          )**3.0/ ( ( 1.0d0-por(jx,jy,jz) )**2.0 ) * &
                                           ( 1.0d0-porOld(jx,jy,jz) )**2.0/ ( (porOld(jx,jy,jz)     )**3.0 )  
    if (jx==12) then
        continue
    end if
!!    permx(jx,jy,jz) = perminx(jx,jy,jz)*( por(jx,jy,jz)/porin(jx,jy,jz) )**3  * ( (1.0d0-porin(jx,jy,jz))/(1.0d0-por(jx,jy,jz)) )**2.0d0
!!    permy(jx,jy,jz) = perminy(jx,jy,jz)*( por(jx,jy,jz)/porin(jx,jy,jz) )**3  * ( (1.0d0-porin(jx,jy,jz))/(1.0d0-por(jx,jy,jz)) )**2.0d0 
!!    permz(jx,jy,jz) = perminz(jx,jy,jz)*( por(jx,jy,jz)/porin(jx,jy,jz) )**3  * ( (1.0d0-porin(jx,jy,jz))/(1.0d0-por(jx,jy,jz)) )**2.0d0 
  END DO
END DO
end do


!!DO jy = 1,ny
!!  DO jx = 1,nx
!!    IF (por(jx,jy,jz) > 0.20) THEN
!!      permx(jx,jy,jz) = perminx(jx,jy,jz)*10000000.0
!!      permy(jx,jy,jz) = perminy(jx,jy,jz)*10000000.0
!!      permz(jx,jy,jz) = perminz(jx,jy,jz)*10000000.0
!!    ELSE
!!      permx(jx,jy,jz) = perminx(jx,jy,jz)
!!      permy(jx,jy,jz) = perminy(jx,jy,jz)
!!      permz(jx,jy,jz) = perminz(jx,jy,jz)
!!    END IF!!
!!  END DO
!!END DO

RETURN
END SUBROUTINE porperm

