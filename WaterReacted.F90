!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:11:19
 
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

SUBROUTINE WaterReacted(nx,ny,nz)
USE crunchtype
USE runtime
USE params
USE concentration
USE temperature

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nx
INTEGER(I4B), INTENT(IN)                                        :: ny
INTEGER(I4B), INTENT(IN)                                        :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                    :: jx
INTEGER(I4B)                                                    :: jy
INTEGER(I4B)                                                    :: jz
INTEGER(I4B)                                                    :: ik

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx    
      H2Oreacted(jx,jy,jz) = s(ikh2o,jx,jy,jz)/sn(ikh2o,jx,jy,jz)
      H2Oreacted(jx,jy,jz) = 1.0d0
    END DO
  END DO
END DO

RETURN
END SUBROUTINE WaterReacted




