!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:10:20
 
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

SUBROUTINE velocalc(nx,ny,nz)
USE crunchtype
USE params
USE medium
USE transport
USE temperature, ONLY: ro
USE flow
USE CrunchFunctions

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: nx
INTEGER(I4B), INTENT(IN)                                       :: ny
INTEGER(I4B), INTENT(IN)                                       :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                          :: jx
INTEGER(I4B)                                                          :: jy
INTEGER(I4B)                                                          :: jz

REAL(DP)                                                              :: vv

REAL(DP)                                                              :: RoAveLeft
REAL(DP)                                                              :: RoAveRight

!  ****** PARAMETERS  ****************************

REAL(DP), PARAMETER                                                   :: visc=0.001d0
REAL(DP), PARAMETER                                                   :: ct=9.135E-10
REAL(DP), PARAMETER                                                   :: big=1.0d0
REAL(DP), PARAMETER                                                   :: zero=0.0d0
REAL(DP), PARAMETER                                                   :: grav=9.806d0
REAL(DP)                                                              :: term1
REAL(DP)                                                              :: term2

CHARACTER (LEN=1)                                                     :: Coordinate

vv = secyr/visc  ! Convert to m/yr

!   calculate darcy fluxes

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      IF (jx /= nx) THEN
        Coordinate = 'X'  
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft) 
        qx(jx,jy,jz)= -2.0d0*vv*(pres(jx+1,jy,jz)-pres(jx,jy,jz))*harx(jx,jy,jz) /(dxx(jx)+dxx(jx+1))         &
                     + SignGravity*vv*harx(jx,jy,jz)*RoAveRight*grav*COSD(x_angle)
      END IF
      IF (jy /= ny) THEN
        Coordinate = 'Y'  
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft) 
        term1 = -2.0d0*vv*(pres(jx,jy+1,jz)-pres(jx,jy,jz))*hary(jx,jy,jz) /(dyy(jy)+dyy(jy+1))
        term2 = SignGravity*vv*hary(jx,jy,jz)*RoAveRight*grav*COSD(y_angle)
        qy(jx,jy,jz)= -2.0d0*vv*(pres(jx,jy+1,jz)-pres(jx,jy,jz))*hary(jx,jy,jz) /(dyy(jy)+dyy(jy+1))     &
                     + SignGravity*vv*hary(jx,jy,jz)*RoAveRight*grav*COSD(y_angle)
      END IF
      IF (jz /= nz) THEN
        Coordinate = 'Z'  
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
        qz(jx,jy,jz)= -2.0d0*vv*(pres(jx,jy,jz+1)-pres(jx,jy,jz))*harz(jx,jy,jz) /(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))     &
                     + SignGravity*vv*harz(jx,jy,jz)*RoAveRight*grav*COSD(z_angle)
      END IF
    
    END DO
  END DO
END DO

DO jz = 1,nz
  DO jx = 1,nx
    IF (activecellPressure(jx,0,jz) == 0) THEN
      term1 = -2.0d0*vv*(pres(jx,1,jz)-pres(jx,0,jz))*hary(jx,1,jz) /(dyy(1))
      term2 = SignGravity*vv*hary(jx,1,jz)*ro(jx,1,jz)*grav*COSD(y_angle)
      qy(jx,0,jz)= term1 + term2    
    ELSE 
      qy(jx,0,jz) = 0.0d0
    END IF
    IF (activecellPressure(jx,ny+1,jz) == 0) THEN
      term1 = -2.0d0*vv*(pres(jx,ny+1,jz)-pres(jx,ny,jz))*hary(jx,ny,jz) /(dyy(ny))
      term2 = SignGravity*vv*hary(jx,ny,jz)*ro(jx,ny,jz)*grav*COSD(y_angle)
      qy(jx,ny,jz)= term1 + term2    
    ELSE 
      qy(jx,ny,jz) = 0.0d0
    END IF
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    IF (activecellPressure(0,jy,jz) == 0) THEN
      term1 = -2.0d0*vv*(pres(1,jy,jz)-pres(0,jy,jz))*harx(1,jy,jz) /(dxx(1))
      term2 = SignGravity*vv*harx(1,jy,jz)*ro(1,jy,jz)*grav*COSD(x_angle)
      qx(0,jy,jz)= term1 + term2    
    ELSE 
      qx(0,jy,jz) = 0.0d0
    END IF
    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
      term1 = -2.0d0*vv*(pres(nx+1,jy,jz)-pres(nx,jy,jz))*harx(nx,jy,jz) /(dxx(nx))
      term2 = SignGravity*vv*harx(nx,jy,jz)*ro(nx,jy,jz)*grav*COSD(x_angle)
      qx(nx,jy,jz)= term1 + term2    
    ELSE 
      qx(nx,jy,jz) = 0.0d0
    END IF
  END DO
END DO

DO jy = 1,ny
  DO jx = 1,nx
    IF (activecellPressure(jx,jy,0) == 0) THEN
      term1 = -2.0d0*vv*(pres(jx,jy,1)-pres(jx,jy,0))*harz(jx,jy,1) /(dzz(jx,jy,1))
      term2 = SignGravity*vv*harz(jx,jy,1)*ro(jx,jy,1)*grav*COSD(z_angle)
      qz(jx,jy,0)= term1 + term2    
    ELSE 
      qz(jx,jy,0) = 0.0d0
    END IF
    IF (activecellPressure(jx,jy,nz+1) == 0) THEN
      term1 = -2.0d0*vv*(pres(jx,jy,nz+1)-pres(jx,jy,nz))*harz(jx,jy,nz) /(dzz(jx,jy,nz))
      term2 = SignGravity*vv*harz(jx,jy,nz)*ro(jx,jy,nz)*grav*COSD(z_angle)
      qz(jx,jy,nz)= term1 + term2    
    ELSE 
      qz(jx,jy,nz) = 0.0d0
    END IF
  END DO
END DO

RETURN
END SUBROUTINE velocalc
