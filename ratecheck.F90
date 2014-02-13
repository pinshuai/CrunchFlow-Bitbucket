SUBROUTINE ratecheck(nx,ny,nz,nrct,dt,delt)
USE crunchtype
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                               :: nx
INTEGER(I4B), INTENT(IN)                               :: ny
INTEGER(I4B), INTENT(IN)                               :: nz
INTEGER(I4B), INTENT(IN)                               :: nrct

REAL(DP), INTENT(OUT)                                  :: dt
REAL(DP), INTENT(IN)                                   :: delt

!  Internal variables

INTEGER(I4B)                                           :: jx
INTEGER(I4B)                                           :: jy
INTEGER(I4B)                                           :: jz
INTEGER(I4B)                                           :: k

REAL(DP)                                               :: dtemp1
REAL(DP)                                               :: dtemp2
REAL(DP)                                               :: cutdelt

!  Check the time step against maximum allowed volume fraction change

cutdelt = 0.5*delt
dt     = 10000000.0
dtemp1 = 10000000.0
dtemp2 = 10000000.0

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      DO k = 1,nrct
 
        IF (dppt(k,jx,jy,jz) > 0.0) THEN
          dtemp1 = vpptmax/(volmol(k)*dppt(k,jx,jy,jz))
        ELSE IF (dppt(k,jx,jy,jz) < 0.0) THEN
          dtemp1 = ABS(vdissmax/(volmol(k)*dppt(k,jx,jy,jz)))
        ELSE
          CONTINUE
        END IF

!        IF (dppt(k,jx,jy,jz) < 0.0 .AND. volfx(k,jx,jy,jz) > 0.0001) THEN
!          dtemp2 = ABS(volfx(k,jx,jy,jz)/(volmol(k)*dppt(k,jx,jy,jz)))
!        END IF

!        dtemp2 = MAX(dtemp2,cutdelt)
 
!        dt = MIN(dt,dtemp1,dtemp2)
        dt = MIN(dt,dtemp1)

      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE ratecheck


