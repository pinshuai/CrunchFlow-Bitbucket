SUBROUTINE FindMaxFlow(nx,ny,nz)

USE crunchtype
USE params
USE flow
USE transport

IMPLICIT NONE

!!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                :: nx
INTEGER(I4B), INTENT(IN)                                :: ny
INTEGER(I4B), INTENT(IN)                                :: nz

!! Internal variables and arrays

maxQx = MAXVAL(DABS(qx(0:nx,1:ny,1:nz)))
IF (maxQx > 0.0d0) THEN
  xflow = .TRUE.
ELSE
  xflow = .FALSE.
END IF

maxQy = MAXVAL(DABS(qy(1:nx,0:ny,1:nz)))
IF (maxQy > 0.0d0) THEN
  yflow = .TRUE.
ELSE
  yflow = .FALSE.
END IF

maxQz = MAXVAL(DABS(qz(1:nx,1:ny,0:nz)))
IF (maxQz > 0.0d0) THEN
  zflow = .TRUE.
ELSE
  zflow = .FALSE.
END IF

RETURN
END SUBROUTINE FindMaxFlow

