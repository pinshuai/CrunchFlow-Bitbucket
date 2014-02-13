SUBROUTINE AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
USE crunchtype
USE temperature

IMPLICIT NONE

!! External arrays and variables

CHARACTER (LEN=1), INTENT(IN)                       :: Coordinate
INTEGER(I4B), INTENT(IN)                            :: jx
INTEGER(I4B), INTENT(IN)                            :: jy
INTEGER(I4B), INTENT(IN)                            :: jz
REAL(DP), INTENT(OUT)                               :: RoAveRight
REAL(DP), INTENT(OUT)                               :: RoAveLeft

!! Internal arrays and variables

RoAveRight = 0.00
RoAveLeft = 0.00

IF (Coordinate == 'X') THEN
    RoAveRight = (ro(jx+1,jy,jz) + ro(jx,jy,jz))/2.0d0
    RoAveLeft = (ro(jx-1,jy,jz) + ro(jx,jy,jz))/2.0d0
END IF

IF (Coordinate == 'Y') THEN
    RoAveRight = (ro(jx,jy+1,jz) + ro(jx,jy,jz))/2.0d0
    RoAveLeft = (ro(jx,jy-1,jz) + ro(jx,jy,jz))/2.0d0
END IF

IF (Coordinate == 'Z') THEN
    RoAveRight = (ro(jx,jy,jz+1) + ro(jx,jy,jz))/2.0d0
    RoAveLeft = (ro(jx,jy,jz-1) + ro(jx,jy,jz))/2.0d0
END IF

IF (RoAveLeft == 0.0 .OR. RoAveRight == 0.0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Average density = 0'
  WRITE(*,*)
  STOP
END IF

RETURN
END SUBROUTINE AverageRo



