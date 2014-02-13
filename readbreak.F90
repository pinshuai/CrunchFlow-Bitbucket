SUBROUTINE readbreak
USE crunchtype
USE params
 
IMPLICIT NONE

CHARACTER (LEN=mls)                                        :: dummy1

100 READ(18,'(a)',END=300) dummy1

IF (dummy1(1:1) == '+') THEN
  RETURN
END IF

GO TO 100

300 RETURN
END SUBROUTINE readbreak
