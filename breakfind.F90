SUBROUTINE BreakFind
USE params

IMPLICIT NONE

CHARACTER (LEN=mls)                     :: dummy1

main:  DO
  READ(18,'(a)') dummy1
  IF (dummy1(1:1) == '+') THEN
    EXIT main
  ELSE IF (dummy1 == 'stop.') THEN
    EXIT main
  ELSE
    CONTINUE
  END IF
END DO main

RETURN
END SUBROUTINE BreakFind
