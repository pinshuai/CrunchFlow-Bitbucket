SUBROUTINE BreakFindPlusCount(iunit5,ncount)
USE params

IMPLICIT NONE

!!  External arrays and variables
INTEGER, INTENT(IN)                     :: iunit5
INTEGER, INTENT(OUT)                    :: ncount

!!  Internal arrays and variables
CHARACTER (LEN=mls)                     :: dummy1

main:  DO
  READ(iunit5,'(a)') dummy1
  ncount = ncount + 1
  IF (dummy1(1:1) == '+') THEN
    ncount = ncount - 1
    EXIT main
  ELSE IF (dummy1 == 'stop.') THEN
    WRITE(*,*)
    WRITE(*,*) ' Looking for "+" in first character position as marker for end of list'
    WRITE(*,*)
    READ(*,*)   
    STOP
  ELSE
    CONTINUE
  END IF
END DO main

RETURN
END SUBROUTINE BreakFindPlusCount
