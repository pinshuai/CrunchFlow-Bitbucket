SUBROUTINE find_string(iunit,dumstring,ifind)
USE crunchtype
USE params
 
IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                      :: iunit
CHARACTER (LEN=mls), INTENT(IN)               :: dumstring
INTEGER(I4B), INTENT(OUT)                     :: ifind

!  Internal variables

CHARACTER (LEN=mls)                           :: dummy1

ifind = 0

main:  DO
  READ(iunit,'(a)') dummy1
  IF (dummy1 == dumstring) THEN
    ifind = 1
    EXIT main
  ELSE IF (dummy1 == 'stop.') THEN
    WRITE(*,*) 
    WRITE(*,*) ' String not found: ', dumstring
    WRITE(*,*)
    READ(*,*)
    STOP
  ELSE
    CONTINUE
  END IF
END DO main

RETURN
END SUBROUTINE find_string

