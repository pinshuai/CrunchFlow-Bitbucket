SUBROUTINE FindFirstStringWithinBlock(iunit,dumstring,ifind)
USE crunchtype
USE strings
USE params
 
IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                    :: iunit
CHARACTER (LEN=mls), INTENT(IN)                             :: dumstring
INTEGER(I4B), INTENT(OUT)                                   :: ifind

!  Internal variables

CHARACTER (LEN=mls)                                         :: dummy1
INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs

ifind = 0

main:  DO
  READ(iunit,'(a)') dummy1
  id = 1
  iff = mls
  CALL sschaine(dummy1,id,iff,ssch,ids,ls)
  CALL squeeze(ssch,ls)
  IF (ssch == dumstring) THEN
    ifind = 1
    EXIT main
  ELSE IF (dummy1(1:1) == '+') THEN
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
END SUBROUTINE FindFirstStringWithinBlock

