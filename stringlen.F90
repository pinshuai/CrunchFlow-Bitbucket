SUBROUTINE stringlen(substring,ls)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:08:46
 
USE crunchtype
USE params
USE strings

IMPLICIT NONE

CHARACTER (LEN=mls), INTENT(IN)                      :: substring
INTEGER(I4B), INTENT(OUT)                            :: ls

INTEGER(I4B)                                         :: i

ls=0
DO i=1,mls
  ls = i
  IF(substring(i:i) == ' ') THEN
    ls=ls -1
    GO TO 100
  END IF
END DO

100 RETURN
END SUBROUTINE stringlen

