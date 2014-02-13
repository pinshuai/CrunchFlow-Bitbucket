SUBROUTINE CountSubstrings(nout,StringCount)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
INTEGER(I4B), INTENT(OUT)                                       :: StringCount

!  Internal variables and arrays

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: lchar


REWIND nout

StringCount = 0

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
lzs=ls
StringCount = StringCount + 1
300 RETURN

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  StringCount = StringCount + 1
  GO TO 200
ELSE
  CONTINUE
END IF

RETURN
END SUBROUTINE CountSubstrings
