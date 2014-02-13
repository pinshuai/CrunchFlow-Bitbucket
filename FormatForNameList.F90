SUBROUTINE FormatForNameList(nout,NameListFormat)
USE crunchtype
USE params, ONLY: mls
USE strings, ONLY: zone, ssch

IMPLICIT NONE

!! External arrays and variables

INTEGER(I4B), INTENT(IN)                              :: nout
LOGICAL(LGT), INTENT(OUT)                             :: NameListFormat

!! Internal arrays and variables

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls

NAMELIST / BlockReadFormat / NameListFormat

NameListFormat = .FALSE.

READ(nout,'(a)') zone

id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)

IF (ssch == '&BlockReadFormat' .OR. ssch == '&BLOCKREADFORMAT' .OR. ssch == '&blockreadformat') THEN
  BACKSPACE nout
  READ(nout, NML = BlockReadFormat) 
ELSE
  BACKSPACE nout
  RETURN
END IF

RETURN
END SUBROUTINE FormatForNameList