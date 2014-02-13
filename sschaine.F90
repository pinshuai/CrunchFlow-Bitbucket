SUBROUTINE sschaine(zone,id,iff,ssch,ids,ls)
USE crunchtype
USE params
 
IMPLICIT NONE

!  External variables and arrays

CHARACTER (LEN=mls), INTENT(OUT)                           :: ssch
CHARACTER (LEN=mls), INTENT(IN)                            :: zone

INTEGER(I4B), INTENT(IN)                                   :: id
INTEGER(I4B), INTENT(IN)                                   :: iff

INTEGER(I4B), INTENT(OUT)                                  :: ids
INTEGER(I4B), INTENT(OUT)                                  :: ls

!  Internal variables and arrays

CHARACTER (LEN=1)                                          :: delim

INTEGER(I4B)                                               :: is
INTEGER(I4B)                                               :: i

ls=0
ssch=' '
DO i=id,iff
  IF(zone(i:i) /= ' ') THEN
    is=i
    ls=1
    ids=i
    delim=' '
    DO WHILE(zone(is:is) /= delim)
      IF(zone(i:i) == '"'.AND.is > i) delim='"'
      ssch(ls:ls)=zone(is:is)
      ls=ls+1
      is=is+1
    END DO
    IF(delim /= '"') ls=ls-1
    RETURN
  END IF
END DO

RETURN
END SUBROUTINE sschaine

