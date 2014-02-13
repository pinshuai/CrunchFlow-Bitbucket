SUBROUTINE majuscules(zone,nbc)
USE crunchtype
USE params
 
IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: nbc

CHARACTER (LEN=mls), INTENT(IN OUT)                   :: zone

!  Internal variables

INTEGER(I4B)                                               :: ic
INTEGER(I4B)                                               :: i

DO i=1,nbc
  ic=ICHAR(zone(i:i))
!         if(ic.gt.96.and.ic.lt.123) then
!           zone(i:i)=char(ic-32)
!         endif
!  Convert to lower case
  IF(ic > 64.AND.ic < 91) THEN
    zone(i:i)=CHAR(ic+32)
  END IF
END DO
END SUBROUTINE majuscules


