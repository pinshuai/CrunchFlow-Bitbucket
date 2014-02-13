SUBROUTINE GetPrimarySpeciesNumber (ncomp,dumstring,ilabel)
USE CrunchType
USE params, ONLY: mls
USE concentration,ONLY: ulab

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)               :: ncomp
CHARACTER (LEN=mls), INTENT(IN)        :: dumstring
INTEGER(I4B), INTENT(OUT)              :: ilabel
INTEGER(I4B)                           :: i
INTEGER(I4B)                           :: ls
LOGICAL (LGT)                          :: SpeciesFound

SpeciesFound = .FALSE.
DO i = 1,ncomp
  IF (dumstring == ulab(i)) THEN
    SpeciesFound = .TRUE.
    ilabel = i
  END IF
END DO
IF (.NOT. SpeciesFound) THEN
  CALL stringlen(dumstring,ls)
  WRITE(*,*) 
  WRITE(*,*) ' Primary species not found in list'
  WRITE(*,*) ' Looking for: ', dumstring(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

END SUBROUTINE GetPrimarySpeciesNumber