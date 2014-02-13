SUBROUTINE PrimarySpeciesCheck(ncomp,npar,DummyStringArray)
USE CrunchType
USE concentration

IMPLICIT NONE

!! External arrays and variables

INTEGER(I4B), INTENT(IN)                                        :: ncomp
INTEGER(I4B), INTENT(IN)                                        :: npar

CHARACTER (LEN=mls), DIMENSION(:)                               :: DummyStringArray

!! Internal arrays and variables

LOGICAL(LGT)                                :: PrimarySpeciesFound

CHARACTER (LEN=mls)                         :: dumstring

INTEGER(I4B)                                :: i
INTEGER(I4B)                                :: ndum
INTEGER(I4B)                                :: ls

DO ndum = 1,npar
  dumstring = DummyStringArray(ndum)
  CALL stringlen(dumstring,ls)
  DO i = 1,ncomp
    IF (dumstring == ulab(i)) THEN
      PrimarySpeciesFound = .TRUE.
    END IF
  END DO
  
  IF (PrimarySpeciesFound) THEN
    CONTINUE
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Primary species not found in list: ', dumstring(1:ls)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

END DO

RETURN
END SUBROUTINE PrimarySpeciesCheck

