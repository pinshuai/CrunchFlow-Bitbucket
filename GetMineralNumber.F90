SUBROUTINE GetMineralNumber (nrct,dumstring,klabel)
USE CrunchType
USE params, ONLY: mls
USE mineral,ONLY: umin

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)               :: nrct
CHARACTER (LEN=mls), INTENT(IN)        :: dumstring
INTEGER(I4B), INTENT(OUT)              :: klabel
INTEGER(I4B)                           :: k
INTEGER(I4B)                           :: ls
LOGICAL (LGT)                          :: MineralFound

MineralFound = .FALSE.
DO k = 1,nrct
  IF (dumstring == umin(k)) THEN
    MineralFound = .TRUE.
    klabel = k
  END IF
END DO
IF (.NOT. MineralFound) THEN
  CALL stringlen(dumstring,ls)
  WRITE(*,*) 
  WRITE(*,*) ' Mineral not found in list'
  WRITE(*,*) ' Looking for: ', dumstring(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

END SUBROUTINE GetMineralNumber