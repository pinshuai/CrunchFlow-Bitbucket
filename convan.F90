SUBROUTINE convan(ssch,ls,res)
USE crunchtype
USE params

IMPLICIT NONE

!  External arrays and variables

CHARACTER (LEN=mls), INTENT(IN OUT)                    :: ssch
CHARACTER (LEN=1), INTENT(OUT)                         :: res
INTEGER(I4B), INTENT(IN OUT)                           :: ls

!  Internal arrays and variables

INTEGER(I4B)                                           :: i
INTEGER(I4B)                                           :: ic

!  verification si ssch est numerique ou alphabetique
res='n'
CALL majuscules(ssch,ls)
DO i=1,ls
  ic=ICHAR(ssch(i:i))
  IF((ic >= 97.AND.ic <= 99).OR.(ic >= 102.AND.ic <= 122)) THEN
    res='a'
    GO TO 100
  END IF
END DO
100   CONTINUE

!   recherche correspondance nom/indice et stockage de l'indice
IF(res == 'n') RETURN
END SUBROUTINE convan

