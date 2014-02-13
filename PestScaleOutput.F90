SUBROUTINE PestScaleOutput(nout)
USE crunchtype
USE CrunchFunctions
USE params
USE strings
USE concentration

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: idum

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  IF (ssch == 'PestScaleOutput' .OR. ssch == 'pestscaleoutput') THEN
    GO TO 200
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  DO idum = 1,nplot
    IF (ssch == ulab(iplot(idum))) THEN
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL stringtype(ssch,lzs,res)
        IF (res /= 'n') THEN
          WRITE(*,*)
          WRITE(*,*) '  Species name in PestScaleOutput keyword should be followed by a numeric value'
          WRITE(*,*) '       ABORTING RUN  '
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        PestScale(idum)= DNUM(ssch)
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' No numerical value for PestScaleOutput given'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    END IF
  END DO
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No information given after PestScaleOutput keyword'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

GOTO 100

300 RETURN

END SUBROUTINE PestScaleOutput