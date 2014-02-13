SUBROUTINE ReadSaturation(nout,nco,SaturationTemp,SaturationFound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nco
REAL(DP), INTENT(OUT)                                       :: SaturationTemp
LOGICAL(LGT), INTENT(IN OUT)                                :: SaturationFound

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: lcond

CHARACTER (LEN=mls)                                         :: dumstring

dumstring = condlabel(nco)
CALL stringlen(dumstring,lcond)

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'set_saturation') THEN
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
  CALL convan(ssch,lzs,res)
  IF (res /= 'n') THEN
    WRITE(*,*)
    WRITE(*,*) ' "Set_saturation" should be followed by a numeric value'
    WRITE(*,*) '   In condition: ', dumstring(1:lcond)
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  SaturationFound = .true.
  SaturationTemp = DNUM(ssch)
  RETURN
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No numerical value for liquid saturation given '
  WRITE(*,*) '   In condition: ', dumstring(1:lcond)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

300 RETURN

END SUBROUTINE ReadSaturation
