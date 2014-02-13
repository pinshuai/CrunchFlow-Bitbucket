SUBROUTINE PestExchange(nout,ncomp,nchem)
USE CrunchType
USE CrunchFunctions
USE params
USE RunTime, ONLY: CreatePestInstructionFile,PestExchangeOutputFile
USE strings
USE concentration

IMPLICIT NONE

INTERFACE
  SUBROUTINE PrimarySpeciesCheck(ncomp,npar,DummyStringArray)
  USE CrunchType
  USE concentration
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                        :: ncomp
  INTEGER(I4B), INTENT(IN)                                        :: npar
  CHARACTER (LEN=mls), DIMENSION(:)                               :: DummyStringArray
  END SUBROUTINE PrimarySpeciesCheck
END INTERFACE

!! External arrays and variables

INTEGER(I4B), INTENT(IN)                                      :: nout
INTEGER(I4B), INTENT(IN)                                      :: ncomp
INTEGER(I4B), INTENT(IN)                                      :: nchem

!! Internal arrays and variables
INTEGER(I4B)                                                  :: id
INTEGER(I4B)                                                  :: iff
INTEGER(I4B)                                                  :: ids
INTEGER(I4B)                                                  :: ls
INTEGER(I4B)                                                  :: lzs
INTEGER(I4B)                                                  :: lchar
INTEGER(I4B)                                                  :: l_string
INTEGER(I4B)                                                  :: lenarray
INTEGER(I4B)                                                  :: npar
INTEGER(I4B)                                                  :: nco
INTEGER(I4B)                                                  :: nl

CHARACTER (LEN=mls)                                           :: dumstring
CHARACTER (LEN=mls)                                           :: dumstring2
CHARACTER (LEN=mls)                                           :: dumstring3
CHARACTER (LEN=mls)                                           :: OpenBracket
CHARACTER (LEN=mls)                                           :: CloseBracket
CHARACTER (LEN=mls)                                           :: CharLengthString
CHARACTER (LEN=mls)                                           :: parchar
CHARACTER (LEN=mls)                                           :: section
CHARACTER (LEN=mls)                                           :: parfind

LOGICAL(LGT)                                                  :: continuation
LOGICAL(LGT)                                                  :: ext
LOGICAL(LGT)                                                  :: FileOpen 

continuation = .FALSE.
UnitPestExchange = 51

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
lzs=ls
CALL convan(ssch,lzs,res)
IF (ssch == 'exchange') THEN
  parfind = parchar
  lchar = ls
  GO TO 200
ELSE
  GO TO 100
END IF
300 RETURN

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)                   !!  Now, look for units for exchange
l_string = ls
IF(ls /= 0) THEN
  npar = 0
  lzs=ls  
  IF (ssch == 'mol/g' .OR. ssch == 'mole/g' .OR. ssch == 'moles/g') THEN
    PestExchangeUnits = 'mol/g'
  ELSE IF (ssch == 'equiv/g' .OR. ssch == 'equivalents/g' .OR. ssch == 'equivalent/g' .OR. ssch == 'eq/g') THEN
    PestExchangeUnits = 'equiv/g'
  ELSE IF (ssch == 'milliequiv/g' .OR. ssch == 'milliequivalents/g' .OR. ssch == 'milliequivalent/g' .OR.  &
                  ssch == 'mequiv/g' .OR. ssch == 'mequivalents/g' .OR. ssch == 'mequivalent/g' .OR. ssch == 'meq/g') THEN
    PestExchangeUnits = 'mequiv/g'
  ELSE IF (ssch == 'millimol/g' .OR. ssch == 'millimol/g' .OR. ssch == 'millimoles/g' .OR.  &
                  ssch == 'mmol/g' .OR. ssch == 'mmole/g' .OR. ssch == 'mmoles/g') THEN
    PestExchangeUnits = 'mmol/g'
  ELSE IF (ssch == 'microequiv/g' .OR. ssch == 'microequivalents/g' .OR. ssch == 'microequivalent/g' .OR.  &
                  ssch == 'uequiv/g' .OR. ssch == 'uequivalents/g' .OR. ssch == 'uequivalent/g' .OR. ssch == 'ueq/g') THEN
    PestExchangeUnits = 'uequiv/g'
  ELSE IF (ssch == 'micromol/g' .OR. ssch == 'micromole/g' .OR. ssch == 'micromoles/g' .OR. &
             ssch == 'umol/g' .OR. ssch == 'umole/g' .OR. ssch == 'umoles/g') THEN
    PestExchangeUnits = 'umol/g'
  ELSE IF (ssch == 'log' .OR. ssch == 'logs' .OR. ssch == 'logarithm' .OR. ssch == 'logarithms') THEN
    PestExchangeUnits = 'logequivalents'
  ELSE IF (ssch == 'logequivalents' .OR. ssch == 'logequivalent') THEN
    PestExchangeUnits = 'logequivalents'
  ELSE IF (ssch == 'logmoles' .OR. ssch == 'logmole') THEN
    PestExchangeUnits = 'logmoles'
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Pest exchange units not recognized in PEST block'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  500 CONTINUE
  IF (continuation) THEN
    id = 1
    continuation = .FALSE.
  ELSE
    id = ids + ls
  END IF

  CALL sschaine(zone,id,iff,ssch,ids,ls)

  IF (ssch == '&') THEN
    READ(nout,'(a)') zone
    continuation = .TRUE.
    GOTO 500
  END IF

  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    npar = npar + 1
    stringarray(npar) = ssch
    GO TO 500
  ELSE
    CONTINUE
  END IF

  lenarray = npar

  IF (npar == 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' No species specified after "exchange" keyword in PEST block'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  NPestExchange = npar
  IF (ALLOCATED(PestExchangeList)) THEN
    DEALLOCATE(PestExchangeList)
    ALLOCATE(PestExchangeList(NPestExchange))
  ELSE
    ALLOCATE(PestExchangeList(NPestExchange))
  END IF
  PestExchangeList(1:npar) = StringArray(1:NPestExchange)
    
  CALL PrimarySpeciesCheck(ncomp,NPestExchange,PestExchangeList)

  IF (CreatePestInstructionFile) THEN
    OPEN(UNIT=UnitPestExchange,FILE='PestExchange.ins',STATUS='unknown')         !! Create the Pest instruction (format) file for the output in PestExchange.out

    dumstring2 = 'l1'
    CALL squeeze(dumstring2,ls)
    OpenBracket = ' ['
    CloseBracket = ']'
    CharLengthString = '1:12'

    WRITE(UnitPestExchange,501)
    DO nco = 1,nchem
      dumstring = condlabel(nco)
      CALL stringlen(dumstring,ls)
      DO nl = 1,NPestExchange
        dumstring3 = IntegerToCharacter(nl)
        CALL squeeze(dumstring3,lchar)
 !!       WRITE(UnitPestExchange,502) dumstring2, OpenBracket(1:2), dumstring(1:ls), dumstring3(1:lchar), CloseBracket(1:1), CharLengthString(1:4)
       WRITE(UnitPestExchange,502) dumstring2, OpenBracket, dumstring, dumstring3, CloseBracket, CharLengthString
      END DO
    END DO
    CLOSE(UNIT=UnitPestExchange,STATUS='keep')
  END IF

501 FORMAT('pif @')
502 FORMAT(a2,a2,a<ls>,a<lchar>,a1,a4)

  INQUIRE(FILE=PestExchangeOutputFile,OPENED=FileOpen)
  IF (FileOpen) THEN
    CLOSE(UNIT=UnitPestExchange,STATUS='delete')
    OPEN(UNIT=UnitPestExchange,FILE=PestExchangeOutputFile,STATUS='new')
  ELSE
    OPEN(UNIT=UnitPestExchange,FILE=PestExchangeOutputFile,STATUS='unknown')
  END IF

ELSE
  WRITE(*,*)
  WRITE(*,*) ' No trailing string found after'
  WRITE(*,*) parchar
  WRITE(*,*) ' In section ',section
  WRITE(*,*)
  parfind = ' '
END IF

RETURN
END SUBROUTINE PestExchange
