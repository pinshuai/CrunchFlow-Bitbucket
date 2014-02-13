!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:03:46
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE read_ionexchangeMIN(nout,ix,isolution,ncomp,nexchange,speciesfound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ix
INTEGER(I4B), INTENT(IN)                                    :: isolution
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nexchange
LOGICAL(LGT), INTENT(OUT)                                   :: speciesfound

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lsave

!!   This routine reads ion exchange concentrations where assumed to be on a specific mineral
!!   NOTE:  Called by find_condition.F90 to read CEC etc from individual geochemical conditions

icec(ix) = 0                     !! Flag to switch between specification of totexch (eq/kgw), to CEC plus solid:solution ratio (icec=1)
speciesfound = .false.
REWIND nout

100 READ(nout,'(a)',END=300) zone
nlen1 = LEN(zone)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL stringtype(ssch,lzs,res)
  IF (res /= 'a') THEN
    WRITE(*,*)
    WRITE(*,*) ' Geochemical input should start with a string'
    WRITE(*,*) ssch
    WRITE(*,5050) isolution
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF
IF (ssch == namexc(ix)) THEN
  speciesfound = .true.     !! Exchange species found
ELSE
  GO TO 100                 !! Cycles back to read next line in geochemical condition
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  lsave = lzs
  IF (res == 'n') THEN
!  If a number, assume it is the total exchange concentration (eq/kgw)
    totexch(ix,isolution) = DNUM(ssch)
    IF (totexch(ix,isolution) == 0.0) THEN
      totexch(ix,isolution) = 1.e-30
    END IF
    icec(ix) = 0                    !! For this option, will need to calculate CEC in find_condition.F90 based on this information
  ELSE     !  An ascii string, so...
    IF (ssch == '-cec') THEN
      totexch(ix,isolution) = 0.0
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        lsave = lzs
        IF (res == 'n') THEN
          cec(ix,isolution) = DNUM(ssch)
          IF (cec(ix,isolution) == 0.0) THEN
            cec(ix,isolution) = 1.e-30
          END IF
          icec(ix) = 1
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' A numerical value for the CEC should follow '
          WRITE(*,*) '   the parameter flag "-cec" '
          WRITE(*,*) ssch(1:lsave)
          WRITE(*,*) ' For exchanger ',namexc(ix)
          WRITE(*,5050) isolution
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' A numerical value for the CEC should follow "-cec" '
        WRITE(*,*) ssch(1:lsave)
        WRITE(*,*) ' For exchanger ',namexc(ix)
        WRITE(*,5050) isolution
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF 
    ELSE    !  Leading string is not '-cec'
      WRITE(*,*)
      WRITE(*,*) ' Either a number for the total exchange capacity (mol/kgw)'
      WRITE(*,*) '   or the parameter "-cec" are necessary'
      WRITE(*,*) ' For exchanger',namexc(ix)
      WRITE(*,5050) isolution
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Looking for total exchanger concentration'
  WRITE(*,*) ' For exchanger',namexc(ix)
  WRITE(*,5050) isolution
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

5050 FORMAT(1X,'Condition number ',i2,' in input file')

300  RETURN
END SUBROUTINE read_ionexchangeMIN
