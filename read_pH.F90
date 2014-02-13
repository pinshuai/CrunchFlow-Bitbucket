!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:06:48
 
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

SUBROUTINE read_ph(nout,ph,guessph,i,isolution,constraint,nrct,phfound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
REAL(DP), DIMENSION(:), INTENT(IN OUT)                      :: ph
REAL(DP), DIMENSION(:), INTENT(IN OUT)                      :: guessph
INTEGER(I4B), INTENT(IN)                                    :: i
INTEGER(I4B), INTENT(IN)                                    :: isolution
CHARACTER (LEN=mls), DIMENSION(:,:), INTENT(IN OUT)         :: constraint
INTEGER(I4B), INTENT(IN)                                    :: nrct
LOGICAL(LGT), INTENT(IN OUT)                                :: phfound

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: tempstring

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: k


REWIND nout

100 READ(nout,'(a)',END=300) zone
nlen1 = LEN(zone)
!      call majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
!        call convan(ssch,lzs,res)
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
IF (ssch == 'pH' .OR. ssch == 'ph' .OR. ssch == 'PH') THEN
  phfound = .true.
ELSE
  GO TO 100
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
!        call convan(ssch,lzs,res)
  CALL stringtype(ssch,lzs,res)
  IF (res == 'n') THEN
    itype(i,isolution) = 7
!  Read the pH value
    ph(isolution) = DNUM(ssch)
  ELSE
!  An ascii string, so check for charge balance, then a mineral name
    IF (ssch == 'charge' .OR. ssch == 'Charge' .OR. ssch == 'CHARGE') THEN
      itype(i,isolution) = 2
      GO TO 200
    END IF
    DO k = 1,nrct
      tempstring = umin(k)
      IF (ssch == tempstring) THEN
        constraint(i,isolution) = umin(k)
        itype(i,isolution) = 3
        GO TO 200
      END IF
    END DO
    WRITE(*,*)
    WRITE(*,*) ' Mineral constraint for pH not found'
    WRITE(*,5050) isolution
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF
!  Now, check to see if there is an optional guess for pH
200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
!  Read the guess for the pH value
    guessph(isolution) = DNUM(ssch)
  ELSE
!  Anything else, ignore it
    CONTINUE
  END IF
  
END IF

5050 FORMAT(1X,'Condition number ',i2,' in input file')

300  RETURN
END SUBROUTINE read_ph
