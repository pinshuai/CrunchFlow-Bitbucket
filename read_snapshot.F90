!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:05:49
 
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

SUBROUTINE read_snapshot(nout,lchar,parchar,parchar2,parfind,realmult,lenarray,section)
USE crunchtype
USE CrunchFunctions
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: lchar
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar2
CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind
INTEGER(I4B), INTENT(OUT)                                   :: lenarray

REAl(DP), DIMENSION(:), INTENT(IN OUT)                      :: realmult

CHARACTER (LEN=mls), INTENT(IN)                             :: section

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: npar

REAL(DP)                                                    :: realjunk

LOGICAL(LGT)                                                :: continuation

continuation = .FALSE.

REWIND nout

npar = 0

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
ELSE
  GOTO 100
END IF
IF (ssch == parchar .OR. ssch == parchar2) THEN
  parfind = parchar
  lchar = ls
  GO TO 200
ELSE
  GO TO 100
END IF
300 RETURN

200 CONTINUE

IF (continuation) THEN
  id = 1
  continuation = .false.
ELSE
  id = ids + ls
END IF


CALL sschaine(zone,id,iff,ssch,ids,ls)

IF (ssch == '&') THEN
  READ(nout,'(a)') zone
  continuation = .TRUE.
  GOTO 200
END IF

IF(ls /= 0) THEN
  npar = npar + 1
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res /= 'n') THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameter should be followed by a numeric value'
    WRITE(*,*) ' In section ',section
    WRITE(*,*) ' Following parameter ',parchar
    WRITE(*,*) ' Aborting run'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  realmult(npar) = DNUM(ssch)
  lenarray = npar
  GO TO 200
ELSE
  npar = npar - 1
  realjunk = -500.0
END IF

RETURN
END SUBROUTINE read_snapshot
