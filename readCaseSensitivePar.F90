!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:06:27
 
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

SUBROUTINE readCaseSensitivePar(nout,lchar,parchar,parfind,realjunk,section)
USE crunchtype
USE CrunchFunctions
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
INTEGER(I4B), INTENT(OUT)                                       :: lchar

CHARACTER (LEN=mls), INTENT(IN)                                 :: section
CHARACTER (LEN=mls), INTENT(IN)                                 :: parchar
CHARACTER (LEN=mls), INTENT(IN OUT)                             :: parfind
REAL(DP), INTENT(IN OUT)                                        :: realjunk

!  Internal variables and arrays

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: lensection
INTEGER(I4B)                                                    :: lenparchar

REWIND nout

CALL stringlen(section,lensection)
CALL stringlen(parchar,lenparchar)

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL stringtype(ssch,lzs,res)
  IF (res /= 'a') THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameters should start with an ASCII string'
    WRITE(*,*) '   In block: ',section(1:lensection)
    WRITE(*,*) '   Looking for string: ',parchar(1:lenparchar)
    WRITE(*,*) ' ABORTING RUN'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF
IF (ssch == parchar) THEN
  parfind = parchar
  lchar = ls
  GO TO 200
ELSE
  GO TO 100
END IF
300 RETURN

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL stringtype(ssch,lzs,res)
  IF (res /= 'n') THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameter should be followed by a numeric value'
    WRITE(*,*) '   In section: ',section(1:lensection)
    WRITE(*,*) '   Following parameter: ',parchar(1:lenparchar)
    WRITE(*,*) ' ABORTING RUN'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  realjunk = DNUM(ssch)
  RETURN
ELSE
  realjunk = -500.0
END IF

RETURN
END SUBROUTINE readCaseSensitivePar
