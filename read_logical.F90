!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:05:14
 
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

SUBROUTINE read_logical(nout,lchar,parchar,parfind, parlogic)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
INTEGER(I4B), INTENT(OUT)                                       :: lchar

CHARACTER (LEN=mls), INTENT(IN)                                 :: parchar
CHARACTER (LEN=mls), INTENT(IN OUT)                             :: parfind
LOGICAL(LGT), INTENT(IN OUT)                                    :: parlogic

!  Internal variables and arrays

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: iflaglogic

CHARACTER (LEN=mls)                                             :: LowerCase
CHARACTER (LEN=mls)                                             :: LowerCase_ssch
INTEGER(I4B)                                                    :: lenparchar


CALL stringlen(parchar,lenparchar)

LowerCase = parchar
call majuscules(LowerCase,lenparchar)

REWIND nout

iflaglogic = 0

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL stringtype(ssch,lzs,res)
END IF

LowerCase_ssch = ssch
call majuscules(LowerCase_ssch,lzs)
IF (LowerCase_ssch == LowerCase) THEN
!!IF (ssch == parchar .OR. ssch == LowerCase) THEN
!        parfind(1:ls) = parchar(1:ls)
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
  CALL convan(ssch,lzs,res)
  IF (res /= 'a') THEN
    WRITE(*,*) ' Initial character string should be followed by a LOGICAL string'
    READ(*,*)
    STOP
  END IF
  IF (ssch == 'true'.OR.ssch == 'yes'.OR.ssch == 'on') THEN
    parlogic = .true.
  ELSE IF (ssch == 'false'.OR.ssch == 'no'.OR.ssch == 'off') THEN
    parlogic = .false.
  ELSE
    WRITE(*,400) zone(1:lchar)
    WRITE(*,*) ' Should be followed by a true or false or yes or no'
    READ(*,*)
    STOP
  END IF
  RETURN
ELSE
  parfind =  ' '
END IF

!  400 format(1x,'Parameter ',a<lchar>,' is a logical')
400 FORMAT(1X,'Parameter ',a35,' is a logical')

RETURN
END SUBROUTINE read_logical
