!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:03:00
 
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

SUBROUTINE read_erosion(nout,erodex,erodey)
USE crunchtype
USE CrunchFunctions
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
REAL(DP), INTENT(OUT)                                       :: erodex
REAL(DP), INTENT(OUT)                                       :: erodey

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs

REWIND nout

10  READ(nout,'(a)',END=1000) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF (ls == 0) THEN
  GO TO 10
END IF

lzs=ls
CALL convan(ssch,lzs,res)

!  Check to see if initial substring is "erode_x"

IF (ssch == 'erode_x') THEN
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls == 0) THEN
    GO TO 10
  END IF
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    erodex = DNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' Cant interpret string following "erode_x" '
    WRITE(*,*) ' A numerical value should follow'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE IF (ssch == 'erode_y') THEN
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls == 0) THEN
    GO TO 10
  END IF
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    erodey = DNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' Cant interpret string following "erode_y" '
    WRITE(*,*) ' A numerical value should follow'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  GO TO 10
END IF

GO TO 10

1000 RETURN
END SUBROUTINE read_erosion
