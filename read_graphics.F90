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

SUBROUTINE read_graphics(nout)
USE crunchtype
USE params
USE runtime
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout

!  Internal variables and arrays

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'graphics' .OR. ssch == 'graphic') THEN
    GO TO 200
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF
300 RETURN

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'xtool') THEN
    xtool = .TRUE.
  ELSE IF (ssch == 'tecplot') THEN
    tecplot = .TRUE.
  ELSE IF (ssch == 'originlab') THEN
    originlab = .TRUE.
  ELSE IF (ssch == 'xmgr') THEN
    xmgr = .TRUE.
  ELSE IF (ssch == 'kaleidagraph') THEN
    kaleidagraph = .TRUE.
  ELSE IF (ssch == 'nview') THEN
    nview = .TRUE.
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Dont recognize this graphics option: ', ssch(1:ls)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No graphics option given in RUNTIME block'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

!!  Look for another graphics specification

GO TO 100

RETURN
END SUBROUTINE read_graphics
