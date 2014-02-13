!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:07:17
 
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

SUBROUTINE read_gravity(nout)
USE crunchtype
USE CrunchFunctions
USE params
USE flow
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs

CHARACTER (LEN=mls)                                         :: dumstring

x_angle = 90.0
y_angle = 90.0
z_angle = 90.0
SignGravity = 1.000

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'gravity') THEN
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
    WRITE(*,*) ' "Gravity" should be followed by a numeric value giving the angle of the gravity vector relative to X coordinate'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  x_angle = DNUM(ssch)
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No numerical value for X component of gravity given '
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res /= 'n') THEN
    WRITE(*,*)
    WRITE(*,*) ' "Gravity" should be followed by a numeric value giving the angle of the gravity vector relative to Y coordinate'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  y_angle = DNUM(ssch)
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No numerical value for Y component of gravity given '
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res /= 'n') THEN
    WRITE(*,*)
    WRITE(*,*) ' "Gravity" should be followed by a numeric value giving the angle of the gravity vector relative to Z coordinate'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  z_angle = DNUM(ssch)
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No numerical value for Z component of gravity given '
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'down') THEN
    SignGravity = 1.00
  ELSE IF (ssch == 'up') THEN
    SignGravity = -1.00
  ELSE
    SignGravity = 1.00
  END IF
ELSE
  SignGravity = 1.00
END IF

300 RETURN

END SUBROUTINE read_gravity
