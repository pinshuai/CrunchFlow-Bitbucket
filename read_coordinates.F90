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

SUBROUTINE read_coordinates(nout)
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
  IF (ssch == 'coordinates' .OR. ssch == 'coordinate') THEN
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
  IF (ssch == 'rectangular') THEN
    cylindrical = .FALSE.
    rectangular = .TRUE.
    spherical = .FALSE.
  ELSE IF (ssch == 'cylindrical') THEN
    cylindrical = .TRUE.
    rectangular = .FALSE.
    spherical = .FALSE.
  ELSE IF (ssch == 'spherical') THEN
    spherical = .TRUE.
    rectangular = .FALSE.
    cylindrical = .FALSE.
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Coordinates option not recognized: ',ssch(1:ls)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE   
  WRITE(*,*)
  WRITE(*,*) ' No coordinates option given'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

RETURN
END SUBROUTINE read_coordinates
