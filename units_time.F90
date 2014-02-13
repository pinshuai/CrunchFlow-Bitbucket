!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:10:17
 
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

SUBROUTINE units_time(nout,section,time_scale)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
CHARACTER (LEN=mls), INTENT(IN)                                 :: section
REAL(DP), INTENT(OUT)                                           :: time_scale

!  Internal variables and arrays

CHARACTER (LEN=mls)                                             :: parchar
CHARACTER (LEN=mls)                                             :: parfind
CHARACTER (LEN=mls)                                             :: dumstring
CHARACTER (LEN=mls)                                             :: time_units

INTEGER(I4B)                                                    :: lchar

time_scale = 1.0d0
parchar = 'time_units'
parfind = ' '
time_units = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN  ! Parameter "time_units" not found
  time_units = 'years'             ! Use default
ELSE
  time_units = dumstring
  
!  Check to see that the units are recognized
  IF (time_units == 'year' .OR. time_units == 'years') THEN
    time_units = 'years'
  END IF
  IF (time_units == 'day' .OR. time_units == 'days') THEN
    time_units = 'days'
  END IF
  IF (time_units == 'hour' .OR. time_units == 'hours') THEN
    time_units = 'hours'
  END IF
  IF (time_units == 'minute' .OR. time_units == 'minutes') THEN
    time_units = 'minutes'
  END IF
  IF (time_units == 'second' .OR. time_units == 'seconds') THEN
    time_units = 'seconds'
  END IF
  
  IF (time_units == 'years') THEN
    time_scale = 1.0d0
  ELSE IF (time_units == 'days') THEN
    time_scale = 1.0d0/(365.0d0)
  ELSE IF (time_units == 'hours') THEN
    time_scale = 1.0d0/(365.0d0*24.0d0)
  ELSE IF (time_units == 'minutes') THEN
    time_scale = 1.0d0/(365.0d0*24.0d0*60.0d0)
  ELSE IF (time_units == 'seconds') THEN
    time_scale = 1.0d0/(365.0d0*24.0d0*60.0d0*60.0d0)
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Time units not recognized'
    WRITE(*,*) ' Section = ',section
    WRITE(*,*) ' Using "years" as time unit'
    WRITE(*,*)
    time_scale = 1.0d0
    READ(*,*)
    STOP
  END IF
END IF

RETURN
END SUBROUTINE units_time
