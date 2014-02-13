!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:10:14
 
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

SUBROUTINE units_distance(nout,section,dist_scale)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
CHARACTER (LEN=mls), INTENT(IN)                                 :: section
REAL(DP), INTENT(OUT)                                           :: dist_scale

!  Internal variables and arrays

CHARACTER (LEN=mls)                                             :: parchar
CHARACTER (LEN=mls)                                             :: parfind
CHARACTER (LEN=mls)                                             :: dumstring
CHARACTER (LEN=mls)                                             :: distance_units

INTEGER(I4B)                                                    :: lchar

!  Search for time units (default is years)

dist_scale = 1.0
parchar = 'distance_units'
parfind = ' '
distance_units = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN  ! Parameter "distance_units" not found
  distance_units = 'meters'             ! Use default
ELSE
  distance_units = dumstring
  
!  Check to see that the units are recognized
  IF (distance_units == 'meter' .OR. distance_units == 'meters') THEN
    distance_units = 'meters'
  END IF
  IF (distance_units == 'centimeter' .OR. distance_units == 'centimeters') THEN
    distance_units = 'centimeters'
  END IF
  IF (distance_units == 'millimeter' .OR. distance_units == 'millimeters') THEN
    distance_units = 'millimeters'
  END IF
  IF (distance_units == 'kilometer' .OR. distance_units == 'kilometers') THEN
    distance_units = 'kilometers'
  END IF
  IF (distance_units == 'micrometer' .OR. distance_units == 'micrometers') THEN
    distance_units = 'micrometers'
  END IF
  IF (distance_units == 'micron' .OR. distance_units == 'microns') THEN
    distance_units = 'micrometers'
  END IF
  IF (distance_units == 'nanometer' .OR. distance_units == 'nanometers') THEN
    distance_units = 'nanometers'
  END IF
  IF (distance_units == 'm') THEN
    distance_units = 'meters'
  END IF
  IF (distance_units == 'cm') THEN
    distance_units = 'centimeters'
  END IF
  IF (distance_units == 'mm') THEN
    distance_units = 'millimeters'
  END IF
  IF (distance_units == 'km') THEN
    distance_units = 'kilometers'
  END IF
  IF (distance_units == 'um') THEN
    distance_units = 'micrometers'
  END IF
  IF (distance_units == 'nm') THEN
    distance_units = 'nanometers'
  END IF
  
  IF (distance_units == 'meters') THEN
    dist_scale = 1.0
  ELSE IF (distance_units == 'centimeters') THEN
    dist_scale = 100.0
  ELSE IF (distance_units == 'millimeters') THEN
    dist_scale = 1000.0
  ELSE IF (distance_units == 'kilometers') THEN
    dist_scale = 0.001
  ELSE IF (distance_units == 'micrometers') THEN
    dist_scale = 1.0E06
  ELSE IF (distance_units == 'nanometers') THEN
    dist_scale = 1.0E09
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Distance units not recognized'
    WRITE(*,*) ' Section = ',section
    WRITE(*,*) ' Using "meters" as distance unit'
    WRITE(*,*)
    dist_scale = 1.0
    READ(*,*)
    STOP
  END IF
END IF

RETURN
END SUBROUTINE units_distance
