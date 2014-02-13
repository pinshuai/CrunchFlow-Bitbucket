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

SUBROUTINE units_concentration(nout,unitsflag,nchem,labeltemp,lcond)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
INTEGER(I4B), DIMENSION(:), INTENT(INOUT)                       :: unitsflag
INTEGER(I4B), INTENT(IN)                                        :: nchem
CHARACTER (LEN=mls), INTENT(IN)                                 :: labeltemp
INTEGER(I4B), INTENT(IN)                                        :: lcond 

!  Internal variables and arrays

CHARACTER (LEN=mls)                                             :: parchar
CHARACTER (LEN=mls)                                             :: parfind
CHARACTER (LEN=mls)                                             :: dumstring
CHARACTER (LEN=mls)                                             :: concentration_units

INTEGER(I4B)                                                    :: lchar
INTEGER(I4B)                                                    :: ls

parchar = 'units'
parfind = ' '
concentration_units = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,labeltemp)
IF (parfind == ' ') THEN  ! Parameter "units" not found
  concentration_units = 'mol/kg'             ! Use default
ELSE
  concentration_units = dumstring
  
  CALL stringlen(concentration_units,ls)

!  Check to see that the units are recognized

  IF (concentration_units == 'mol/kg' .OR. concentration_units == 'mol/kgw' .OR.        &
          concentration_units == 'mole/kg' .OR. concentration_units == 'mole/kgw') THEN
    unitsflag(nchem) = 1
  ELSE IF (concentration_units == 'ppm') THEN
    unitsflag(nchem) = 2
  ELSE IF (concentration_units == 'mmol/kg') THEN
    unitsflag(nchem) = 3
  ELSE IF (concentration_units == 'umol/kg') THEN
    unitsflag(nchem) = 4
  ELSE IF (concentration_units == 'micromol/kg') THEN
    unitsflag(nchem) = 4
  ELSE IF (concentration_units == 'molar' .OR. concentration_units == 'molarity' .OR.    &
       concentration_units == 'mol/L') THEN
    unitsflag(nchem) = 5
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Concentration units not recognized'
    WRITE(*,*) ' Units: ',concentration_units(1:ls)
    WRITE(*,*) ' In condition: ', labeltemp(1:lcond)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

RETURN
END SUBROUTINE units_concentration
