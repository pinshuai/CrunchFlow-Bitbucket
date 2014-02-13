!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:09:27
 
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

SUBROUTINE totexchange_init(ncomp,nexchange,nexch_sec,nco)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nexchange
INTEGER(I4B), INTENT(IN)                                    :: nexch_sec
INTEGER(I4B), INTENT(IN)                                    :: nco

!  Internal variables and arrays

INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: ix
INTEGER(I4B)                                                :: nex

REAL(DP)                                                    :: sum

DO i = 1,ncomp
  sum=0.0
  IF (equilibrate(i,nco)) THEN
    DO nex = 1,nexch_sec
      sum = sum + muexc(nex,i)*spextmp10(nex+nexchange)
    END DO
  END IF
  sexch(i) = sum
END DO

DO ix = 1,nexchange
  sum = 0.0
  DO nex = 1,nexch_sec
    sum = sum + muexc(nex,ix+ncomp)*spextmp10(nex+nexchange)
  END DO
  totextmp(ix) = sum
END DO

RETURN
END SUBROUTINE totexchange_init
!  **************************************************************
