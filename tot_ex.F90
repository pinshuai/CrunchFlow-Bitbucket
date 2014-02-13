!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:09:10
 
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

SUBROUTINE tot_ex(ncomp,nexchange,nexch_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                       :: ncomp
INTEGER(I4B), INTENT(IN)                                       :: nexchange
INTEGER(I4B), INTENT(IN)                                       :: nexch_sec
INTEGER(I4B), INTENT(IN)                                       :: jx
INTEGER(I4B), INTENT(IN)                                       :: jy
INTEGER(I4B), INTENT(IN)                                       :: jz

!  Internal variables and arrays

INTEGER(I4B)                                                   :: ix
INTEGER(I4B)                                                   :: nex

REAL(DP)                                                       :: sum


DO ix = 1,nexchange
  sum = 0.0
  DO nex = 1,nexch_sec
    sum = sum + muexc(nex,ix+ncomp)*spex10(nex+nexchange,jx,jy,jz)
  END DO
  totex(ix) = sum
END DO

RETURN
END SUBROUTINE tot_ex
!  **************************************************************
