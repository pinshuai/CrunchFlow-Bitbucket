!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:56:23
 
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

SUBROUTINE exactivity_init(ncomp,nexchange,nexch_sec,nco)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nexchange
INTEGER(I4B), INTENT(IN)                                   :: nexch_sec
INTEGER(I4B), INTENT(IN)                                   :: nco

!  Internal variables

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: nex

REAL(DP)                                                   :: activity

IF (iexc == 1) THEN        ! Gaines-Thomas convention
  DO ix = 1,nexchange
    sumactivity(ix) = 0.0
  END DO
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    activity = spextmp10(nex+nexchange)*muexc(nex,ix+ncomp)/totexch(ix,nco)
    aexch(nex) = activity
    sumactivity(ix) = sumactivity(ix) + activity
  END DO
ELSE IF (iexc == 2) THEN   ! Vanselow convention
  DO ix = 1,nexchange
    tec(ix) = 0.0
  END DO
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    tec(ix) = tec(ix) + spextmp10(nex+nexchange)
  END DO
  DO ix = 1,nexchange
    sumactivity(ix) = 0.0
  END DO
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    aexch(nex) = spextmp10(nex+nexchange)/tec(ix)
    sumactivity(ix) = sumactivity(ix) + aexch(nex)
  END DO
ELSE IF (iexc == 3) THEN   ! Gapon convention
ELSE
  CONTINUE
END IF

RETURN
END SUBROUTINE exactivity_init
