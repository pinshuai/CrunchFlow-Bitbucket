!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:56:08
 
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

SUBROUTINE ex_activity(ncomp,nexchange,nexch_sec,jx,jy,jz)
USE crunchtype
USE concentration
USE params
USE transport
USE medium
USE temperature

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nexchange
INTEGER(I4B), INTENT(IN)                                   :: nexch_sec
INTEGER(I4B), INTENT(IN)                                   :: jx
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz

!  Internal variables

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: nex

REAL(DP)                                                   :: activity
REAL(DP)                                                   :: PorSatRo

PorSatRo = por(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz)

IF (iexc == 1) THEN        ! Gaines-Thomas convention
  DO ix = 1,nexchange
    sumactivity(ix) = 0.0
  END DO
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
!CIS    activity = spex10(nex+nexchange,jx,jy,jz)*PorSatRo*muexc(nex,ix+ncomp)/  &
    activity = spex10(nex+nexchange,jx,jy,jz)*muexc(nex,ix+ncomp)/  &
        exchangesites(ix,jx,jy,jz)
    aexch(nex) = activity
    sumactivity(ix) = sumactivity(ix) + activity
  END DO
ELSE IF (iexc == 2) THEN   ! Vanselow convention
  DO ix = 1,nexchange
    tec(ix) = 0.0
  END DO
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    tec(ix) = tec(ix) + spex10(nex+nexchange,jx,jy,jz)
  END DO
  DO ix = 1,nexchange
    sumactivity(ix) = 0.0
  END DO
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    aexch(nex) = spex10(nex+nexchange,jx,jy,jz)/tec(ix)
    sumactivity(ix) = sumactivity(ix) + aexch(nex)
  END DO
ELSE IF (iexc == 3) THEN   ! Gapon convention
ELSE
  CONTINUE
END IF

RETURN
END SUBROUTINE ex_activity
