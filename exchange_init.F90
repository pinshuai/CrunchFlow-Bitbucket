!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:56:31
 
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

SUBROUTINE exchange_init(ncomp,nspec,nexchange,nexch_sec,nco)
USE crunchtype
USE concentration
USE params

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nspec
INTEGER(I4B), INTENT(IN)                                   :: nexchange
INTEGER(I4B), INTENT(IN)                                   :: nexch_sec
INTEGER(I4B), INTENT(IN)                                   :: nco

!  Internal variables

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: nex
INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: ik

REAL(DP)                                                   :: activity
REAL(DP)                                                   :: sum
REAL(DP)                                                   :: sion_tmp

sum = 0.0D0
DO ik = 1,ncomp+nspec
  sum = sum + sptmp10(ik)*chg(ik)*chg(ik)
END DO
sion_tmp = 0.50D0*sum

IF (iexc == 1) THEN        ! Gaines-Thomas convention
  DO ix = 1,nexchange
    sumactivity(ix) = 0.0
  END DO
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    sum = 0.0
    DO i = 1,ncomp
      sum = sum + muexc(nex,i)*(sptmp(i)+gamtmp(i))
    END DO
    sum = sum + muexc(nex,ix+ncomp)*spextmp(ix)
    activity = EXP(-keqexc(nex) + sum + bfit(nex)*sion_tmp )
    aexch(nex) = activity
    spextmp10(nex+nexchange) = activity*totexch(ix,nco)/muexc(nex,ix+ncomp)
    sumactivity(ix) = sumactivity(ix) + aexch(nex)
  END DO
ELSE IF (iexc == 2) THEN   ! Vanselow convention
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    sum = 0.0
    DO i = 1,ncomp
      sum = sum + muexc(nex,i)*(sptmp(i)+gamtmp(i))
    END DO
    sum = sum + muexc(nex,ix+ncomp)*spextmp(ix)
    aexch(nex) = EXP(-keqexc(nex) + sum + bfit(nex)*sion_tmp )
  END DO
  DO ix = 1,nexchange
    wt_aexch(ix) = 0.0
    sumactivity(ix) = 0.0
  END DO
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    wt_aexch(ix) = wt_aexch(ix) + muexc(nex,ix+ncomp)*aexch(nex)
    sumactivity(ix) = sumactivity(ix) + aexch(nex)
  END DO
  DO ix = 1,nexchange
    tec(ix) = totexch(ix,nco)/wt_aexch(ix)
  END DO
  DO nex=1,nexch_sec
    ix = ixlink(nex)
    spextmp10(nex+nexchange) = aexch(nex)*tec(ix)
  END DO
ELSE IF (iexc == 3) THEN   ! Gapon convention
ELSE
  CONTINUE
END IF

RETURN
END SUBROUTINE exchange_init
!  **************************************************************
