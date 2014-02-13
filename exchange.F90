!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:56:28
 
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

SUBROUTINE exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE transport
USE medium
USE temperature

IMPLICIT NONE
!fp! auto_par_loops=0;

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
INTEGER(I4B)                                               :: i


REAL(DP)                                                   :: activity
REAL(DP)                                                   :: sum
REAL(DP)                                                   :: check
REAL(DP)                                                   :: sumtemp
REAL(DP)                                                   :: exchangetemp

IF (iexc == 1) THEN        ! Gaines-Thomas convention
  sumactivity = 0.0
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    exchangetemp = exchangesites(ix,jx,jy,jz)

    IF (ikh2o /= 0) THEN
      sum = 0.0
      DO i = 1,ncomp
        IF (ulab(i) == 'H2O') THEN
          sum = sum + muexc(nex,i)*(gam(i,jx,jy,jz))
        ELSE
          sum = sum + muexc(nex,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
        END IF
      END DO

    ELSE

      sum = 0.0
      DO i = 1,ncomp
        sum = sum + muexc(nex,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
      END DO

    END IF

    sum = sum + muexc(nex,ix+ncomp)*spex(ix,jx,jy,jz)
    activity = EXP(-keqexc(nex) + sum + bfit(nex)*sion(jx,jy,jz) )
    aexch(nex) = activity
    spex10(nex+nexchange,jx,jy,jz) = activity*exchangetemp/(muexc(nex,ix+ncomp))  
    sumactivity(ix) = sumactivity(ix) + aexch(nex)
    spex(nex+nexchange,jx,jy,jz) = LOG(spex10(nex+nexchange,jx,jy,jz))
  END DO
ELSE IF (iexc == 2) THEN   ! Vanselow convention
  DO nex = 1,nexch_sec
    ix = ixlink(nex)

    IF (ikh2o /= 0) THEN
      sum = 0.0
      DO i = 1,ncomp
        IF (ulab(i) == 'H2O') THEN
          sum = sum + muexc(nex,i)*(gam(i,jx,jy,jz))
        ELSE
          sum = sum + muexc(nex,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
        END IF
      END DO

    ELSE

      sum = 0.0
      DO i = 1,ncomp
        sum = sum + muexc(nex,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
      END DO

    END IF
    sum = sum + muexc(nex,ix+ncomp)*spex(ix,jx,jy,jz)
    aexch(nex) = EXP(-keqexc(nex) + sum + bfit(nex)*sion(jx,jy,jz))
  END DO
  wt_aexch = 0.0
  sumactivity = 0.0
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    wt_aexch(ix) = wt_aexch(ix) + muexc(nex,ix+ncomp)*aexch(nex)
    sumactivity(ix) = sumactivity(ix) + aexch(nex)
  END DO
  DO ix = 1,nexchange
    tec(ix) = exchangesites(ix,jx,jy,jz)/wt_aexch(ix)
  END DO
  DO nex=1,nexch_sec
    ix = ixlink(nex)
    spex10(nex+nexchange,jx,jy,jz) = aexch(nex)*tec(ix)
    spex(nex+nexchange,jx,jy,jz) = LOG(spex10(nex+nexchange,jx,jy,jz))
  END DO
ELSE IF (iexc == 3) THEN   ! Gapon convention--same as Gaines-Thomas, but with different reaction stoichiometry
ELSE
  CONTINUE
END IF

RETURN
END SUBROUTINE exchange
!  **************************************************************
