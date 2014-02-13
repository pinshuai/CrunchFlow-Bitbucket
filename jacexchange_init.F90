
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:59:26
 
!******************        GIMRT98     ************************
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

SUBROUTINE jacexchange_init(ncomp,nexchange,nexch_sec,neqn,nco)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                           :: neqn
INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nexch_sec
INTEGER(I4B), INTENT(IN)                           :: nco

!  Internal variables and arrays

REAL(DP)                                           :: sum
REAL(DP)                                           :: sum1
REAL(DP)                                           :: sum2

INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: i2
INTEGER(I4B)                                       :: nex
INTEGER(I4B)                                       :: ix
INTEGER(I4B)                                       :: ixcheck

fexch = 0.0

IF (iexc == 1 .OR. iexc == 3) THEN
  
  DO i = 1,ncomp
    IF (itype(i,nco) == 1) THEN
      DO i2 = 1,ncomp+nexchange
        sum = 0.0
        IF (equilibrate(i,nco)) THEN
          DO nex = 1,nexch_sec
            sum = sum + muexc(nex,i)*muexc(nex,i2)* spextmp10(nexchange+nex)
          END DO
        END IF
        fexch(i,i2) = sum
      END DO
    ELSE
      DO i2 = 1,ncomp+nexchange
        fexch(i,i2) = 0.0
      END DO
    END IF
  END DO
  
!  Part of "fexch" referring to total exchange concentrations (not aqueous components)
  
  DO ix = 1,nexchange
    DO i2 = 1,ncomp+nexchange
      sum2 = 0.0
      DO nex = 1,nexch_sec
        ixcheck = ixlink(nex)
        IF (ixcheck == ix) THEN
          sum2 = sum2 + muexc(nex,i2)*aexch(nex)
        END IF
        sum2 = sum2 + muexc(nex,i2)*aexch(nex)
      END DO
      fexch(ix+ncomp,i2) = sum2
    END DO
  END DO
  
ELSE    !                  Vanselow convention
  
  DO ix = 1,nexchange
    DO i2 = 1,ncomp+nexchange
      sum1 = 0.0
      sum2 = 0.0
      DO nex = 1,nexch_sec
        ixcheck = ixlink(nex)
        sum1 = sum1 + muexc(nex,ix+ncomp)*muexc(nex,i2)* aexch(nex)
        IF (ixcheck == ix) THEN
          sum2 = sum2 + muexc(nex,i2)*aexch(nex)
        END IF
      END DO
      fweight(ix+ncomp,i2) = sum1
      fexch(ix+ncomp,i2) = sum2
    END DO
  END DO
  
  
  DO i = 1,ncomp
    IF (itype(i,nco) == 1) THEN
      DO i2 = 1,ncomp+nexchange
        sum1 = 0.0
        sum2 = 0.0
        IF (equilibrate(i,nco)) THEN
          DO nex = 1,nexch_sec
            ix = ixlink(nex)
            sum1 = sum1 + muexc(nex,i)*muexc(nex,i2)* aexch(nex)*tec(ix)
            sum2 = sum2 - totexch(ix,nco)*muexc(nex,i)*aexch(nex)/  &
                (wt_aexch(ix)*wt_aexch(ix))*fweight(ix+ncomp,i2)
          END DO
        END IF
        fexch(i,i2) = sum1 + sum2
      END DO
    ELSE
      DO i2 = 1,ncomp+nexchange
        fexch(i,i2) = 0.0
      END DO
    END IF
  END DO
  
END IF

RETURN
END SUBROUTINE jacexchange_init
!***********************************************************
